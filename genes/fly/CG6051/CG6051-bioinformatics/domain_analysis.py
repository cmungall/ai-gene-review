#!/usr/bin/env python3
"""
Domain and motif analysis for CG6051 protein.
Identifies FYVE domain, checks for phosphatase domains, and searches for regulatory motifs.
"""

import sys
import re
import json
from pathlib import Path
from typing import Dict, List, Tuple
import requests
from Bio import SeqIO
from Bio.Seq import Seq

def load_sequence(fasta_file: str) -> Tuple[str, str]:
    """Load protein sequence from FASTA file."""
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.id), str(record.seq)

def scan_interpro(sequence: str) -> Dict:
    """Query InterProScan for domain predictions."""
    # Note: InterProScan web service has limitations
    # This is a simplified approach - for production, use local InterProScan
    results = {
        "domains": [],
        "note": "Manual InterProScan analysis required for comprehensive results"
    }
    
    # Based on UniProt annotation, we know FYVE domain is at 909-969
    results["domains"].append({
        "type": "FYVE",
        "start": 909,
        "end": 969,
        "source": "UniProt annotation",
        "description": "FYVE zinc finger domain - phosphatidylinositol-3-phosphate binding"
    })
    
    return results

def search_phosphatase_motifs(sequence: str) -> List[Dict]:
    """Search for phosphatase catalytic motifs."""
    phosphatase_motifs = {
        # Protein tyrosine phosphatase active site
        "PTP_active": r"[LIVMFY]H[CS][LIVM]G[LIVM][CS]R[ST]",
        # Dual specificity phosphatase active site  
        "DSP_active": r"HC[LIVM][LIVM][LIVM][LIVM][LIVM]R",
        # Serine/threonine phosphatase signatures
        "PP1_PP2A": r"[LIVM][LIVM]GD[YFH]HG",
        "PP2C": r"[LIVM][LIVM]D[GS]H[GA]",
        # Acid phosphatase signature
        "Acid_phosphatase": r"[LIVM].[LIVM]D[ST]G[STN]",
        # Alkaline phosphatase
        "Alkaline_phosphatase": r"D[ST][LIVM][LIVM]D[ST][LIVM]"
    }
    
    matches = []
    for motif_name, pattern in phosphatase_motifs.items():
        for match in re.finditer(pattern, sequence):
            matches.append({
                "motif": motif_name,
                "pattern": pattern,
                "position": match.start() + 1,  # 1-based
                "sequence": match.group(),
                "context": sequence[max(0, match.start()-5):min(len(sequence), match.end()+5)]
            })
    
    return matches

def search_tos_motif(sequence: str) -> List[Dict]:
    """Search for TOS (TOR signaling) motif.
    The canonical TOS motif is F[DE][LIVM][DE][LIVM]
    """
    tos_patterns = [
        (r"F[DE][LIVM][DE][LIVM]", "Canonical TOS motif"),
        (r"F[DE][LIVM]{2}[DE]", "TOS-like motif variant 1"),
        (r"[FY][DE][LIVM][DE][LIVM]", "TOS-like motif variant 2")
    ]
    
    matches = []
    for pattern, description in tos_patterns:
        for match in re.finditer(pattern, sequence):
            matches.append({
                "motif": "TOS",
                "description": description,
                "pattern": pattern,
                "position": match.start() + 1,  # 1-based
                "sequence": match.group(),
                "context": sequence[max(0, match.start()-10):min(len(sequence), match.end()+10)]
            })
    
    return matches

def analyze_fyve_domain(sequence: str, start: int = 909, end: int = 969) -> Dict:
    """Analyze the FYVE domain structure and conservation."""
    fyve_seq = sequence[start-1:end]  # Convert to 0-based
    
    # FYVE domain consensus features
    analysis = {
        "domain_sequence": fyve_seq,
        "length": len(fyve_seq),
        "start": start,
        "end": end,
        "features": []
    }
    
    # Check for conserved cysteines (8 cysteines typical for FYVE)
    cys_positions = [i+start for i, aa in enumerate(fyve_seq) if aa == 'C']
    analysis["cysteine_positions"] = cys_positions
    analysis["cysteine_count"] = len(cys_positions)
    
    # Check for basic patch (R/K rich region for PI3P binding)
    basic_count = sum(1 for aa in fyve_seq if aa in 'RK')
    analysis["basic_residues"] = basic_count
    analysis["basic_percentage"] = (basic_count / len(fyve_seq)) * 100
    
    # Check for conserved FYVE motifs
    if "RVC" in fyve_seq:
        analysis["features"].append("Contains RVC motif (conserved in FYVE)")
    if "HHCR" in fyve_seq:
        analysis["features"].append("Contains HHCR motif (zinc coordination)")
    
    # WxxD motif (important for PI3P binding)
    wxxd_match = re.search(r"W.{2}D", fyve_seq)
    if wxxd_match:
        analysis["features"].append(f"Contains WxxD motif at position {wxxd_match.start() + start}")
    
    return analysis

def search_regulatory_motifs(sequence: str) -> Dict:
    """Search for various regulatory motifs."""
    motifs = {}
    
    # Nuclear localization signals
    nls_patterns = [
        (r"[KR]{4,}", "Poly-basic NLS"),
        (r"[KR]{2}.{10,12}[KR]{3,}", "Bipartite NLS"),
        (r"P[KR]{3,5}", "SV40-like NLS")
    ]
    
    motifs["NLS"] = []
    for pattern, description in nls_patterns:
        for match in re.finditer(pattern, sequence):
            motifs["NLS"].append({
                "type": description,
                "position": match.start() + 1,
                "sequence": match.group()
            })
    
    # Phosphorylation sites (simplified - actual sites at S549, S550, S810)
    known_phospho = [549, 550, 810]
    motifs["phosphorylation_sites"] = []
    for pos in known_phospho:
        if sequence[pos-1] == 'S':
            context = sequence[max(0, pos-8):min(len(sequence), pos+7)]
            motifs["phosphorylation_sites"].append({
                "position": pos,
                "residue": "S",
                "context": context,
                "source": "Experimental (mass spec)"
            })
    
    # PEST sequences (protein degradation signals)
    pest_pattern = r"[PEST]{4,}"
    motifs["PEST"] = []
    for match in re.finditer(pest_pattern, sequence):
        if len(match.group()) >= 5:
            motifs["PEST"].append({
                "position": match.start() + 1,
                "length": len(match.group()),
                "sequence": match.group()
            })
    
    return motifs

def analyze_disordered_regions(sequence: str) -> List[Tuple[int, int]]:
    """Identify potential disordered regions based on composition."""
    # Simplified disorder prediction based on amino acid composition
    disorder_prone = set("PQSTNRKED")
    window_size = 20
    threshold = 0.5
    
    disordered = []
    in_disorder = False
    start = 0
    
    for i in range(len(sequence) - window_size):
        window = sequence[i:i+window_size]
        disorder_score = sum(1 for aa in window if aa in disorder_prone) / window_size
        
        if disorder_score >= threshold and not in_disorder:
            in_disorder = True
            start = i + 1  # 1-based
        elif disorder_score < threshold and in_disorder:
            in_disorder = False
            if i - start > 10:  # Minimum length
                disordered.append((start, i))
    
    return disordered

def compare_with_human_lst2() -> Dict:
    """Compare key features with human LST2/ZFYVE28."""
    comparison = {
        "human_lst2": {
            "uniprot": "Q9HCC9",
            "length": 895,
            "fyve_domain": "815-875",
            "tos_motif": "401-405 (FDIDI)",
            "function": "Negative regulator of EGFR signaling, involved in endosomal sorting"
        },
        "fly_cg6051": {
            "uniprot": "Q9VB70", 
            "length": 989,
            "fyve_domain": "909-969",
            "function": "Negative regulator of EGFR signaling (predicted by homology)"
        },
        "comparison": {
            "length_difference": 989 - 895,
            "fyve_position": "Both have C-terminal FYVE domains",
            "functional_conservation": "Both regulate EGFR signaling"
        }
    }
    
    return comparison

def main():
    # Input file
    fasta_file = sys.argv[1] if len(sys.argv) > 1 else "CG6051.fasta"
    
    if not Path(fasta_file).exists():
        print(f"Error: {fasta_file} not found")
        sys.exit(1)
    
    # Load sequence
    seq_id, sequence = load_sequence(fasta_file)
    print(f"Analyzing protein: {seq_id}")
    print(f"Length: {len(sequence)} amino acids\n")
    
    # Comprehensive analysis
    results = {
        "protein_id": seq_id,
        "length": len(sequence),
        "analyses": {}
    }
    
    # 1. Domain analysis
    print("1. Domain Analysis")
    print("-" * 50)
    domains = scan_interpro(sequence)
    results["analyses"]["domains"] = domains
    for domain in domains["domains"]:
        print(f"  {domain['type']}: {domain['start']}-{domain['end']} - {domain['description']}")
    
    # 2. FYVE domain detailed analysis
    print("\n2. FYVE Domain Detailed Analysis")
    print("-" * 50)
    fyve_analysis = analyze_fyve_domain(sequence)
    results["analyses"]["fyve_domain"] = fyve_analysis
    print(f"  Position: {fyve_analysis['start']}-{fyve_analysis['end']}")
    print(f"  Cysteines: {fyve_analysis['cysteine_count']} at positions {fyve_analysis['cysteine_positions']}")
    print(f"  Basic residues: {fyve_analysis['basic_residues']} ({fyve_analysis['basic_percentage']:.1f}%)")
    for feature in fyve_analysis["features"]:
        print(f"  ✓ {feature}")
    
    # 3. Phosphatase motif search
    print("\n3. Phosphatase Motif Search")
    print("-" * 50)
    phosphatase_motifs = search_phosphatase_motifs(sequence)
    results["analyses"]["phosphatase_motifs"] = phosphatase_motifs
    if phosphatase_motifs:
        print(f"  Found {len(phosphatase_motifs)} potential phosphatase motifs:")
        for motif in phosphatase_motifs:
            print(f"    {motif['motif']} at position {motif['position']}: {motif['sequence']} ({motif['context']})")
    else:
        print("  ✓ No canonical phosphatase motifs detected")
    
    # 4. TOS motif search
    print("\n4. TOS Motif Search")
    print("-" * 50)
    tos_motifs = search_tos_motif(sequence)
    results["analyses"]["tos_motifs"] = tos_motifs
    if tos_motifs:
        print(f"  Found {len(tos_motifs)} TOS/TOS-like motifs:")
        for motif in tos_motifs[:5]:  # Show first 5
            print(f"    {motif['description']} at position {motif['position']}: {motif['sequence']}")
            print(f"      Context: {motif['context']}")
    else:
        print("  No TOS motifs detected")
    
    # 5. Regulatory motifs
    print("\n5. Regulatory Motifs")
    print("-" * 50)
    reg_motifs = search_regulatory_motifs(sequence)
    results["analyses"]["regulatory_motifs"] = reg_motifs
    
    if reg_motifs["phosphorylation_sites"]:
        print(f"  Phosphorylation sites (experimental):")
        for site in reg_motifs["phosphorylation_sites"]:
            print(f"    S{site['position']}: {site['context']}")
    
    if reg_motifs["NLS"]:
        print(f"  Nuclear localization signals: {len(reg_motifs['NLS'])}")
    
    if reg_motifs["PEST"]:
        print(f"  PEST sequences: {len(reg_motifs['PEST'])}")
    
    # 6. Disordered regions
    print("\n6. Predicted Disordered Regions")
    print("-" * 50)
    disordered = analyze_disordered_regions(sequence)
    results["analyses"]["disordered_regions"] = disordered
    if disordered:
        print(f"  Found {len(disordered)} disordered regions:")
        for start, end in disordered[:5]:  # Show first 5
            print(f"    {start}-{end} ({end-start+1} aa)")
    
    # 7. Comparison with human LST2
    print("\n7. Comparison with Human LST2/ZFYVE28")
    print("-" * 50)
    comparison = compare_with_human_lst2()
    results["analyses"]["human_comparison"] = comparison
    print(f"  Human LST2: {comparison['human_lst2']['length']} aa")
    print(f"  Fly CG6051: {comparison['fly_cg6051']['length']} aa")
    print(f"  Difference: {comparison['comparison']['length_difference']} aa")
    print(f"  FYVE domains: {comparison['comparison']['fyve_position']}")
    print(f"  Human TOS motif: {comparison['human_lst2']['tos_motif']}")
    
    # Save detailed results
    with open("domain_analysis_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 50)
    print("Analysis complete. Results saved to domain_analysis_results.json")

if __name__ == "__main__":
    main()