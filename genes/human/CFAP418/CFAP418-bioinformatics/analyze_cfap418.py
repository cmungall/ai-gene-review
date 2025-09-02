#!/usr/bin/env python3
"""
Bioinformatics analysis of CFAP418 (Q6ZT21)
Focus: Domain architecture, coiled-coil predictions, and ciliary localization signals
"""

import requests
from Bio import SeqIO
from io import StringIO
import json
import re

def fetch_uniprot_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        record = SeqIO.read(StringIO(response.text), "fasta")
        return str(record.seq)
    return None

def fetch_uniprot_features(uniprot_id):
    """Fetch protein features from UniProt JSON"""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return None

def analyze_coiled_coils(sequence):
    """Simple coiled-coil prediction using heptad repeat patterns"""
    print("\n=== Coiled-Coil Analysis ===")
    
    # Simplified coiled-coil detection based on hydrophobic patterns
    # Real analysis would use COILS, Paircoil2, or similar
    
    # Look for regions enriched in L, I, V, M at positions a and d of heptad
    window_size = 28  # 4 heptads
    threshold = 0.35  # threshold for hydrophobic content
    
    hydrophobic = set('LIVMF')
    potential_cc = []
    
    for i in range(0, len(sequence) - window_size, 7):
        window = sequence[i:i+window_size]
        # Check positions a (0, 7, 14, 21) and d (3, 10, 17, 24)
        heptad_hydrophobic = 0
        for j in [0, 3, 7, 10, 14, 17, 21, 24]:
            if j < len(window) and window[j] in hydrophobic:
                heptad_hydrophobic += 1
        
        ratio = heptad_hydrophobic / 8.0
        if ratio >= threshold:
            potential_cc.append((i, i+window_size, ratio))
    
    if potential_cc:
        print(f"Potential coiled-coil regions detected: {len(potential_cc)}")
        # Merge overlapping regions
        merged = []
        for start, end, score in potential_cc:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]), max(score, merged[-1][2]))
            else:
                merged.append((start, end, score))
        
        for start, end, score in merged[:5]:  # Show top 5
            print(f"  Position {start}-{end} (score: {score:.2f})")
            print(f"    Sequence: {sequence[start:start+20]}...")
    else:
        print("No strong coiled-coil regions detected with simple analysis")
    
    return len(potential_cc) > 0

def analyze_domains_interpro(sequence):
    """Analyze domains using sequence features"""
    print("\n=== Domain Architecture Analysis ===")
    
    # Check for WD40 repeats (simplified pattern matching)
    wd40_pattern = r'[WF].{5,7}[LIVMF].{10,15}[DEN]'
    wd40_matches = list(re.finditer(wd40_pattern, sequence))
    
    if wd40_matches:
        print(f"Potential WD40-like repeats: {len(wd40_matches)}")
        for i, match in enumerate(wd40_matches[:3], 1):
            print(f"  Repeat {i}: Position {match.start()}-{match.end()}")
    
    # Check for tetratricopeptide repeats (TPR) - simplified
    tpr_pattern = r'[ALP].{7}[ALP].{7}[FYL]'
    tpr_matches = list(re.finditer(tpr_pattern, sequence))
    
    if tpr_matches:
        print(f"Potential TPR-like motifs: {len(tpr_matches)}")
    
    # Check for EF-hand calcium binding motifs
    ef_hand_pattern = r'D.{3}[DN].{3}[DN]'
    ef_matches = list(re.finditer(ef_hand_pattern, sequence))
    
    if ef_matches:
        print(f"Potential EF-hand-like motifs: {len(ef_matches)}")
    
    print("\nNote: These are simplified pattern matches.")
    print("For accurate domain prediction, use InterProScan or Pfam HMMs")

def analyze_ciliary_targeting(sequence):
    """Check for ciliary targeting signals"""
    print("\n=== Ciliary Targeting Signal Analysis ===")
    
    # Check for RVxP ciliary targeting signal
    rvxp_pattern = r'RV.P'
    rvxp_matches = list(re.finditer(rvxp_pattern, sequence))
    
    if rvxp_matches:
        print(f"RVxP ciliary targeting signal found: {len(rvxp_matches)} instance(s)")
        for match in rvxp_matches:
            pos = match.start()
            print(f"  Position {pos}: {sequence[pos:pos+4]}")
    
    # Check for VxPx ciliary targeting signal  
    vxpx_pattern = r'V.P.'
    vxpx_matches = list(re.finditer(vxpx_pattern, sequence))
    
    if vxpx_matches:
        print(f"VxPx-like motifs found: {len(vxpx_matches)} instance(s)")
    
    # Check for polyglutamylation sites (E-rich regions)
    e_rich_pattern = r'E{3,}'
    e_rich = list(re.finditer(e_rich_pattern, sequence))
    
    if e_rich:
        print(f"Polyglutamylation target sites (E-rich): {len(e_rich)}")
        for match in e_rich[:3]:
            print(f"  Position {match.start()}: {match.group()}")
    
    if not (rvxp_matches or vxpx_matches):
        print("No canonical ciliary targeting signals detected")

def analyze_protein_properties(sequence):
    """Analyze basic protein properties"""
    print("\n=== Protein Properties ===")
    
    length = len(sequence)
    print(f"Length: {length} amino acids")
    
    # Calculate molecular weight (approximate)
    mw_table = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    mw = sum(mw_table.get(aa, 110) for aa in sequence) / 1000
    print(f"Molecular weight: ~{mw:.1f} kDa")
    
    # Calculate pI (simplified)
    basic = sum(1 for aa in sequence if aa in 'RKH')
    acidic = sum(1 for aa in sequence if aa in 'DE')
    charge_ratio = basic / (acidic + 1)
    
    print(f"Basic residues (R,K,H): {basic}")
    print(f"Acidic residues (D,E): {acidic}")
    
    if charge_ratio > 1.1:
        print("Theoretical pI: Basic (>7.0)")
    elif charge_ratio < 0.9:
        print("Theoretical pI: Acidic (<7.0)")
    else:
        print("Theoretical pI: Near neutral (~7.0)")
    
    # Composition analysis
    proline_count = sequence.count('P')
    glycine_count = sequence.count('G')
    serine_count = sequence.count('S')
    
    print(f"\nComposition features:")
    print(f"Proline content: {proline_count/length*100:.1f}%")
    print(f"Glycine content: {glycine_count/length*100:.1f}%")
    print(f"Serine content: {serine_count/length*100:.1f}%")

def check_uniprot_annotations(uniprot_data):
    """Extract relevant annotations from UniProt"""
    print("\n=== UniProt Annotations ===")
    
    if not uniprot_data:
        print("Unable to fetch UniProt data")
        return
    
    # Check for subcellular location
    if 'comments' in uniprot_data:
        for comment in uniprot_data['comments']:
            if comment.get('type') == 'SUBCELLULAR_LOCATION':
                print("Subcellular location annotation found:")
                if 'locations' in comment:
                    for loc in comment['locations']:
                        if 'location' in loc:
                            print(f"  - {loc['location'].get('value', 'Unknown')}")
    
    # Check for features
    if 'features' in uniprot_data:
        feature_types = {}
        for feature in uniprot_data['features']:
            ftype = feature.get('type', 'Unknown')
            feature_types[ftype] = feature_types.get(ftype, 0) + 1
        
        print("\nFeature summary:")
        for ftype, count in sorted(feature_types.items()):
            if ftype in ['Domain', 'Region', 'Coiled coil', 'Compositional bias', 'Modified residue']:
                print(f"  {ftype}: {count}")

def main():
    print("=" * 60)
    print("CFAP418 Bioinformatics Analysis")
    print("UniProt: Q6ZT21 (Human)")
    print("=" * 60)
    
    # Fetch sequence
    uniprot_id = "Q6ZT21"
    sequence = fetch_uniprot_sequence(uniprot_id)
    if not sequence:
        print("Failed to fetch sequence from UniProt")
        return
    
    print(f"\nSequence retrieved: {len(sequence)} amino acids")
    
    # Fetch UniProt annotations
    uniprot_data = fetch_uniprot_features(uniprot_id)
    
    # Perform analyses
    analyze_protein_properties(sequence)
    check_uniprot_annotations(uniprot_data)
    has_cc = analyze_coiled_coils(sequence)
    analyze_domains_interpro(sequence)
    analyze_ciliary_targeting(sequence)
    
    # Conclusions
    print("\n" + "=" * 60)
    print("CONCLUSIONS")
    print("=" * 60)
    
    print("\n1. PROTEIN CHARACTERISTICS:")
    print(f"   - Large protein ({len(sequence)} aa)")
    print("   - Contains potential coiled-coil regions" if has_cc else "   - Limited coiled-coil potential")
    
    print("\n2. FUNCTIONAL PREDICTIONS:")
    print("   - Likely ciliary/flagellar protein (CFAP family)")
    print("   - May contain protein-protein interaction domains")
    print("   - Potential structural role based on size and composition")
    
    print("\n3. LIMITATIONS:")
    print("   - Pattern-based predictions are approximate")
    print("   - Full domain analysis requires HMM-based tools (InterProScan)")
    print("   - Experimental validation needed for functional assignments")
    
    print("\n4. RECOMMENDATIONS:")
    print("   - Run InterProScan for comprehensive domain analysis")
    print("   - Use COILS/Paircoil2 for accurate coiled-coil prediction")
    print("   - Check ciliary proteomics databases for validation")

if __name__ == "__main__":
    main()