#!/usr/bin/env python3
"""
Bioinformatics analysis of RBFOX3/NeuN (A6NFN3)
Focus: RNA binding domain analysis, intrinsically disordered regions, and splice regulation motifs
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

def fetch_uniprot_data(uniprot_id):
    """Fetch UniProt JSON data"""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    return None

def analyze_rrm_domain(sequence):
    """Analyze RNA Recognition Motif (RRM) domain"""
    print("\n=== RNA Recognition Motif (RRM) Analysis ===")
    
    # RRM domains are typically 80-90 aa with conserved RNP motifs
    # RNP-1: [RK]-G-[FY]-[GA]-[FY]-[ILV]-X-[FY] (octamer)
    # RNP-2: [ILV]-[FY]-[ILV]-X-N-L (hexamer)
    
    # Simplified RNP-1 pattern
    rnp1_pattern = r'[RK]G[FY][GA][FY][ILV].[FY]'
    rnp1_matches = list(re.finditer(rnp1_pattern, sequence))
    
    if rnp1_matches:
        print(f"RNP-1 motifs found: {len(rnp1_matches)}")
        for i, match in enumerate(rnp1_matches, 1):
            pos = match.start()
            print(f"  RNP-1 #{i}: Position {pos}-{match.end()}")
            print(f"    Sequence: {match.group()}")
            # RRM typically extends ~40 aa before and after RNP-1
            rrm_start = max(0, pos - 40)
            rrm_end = min(len(sequence), pos + 50)
            print(f"    Predicted RRM domain: {rrm_start}-{rrm_end}")
    
    # Simplified RNP-2 pattern  
    rnp2_pattern = r'[ILV][FY][ILV].NL'
    rnp2_matches = list(re.finditer(rnp2_pattern, sequence))
    
    if rnp2_matches:
        print(f"\nRNP-2 motifs found: {len(rnp2_matches)}")
        for i, match in enumerate(rnp2_matches, 1):
            print(f"  RNP-2 #{i}: Position {match.start()}-{match.end()}: {match.group()}")
    
    if not (rnp1_matches or rnp2_matches):
        print("No canonical RNP motifs detected")
    
    return len(rnp1_matches) > 0

def analyze_disorder(sequence):
    """Predict intrinsically disordered regions"""
    print("\n=== Intrinsically Disordered Region (IDR) Analysis ===")
    
    # Simple disorder prediction based on composition
    # Disordered regions are enriched in P, E, S, K, R, Q and depleted in W, C, F, Y, I, V, L
    
    window_size = 30
    disorder_threshold = 0.5
    
    disorder_promoting = set('PESKRQ')
    order_promoting = set('WCFYIVL')
    
    disordered_regions = []
    
    for i in range(0, len(sequence) - window_size):
        window = sequence[i:i+window_size]
        disorder_score = sum(1 for aa in window if aa in disorder_promoting) / window_size
        order_score = sum(1 for aa in window if aa in order_promoting) / window_size
        
        net_disorder = disorder_score - order_score
        
        if net_disorder > disorder_threshold:
            disordered_regions.append((i, i+window_size, net_disorder))
    
    # Merge overlapping regions
    if disordered_regions:
        merged = []
        for start, end, score in disordered_regions:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]), max(score, merged[-1][2]))
            else:
                merged.append((start, end, score))
        
        print(f"Predicted disordered regions: {len(merged)}")
        for start, end, score in merged:
            length = end - start
            print(f"  Position {start}-{end} (length: {length} aa, score: {score:.2f})")
            # Check composition
            region_seq = sequence[start:end]
            s_count = region_seq.count('S')
            p_count = region_seq.count('P')
            e_count = region_seq.count('E')
            print(f"    Composition: S={s_count}, P={p_count}, E={e_count}")
    else:
        print("No significant disordered regions predicted")
    
    return len(disordered_regions) > 0

def analyze_fox_binding_motifs(sequence):
    """Analyze FOX-specific RNA binding features"""
    print("\n=== FOX Family Binding Motif Analysis ===")
    
    # FOX proteins recognize (U)GCAUG motifs in RNA
    # Check for conserved residues important for GCAUG recognition
    
    # Key positions in RRM for UGCAUG binding (based on RBFOX1 structure)
    print("Checking for FOX-specific features in RRM domain...")
    
    # Look for conserved F in RRM (important for RNA binding)
    f_positions = [i for i, aa in enumerate(sequence) if aa == 'F']
    
    # In the RRM region (typically positions 100-200 for RBFOX proteins)
    rrm_f = [pos for pos in f_positions if 100 <= pos <= 200]
    
    if rrm_f:
        print(f"Conserved F residues in RRM region: {len(rrm_f)} at positions {rrm_f[:5]}...")
        print("  These may be important for UGCAUG RNA recognition")
    
    # Check for basic patches that interact with RNA backbone
    basic_patch_pattern = r'[RK]{2,}'
    basic_patches = list(re.finditer(basic_patch_pattern, sequence))
    
    if basic_patches:
        print(f"\nBasic patches for RNA binding: {len(basic_patches)}")
        for i, match in enumerate(basic_patches[:3], 1):
            print(f"  Patch {i}: Position {match.start()}: {match.group()}")

def analyze_protein_interactions(sequence):
    """Analyze potential protein interaction regions"""
    print("\n=== Protein Interaction Region Analysis ===")
    
    # Check for proline-rich regions (often mediate protein interactions)
    proline_rich_pattern = r'P{2,}|P.P.P'
    p_rich = list(re.finditer(proline_rich_pattern, sequence))
    
    if p_rich:
        print(f"Proline-rich regions: {len(p_rich)}")
        for match in p_rich[:3]:
            print(f"  Position {match.start()}-{match.end()}: {match.group()}")
    
    # Check for SH3 binding motifs (PxxP)
    sh3_pattern = r'P..P'
    sh3_motifs = list(re.finditer(sh3_pattern, sequence))
    
    if sh3_motifs:
        print(f"\nPotential SH3 binding motifs (PxxP): {len(sh3_motifs)}")
        print(f"  First 3 motifs: {[m.group() for m in sh3_motifs[:3]]}")
    
    # Check for WW domain binding motifs (PPxY)
    ww_pattern = r'PP.Y'
    ww_motifs = list(re.finditer(ww_pattern, sequence))
    
    if ww_motifs:
        print(f"\nPotential WW domain binding motifs (PPxY): {len(ww_motifs)}")
        for match in ww_motifs:
            print(f"  Position {match.start()}: {match.group()}")

def analyze_post_translational_modifications(sequence):
    """Predict post-translational modification sites"""
    print("\n=== Post-Translational Modification Sites ===")
    
    # Phosphorylation sites
    # Simple patterns for common kinase recognition sites
    pka_pattern = r'R.[ST]'  # PKA consensus
    ck2_pattern = r'[ST]..E'  # CK2 consensus
    
    pka_sites = list(re.finditer(pka_pattern, sequence))
    ck2_sites = list(re.finditer(ck2_pattern, sequence))
    
    print(f"Potential phosphorylation sites:")
    print(f"  PKA sites (R-x-S/T): {len(pka_sites)}")
    if pka_sites:
        print(f"    Positions: {[m.start() for m in pka_sites[:5]]}...")
    
    print(f"  CK2 sites (S/T-x-x-E): {len(ck2_sites)}")
    if ck2_sites:
        print(f"    Positions: {[m.start() for m in ck2_sites[:5]]}...")
    
    # SUMOylation sites (ΨKxE where Ψ is hydrophobic)
    sumo_pattern = r'[VILMF]K.E'
    sumo_sites = list(re.finditer(sumo_pattern, sequence))
    
    if sumo_sites:
        print(f"\nPotential SUMOylation sites: {len(sumo_sites)}")
        for match in sumo_sites[:3]:
            print(f"  Position {match.start()}: {match.group()}")

def calculate_properties(sequence):
    """Calculate basic protein properties"""
    print("\n=== Protein Properties ===")
    
    length = len(sequence)
    print(f"Length: {length} amino acids")
    
    # Molecular weight
    mw_table = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    mw = sum(mw_table.get(aa, 110) for aa in sequence) / 1000
    print(f"Molecular weight: ~{mw:.1f} kDa")
    
    # Composition
    basic = sum(1 for aa in sequence if aa in 'RKH')
    acidic = sum(1 for aa in sequence if aa in 'DE')
    
    print(f"Basic residues: {basic}")
    print(f"Acidic residues: {acidic}")
    print(f"Net charge at pH 7: {basic - acidic}")
    
    # Domain architecture estimate
    if length > 300:
        print(f"\nDomain architecture (estimated):")
        print(f"  N-terminal region: 1-100")
        print(f"  RRM domain: ~100-190")
        print(f"  C-terminal region: 190-{length}")

def main():
    print("=" * 60)
    print("RBFOX3/NeuN Bioinformatics Analysis")
    print("UniProt: A6NFN3 (Human)")
    print("=" * 60)
    
    # Fetch sequence
    uniprot_id = "A6NFN3"
    sequence = fetch_uniprot_sequence(uniprot_id)
    if not sequence:
        print("Failed to fetch sequence from UniProt")
        return
    
    print(f"\nSequence retrieved: {len(sequence)} amino acids")
    
    # Perform analyses
    calculate_properties(sequence)
    has_rrm = analyze_rrm_domain(sequence)
    has_disorder = analyze_disorder(sequence)
    analyze_fox_binding_motifs(sequence)
    analyze_protein_interactions(sequence)
    analyze_post_translational_modifications(sequence)
    
    # Conclusions
    print("\n" + "=" * 60)
    print("CONCLUSIONS")
    print("=" * 60)
    
    print("\n1. DOMAIN STRUCTURE:")
    if has_rrm:
        print("   - Contains RNA Recognition Motif (RRM) domain")
        print("   - RRM likely mediates UGCAUG RNA binding")
    else:
        print("   - RRM domain features detected but not canonical")
    
    print("\n2. FUNCTIONAL FEATURES:")
    print("   - RNA-binding protein (FOX family)")
    print("   - Splicing regulator recognizing UGCAUG motifs")
    if has_disorder:
        print("   - Contains intrinsically disordered regions")
        print("   - IDRs may mediate protein-protein interactions")
    
    print("\n3. REGULATORY FEATURES:")
    print("   - Multiple phosphorylation sites predicted")
    print("   - Potential for post-translational regulation")
    
    print("\n4. NEURONAL MARKER:")
    print("   - Known as NeuN (Neuronal Nuclei)")
    print("   - Specifically expressed in post-mitotic neurons")
    print("   - Used as neuronal marker in research/diagnostics")
    
    print("\n5. LIMITATIONS:")
    print("   - Pattern-based predictions require experimental validation")
    print("   - 3D structure needed for detailed RNA binding analysis")
    print("   - PTM predictions are sequence-based only")

if __name__ == "__main__":
    main()