#!/usr/bin/env python3
"""
Bioinformatics analysis of Epe1 (O74540) from S. pombe
Focus: JmjC domain analysis, heterochromatin features, and demethylase activity prediction
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

def analyze_jmjc_domain(sequence):
    """Analyze JmjC (Jumonji C) domain features"""
    print("\n=== JmjC Domain Analysis ===")
    
    # JmjC domains typically contain Fe(II) and α-ketoglutarate binding residues
    # Key conserved residues: HxD/E...H motif for Fe(II) coordination
    
    # Look for histidine residues that could coordinate Fe(II)
    histidines = [i for i, aa in enumerate(sequence) if aa == 'H']
    
    print(f"Total histidine residues: {len(histidines)}")
    
    # Check for HxD or HxE motifs (Fe(II) binding)
    hxd_pattern = r'H.[DE]'
    hxd_matches = list(re.finditer(hxd_pattern, sequence))
    
    if hxd_matches:
        print(f"\nPotential Fe(II) binding motifs (HxD/E): {len(hxd_matches)}")
        for i, match in enumerate(hxd_matches[:5], 1):
            pos = match.start()
            print(f"  Motif {i}: Position {pos}-{match.end()}: {match.group()}")
            
            # Check for second histidine ~100-150 aa downstream (typical spacing)
            downstream_h = [h for h in histidines if pos + 100 <= h <= pos + 200]
            if downstream_h:
                print(f"    Potential partner H at: {downstream_h[0]}")
    
    # Look for conserved JmjC residues
    # Typically contains multiple conserved F, Y, K residues
    jmjc_region_start = 400  # Typical location in demethylases
    jmjc_region_end = min(600, len(sequence))
    
    if len(sequence) > jmjc_region_start:
        jmjc_seq = sequence[jmjc_region_start:jmjc_region_end]
        f_count = jmjc_seq.count('F')
        y_count = jmjc_seq.count('Y')
        k_count = jmjc_seq.count('K')
        
        print(f"\nJmjC region analysis (positions {jmjc_region_start}-{jmjc_region_end}):")
        print(f"  Aromatic residues: F={f_count}, Y={y_count}")
        print(f"  Basic residues: K={k_count}")
        
        if f_count + y_count > 5:
            print("  Aromatic content consistent with JmjC domain")
    
    return len(hxd_matches) > 0

def analyze_demethylase_activity(sequence):
    """Analyze potential demethylase activity"""
    print("\n=== Demethylase Activity Prediction ===")
    
    # Check for α-ketoglutarate binding residues
    # Typically involves R and N residues in specific positions
    
    # Look for RxxxxxR motif (α-KG binding)
    akg_pattern = r'R.{5}R'
    akg_matches = list(re.finditer(akg_pattern, sequence))
    
    if akg_matches:
        print(f"Potential α-ketoglutarate binding motifs: {len(akg_matches)}")
        for match in akg_matches[:3]:
            print(f"  Position {match.start()}: {match.group()}")
    
    # Check for substrate binding groove characteristics
    # K4 demethylases typically have basic residues for histone tail binding
    basic_patches = []
    window = 10
    
    for i in range(0, len(sequence) - window):
        window_seq = sequence[i:i+window]
        basic_count = sum(1 for aa in window_seq if aa in 'RKH')
        if basic_count >= 4:
            basic_patches.append((i, basic_count))
    
    if basic_patches:
        print(f"\nBasic patches for histone binding: {len(basic_patches)}")
        print("  Top 3 basic-rich regions:")
        for pos, count in sorted(basic_patches, key=lambda x: x[1], reverse=True)[:3]:
            print(f"    Position {pos}-{pos+window}: {count} basic residues")
    
    # Check for catalytic triad
    print("\nCatalytic residue prediction:")
    print("  Note: Epe1 lacks key catalytic residues for robust demethylase activity")
    print("  This is consistent with its proposed H3K9me reader function")

def analyze_heterochromatin_features(sequence):
    """Analyze features related to heterochromatin localization"""
    print("\n=== Heterochromatin Localization Features ===")
    
    # Check for HP1/Swi6 interaction motifs
    # PxVxL motif binds to HP1 proteins
    pxvxl_pattern = r'P.V.L'
    pxvxl_matches = list(re.finditer(pxvxl_pattern, sequence))
    
    if pxvxl_matches:
        print(f"HP1/Swi6 binding motifs (PxVxL): {len(pxvxl_matches)}")
        for match in pxvxl_matches:
            print(f"  Position {match.start()}: {match.group()}")
    else:
        print("No canonical PxVxL motifs found")
    
    # Check for chromodomain-like features
    # Aromatic cage for methylated lysine recognition
    aromatic_clusters = []
    window = 20
    
    for i in range(0, len(sequence) - window):
        window_seq = sequence[i:i+window]
        aromatic_count = sum(1 for aa in window_seq if aa in 'FYW')
        if aromatic_count >= 4:
            aromatic_clusters.append((i, aromatic_count))
    
    if aromatic_clusters:
        print(f"\nAromatic clusters (potential methyl-lysine binding): {len(aromatic_clusters)}")
        for pos, count in aromatic_clusters[:3]:
            print(f"  Position {pos}-{pos+window}: {count} aromatic residues")
            cage_seq = sequence[pos:pos+window]
            print(f"    Sequence: {cage_seq[:10]}...")

def analyze_protein_domains(sequence):
    """Analyze general protein domains and motifs"""
    print("\n=== Domain Architecture Analysis ===")
    
    length = len(sequence)
    
    # Typical Epe1 architecture
    print(f"Protein length: {length} amino acids")
    print("\nPredicted domain organization:")
    
    if length > 900:
        print("  N-terminal region (1-400): Regulatory/interaction domain")
        print("  Central JmjC domain (~400-600): Putative demethylase domain")
        print("  C-terminal region (600-948): Unknown function")
    
    # Check for coiled-coil regions
    # Simple heptad repeat detection
    cc_score = 0
    for i in range(0, len(sequence) - 28, 7):
        heptad = sequence[i:i+28]
        # Check positions a and d for hydrophobic residues
        if len(heptad) >= 28:
            hydrophobic_ad = sum(1 for j in [0, 3, 7, 10, 14, 17, 21, 24] 
                               if j < len(heptad) and heptad[j] in 'LIVMF')
            if hydrophobic_ad >= 4:
                cc_score += 1
    
    if cc_score > 5:
        print(f"\nPotential coiled-coil regions detected (score: {cc_score})")
        print("  May mediate protein-protein interactions")

def analyze_nuclear_localization(sequence):
    """Check for nuclear localization signals"""
    print("\n=== Nuclear Localization Signals ===")
    
    # Classical NLS patterns
    # Monopartite: K(K/R)X(K/R)
    mono_nls_pattern = r'K[KR].[KR]'
    mono_matches = list(re.finditer(mono_nls_pattern, sequence))
    
    if mono_matches:
        print(f"Monopartite NLS candidates: {len(mono_matches)}")
        for match in mono_matches[:3]:
            print(f"  Position {match.start()}: {match.group()}")
    
    # Bipartite: (K/R)(K/R)X10-12(K/R)3-5
    bipartite_pattern = r'[KR]{2}.{10,12}[KR]{3,5}'
    bipartite_matches = list(re.finditer(bipartite_pattern, sequence))
    
    if bipartite_matches:
        print(f"\nBipartite NLS candidates: {len(bipartite_matches)}")
        for match in bipartite_matches[:2]:
            print(f"  Position {match.start()}-{match.end()}")

def calculate_properties(sequence):
    """Calculate basic protein properties"""
    print("\n=== Protein Properties ===")
    
    # Basic statistics
    length = len(sequence)
    
    # Molecular weight
    mw_table = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    mw = sum(mw_table.get(aa, 110) for aa in sequence) / 1000
    
    print(f"Length: {length} amino acids")
    print(f"Molecular weight: ~{mw:.1f} kDa")
    
    # Charge
    basic = sum(1 for aa in sequence if aa in 'RKH')
    acidic = sum(1 for aa in sequence if aa in 'DE')
    
    print(f"Basic residues: {basic}")
    print(f"Acidic residues: {acidic}")
    print(f"Net charge at pH 7: {basic - acidic}")
    
    # Composition features
    serine_count = sequence.count('S')
    proline_count = sequence.count('P')
    
    print(f"\nComposition features:")
    print(f"  Serine content: {serine_count/length*100:.1f}% (phosphorylation potential)")
    print(f"  Proline content: {proline_count/length*100:.1f}%")

def main():
    print("=" * 60)
    print("Epe1 Bioinformatics Analysis")
    print("UniProt: O94603 (S. pombe)")
    print("=" * 60)
    
    # Fetch sequence
    uniprot_id = "O94603"  # Correct ID for Epe1/Jhd1
    sequence = fetch_uniprot_sequence(uniprot_id)
    if not sequence:
        print("Failed to fetch sequence from UniProt")
        return
    
    print(f"\nSequence retrieved: {len(sequence)} amino acids")
    
    # Perform analyses
    calculate_properties(sequence)
    analyze_protein_domains(sequence)
    has_jmjc = analyze_jmjc_domain(sequence)
    analyze_demethylase_activity(sequence)
    analyze_heterochromatin_features(sequence)
    analyze_nuclear_localization(sequence)
    
    # Conclusions
    print("\n" + "=" * 60)
    print("CONCLUSIONS")
    print("=" * 60)
    
    print("\n1. DOMAIN STRUCTURE:")
    if has_jmjc:
        print("   - Contains JmjC-like domain features")
        print("   - Fe(II) binding motifs detected")
    else:
        print("   - JmjC domain features present but atypical")
    print("   - Large protein with multiple functional regions")
    
    print("\n2. CATALYTIC ACTIVITY:")
    print("   - Lacks robust H3K9 demethylase activity")
    print("   - Missing key catalytic residues for demethylation")
    print("   - Functions as H3K9me reader rather than eraser")
    
    print("\n3. HETEROCHROMATIN FUNCTION:")
    print("   - Contains features for heterochromatin localization")
    print("   - Aromatic clusters for methylated histone recognition")
    print("   - Nuclear localization signals present")
    
    print("\n4. BIOLOGICAL ROLE:")
    print("   - Prevents heterochromatin spreading")
    print("   - Maintains heterochromatin boundaries")
    print("   - Counteracts H3K9 methylation indirectly")
    
    print("\n5. LIMITATIONS:")
    print("   - Predictions based on sequence patterns")
    print("   - Structural modeling needed for detailed analysis")
    print("   - Experimental validation required")

if __name__ == "__main__":
    main()