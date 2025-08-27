#!/usr/bin/env python3
"""
RBFOX3 Protein Sequence Analysis

This script performs comprehensive sequence analysis of RBFOX3 including:
- Domain architecture analysis
- Conservation analysis across species
- RNA-binding motif prediction
- Secondary structure prediction
"""

from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from pathlib import Path
import json

# RBFOX3 human sequence from UniProt (A6NFN3)
RBFOX3_SEQUENCE = """MAQPYPPAQYPPPPQNGIPAYAPPPPHPTQDYSGQTPVPTEHGMTLYTP
AQTHPEQPGSEASTQPIAGTQTVPQTDEAAQTDSQPLHPSDPTEKQQPKRLHVSNIPFRFRDPDLRQMFG
QFGKILDVEIIFNERGSKGFGFVTFETSSDADRAREKLNGTIVEGRKIEVNNATARVMTNKKTGNPYTNG
WKLNPVVGAVYGPEFYAVTGFPYPTTGTAVAYRGAHLRGRGRRAVYNTFRAPPPPPIPTYGGAVVYQDGF
YGAEIYGGYAAYRYAQPAAAAAYSDSYGRVYAAADPYHHTTIGPAATYSIGTTM""".replace('\n', '')

def analyze_basic_properties():
    """Analyze basic biochemical properties of RBFOX3"""
    print("=== RBFOX3 Basic Sequence Properties ===")
    
    seq = Seq(RBFOX3_SEQUENCE)
    protein_analysis = ProteinAnalysis(str(seq))
    
    print(f"Sequence length: {len(seq)} amino acids")
    print(f"Molecular weight: {protein_analysis.molecular_weight():.1f} Da")
    print(f"Theoretical pI: {protein_analysis.isoelectric_point():.2f}")
    print(f"Instability index: {protein_analysis.instability_index():.2f}")
    print(f"Gravy score: {protein_analysis.gravy():.3f}")
    
    # Amino acid composition
    aa_percent = protein_analysis.get_amino_acids_percent()
    print("\nAmino acid composition (top 5):")
    sorted_aa = sorted(aa_percent.items(), key=lambda x: x[1], reverse=True)
    for aa, percent in sorted_aa[:5]:
        print(f"  {aa}: {percent:.1%}")
    
    return protein_analysis

def identify_domains():
    """Identify known domains in RBFOX3"""
    print("\n=== RBFOX3 Domain Architecture ===")
    
    # Based on UniProt annotation
    domains = {
        "N-terminal disordered region": (1, 104),
        "RNA Recognition Motif (RRM)": (100, 175), 
        "Fox-1 C-terminal domain": (176, 312),
        "Pro-rich region": (1, 29),
        "Polar residue region": (49, 87)
    }
    
    print("Known domains and regions:")
    for domain_name, (start, end) in domains.items():
        length = end - start + 1
        print(f"  {domain_name}: {start}-{end} ({length} aa)")
    
    # Analyze RRM domain sequence
    rrm_seq = RBFOX3_SEQUENCE[99:175]  # 0-indexed
    print(f"\nRRM domain sequence (aa 100-175):")
    print(f"  {rrm_seq}")
    
    # RNA interaction sites from UniProt
    rna_sites = [100, 108, 109, 133, 138, 142, 166, 176]
    print(f"\nRNA interaction sites: {rna_sites}")
    
    return domains

def analyze_sequence_features():
    """Analyze sequence features relevant to RNA binding"""
    print("\n=== RBFOX3 Sequence Features Analysis ===")
    
    seq = RBFOX3_SEQUENCE
    
    # Look for RNA-binding motifs
    rnp1_pattern = r'[RK].{2}[FY].{2}[FY]'  # RNP1 motif pattern
    rnp2_pattern = r'[LIV].{2}[FY].{2}[RKH]'  # RNP2 motif pattern
    
    rnp1_matches = [(m.start()+1, m.end(), m.group()) for m in re.finditer(rnp1_pattern, seq)]
    rnp2_matches = [(m.start()+1, m.end(), m.group()) for m in re.finditer(rnp2_pattern, seq)]
    
    print("RNP motif analysis:")
    print(f"  RNP1-like motifs: {rnp1_matches}")
    print(f"  RNP2-like motifs: {rnp2_matches}")
    
    # Analyze charge distribution
    positive_residues = seq.count('K') + seq.count('R') + seq.count('H')
    negative_residues = seq.count('D') + seq.count('E')
    net_charge = positive_residues - negative_residues
    
    print(f"\nCharge analysis:")
    print(f"  Positive residues (K,R,H): {positive_residues}")
    print(f"  Negative residues (D,E): {negative_residues}")
    print(f"  Net charge: {net_charge:+d}")
    
    # Proline content (relevant for disordered regions)
    proline_count = seq.count('P')
    proline_percent = (proline_count / len(seq)) * 100
    print(f"\nProline content: {proline_count} ({proline_percent:.1f}%)")
    
    return {
        'rnp1_motifs': rnp1_matches,
        'rnp2_motifs': rnp2_matches,
        'net_charge': net_charge,
        'proline_percent': proline_percent
    }

def fetch_rbfox_family_sequences():
    """Fetch sequences of RBFOX family members for comparative analysis"""
    print("\n=== Fetching RBFOX Family Sequences ===")
    
    # UniProt IDs for human RBFOX proteins
    rbfox_proteins = {
        'RBFOX1': 'Q9NWB1',
        'RBFOX2': 'O43251', 
        'RBFOX3': 'A6NFN3'
    }
    
    sequences = {}
    
    for protein, uniprot_id in rbfox_proteins.items():
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                fasta_content = response.text
                # Parse the FASTA content
                lines = fasta_content.strip().split('\n')
                if len(lines) >= 2:
                    sequence = ''.join(lines[1:])
                    sequences[protein] = sequence
                    print(f"  Retrieved {protein} ({uniprot_id}): {len(sequence)} aa")
                else:
                    print(f"  Failed to parse {protein} sequence")
            else:
                print(f"  Failed to fetch {protein}: HTTP {response.status_code}")
                
        except requests.RequestException as e:
            print(f"  Error fetching {protein}: {e}")
    
    return sequences

def analyze_family_conservation(sequences):
    """Analyze conservation across RBFOX family members"""
    print("\n=== RBFOX Family Conservation Analysis ===")
    
    if len(sequences) < 2:
        print("Not enough sequences for conservation analysis")
        return
    
    # Compare sequence lengths
    print("Sequence lengths:")
    for protein, seq in sequences.items():
        print(f"  {protein}: {len(seq)} aa")
    
    # Focus on RRM domain conservation (approximate positions)
    rrm_start, rrm_end = 80, 160  # Approximate RRM region
    
    print(f"\nRRM domain comparison (approximate positions {rrm_start}-{rrm_end}):")
    for protein, seq in sequences.items():
        if len(seq) > rrm_end:
            rrm_region = seq[rrm_start:rrm_end]
            print(f"  {protein}: {rrm_region}")
    
    # Calculate pairwise identity in RRM region
    if 'RBFOX3' in sequences and 'RBFOX1' in sequences:
        rbfox3_rrm = sequences['RBFOX3'][rrm_start:rrm_end] if len(sequences['RBFOX3']) > rrm_end else sequences['RBFOX3']
        rbfox1_rrm = sequences['RBFOX1'][rrm_start:rrm_end] if len(sequences['RBFOX1']) > rrm_end else sequences['RBFOX1']
        
        # Simple identity calculation
        min_len = min(len(rbfox3_rrm), len(rbfox1_rrm))
        identical = sum(1 for i in range(min_len) if rbfox3_rrm[i] == rbfox1_rrm[i])
        identity = (identical / min_len) * 100
        
        print(f"\nRBFOX3 vs RBFOX1 RRM identity: {identity:.1f}%")

def create_domain_visualization():
    """Create a visualization of RBFOX3 domain architecture"""
    print("\n=== Creating Domain Visualization ===")
    
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Domain definitions
    domains = [
        ("N-term disordered", 1, 104, "lightgray"),
        ("RRM domain", 100, 175, "red"),
        ("Fox-1 C-domain", 176, 312, "blue")
    ]
    
    # Draw protein length
    protein_length = 312
    ax.plot([1, protein_length], [1, 1], 'k-', linewidth=8, alpha=0.3)
    
    # Draw domains
    for name, start, end, color in domains:
        ax.plot([start, end], [1, 1], color=color, linewidth=12, alpha=0.8)
        # Add domain label
        mid_point = (start + end) / 2
        ax.text(mid_point, 1.2, name, ha='center', va='bottom', fontsize=10, rotation=0)
    
    # Mark RNA interaction sites
    rna_sites = [100, 108, 109, 133, 138, 142, 166, 176]
    for site in rna_sites:
        ax.plot(site, 1, 'ko', markersize=6)
    
    ax.set_xlim(0, 320)
    ax.set_ylim(0.5, 1.8)
    ax.set_xlabel('Amino Acid Position')
    ax.set_title('RBFOX3 Domain Architecture')
    ax.set_yticks([])
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], color='lightgray', lw=6, label='Disordered region'),
        plt.Line2D([0], [0], color='red', lw=6, label='RRM domain'),
        plt.Line2D([0], [0], color='blue', lw=6, label='Fox-1 C-domain'),
        plt.Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=6, label='RNA interaction sites')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig('rbfox3_domains.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Domain visualization saved as 'rbfox3_domains.png'")

def main():
    """Main analysis function"""
    print("RBFOX3 Comprehensive Sequence Analysis")
    print("=" * 50)
    
    # Basic properties
    protein_analysis = analyze_basic_properties()
    
    # Domain identification
    domains = identify_domains()
    
    # Sequence features
    features = analyze_sequence_features()
    
    # Family comparison
    family_sequences = fetch_rbfox_family_sequences()
    if family_sequences:
        analyze_family_conservation(family_sequences)
    
    # Create visualization
    create_domain_visualization()
    
    # Save results summary
    results = {
        'sequence_length': len(RBFOX3_SEQUENCE),
        'molecular_weight': protein_analysis.molecular_weight(),
        'isoelectric_point': protein_analysis.isoelectric_point(),
        'domains': domains,
        'sequence_features': features,
        'family_sequences_retrieved': list(family_sequences.keys()) if family_sequences else []
    }
    
    with open('rbfox3_analysis_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print("\n" + "=" * 50)
    print("Analysis complete. Results saved to 'rbfox3_analysis_results.json'")

if __name__ == "__main__":
    main()