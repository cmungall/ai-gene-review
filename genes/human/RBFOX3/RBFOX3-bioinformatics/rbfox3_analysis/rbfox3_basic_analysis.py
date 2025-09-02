#!/usr/bin/env python3
"""
RBFOX3 Basic Sequence Analysis

This script performs fundamental sequence analysis of RBFOX3
"""

from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import re
import json
from pathlib import Path

def load_sequence_from_fasta(fasta_path):
    """Load sequence from FASTA file"""
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    # Skip header, join sequence lines
    sequence = ''.join(line.strip() for line in lines[1:])
    return sequence

# Load RBFOX3 sequence from FASTA file
fasta_path = Path(__file__).parent.parent.parent / "RBFOX3.fasta"
if fasta_path.exists():
    RBFOX3_SEQUENCE = load_sequence_from_fasta(fasta_path)
else:
    raise FileNotFoundError(f"FASTA file not found at {fasta_path}. Please run: python scripts/uniprot2fasta.py genes/human/RBFOX3/RBFOX3-uniprot.txt")

def analyze_basic_properties():
    """Analyze basic biochemical properties of RBFOX3"""
    print("=== RBFOX3 Basic Sequence Properties ===")
    
    seq = Seq(RBFOX3_SEQUENCE)
    protein_analysis = ProteinAnalysis(str(seq))
    
    print(f"Sequence length: {len(seq)} amino acids")
    print(f"Molecular weight: {protein_analysis.molecular_weight():.1f} Da")
    print(f"Theoretical pI: {protein_analysis.isoelectric_point():.2f}")
    
    # Amino acid composition
    aa_percent = protein_analysis.amino_acids_percent
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
    
    return domains

def analyze_sequence_features():
    """Analyze sequence features relevant to RNA binding"""
    print("\n=== RBFOX3 Sequence Features Analysis ===")
    
    seq = RBFOX3_SEQUENCE
    
    # Look for RNA-binding motifs
    rnp1_pattern = r'[RK].{2}[FY].{2}[FY]'
    rnp2_pattern = r'[LIV].{2}[FY].{2}[RKH]'
    
    rnp1_matches = [(m.start()+1, m.end(), m.group()) for m in re.finditer(rnp1_pattern, seq)]
    rnp2_matches = [(m.start()+1, m.end(), m.group()) for m in re.finditer(rnp2_pattern, seq)]
    
    print("RNP motif analysis:")
    print(f"  RNP1-like motifs: {rnp1_matches}")
    print(f"  RNP2-like motifs: {rnp2_matches}")
    
    # Charge analysis
    positive_residues = seq.count('K') + seq.count('R') + seq.count('H')
    negative_residues = seq.count('D') + seq.count('E')
    net_charge = positive_residues - negative_residues
    
    print(f"\nCharge analysis:")
    print(f"  Positive residues (K,R,H): {positive_residues}")
    print(f"  Negative residues (D,E): {negative_residues}")
    print(f"  Net charge: {net_charge:+d}")
    
    # Proline content
    proline_count = seq.count('P')
    proline_percent = (proline_count / len(seq)) * 100
    print(f"\nProline content: {proline_count} ({proline_percent:.1f}%)")
    
    return {
        'rnp1_motifs': rnp1_matches,
        'rnp2_motifs': rnp2_matches,
        'net_charge': net_charge,
        'proline_percent': proline_percent
    }

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
    
    # Draw protein backbone
    protein_length = 312
    ax.plot([1, protein_length], [1, 1], 'k-', linewidth=8, alpha=0.3)
    
    # Draw domains
    for name, start, end, color in domains:
        ax.plot([start, end], [1, 1], color=color, linewidth=12, alpha=0.8)
        mid_point = (start + end) / 2
        ax.text(mid_point, 1.2, name, ha='center', va='bottom', fontsize=10)
    
    # Mark RNA interaction sites from UniProt
    rna_sites = [100, 108, 109, 133, 138, 142, 166, 176]
    for site in rna_sites:
        ax.plot(site, 1, 'ko', markersize=6)
    
    ax.set_xlim(0, 320)
    ax.set_ylim(0.5, 1.8)
    ax.set_xlabel('Amino Acid Position')
    ax.set_title('RBFOX3 Domain Architecture')
    ax.set_yticks([])
    
    # Legend
    legend_elements = [
        plt.Line2D([0], [0], color='lightgray', lw=6, label='Disordered region'),
        plt.Line2D([0], [0], color='red', lw=6, label='RRM domain'),
        plt.Line2D([0], [0], color='blue', lw=6, label='Fox-1 C-domain'),
        plt.Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=6, label='RNA interaction sites')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig('rbfox3_domains.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Domain visualization saved as 'rbfox3_domains.png'")

def main():
    """Main analysis function"""
    print("RBFOX3 Basic Sequence Analysis")
    print("=" * 40)
    
    # Basic properties
    protein_analysis = analyze_basic_properties()
    
    # Domain identification
    domains = identify_domains()
    
    # Sequence features
    features = analyze_sequence_features()
    
    # Create visualization
    create_domain_visualization()
    
    # Save results
    results = {
        'sequence_length': len(RBFOX3_SEQUENCE),
        'molecular_weight': float(protein_analysis.molecular_weight()),
        'isoelectric_point': float(protein_analysis.isoelectric_point()),
        'domains': {k: list(v) for k, v in domains.items()},
        'sequence_features': features
    }
    
    with open('rbfox3_basic_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 40)
    print("Basic analysis complete. Results saved to 'rbfox3_basic_results.json'")

if __name__ == "__main__":
    main()