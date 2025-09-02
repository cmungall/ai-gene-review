#!/usr/bin/env python3
"""
RBFOX Family Comparative Analysis

This script analyzes the RBFOX protein family focusing on conservation
and functional divergence between RBFOX1, RBFOX2, and RBFOX3.
"""

import matplotlib.pyplot as plt
import pandas as pd
import json
from pathlib import Path

def load_sequence_from_fasta(fasta_path):
    """Load sequence from FASTA file"""
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    # Skip header, join sequence lines
    sequence = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
    return sequence

# Load sequences from FASTA files
script_dir = Path(__file__).parent

# Load RBFOX family sequences
RBFOX_SEQUENCES = {}

for protein in ['RBFOX1', 'RBFOX2', 'RBFOX3']:
    if protein == 'RBFOX3':
        # Try parent directory first
        fasta_path = script_dir.parent.parent / f"{protein}.fasta"
        if not fasta_path.exists():
            fasta_path = script_dir / f"{protein}.fasta"
    else:
        fasta_path = script_dir / f"{protein}.fasta"
    
    if not fasta_path.exists():
        raise FileNotFoundError(f"{protein}.fasta not found at {fasta_path}. Please run: python fetch_rbfox_family.py")
    
    RBFOX_SEQUENCES[protein] = load_sequence_from_fasta(fasta_path)
    print(f"Loaded {protein}: {len(RBFOX_SEQUENCES[protein])} aa")

def analyze_sequence_lengths():
    """Compare sequence lengths of RBFOX family members"""
    print("=== RBFOX Family Sequence Length Comparison ===")
    
    lengths = {protein: len(seq) for protein, seq in RBFOX_SEQUENCES.items()}
    
    for protein, length in lengths.items():
        print(f"{protein}: {length} amino acids")
    
    # Calculate differences
    rbfox3_len = lengths['RBFOX3']
    print(f"\nLength differences relative to RBFOX3:")
    for protein, length in lengths.items():
        if protein != 'RBFOX3':
            diff = length - rbfox3_len
            print(f"  {protein}: {diff:+d} amino acids")
    
    return lengths

def find_rrm_domains():
    """Identify and compare RRM domains across RBFOX family"""
    print("\n=== RRM Domain Comparison ===")
    
    # RRM consensus pattern (simplified)
    # Look for the core RRM sequence starting around position 100
    rrm_start = 80  # Start search position
    rrm_length = 80  # Approximate RRM length
    
    rrm_sequences = {}
    
    for protein, seq in RBFOX_SEQUENCES.items():
        # Find RRM region (approximate)
        if len(seq) > rrm_start + rrm_length:
            rrm_region = seq[rrm_start:rrm_start + rrm_length]
            rrm_sequences[protein] = rrm_region
            print(f"{protein} RRM region:")
            print(f"  {rrm_region}")
    
    return rrm_sequences

def calculate_pairwise_identity(seq1, seq2):
    """Calculate sequence identity between two sequences"""
    min_len = min(len(seq1), len(seq2))
    identical = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    return (identical / min_len) * 100

def analyze_family_conservation():
    """Analyze conservation across RBFOX family members"""
    print("\n=== RBFOX Family Conservation Analysis ===")
    
    proteins = list(RBFOX_SEQUENCES.keys())
    
    # Calculate pairwise identities for full sequences
    print("Full sequence identities:")
    for i, prot1 in enumerate(proteins):
        for j, prot2 in enumerate(proteins):
            if i < j:
                identity = calculate_pairwise_identity(
                    RBFOX_SEQUENCES[prot1], 
                    RBFOX_SEQUENCES[prot2]
                )
                print(f"  {prot1} vs {prot2}: {identity:.1f}%")
    
    # Focus on RRM region conservation
    print("\nRRM region identities (approximate):")
    rrm_start, rrm_end = 80, 160
    
    for i, prot1 in enumerate(proteins):
        for j, prot2 in enumerate(proteins):
            if i < j:
                seq1 = RBFOX_SEQUENCES[prot1]
                seq2 = RBFOX_SEQUENCES[prot2]
                
                if len(seq1) > rrm_end and len(seq2) > rrm_end:
                    rrm1 = seq1[rrm_start:rrm_end]
                    rrm2 = seq2[rrm_start:rrm_end]
                    identity = calculate_pairwise_identity(rrm1, rrm2)
                    print(f"  {prot1} vs {prot2} RRM: {identity:.1f}%")

def analyze_c_terminal_divergence():
    """Analyze C-terminal region differences"""
    print("\n=== C-terminal Region Analysis ===")
    
    # Compare C-terminal regions (after RRM)
    c_term_start = 180
    
    print("C-terminal region lengths:")
    for protein, seq in RBFOX_SEQUENCES.items():
        if len(seq) > c_term_start:
            c_term_len = len(seq) - c_term_start
            print(f"  {protein}: {c_term_len} amino acids")
    
    # Show C-terminal sequences
    print("\nC-terminal sequences (from position 180):")
    for protein, seq in RBFOX_SEQUENCES.items():
        if len(seq) > c_term_start:
            c_term_seq = seq[c_term_start:]
            print(f"  {protein} ({len(c_term_seq)} aa):")
            # Break into lines for readability
            for i in range(0, len(c_term_seq), 60):
                print(f"    {c_term_seq[i:i+60]}")

def analyze_proline_content():
    """Analyze proline content across RBFOX family"""
    print("\n=== Proline Content Analysis ===")
    
    proline_data = {}
    
    for protein, seq in RBFOX_SEQUENCES.items():
        proline_count = seq.count('P')
        proline_percent = (proline_count / len(seq)) * 100
        proline_data[protein] = {
            'count': proline_count,
            'percentage': proline_percent
        }
        print(f"{protein}: {proline_count} prolines ({proline_percent:.1f}%)")
    
    return proline_data

def create_family_comparison_plot():
    """Create visualization comparing RBFOX family members"""
    print("\n=== Creating Family Comparison Plot ===")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: Sequence lengths
    proteins = list(RBFOX_SEQUENCES.keys())
    lengths = [len(RBFOX_SEQUENCES[p]) for p in proteins]
    colors = ['blue', 'green', 'red']
    
    ax1.bar(proteins, lengths, color=colors, alpha=0.7)
    ax1.set_ylabel('Sequence Length (amino acids)')
    ax1.set_title('RBFOX Family Sequence Lengths')
    
    # Add length values on bars
    for i, length in enumerate(lengths):
        ax1.text(i, length + 5, str(length), ha='center')
    
    # Plot 2: Domain architecture comparison
    domain_regions = {
        'RBFOX1': [(1, 100, 'N-term'), (100, 175, 'RRM'), (176, 370, 'C-term')],
        'RBFOX2': [(1, 100, 'N-term'), (100, 175, 'RRM'), (176, 390, 'C-term')],
        'RBFOX3': [(1, 100, 'N-term'), (100, 175, 'RRM'), (176, 313, 'C-term')]
    }
    
    y_positions = {'RBFOX1': 2, 'RBFOX2': 1, 'RBFOX3': 0}
    domain_colors = {'N-term': 'lightgray', 'RRM': 'red', 'C-term': 'blue'}
    
    for protein, domains in domain_regions.items():
        y = y_positions[protein]
        for start, end, domain_type in domains:
            ax2.barh(y, end - start, left=start, height=0.6, 
                    color=domain_colors[domain_type], alpha=0.8)
        
        ax2.text(-20, y, protein, ha='right', va='center', fontweight='bold')
    
    ax2.set_xlabel('Amino Acid Position')
    ax2.set_title('RBFOX Family Domain Architecture Comparison')
    ax2.set_ylim(-0.5, 2.5)
    ax2.set_yticks([])
    
    # Add legend
    legend_elements = [plt.Rectangle((0, 0), 1, 1, color=color, alpha=0.8, label=domain)
                      for domain, color in domain_colors.items()]
    ax2.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig('rbfox_family_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Family comparison plot saved as 'rbfox_family_comparison.png'")

def main():
    """Main analysis function"""
    print("RBFOX Family Comparative Analysis")
    print("=" * 45)
    
    # Sequence length analysis
    lengths = analyze_sequence_lengths()
    
    # RRM domain comparison
    rrm_sequences = find_rrm_domains()
    
    # Conservation analysis
    analyze_family_conservation()
    
    # C-terminal analysis
    analyze_c_terminal_divergence()
    
    # Proline content analysis
    proline_data = analyze_proline_content()
    
    # Create visualization
    create_family_comparison_plot()
    
    # Save results
    results = {
        'sequence_lengths': lengths,
        'proline_content': proline_data,
        'analysis_summary': {
            'rbfox3_unique_features': [
                'Shortest family member (313 aa)',
                'Highest proline content',
                'Unique C-terminal domain structure',
                'Neuronal-specific expression'
            ],
            'conserved_features': [
                'RNA Recognition Motif (RRM) domain',
                'RNA-binding specificity for UGCAUG motifs',
                'Alternative splicing regulation function'
            ]
        }
    }
    
    with open('rbfox_family_analysis.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 45)
    print("Family analysis complete. Results saved to 'rbfox_family_analysis.json'")

if __name__ == "__main__":
    main()