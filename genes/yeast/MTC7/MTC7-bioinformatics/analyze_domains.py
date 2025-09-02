#!/usr/bin/env python3
"""
Analyze functional domains and motifs in MTC7 protein.

This script searches for known domains, motifs, and sequence patterns
including the C-terminal polybasic tract.
"""

import json
import re
import click
from Bio import SeqIO
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def find_polybasic_regions(sequence, min_basic=4, window=10):
    """
    Find polybasic regions in the sequence.
    
    Args:
        sequence: Protein sequence
        min_basic: Minimum number of basic residues (K/R) in window
        window: Window size for searching
    """
    basic_aa = set('KR')
    polybasic_regions = []
    
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        basic_count = sum(1 for aa in window_seq if aa in basic_aa)
        
        if basic_count >= min_basic:
            # Extend to find the full polybasic region
            start = i
            end = i + window
            
            # Extend left
            while start > 0 and sequence[start - 1] in basic_aa:
                start -= 1
            
            # Extend right
            while end < len(sequence) and sequence[end] in basic_aa:
                end += 1
            
            # Check if this overlaps with previously found regions
            if not polybasic_regions or start > polybasic_regions[-1][1]:
                polybasic_regions.append((start + 1, end))  # 1-indexed
    
    return polybasic_regions


def analyze_sequence_patterns(sequence):
    """
    Search for common sequence patterns and motifs.
    """
    patterns = {
        'Nuclear Localization Signal (NLS)': [
            r'[KR]{2}[A-Z]{0,2}[KR]{2}',  # Bipartite NLS
            r'[KR]{4,}',  # Classical monopartite NLS
            r'P[KR]{3,5}',  # Pat7 NLS
        ],
        'ER retention signal': [
            r'[KR]DEL$',  # KDEL/RDEL at C-terminus
            r'KK[A-Z]{0,2}$',  # Di-lysine motif at C-terminus
        ],
        'Phosphorylation sites': [
            r'[ST][A-Z][KR]',  # PKA/PKC consensus
            r'[ST]P',  # Proline-directed kinases
        ],
        'Myristoylation': [
            r'^MG[A-Z]{3}[ST]',  # N-terminal myristoylation
        ],
        'Palmitoylation': [
            r'C[A-Z]{2,4}C',  # Cysteine-rich regions
        ]
    }
    
    found_patterns = {}
    
    for pattern_name, pattern_list in patterns.items():
        matches = []
        for pattern in pattern_list:
            for match in re.finditer(pattern, sequence):
                matches.append({
                    'position': match.start() + 1,  # 1-indexed
                    'sequence': match.group(),
                    'pattern': pattern
                })
        
        if matches:
            found_patterns[pattern_name] = matches
    
    return found_patterns


def analyze_compositional_bias(sequence):
    """
    Analyze compositional biases in the sequence.
    """
    # Count amino acids
    aa_counts = {}
    for aa in sequence:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1
    
    total = len(sequence)
    
    # Calculate percentages
    aa_percentages = {aa: (count / total) * 100 for aa, count in aa_counts.items()}
    
    # Group by properties
    groups = {
        'Hydrophobic': 'AILMFWV',
        'Aromatic': 'FWY',
        'Positive': 'KRH',
        'Negative': 'DE',
        'Polar': 'NQST',
        'Small': 'AGST',
        'Tiny': 'AGS'
    }
    
    group_percentages = {}
    for group_name, aas in groups.items():
        percentage = sum(aa_percentages.get(aa, 0) for aa in aas)
        group_percentages[group_name] = percentage
    
    return aa_percentages, group_percentages


def analyze_regions(sequence):
    """
    Analyze different regions of the protein.
    """
    # Define regions based on known TM domains
    regions = {
        'N-terminal (1-12)': sequence[0:12],
        'TM1 (13-33)': sequence[12:33],
        'Loop1 (34-41)': sequence[33:41],
        'TM2 (42-62)': sequence[41:62],
        'C-terminal (63-end)': sequence[62:]
    }
    
    region_analysis = {}
    
    for region_name, region_seq in regions.items():
        aa_comp, group_comp = analyze_compositional_bias(region_seq)
        
        # Find polybasic stretches in region
        polybasic = find_polybasic_regions(region_seq, min_basic=3, window=6)
        
        region_analysis[region_name] = {
            'length': len(region_seq),
            'sequence': region_seq,
            'group_composition': group_comp,
            'polybasic_regions': polybasic,
            'lysine_count': region_seq.count('K'),
            'arginine_count': region_seq.count('R')
        }
    
    return region_analysis


def plot_charge_distribution(sequence, output_prefix='domain_analysis', protein_name='Protein'):
    """
    Plot the charge distribution along the protein sequence.
    """
    charges = {'K': 1, 'R': 1, 'H': 0.5, 'D': -1, 'E': -1}
    
    charge_values = []
    for aa in sequence:
        charge_values.append(charges.get(aa, 0))
    
    # Calculate running average with window
    window = 5
    running_avg = []
    positions = []
    
    for i in range(len(sequence) - window + 1):
        avg_charge = sum(charge_values[i:i+window]) / window
        running_avg.append(avg_charge)
        positions.append(i + window // 2 + 1)
    
    plt.figure(figsize=(12, 4))
    plt.plot(positions, running_avg, 'b-', linewidth=1.5)
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Mark TM regions
    tm_regions = [(13, 33), (42, 62)]
    for i, (start, end) in enumerate(tm_regions, 1):
        plt.axvspan(start, end, alpha=0.2, color='gray', label=f'TM{i}' if i == 1 else '')
    
    # Mark C-terminal region
    plt.axvspan(63, len(sequence), alpha=0.2, color='yellow', label='C-terminal')
    
    plt.xlabel('Position')
    plt.ylabel('Average Charge (window=5)')
    plt.title(f'{protein_name} Charge Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plot_file = f'{output_prefix}_charge_distribution.png'
    plt.savefig(plot_file, dpi=150)
    plt.close()
    return plot_file


@click.command()
@click.argument('fasta_file', type=click.Path(exists=True))
@click.option('--output-prefix', default='domain_analysis', help='Prefix for output files')
def main(fasta_file, output_prefix):
    """
    Analyze domains and motifs in a protein sequence from FASTA file.
    
    FASTA_FILE: Path to the input FASTA file containing protein sequence
    """
    # Read the sequence
    with open(fasta_file, 'r') as f:
        record = next(SeqIO.parse(f, "fasta"))
    
    sequence = str(record.seq)
    protein_name = record.id.split('|')[-1] if '|' in record.id else record.id
    
    print(f"Analyzing {protein_name} domains and motifs")
    print("=" * 60)
    
    # Find polybasic regions
    polybasic_regions = find_polybasic_regions(sequence, min_basic=4, window=8)
    
    print("\nPolybasic Regions:")
    print("-" * 40)
    for start, end in polybasic_regions:
        region_seq = sequence[start-1:end]
        k_count = region_seq.count('K')
        r_count = region_seq.count('R')
        print(f"  Position {start}-{end}: {region_seq}")
        print(f"    Lysines: {k_count}, Arginines: {r_count}")
    
    # Analyze C-terminal specifically
    c_terminal = sequence[62:]  # After TM2
    print(f"\nC-terminal Analysis (positions 63-{len(sequence)}):")
    print("-" * 40)
    print(f"  Sequence: {c_terminal}")
    print(f"  Length: {len(c_terminal)} aa")
    print(f"  Lysine count: {c_terminal.count('K')}")
    print(f"  Arginine count: {c_terminal.count('R')}")
    print(f"  Basic residue percentage: {((c_terminal.count('K') + c_terminal.count('R')) / len(c_terminal)) * 100:.1f}%")
    
    # Search for sequence patterns
    patterns = analyze_sequence_patterns(sequence)
    
    print("\nSequence Patterns and Motifs:")
    print("-" * 40)
    for pattern_type, matches in patterns.items():
        print(f"\n{pattern_type}:")
        for match in matches:
            print(f"  Position {match['position']}: {match['sequence']}")
    
    # Analyze regional composition
    regions = analyze_regions(sequence)
    
    print("\nRegional Analysis:")
    print("-" * 40)
    for region_name, region_data in regions.items():
        print(f"\n{region_name}:")
        print(f"  Length: {region_data['length']} aa")
        if region_data['lysine_count'] > 0 or region_data['arginine_count'] > 0:
            print(f"  Basic residues: K={region_data['lysine_count']}, R={region_data['arginine_count']}")
        if region_data['polybasic_regions']:
            print(f"  Contains polybasic region(s)")
    
    # Analyze overall composition
    aa_comp, group_comp = analyze_compositional_bias(sequence)
    
    print("\nOverall Amino Acid Composition by Groups:")
    print("-" * 40)
    for group, percentage in sorted(group_comp.items(), key=lambda x: x[1], reverse=True):
        print(f"  {group}: {percentage:.1f}%")
    
    # Plot charge distribution
    plot_file = plot_charge_distribution(sequence, output_prefix, protein_name)
    print(f"\nCharge distribution plot saved as '{plot_file}'")
    
    # Save results
    results = {
        'protein': protein_name,
        'protein_length': len(sequence),
        'polybasic_regions': polybasic_regions,
        'c_terminal_analysis': {
            'sequence': c_terminal,
            'length': len(c_terminal),
            'lysine_count': c_terminal.count('K'),
            'arginine_count': c_terminal.count('R'),
            'basic_percentage': ((c_terminal.count('K') + c_terminal.count('R')) / len(c_terminal)) * 100
        },
        'sequence_patterns': {k: [m['sequence'] + f" (pos {m['position']})" for m in v] 
                              for k, v in patterns.items()},
        'regional_analysis': {k: {
            'length': v['length'],
            'lysine_count': v['lysine_count'],
            'arginine_count': v['arginine_count']
        } for k, v in regions.items()},
        'composition_groups': group_comp
    }
    
    json_file = f'{output_prefix}_results.json'
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to '{json_file}'")


if __name__ == "__main__":
    main()