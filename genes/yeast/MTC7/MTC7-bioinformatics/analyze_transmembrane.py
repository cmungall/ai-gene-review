#!/usr/bin/env python3
"""
Analyze transmembrane regions of a protein.

This script analyzes the hydrophobicity and transmembrane regions of a protein
using various methods including Kyte-Doolittle hydropathy analysis.
"""

import json
import click
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from pathlib import Path


def kyte_doolittle_hydropathy(sequence, window_size=19):
    """
    Calculate Kyte-Doolittle hydropathy values.
    
    Uses the Kyte-Doolittle scale for amino acid hydrophobicity.
    Window size of 19 is standard for transmembrane helix prediction.
    """
    # Kyte-Doolittle hydrophobicity scale
    kd_scale = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
        'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
        'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
    
    hydropathy_values = []
    positions = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        hydropathy = sum(kd_scale.get(aa, 0) for aa in window) / window_size
        hydropathy_values.append(hydropathy)
        positions.append(i + window_size // 2 + 1)  # Center of window, 1-indexed
    
    return positions, hydropathy_values


def predict_tm_regions(positions, hydropathy_values, threshold=1.6, min_length=17):
    """
    Predict transmembrane regions based on hydropathy values.
    
    Args:
        positions: Position indices
        hydropathy_values: Hydropathy scores
        threshold: Minimum hydropathy score for TM region (default 1.6)
        min_length: Minimum length for a TM helix (default 17 amino acids)
    """
    tm_regions = []
    in_tm = False
    start = None
    
    for i, (pos, hydro) in enumerate(zip(positions, hydropathy_values)):
        if hydro >= threshold and not in_tm:
            in_tm = True
            start = pos
        elif hydro < threshold and in_tm:
            in_tm = False
            if pos - start >= min_length:
                tm_regions.append((start, pos - 1))
    
    # Check if we're still in a TM region at the end
    if in_tm and positions[-1] - start >= min_length:
        tm_regions.append((start, positions[-1]))
    
    return tm_regions


def analyze_tm_composition(sequence, tm_regions):
    """
    Analyze amino acid composition of transmembrane regions.
    """
    hydrophobic_aa = set('AILMFWV')
    aromatic_aa = set('FWY')
    charged_aa = set('DEKR')
    
    results = []
    
    for start, end in tm_regions:
        tm_seq = sequence[start-1:end]  # Convert to 0-indexed
        
        hydrophobic_count = sum(1 for aa in tm_seq if aa in hydrophobic_aa)
        aromatic_count = sum(1 for aa in tm_seq if aa in aromatic_aa)
        charged_count = sum(1 for aa in tm_seq if aa in charged_aa)
        
        results.append({
            'region': f'{start}-{end}',
            'sequence': tm_seq,
            'length': len(tm_seq),
            'hydrophobic_percent': (hydrophobic_count / len(tm_seq)) * 100,
            'aromatic_percent': (aromatic_count / len(tm_seq)) * 100,
            'charged_percent': (charged_count / len(tm_seq)) * 100
        })
    
    return results


@click.command()
@click.argument('fasta_file', type=click.Path(exists=True))
@click.option('--output-prefix', default='tm_analysis', help='Prefix for output files')
@click.option('--known-tm', multiple=True, help='Known TM regions in format start-end (e.g., 13-33)')
def main(fasta_file, output_prefix, known_tm):
    """
    Analyze transmembrane regions in a protein sequence from FASTA file.
    
    FASTA_FILE: Path to the input FASTA file containing protein sequence
    """
    # Read the sequence
    with open(fasta_file, 'r') as f:
        record = next(SeqIO.parse(f, "fasta"))
    
    sequence = str(record.seq)
    protein_name = record.id.split('|')[-1] if '|' in record.id else record.id
    
    print(f"Analyzing {protein_name} ({len(sequence)} amino acids)")
    print("=" * 60)
    
    # Calculate hydropathy
    positions, hydropathy_values = kyte_doolittle_hydropathy(sequence)
    
    # Predict TM regions
    predicted_tm = predict_tm_regions(positions, hydropathy_values)
    
    print("\nPredicted Transmembrane Regions (Kyte-Doolittle):")
    print("-" * 40)
    for start, end in predicted_tm:
        print(f"  TM{predicted_tm.index((start, end)) + 1}: {start}-{end} ({end - start + 1} aa)")
        print(f"       Sequence: {sequence[start-1:end]}")
    
    # Parse known TM regions if provided
    known_tm_list = []
    if known_tm:
        for tm in known_tm:
            start, end = map(int, tm.split('-'))
            known_tm_list.append((start, end))
        
        print("\nKnown Transmembrane Regions:")
        print("-" * 40)
        for i, (start, end) in enumerate(known_tm_list, 1):
            print(f"  TM{i}: {start}-{end} ({end - start + 1} aa)")
            print(f"       Sequence: {sequence[start-1:end]}")
    
    # Analyze composition of known TM regions if provided
    tm_composition = analyze_tm_composition(sequence, known_tm_list) if known_tm_list else []
    
    print("\nTransmembrane Region Composition Analysis:")
    print("-" * 40)
    for tm_info in tm_composition:
        print(f"\nTM Region {tm_info['region']}:")
        print(f"  Length: {tm_info['length']} aa")
        print(f"  Hydrophobic: {tm_info['hydrophobic_percent']:.1f}%")
        print(f"  Aromatic: {tm_info['aromatic_percent']:.1f}%")
        print(f"  Charged: {tm_info['charged_percent']:.1f}%")
    
    # Create hydropathy plot
    plt.figure(figsize=(12, 6))
    plt.plot(positions, hydropathy_values, 'b-', linewidth=1.5, label='Hydropathy')
    plt.axhline(y=1.6, color='r', linestyle='--', alpha=0.5, label='TM threshold')
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Mark known TM regions if provided
    if known_tm_list:
        for i, (start, end) in enumerate(known_tm_list, 1):
            plt.axvspan(start, end, alpha=0.3, color='yellow', label=f'Known TM{i}' if i == 1 else '')
    
    # Mark predicted TM regions
    for i, (start, end) in enumerate(predicted_tm, 1):
        plt.axvspan(start, end, alpha=0.2, color='green', label='Predicted TM' if i == 1 else '')
    
    plt.xlabel('Position')
    plt.ylabel('Hydropathy Score')
    plt.title(f'{protein_name} Hydropathy Plot (Kyte-Doolittle)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot with output prefix
    plot_file = f'{output_prefix}_hydropathy.png'
    plt.savefig(plot_file, dpi=150)
    plt.close()
    
    print(f"\nHydropathy plot saved as '{plot_file}'")
    
    # Save results to JSON
    results = {
        'protein': protein_name,
        'protein_length': len(sequence),
        'known_tm_regions': known_tm_list,
        'predicted_tm_regions': predicted_tm,
        'tm_composition': tm_composition,
        'analysis_method': 'Kyte-Doolittle hydropathy (window=19)'
    }
    
    json_file = f'{output_prefix}_results.json'
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to '{json_file}'")


if __name__ == "__main__":
    main()