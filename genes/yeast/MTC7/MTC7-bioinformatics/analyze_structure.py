#!/usr/bin/env python3
"""
Analyze structural features of MTC7 protein.

This script analyzes secondary structure predictions, disorder regions,
and checks for AlphaFold structures.
"""

import json
import click
import requests
from Bio import SeqIO
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def predict_secondary_structure(sequence):
    """
    Predict secondary structure using simple Chou-Fasman parameters.
    
    This is a simplified implementation for demonstration.
    For production use, consider tools like PSIPRED or JPred.
    """
    # Chou-Fasman propensities (simplified)
    helix_propensity = {
        'E': 1.51, 'A': 1.42, 'L': 1.21, 'M': 1.45, 'Q': 1.11,
        'K': 1.16, 'F': 1.13, 'H': 1.00, 'V': 1.06, 'I': 1.08,
        'Y': 0.69, 'C': 0.70, 'W': 1.08, 'R': 0.98, 'T': 0.83,
        'G': 0.57, 'N': 0.67, 'D': 1.01, 'S': 0.77, 'P': 0.57
    }
    
    sheet_propensity = {
        'M': 1.05, 'V': 1.70, 'I': 1.60, 'C': 1.19, 'Y': 1.47,
        'F': 1.38, 'Q': 1.10, 'W': 1.37, 'L': 1.30, 'T': 1.19,
        'A': 0.83, 'R': 0.93, 'G': 0.75, 'D': 0.54, 'K': 0.74,
        'S': 0.75, 'H': 0.87, 'N': 0.89, 'P': 0.55, 'E': 0.37
    }
    
    window = 6
    helix_scores = []
    sheet_scores = []
    
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        
        helix_score = sum(helix_propensity.get(aa, 1.0) for aa in window_seq) / window
        sheet_score = sum(sheet_propensity.get(aa, 1.0) for aa in window_seq) / window
        
        helix_scores.append(helix_score)
        sheet_scores.append(sheet_score)
    
    # Predict structure based on scores
    predictions = []
    for h, s in zip(helix_scores, sheet_scores):
        if h > 1.03 and h > s:
            predictions.append('H')  # Helix
        elif s > 1.05 and s > h:
            predictions.append('E')  # Beta sheet
        else:
            predictions.append('C')  # Coil/loop
    
    return predictions, helix_scores, sheet_scores


def predict_disorder(sequence):
    """
    Predict disordered regions using a simple composition-based method.
    
    Based on the observation that disordered regions are enriched in
    certain amino acids (P, E, S, K, R, Q, D) and depleted in others (W, Y, F, I, L, V).
    """
    disorder_promoting = set('PESKRQD')
    order_promoting = set('WYFILVMCN')
    
    window = 15
    disorder_scores = []
    
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        
        disorder_count = sum(1 for aa in window_seq if aa in disorder_promoting)
        order_count = sum(1 for aa in window_seq if aa in order_promoting)
        
        # Simple disorder score
        disorder_score = (disorder_count - order_count) / window
        disorder_scores.append(disorder_score)
    
    return disorder_scores


def check_alphafold_structure(uniprot_id):
    """
    Check if AlphaFold structure is available for this protein.
    """
    # Extract just the accession from the full ID
    accession = uniprot_id.split('|')[1] if '|' in uniprot_id else uniprot_id
    
    base_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    
    try:
        response = requests.get(base_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data:
                entry = data[0] if isinstance(data, list) else data
                return {
                    'available': True,
                    'accession': entry.get('entryId', accession),
                    'url': f"https://alphafold.ebi.ac.uk/entry/{accession}",
                    'confidence_version': entry.get('confidenceVersion', 'Unknown'),
                    'model_url': entry.get('pdbUrl', '')
                }
        return {'available': False, 'accession': accession}
    except Exception as e:
        print(f"Error checking AlphaFold: {e}")
        return {'available': False, 'accession': accession, 'error': str(e)}


def analyze_tm_helix_properties(sequence, tm_regions):
    """
    Analyze properties specific to transmembrane helices.
    """
    tm_properties = []
    
    for start, end in tm_regions:
        tm_seq = sequence[start-1:end]
        
        # Calculate hydrophobic moment (simplified)
        hydrophobic_values = {
            'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90,
            'C': 0.29, 'Q': -0.85, 'E': -0.74, 'G': 0.48,
            'H': -0.40, 'I': 1.38, 'L': 1.06, 'K': -1.50,
            'M': 0.64, 'F': 1.19, 'P': 0.12, 'S': -0.18,
            'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
        }
        
        # Calculate average hydrophobicity
        avg_hydrophobicity = sum(hydrophobic_values.get(aa, 0) for aa in tm_seq) / len(tm_seq)
        
        # Check for aromatic residues at interfaces (common in TM helices)
        aromatic_at_start = tm_seq[:3].count('F') + tm_seq[:3].count('W') + tm_seq[:3].count('Y')
        aromatic_at_end = tm_seq[-3:].count('F') + tm_seq[-3:].count('W') + tm_seq[-3:].count('Y')
        
        # Check for positive charges flanking TM (positive-inside rule)
        if start > 5:
            n_flank = sequence[start-6:start-1]
            n_positive = n_flank.count('K') + n_flank.count('R')
        else:
            n_positive = 0
        
        if end < len(sequence) - 5:
            c_flank = sequence[end:end+5]
            c_positive = c_flank.count('K') + c_flank.count('R')
        else:
            c_positive = 0
        
        tm_properties.append({
            'region': f'{start}-{end}',
            'sequence': tm_seq,
            'avg_hydrophobicity': avg_hydrophobicity,
            'aromatic_at_interfaces': aromatic_at_start + aromatic_at_end,
            'n_flank_positive': n_positive,
            'c_flank_positive': c_positive,
            'follows_positive_inside': n_positive > 0 or c_positive > 0
        })
    
    return tm_properties


def plot_structural_features(sequence, ss_predictions, disorder_scores, output_prefix='structure_analysis', protein_name='Protein'):
    """
    Create a comprehensive plot of structural features.
    """
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    
    positions = list(range(1, len(ss_predictions) + 1))
    
    # Plot 1: Secondary structure propensity
    ax1 = axes[0]
    ss_numeric = {'H': 1, 'E': 0, 'C': -1}
    ss_values = [ss_numeric[p] for p in ss_predictions]
    
    ax1.plot(positions, ss_values, 'b-', linewidth=1)
    ax1.fill_between(positions, 0, ss_values, where=[v > 0 for v in ss_values], 
                     alpha=0.3, color='red', label='Helix')
    ax1.fill_between(positions, 0, ss_values, where=[v < 0 for v in ss_values], 
                     alpha=0.3, color='blue', label='Sheet')
    ax1.set_ylabel('Structure Propensity')
    ax1.set_title(f'{protein_name} Secondary Structure Prediction')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Disorder prediction
    ax2 = axes[1]
    disorder_positions = list(range(8, 8 + len(disorder_scores)))  # Centered window
    ax2.plot(disorder_positions, disorder_scores, 'g-', linewidth=1.5)
    ax2.axhline(y=0.2, color='r', linestyle='--', alpha=0.5, label='Disorder threshold')
    ax2.fill_between(disorder_positions, 0.2, disorder_scores, 
                     where=[s > 0.2 for s in disorder_scores],
                     alpha=0.3, color='orange', label='Disordered')
    ax2.set_ylabel('Disorder Score')
    ax2.set_title(f'{protein_name} Disorder Prediction')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Combined view with TM regions
    ax3 = axes[2]
    
    # Mark TM regions
    tm_regions = [(13, 33), (42, 62)]
    for i, (start, end) in enumerate(tm_regions, 1):
        ax3.axvspan(start, end, alpha=0.3, color='cyan', 
                   label=f'TM{i}' if i == 1 else '')
    
    # Mark C-terminal region
    ax3.axvspan(63, len(sequence), alpha=0.2, color='yellow', label='C-terminal')
    
    # Add sequence position reference
    ax3.set_xlim(1, len(sequence))
    ax3.set_ylim(0, 1)
    ax3.set_xlabel('Position')
    ax3.set_ylabel('Regions')
    ax3.set_title(f'{protein_name} Regions Overview')
    ax3.legend()
    
    plt.tight_layout()
    plot_file = f'{output_prefix}_structural_features.png'
    plt.savefig(plot_file, dpi=150)
    plt.close()


@click.command()
@click.argument('fasta_file', type=click.Path(exists=True))
@click.option('--output-prefix', default='structure_analysis', help='Prefix for output files')
@click.option('--uniprot-id', help='UniProt ID for AlphaFold lookup (optional)')
def main(fasta_file, output_prefix, uniprot_id):
    """
    Analyze structural features of a protein sequence from FASTA file.
    
    FASTA_FILE: Path to the input FASTA file containing protein sequence
    """
    # Read the sequence
    with open(fasta_file, 'r') as f:
        record = next(SeqIO.parse(f, "fasta"))
    
    sequence = str(record.seq)
    protein_name = record.id.split('|')[-1] if '|' in record.id else record.id
    
    # Try to extract UniProt ID if not provided
    if not uniprot_id and '|' in record.id:
        parts = record.id.split('|')
        if len(parts) >= 2:
            uniprot_id = parts[1]
    
    print(f"Analyzing {protein_name} structural features")
    print("=" * 60)
    
    # Predict secondary structure
    ss_predictions, helix_scores, sheet_scores = predict_secondary_structure(sequence)
    
    print("\nSecondary Structure Prediction Summary:")
    print("-" * 40)
    helix_count = ss_predictions.count('H')
    sheet_count = ss_predictions.count('E')
    coil_count = ss_predictions.count('C')
    
    print(f"  Helix: {helix_count}/{len(ss_predictions)} ({helix_count/len(ss_predictions)*100:.1f}%)")
    print(f"  Sheet: {sheet_count}/{len(ss_predictions)} ({sheet_count/len(ss_predictions)*100:.1f}%)")
    print(f"  Coil:  {coil_count}/{len(ss_predictions)} ({coil_count/len(ss_predictions)*100:.1f}%)")
    
    # Predict disorder
    disorder_scores = predict_disorder(sequence)
    disordered_regions = []
    in_disorder = False
    start = None
    
    for i, score in enumerate(disorder_scores):
        if score > 0.2 and not in_disorder:
            in_disorder = True
            start = i + 8  # Account for window centering
        elif score <= 0.2 and in_disorder:
            in_disorder = False
            if i + 7 - start > 5:  # Minimum length
                disordered_regions.append((start, i + 7))
    
    print("\nDisorder Prediction:")
    print("-" * 40)
    if disordered_regions:
        for start, end in disordered_regions:
            print(f"  Disordered region: {start}-{end}")
            print(f"    Sequence: {sequence[start-1:end]}")
    else:
        print("  No significant disordered regions predicted")
    
    # Analyze TM helix properties
    tm_regions = [(13, 33), (42, 62)]
    tm_properties = analyze_tm_helix_properties(sequence, tm_regions)
    
    print("\nTransmembrane Helix Analysis:")
    print("-" * 40)
    for tm in tm_properties:
        print(f"\nTM Region {tm['region']}:")
        print(f"  Average hydrophobicity: {tm['avg_hydrophobicity']:.2f}")
        print(f"  Aromatic residues at interfaces: {tm['aromatic_at_interfaces']}")
        print(f"  Positive charges in N-flank: {tm['n_flank_positive']}")
        print(f"  Positive charges in C-flank: {tm['c_flank_positive']}")
        print(f"  Follows positive-inside rule: {tm['follows_positive_inside']}")
    
    # Check AlphaFold
    print("\nChecking AlphaFold Structure Database:")
    print("-" * 40)
    af_result = check_alphafold_structure(uniprot_id) if uniprot_id else {'available': False, 'message': 'No UniProt ID provided'}
    
    if af_result['available']:
        print(f"  AlphaFold structure available!")
        print(f"  URL: {af_result['url']}")
        print(f"  Model URL: {af_result.get('model_url', 'N/A')}")
    else:
        print(f"  No AlphaFold structure found for {af_result['accession']}")
        if 'error' in af_result:
            print(f"  Error: {af_result['error']}")
    
    # Create visualization
    plot_structural_features(sequence, ss_predictions, disorder_scores, output_prefix, protein_name)
    print(f"\nStructural features plot saved as '{output_prefix}_structural_features.png'")
    
    # Save results
    results = {
        'protein': protein_name,
        'protein_length': len(sequence),
        'secondary_structure': {
            'helix_percentage': helix_count/len(ss_predictions)*100,
            'sheet_percentage': sheet_count/len(ss_predictions)*100,
            'coil_percentage': coil_count/len(ss_predictions)*100
        },
        'disordered_regions': disordered_regions,
        'tm_helix_properties': tm_properties,
        'alphafold_structure': af_result
    }
    
    json_file = f'{output_prefix}_results.json'
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to '{json_file}'")


if __name__ == "__main__":
    main()