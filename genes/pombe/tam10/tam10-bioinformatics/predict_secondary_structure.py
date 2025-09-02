#!/usr/bin/env python3
"""
Secondary structure prediction using simple methods.
Uses Chou-Fasman and GOR methods for prediction.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import json

def parse_args():
    parser = argparse.ArgumentParser(description='Predict protein secondary structure')
    parser.add_argument('input_fasta', help='Input FASTA file')
    parser.add_argument('--output_dir', default='results', help='Output directory for results')
    parser.add_argument('--window_size', type=int, default=17, help='Window size for GOR method')
    return parser.parse_args()

# Chou-Fasman parameters
CHOU_FASMAN_ALPHA = {
    'A': 1.42, 'C': 0.70, 'D': 1.01, 'E': 1.51, 'F': 1.13,
    'G': 0.57, 'H': 1.00, 'I': 1.08, 'K': 1.16, 'L': 1.21,
    'M': 1.45, 'N': 0.67, 'P': 0.57, 'Q': 1.11, 'R': 0.98,
    'S': 0.77, 'T': 0.83, 'V': 1.06, 'W': 1.08, 'Y': 0.69
}

CHOU_FASMAN_BETA = {
    'A': 0.83, 'C': 1.19, 'D': 0.54, 'E': 0.37, 'F': 1.38,
    'G': 0.75, 'H': 0.87, 'I': 1.60, 'K': 0.74, 'L': 1.30,
    'M': 1.05, 'N': 0.89, 'P': 0.55, 'Q': 1.10, 'R': 0.93,
    'S': 0.75, 'T': 1.19, 'V': 1.70, 'W': 1.37, 'Y': 1.47
}

CHOU_FASMAN_TURN = {
    'A': 0.66, 'C': 1.19, 'D': 1.46, 'E': 0.74, 'F': 0.60,
    'G': 1.56, 'H': 0.95, 'I': 0.47, 'K': 1.01, 'L': 0.59,
    'M': 0.60, 'N': 1.56, 'P': 1.52, 'Q': 0.98, 'R': 0.95,
    'S': 1.43, 'T': 0.96, 'V': 0.50, 'W': 0.96, 'Y': 1.14
}

def chou_fasman_prediction(sequence):
    """Predict secondary structure using Chou-Fasman method."""
    
    sequence_str = str(sequence)
    length = len(sequence_str)
    
    # Initialize scores
    alpha_scores = []
    beta_scores = []
    turn_scores = []
    
    # Calculate scores for each position using sliding window
    window_size = 6
    for i in range(length):
        # Get window around position i
        start = max(0, i - window_size // 2)
        end = min(length, i + window_size // 2 + 1)
        window = sequence_str[start:end]
        
        # Calculate average propensities
        alpha_score = np.mean([CHOU_FASMAN_ALPHA.get(aa, 1.0) for aa in window])
        beta_score = np.mean([CHOU_FASMAN_BETA.get(aa, 1.0) for aa in window])
        turn_score = np.mean([CHOU_FASMAN_TURN.get(aa, 1.0) for aa in window])
        
        alpha_scores.append(alpha_score)
        beta_scores.append(beta_score)
        turn_scores.append(turn_score)
    
    # Assign secondary structure
    structure = []
    for i in range(length):
        scores = {
            'H': alpha_scores[i],
            'E': beta_scores[i],
            'T': turn_scores[i]
        }
        
        # Assign structure based on highest score above threshold
        if scores['H'] > 1.03 and scores['H'] > scores['E']:
            structure.append('H')
        elif scores['E'] > 1.05 and scores['E'] > scores['H']:
            structure.append('E')
        elif scores['T'] > 1.00:
            structure.append('T')
        else:
            structure.append('C')  # Coil
    
    return {
        'structure': ''.join(structure),
        'alpha_scores': alpha_scores,
        'beta_scores': beta_scores,
        'turn_scores': turn_scores
    }

def gor_method_prediction(sequence, window_size=17):
    """Simplified GOR method for secondary structure prediction."""
    
    sequence_str = str(sequence)
    length = len(sequence_str)
    
    # Simplified GOR parameters (frequencies)
    # These are approximations
    helix_aa = set('AELM')
    sheet_aa = set('VIFY')
    turn_aa = set('NGPS')
    
    structure = []
    helix_scores = []
    sheet_scores = []
    turn_scores = []
    
    for i in range(length):
        # Get window
        start = max(0, i - window_size // 2)
        end = min(length, i + window_size // 2 + 1)
        window = sequence_str[start:end]
        
        # Count amino acid types in window
        helix_count = sum(1 for aa in window if aa in helix_aa)
        sheet_count = sum(1 for aa in window if aa in sheet_aa)
        turn_count = sum(1 for aa in window if aa in turn_aa)
        
        # Normalize by window length
        window_len = len(window)
        helix_score = helix_count / window_len
        sheet_score = sheet_count / window_len
        turn_score = turn_count / window_len
        
        helix_scores.append(helix_score)
        sheet_scores.append(sheet_score)
        turn_scores.append(turn_score)
        
        # Assign structure
        if helix_score > sheet_score and helix_score > turn_score and helix_score > 0.25:
            structure.append('H')
        elif sheet_score > helix_score and sheet_score > turn_score and sheet_score > 0.25:
            structure.append('E')
        elif turn_score > 0.2:
            structure.append('T')
        else:
            structure.append('C')
    
    return {
        'structure': ''.join(structure),
        'helix_scores': helix_scores,
        'sheet_scores': sheet_scores,
        'turn_scores': turn_scores
    }

def consensus_prediction(predictions):
    """Generate consensus secondary structure from multiple predictions."""
    
    structures = [pred['structure'] for pred in predictions]
    consensus = []
    
    for i in range(len(structures[0])):
        votes = [s[i] for s in structures]
        # Most common structure at this position
        from collections import Counter
        count = Counter(votes)
        consensus.append(count.most_common(1)[0][0])
    
    return ''.join(consensus)

def analyze_secondary_structure(structure):
    """Analyze the secondary structure composition."""
    
    from collections import Counter
    
    counts = Counter(structure)
    total = len(structure)
    
    composition = {
        'helix': counts.get('H', 0),
        'sheet': counts.get('E', 0),
        'turn': counts.get('T', 0),
        'coil': counts.get('C', 0),
        'helix_percent': (counts.get('H', 0) / total) * 100,
        'sheet_percent': (counts.get('E', 0) / total) * 100,
        'turn_percent': (counts.get('T', 0) / total) * 100,
        'coil_percent': (counts.get('C', 0) / total) * 100
    }
    
    # Find continuous secondary structure elements
    elements = []
    current_type = structure[0]
    start = 0
    
    for i in range(1, len(structure)):
        if structure[i] != current_type:
            if current_type != 'C':  # Don't record coil regions
                elements.append({
                    'type': current_type,
                    'start': start + 1,  # 1-based
                    'end': i,
                    'length': i - start
                })
            current_type = structure[i]
            start = i
    
    # Don't forget the last element
    if current_type != 'C':
        elements.append({
            'type': current_type,
            'start': start + 1,
            'end': len(structure),
            'length': len(structure) - start
        })
    
    return composition, elements

def plot_secondary_structure(sequence, cf_pred, gor_pred, consensus, output_dir):
    """Generate plots for secondary structure predictions."""
    
    length = len(sequence)
    positions = list(range(1, length + 1))
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 8))
    
    # Plot 1: Chou-Fasman scores
    ax1 = axes[0]
    ax1.plot(positions, cf_pred['alpha_scores'], 'r-', label='Alpha helix', alpha=0.7)
    ax1.plot(positions, cf_pred['beta_scores'], 'b-', label='Beta sheet', alpha=0.7)
    ax1.plot(positions, cf_pred['turn_scores'], 'g-', label='Turn', alpha=0.7)
    ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.3)
    ax1.set_ylabel('Propensity')
    ax1.set_title('Chou-Fasman Secondary Structure Propensities')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: GOR scores
    ax2 = axes[1]
    ax2.plot(positions, gor_pred['helix_scores'], 'r-', label='Helix', alpha=0.7)
    ax2.plot(positions, gor_pred['sheet_scores'], 'b-', label='Sheet', alpha=0.7)
    ax2.plot(positions, gor_pred['turn_scores'], 'g-', label='Turn', alpha=0.7)
    ax2.set_ylabel('Score')
    ax2.set_title('GOR Method Scores')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Structure visualization
    ax3 = axes[2]
    
    # Convert structure to numeric for plotting
    structure_map = {'H': 3, 'E': 2, 'T': 1, 'C': 0}
    consensus_numeric = [structure_map[s] for s in consensus]
    
    # Create color-coded bar plot
    colors = []
    for s in consensus:
        if s == 'H':
            colors.append('red')
        elif s == 'E':
            colors.append('blue')
        elif s == 'T':
            colors.append('green')
        else:
            colors.append('gray')
    
    ax3.bar(positions, [1]*length, color=colors, width=1.0)
    ax3.set_ylabel('Structure')
    ax3.set_xlabel('Position')
    ax3.set_title('Consensus Secondary Structure (Red=Helix, Blue=Sheet, Green=Turn, Gray=Coil)')
    ax3.set_ylim(0, 1.5)
    ax3.set_yticks([])
    
    plt.tight_layout()
    plt.savefig(output_dir / 'secondary_structure.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Create a summary pie chart
    fig, ax = plt.subplots(figsize=(8, 8))
    
    comp, _ = analyze_secondary_structure(consensus)
    labels = ['Helix', 'Sheet', 'Turn', 'Coil']
    sizes = [comp['helix'], comp['sheet'], comp['turn'], comp['coil']]
    colors = ['red', 'blue', 'green', 'gray']
    
    # Filter out zero values
    non_zero = [(l, s, c) for l, s, c in zip(labels, sizes, colors) if s > 0]
    if non_zero:
        labels, sizes, colors = zip(*non_zero)
        ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    
    ax.set_title('Secondary Structure Composition')
    plt.savefig(output_dir / 'ss_composition.png', dpi=150, bbox_inches='tight')
    plt.close()

def main():
    args = parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Read sequence
    record = next(SeqIO.parse(args.input_fasta, "fasta"))
    sequence = record.seq
    
    print(f"Predicting secondary structure for: {record.id}")
    print(f"Sequence length: {len(sequence)} amino acids")
    
    # Perform predictions
    print("\n1. Chou-Fasman prediction...")
    cf_prediction = chou_fasman_prediction(sequence)
    
    print("2. GOR method prediction...")
    gor_prediction = gor_method_prediction(sequence, args.window_size)
    
    print("3. Generating consensus...")
    consensus = consensus_prediction([cf_prediction, gor_prediction])
    
    print("4. Analyzing structure composition...")
    composition, elements = analyze_secondary_structure(consensus)
    
    print("5. Generating plots...")
    plot_secondary_structure(sequence, cf_prediction, gor_prediction, consensus, output_dir)
    
    # Save results
    results = {
        'sequence_id': record.id,
        'sequence_length': len(sequence),
        'chou_fasman_structure': cf_prediction['structure'],
        'gor_structure': gor_prediction['structure'],
        'consensus_structure': consensus,
        'composition': composition,
        'structural_elements': elements
    }
    
    # Save JSON results
    with open(output_dir / 'ss_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save text report
    with open(output_dir / 'ss_report.txt', 'w') as f:
        f.write(f"Secondary Structure Prediction Report\n")
        f.write(f"=====================================\n\n")
        f.write(f"Sequence: {record.id}\n")
        f.write(f"Length: {len(sequence)} aa\n\n")
        
        f.write(f"Composition:\n")
        f.write(f"-----------\n")
        f.write(f"Alpha helix: {composition['helix']} residues ({composition['helix_percent']:.1f}%)\n")
        f.write(f"Beta sheet: {composition['sheet']} residues ({composition['sheet_percent']:.1f}%)\n")
        f.write(f"Turn: {composition['turn']} residues ({composition['turn_percent']:.1f}%)\n")
        f.write(f"Coil/Loop: {composition['coil']} residues ({composition['coil_percent']:.1f}%)\n\n")
        
        f.write(f"Structural Elements:\n")
        f.write(f"-------------------\n")
        if elements:
            for elem in sorted(elements, key=lambda x: x['start']):
                struct_type = {'H': 'Helix', 'E': 'Sheet', 'T': 'Turn'}.get(elem['type'], elem['type'])
                f.write(f"- {struct_type}: positions {elem['start']}-{elem['end']} (length: {elem['length']})\n")
        else:
            f.write("No defined structural elements found\n")
        
        f.write(f"\nConsensus Structure:\n")
        f.write(f"-------------------\n")
        # Write structure in blocks of 60
        for i in range(0, len(consensus), 60):
            f.write(f"{i+1:4d} {consensus[i:i+60]}\n")
    
    print(f"\nSecondary structure prediction complete!")
    print(f"Results saved to {output_dir}/")
    print(f"- ss_results.json: Detailed results")
    print(f"- ss_report.txt: Human-readable report")
    print(f"- secondary_structure.png: Structure predictions")
    print(f"- ss_composition.png: Composition chart")

if __name__ == "__main__":
    main()