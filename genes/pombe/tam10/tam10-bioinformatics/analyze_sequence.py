#!/usr/bin/env python3
"""
Comprehensive sequence analysis for protein sequences.
Accepts FASTA input and performs multiple analyses.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import molecular_weight
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from collections import Counter
import json

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze protein sequence properties')
    parser.add_argument('input_fasta', help='Input FASTA file')
    parser.add_argument('--output_dir', default='results', help='Output directory for results')
    parser.add_argument('--window_size', type=int, default=15, help='Window size for sliding window analyses')
    return parser.parse_args()

def basic_composition_analysis(sequence):
    """Analyze basic sequence composition."""
    
    # Use ProtParam for basic analysis
    analyzed_seq = ProtParam.ProteinAnalysis(str(sequence))
    
    # Basic properties
    results = {
        'length': len(sequence),
        'molecular_weight': analyzed_seq.molecular_weight(),
        'aromaticity': analyzed_seq.aromaticity(),
        'instability_index': analyzed_seq.instability_index(),
        'isoelectric_point': analyzed_seq.isoelectric_point(),
        'gravy': analyzed_seq.gravy(),  # Grand average of hydropathy
    }
    
    # Amino acid composition
    aa_percent = analyzed_seq.amino_acids_percent
    results['amino_acid_percent'] = aa_percent
    
    # Count specific residue types
    sequence_str = str(sequence)
    results['basic_residues'] = sequence_str.count('K') + sequence_str.count('R') + sequence_str.count('H')
    results['acidic_residues'] = sequence_str.count('D') + sequence_str.count('E')
    results['polar_residues'] = sum(sequence_str.count(aa) for aa in 'STNQCY')
    results['hydrophobic_residues'] = sum(sequence_str.count(aa) for aa in 'AVILMFYW')
    results['aromatic_residues'] = sum(sequence_str.count(aa) for aa in 'FWY')
    results['small_residues'] = sum(sequence_str.count(aa) for aa in 'AGSV')
    results['charged_residues'] = results['basic_residues'] + results['acidic_residues']
    
    # Calculate percentages
    length = len(sequence)
    results['percent_basic'] = (results['basic_residues'] / length) * 100
    results['percent_acidic'] = (results['acidic_residues'] / length) * 100
    results['percent_polar'] = (results['polar_residues'] / length) * 100
    results['percent_hydrophobic'] = (results['hydrophobic_residues'] / length) * 100
    results['percent_charged'] = (results['charged_residues'] / length) * 100
    
    # Extinction coefficients
    results['extinction_coefficient_reduced'] = analyzed_seq.molar_extinction_coefficient()[0]
    results['extinction_coefficient_oxidized'] = analyzed_seq.molar_extinction_coefficient()[1]
    
    return results

def analyze_charge_distribution(sequence, window_size=15):
    """Analyze charge distribution along the sequence."""
    
    sequence_str = str(sequence)
    positions = []
    local_charges = []
    local_hydrophobicity = []
    
    # Hydrophobicity scale (Kyte-Doolittle)
    hydrophobicity = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    for i in range(len(sequence_str) - window_size + 1):
        window = sequence_str[i:i + window_size]
        
        # Calculate local charge
        positive = window.count('K') + window.count('R') + window.count('H')
        negative = window.count('D') + window.count('E')
        net_charge = positive - negative
        
        # Calculate local hydrophobicity
        hydro_score = sum(hydrophobicity.get(aa, 0) for aa in window) / len(window)
        
        positions.append(i + window_size // 2)
        local_charges.append(net_charge)
        local_hydrophobicity.append(hydro_score)
    
    return {
        'positions': positions,
        'local_charges': local_charges,
        'local_hydrophobicity': local_hydrophobicity,
        'charge_variance': np.var(local_charges),
        'hydrophobicity_variance': np.var(local_hydrophobicity)
    }

def find_motifs_and_patterns(sequence):
    """Search for common motifs and patterns."""
    
    sequence_str = str(sequence)
    patterns_found = []
    
    # Common motifs to search for
    motifs = {
        'Nuclear Localization Signal (NLS) - Classical': r'[KR]{4,}',
        'NLS - Bipartite': r'[KR]{2}.{10,12}[KR]{3,}',
        'Leucine Zipper': r'L.{6}L.{6}L.{6}L',
        'LXXLL motif': r'L..LL',
        'RGD motif': r'RGD',
        'Proline-rich region': r'P{3,}|P.P.P',
        'Polyglutamine': r'Q{5,}',
        'Serine-rich region': r'S{4,}',
        'Acidic cluster': r'[DE]{4,}',
        'Basic cluster': r'[KR]{4,}',
        'KEN box': r'KEN',
        'D-box': r'R..L',
        'PEST region indicator': r'[PEST]{5,}',
        'Coiled-coil heptad pattern': r'[LIMVF].{5}[LIMVF]',
    }
    
    import re
    for motif_name, pattern in motifs.items():
        matches = list(re.finditer(pattern, sequence_str))
        if matches:
            for match in matches:
                patterns_found.append({
                    'motif': motif_name,
                    'pattern': pattern,
                    'start': match.start() + 1,  # 1-based
                    'end': match.end(),
                    'sequence': match.group()
                })
    
    # Check for low complexity regions
    lc_regions = find_low_complexity_regions(sequence_str)
    if lc_regions:
        patterns_found.extend(lc_regions)
    
    # Check for repetitive sequences
    repeats = find_repeats(sequence_str)
    if repeats:
        patterns_found.extend(repeats)
    
    return patterns_found

def find_low_complexity_regions(sequence, window=20, threshold=0.5):
    """Find low complexity regions in the sequence."""
    
    regions = []
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        
        # Calculate Shannon entropy
        aa_counts = Counter(window_seq)
        total = len(window_seq)
        entropy = -sum((count/total) * np.log2(count/total) for count in aa_counts.values())
        max_entropy = np.log2(min(20, total))  # Maximum possible entropy
        
        normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0
        
        if normalized_entropy < threshold:
            regions.append({
                'motif': 'Low complexity region',
                'start': i + 1,
                'end': i + window,
                'entropy': normalized_entropy,
                'sequence': window_seq
            })
    
    # Merge overlapping regions
    merged = []
    for region in regions:
        if merged and region['start'] <= merged[-1]['end']:
            merged[-1]['end'] = max(merged[-1]['end'], region['end'])
            merged[-1]['sequence'] = sequence[merged[-1]['start']-1:merged[-1]['end']]
        else:
            merged.append(region)
    
    return merged

def find_repeats(sequence, min_repeat_length=3, min_occurrences=3):
    """Find repetitive sequences."""
    
    repeats = []
    found_patterns = set()
    
    for length in range(min_repeat_length, min(10, len(sequence)//min_occurrences + 1)):
        for i in range(len(sequence) - length + 1):
            pattern = sequence[i:i + length]
            
            if pattern in found_patterns:
                continue
                
            # Count non-overlapping occurrences
            count = 0
            pos = 0
            positions = []
            while pos < len(sequence):
                index = sequence.find(pattern, pos)
                if index == -1:
                    break
                count += 1
                positions.append(index + 1)  # 1-based
                pos = index + length
            
            if count >= min_occurrences:
                found_patterns.add(pattern)
                repeats.append({
                    'motif': f'Repeat ({count}x)',
                    'pattern': pattern,
                    'occurrences': count,
                    'positions': positions,
                    'sequence': pattern
                })
    
    return repeats

def predict_coiled_coils(sequence):
    """Simple coiled-coil prediction based on heptad repeats."""
    
    sequence_str = str(sequence)
    window_size = 28  # 4 heptads
    scores = []
    
    # Heptad positions a-g (0-6)
    # Positions a(0) and d(3) are typically hydrophobic
    hydrophobic = set('LIMVFYA')
    
    for i in range(len(sequence_str) - window_size + 1):
        window = sequence_str[i:i + window_size]
        score = 0
        
        for j in range(0, window_size, 7):
            if j < len(window):
                # Check position 'a'
                if window[j] in hydrophobic:
                    score += 1
                # Check position 'd' 
                if j + 3 < len(window) and window[j + 3] in hydrophobic:
                    score += 1
        
        # Normalize score
        max_score = (window_size // 7) * 2
        normalized_score = score / max_score if max_score > 0 else 0
        scores.append((i + window_size // 2, normalized_score))
    
    # Find regions with high scores
    coiled_coil_regions = []
    threshold = 0.5
    in_region = False
    start = 0
    
    for pos, score in scores:
        if score >= threshold and not in_region:
            start = pos
            in_region = True
        elif score < threshold and in_region:
            if pos - start >= 14:  # Minimum 2 heptads
                coiled_coil_regions.append({
                    'start': start - window_size // 2 + 1,
                    'end': pos + window_size // 2,
                    'max_score': max(s for p, s in scores if start <= p <= pos),
                    'length': pos - start + window_size
                })
            in_region = False
    
    return coiled_coil_regions, scores

def analyze_disorder_tendency(sequence):
    """Simple disorder prediction based on composition."""
    
    sequence_str = str(sequence)
    window_size = 15
    disorder_scores = []
    
    # Disorder-promoting residues
    disorder_promoting = set('KRSQENP')
    order_promoting = set('WFYILVM')
    
    for i in range(len(sequence_str) - window_size + 1):
        window = sequence_str[i:i + window_size]
        
        disorder_count = sum(1 for aa in window if aa in disorder_promoting)
        order_count = sum(1 for aa in window if aa in order_promoting)
        
        # Simple disorder score
        disorder_score = (disorder_count - order_count) / window_size
        disorder_scores.append((i + window_size // 2, disorder_score))
    
    # Find disordered regions
    disordered_regions = []
    threshold = 0.2
    in_region = False
    start = 0
    
    for pos, score in disorder_scores:
        if score >= threshold and not in_region:
            start = pos
            in_region = True
        elif score < threshold and in_region:
            if pos - start >= 10:  # Minimum length for disorder
                disordered_regions.append({
                    'start': start - window_size // 2 + 1,
                    'end': pos + window_size // 2,
                    'max_score': max(s for p, s in disorder_scores if start <= p <= pos),
                    'length': pos - start + window_size
                })
            in_region = False
    
    # Check if still in region at the end
    if in_region and disorder_scores:
        last_pos = disorder_scores[-1][0]
        if last_pos - start >= 10:
            disordered_regions.append({
                'start': start - window_size // 2 + 1,
                'end': len(sequence_str),
                'max_score': max(s for p, s in disorder_scores if start <= p),
                'length': len(sequence_str) - start + window_size // 2
            })
    
    return disordered_regions, disorder_scores

def plot_sequence_properties(sequence, charge_data, coiled_coil_scores, disorder_scores, output_dir):
    """Generate plots for sequence properties."""
    
    fig, axes = plt.subplots(4, 1, figsize=(12, 10))
    
    # Plot 1: Charge distribution
    ax1 = axes[0]
    ax1.plot(charge_data['positions'], charge_data['local_charges'], 'b-', alpha=0.7)
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_ylabel('Net Charge')
    ax1.set_title('Local Charge Distribution (15 aa window)')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Hydrophobicity
    ax2 = axes[1]
    ax2.plot(charge_data['positions'], charge_data['local_hydrophobicity'], 'g-', alpha=0.7)
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax2.set_ylabel('Hydrophobicity')
    ax2.set_title('Hydrophobicity Profile (Kyte-Doolittle)')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Coiled-coil propensity
    ax3 = axes[2]
    if coiled_coil_scores:
        positions, scores = zip(*coiled_coil_scores)
        ax3.plot(positions, scores, 'r-', alpha=0.7)
        ax3.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Threshold')
    ax3.set_ylabel('CC Propensity')
    ax3.set_title('Coiled-Coil Propensity')
    ax3.set_ylim(0, 1)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Disorder tendency
    ax4 = axes[3]
    if disorder_scores:
        positions, scores = zip(*disorder_scores)
        ax4.plot(positions, scores, 'm-', alpha=0.7)
        ax4.axhline(y=0.2, color='gray', linestyle='--', alpha=0.5, label='Threshold')
    ax4.set_ylabel('Disorder Score')
    ax4.set_xlabel('Position')
    ax4.set_title('Disorder Tendency')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'sequence_properties.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Create amino acid composition plot
    fig, ax = plt.subplots(figsize=(10, 6))
    sequence_str = str(sequence)
    aa_counts = Counter(sequence_str)
    
    # Group amino acids by properties
    aa_groups = {
        'Hydrophobic': 'AVILMFYW',
        'Polar': 'STNQCY',
        'Basic': 'KRH',
        'Acidic': 'DE',
        'Special': 'GP'
    }
    
    group_counts = {}
    for group, aas in aa_groups.items():
        group_counts[group] = sum(aa_counts.get(aa, 0) for aa in aas)
    
    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99', '#ff99cc']
    plt.pie(group_counts.values(), labels=group_counts.keys(), colors=colors, autopct='%1.1f%%')
    plt.title('Amino Acid Composition by Property')
    plt.savefig(output_dir / 'aa_composition.png', dpi=150, bbox_inches='tight')
    plt.close()

def main():
    args = parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Read the sequence
    record = next(SeqIO.parse(args.input_fasta, "fasta"))
    sequence = record.seq
    
    print(f"Analyzing sequence: {record.id}")
    print(f"Length: {len(sequence)} amino acids")
    
    # Perform analyses
    print("\n1. Basic composition analysis...")
    composition = basic_composition_analysis(sequence)
    
    print("2. Charge and hydrophobicity distribution...")
    charge_data = analyze_charge_distribution(sequence, args.window_size)
    
    print("3. Motif and pattern search...")
    motifs = find_motifs_and_patterns(sequence)
    
    print("4. Coiled-coil prediction...")
    coiled_coils, cc_scores = predict_coiled_coils(sequence)
    
    print("5. Disorder prediction...")
    disorder_regions, disorder_scores = analyze_disorder_tendency(sequence)
    
    print("6. Generating plots...")
    plot_sequence_properties(sequence, charge_data, cc_scores, disorder_scores, output_dir)
    
    # Save results
    results = {
        'sequence_id': record.id,
        'sequence_description': record.description,
        'composition': composition,
        'charge_distribution': {
            'charge_variance': charge_data['charge_variance'],
            'hydrophobicity_variance': charge_data['hydrophobicity_variance']
        },
        'motifs_found': motifs,
        'coiled_coil_regions': coiled_coils,
        'disorder_regions': disorder_regions
    }
    
    # Save JSON results
    with open(output_dir / 'analysis_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    # Generate text report
    with open(output_dir / 'analysis_report.txt', 'w') as f:
        f.write(f"Sequence Analysis Report\n")
        f.write(f"========================\n\n")
        f.write(f"Sequence ID: {record.id}\n")
        f.write(f"Description: {record.description}\n")
        f.write(f"Length: {len(sequence)} amino acids\n\n")
        
        f.write(f"Basic Properties:\n")
        f.write(f"-----------------\n")
        f.write(f"Molecular Weight: {composition['molecular_weight']:.2f} Da\n")
        f.write(f"Isoelectric Point: {composition['isoelectric_point']:.2f}\n")
        f.write(f"Aromaticity: {composition['aromaticity']:.3f}\n")
        f.write(f"Instability Index: {composition['instability_index']:.2f}\n")
        f.write(f"GRAVY: {composition['gravy']:.3f}\n\n")
        
        f.write(f"Residue Composition:\n")
        f.write(f"-------------------\n")
        f.write(f"Basic residues (K,R,H): {composition['basic_residues']} ({composition['percent_basic']:.1f}%)\n")
        f.write(f"Acidic residues (D,E): {composition['acidic_residues']} ({composition['percent_acidic']:.1f}%)\n")
        f.write(f"Polar residues: {composition['polar_residues']} ({composition['percent_polar']:.1f}%)\n")
        f.write(f"Hydrophobic residues: {composition['hydrophobic_residues']} ({composition['percent_hydrophobic']:.1f}%)\n")
        f.write(f"Charged residues: {composition['charged_residues']} ({composition['percent_charged']:.1f}%)\n\n")
        
        f.write(f"Motifs and Patterns Found:\n")
        f.write(f"-------------------------\n")
        if motifs:
            for motif in motifs:
                f.write(f"- {motif['motif']}")
                if 'start' in motif:
                    f.write(f" at position {motif['start']}-{motif['end']}")
                if 'sequence' in motif:
                    f.write(f": {motif['sequence']}")
                f.write("\n")
        else:
            f.write("No significant motifs found\n")
        f.write("\n")
        
        f.write(f"Coiled-Coil Regions:\n")
        f.write(f"-------------------\n")
        if coiled_coils:
            for cc in coiled_coils:
                f.write(f"- Position {cc['start']}-{cc['end']} (score: {cc['max_score']:.2f})\n")
        else:
            f.write("No significant coiled-coil regions detected\n")
        f.write("\n")
        
        f.write(f"Disorder Regions:\n")
        f.write(f"----------------\n")
        if disorder_regions:
            for dr in disorder_regions:
                f.write(f"- Position {dr['start']}-{dr['end']} (score: {dr['max_score']:.2f})\n")
        else:
            f.write("No significant disordered regions detected\n")
    
    print(f"\nAnalysis complete! Results saved to {output_dir}/")
    print(f"- analysis_results.json: Detailed results in JSON format")
    print(f"- analysis_report.txt: Human-readable report")
    print(f"- sequence_properties.png: Visualization of sequence properties")
    print(f"- aa_composition.png: Amino acid composition chart")

if __name__ == "__main__":
    main()