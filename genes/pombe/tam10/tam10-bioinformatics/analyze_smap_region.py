#!/usr/bin/env python3
"""
Analyze potential SMAP domain region in protein sequences.
SMAP (Small acidic protein) domains are typically 50-60 aa regions.
"""

import argparse
from pathlib import Path
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import json

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze potential SMAP domain regions')
    parser.add_argument('input_fasta', help='Input FASTA file')
    parser.add_argument('--output_dir', default='results', help='Output directory for results')
    parser.add_argument('--window_size', type=int, default=60, help='Window size for SMAP search')
    return parser.parse_args()

def find_smap_like_regions(sequence, window_size=60):
    """
    Search for SMAP-like regions based on typical characteristics:
    - Small size (50-60 aa)
    - Often acidic
    - May have specific sequence patterns
    """
    
    sequence_str = str(sequence)
    length = len(sequence_str)
    
    potential_regions = []
    
    # Scan with sliding window
    for i in range(max(0, length - window_size + 1)):
        window = sequence_str[i:i + window_size]
        
        # Calculate properties
        acidic_count = window.count('D') + window.count('E')
        basic_count = window.count('K') + window.count('R') + window.count('H')
        net_charge = basic_count - acidic_count
        
        # Calculate hydrophobicity
        hydrophobic_count = sum(window.count(aa) for aa in 'AVILMFYW')
        hydrophilic_count = sum(window.count(aa) for aa in 'KRHDENQST')
        
        # SMAP domains often have mixed properties
        properties = {
            'start': i + 1,  # 1-based
            'end': i + window_size,
            'sequence': window,
            'acidic_residues': acidic_count,
            'basic_residues': basic_count,
            'net_charge': net_charge,
            'percent_acidic': (acidic_count / len(window)) * 100,
            'percent_basic': (basic_count / len(window)) * 100,
            'hydrophobic_residues': hydrophobic_count,
            'hydrophilic_residues': hydrophilic_count,
            'hydrophobic_percent': (hydrophobic_count / len(window)) * 100
        }
        
        # Score the region for SMAP-like characteristics
        score = 0
        
        # SMAP domains can be acidic
        if properties['percent_acidic'] > 15:
            score += 2
        elif properties['percent_acidic'] > 10:
            score += 1
            
        # But some have mixed charge
        if abs(properties['net_charge']) < 5:
            score += 1
            
        # Check for presence of aromatic residues (often present in SMAP)
        aromatic_count = sum(window.count(aa) for aa in 'FWY')
        if aromatic_count >= 2:
            score += 1
            
        # Check for proline content (sometimes enriched)
        proline_count = window.count('P')
        if proline_count >= 2:
            score += 1
        
        properties['smap_score'] = score
        potential_regions.append(properties)
    
    # Sort by SMAP score
    potential_regions.sort(key=lambda x: x['smap_score'], reverse=True)
    
    return potential_regions

def analyze_domain_boundaries(sequence, candidate_region):
    """
    Analyze potential domain boundaries around a candidate region.
    """
    
    sequence_str = str(sequence)
    start = candidate_region['start'] - 1  # Convert to 0-based
    end = candidate_region['end']
    
    # Extend analysis window
    extended_start = max(0, start - 20)
    extended_end = min(len(sequence_str), end + 20)
    
    extended_seq = sequence_str[extended_start:extended_end]
    
    # Look for potential boundary markers
    boundaries = {
        'n_terminal_extension': sequence_str[extended_start:start] if extended_start < start else '',
        'c_terminal_extension': sequence_str[end:extended_end] if end < extended_end else '',
        'has_n_terminal_disorder': False,
        'has_c_terminal_disorder': False,
        'has_proline_boundary': False,
        'has_charge_transition': False
    }
    
    # Check for disorder-prone regions at boundaries
    disorder_aa = set('KRSQENP')
    
    if boundaries['n_terminal_extension']:
        n_disorder = sum(1 for aa in boundaries['n_terminal_extension'] if aa in disorder_aa)
        boundaries['has_n_terminal_disorder'] = n_disorder / len(boundaries['n_terminal_extension']) > 0.4
        
    if boundaries['c_terminal_extension']:
        c_disorder = sum(1 for aa in boundaries['c_terminal_extension'] if aa in disorder_aa)
        boundaries['has_c_terminal_disorder'] = c_disorder / len(boundaries['c_terminal_extension']) > 0.4
    
    # Check for proline boundaries (common domain separators)
    if start > 0:
        boundaries['has_proline_boundary'] = sequence_str[start-1] == 'P' or sequence_str[start] == 'P'
    if end < len(sequence_str):
        boundaries['has_proline_boundary'] = boundaries['has_proline_boundary'] or sequence_str[end-1] == 'P'
    
    # Check for charge transitions
    if boundaries['n_terminal_extension'] and boundaries['c_terminal_extension']:
        n_charge = sum(1 if aa in 'KRH' else -1 if aa in 'DE' else 0 
                      for aa in boundaries['n_terminal_extension'])
        c_charge = sum(1 if aa in 'KRH' else -1 if aa in 'DE' else 0 
                      for aa in boundaries['c_terminal_extension'])
        boundaries['has_charge_transition'] = abs(n_charge - c_charge) > 5
    
    return boundaries

def compare_with_known_smap(candidate_seq):
    """
    Compare candidate sequence with known SMAP domain characteristics.
    Note: This is a simplified comparison as we don't have a full SMAP profile.
    """
    
    # Known characteristics of SMAP domains from literature
    # These are approximations
    typical_features = {
        'length_range': (45, 65),
        'common_motifs': [
            'D.{2,4}[DE]',  # Acidic clusters
            'W.{2,5}[YF]',  # Aromatic patterns
        ],
        'enriched_aa': 'DEWY',
        'depleted_aa': 'C'
    }
    
    score = 0
    features_found = []
    
    # Check length
    if typical_features['length_range'][0] <= len(candidate_seq) <= typical_features['length_range'][1]:
        score += 1
        features_found.append('appropriate_length')
    
    # Check for enriched amino acids
    enriched_count = sum(candidate_seq.count(aa) for aa in typical_features['enriched_aa'])
    if enriched_count / len(candidate_seq) > 0.15:
        score += 1
        features_found.append('enriched_residues')
    
    # Check for depleted amino acids
    depleted_count = sum(candidate_seq.count(aa) for aa in typical_features['depleted_aa'])
    if depleted_count == 0:
        score += 1
        features_found.append('lacks_cysteine')
    
    # Check for motifs
    import re
    for motif in typical_features['common_motifs']:
        if re.search(motif, candidate_seq):
            score += 1
            features_found.append(f'motif_{motif}')
    
    return score, features_found

def plot_smap_analysis(sequence, regions, output_dir):
    """Generate plots for SMAP region analysis."""
    
    if not regions:
        return
    
    # Take top 3 candidate regions
    top_regions = regions[:min(3, len(regions))]
    
    fig, axes = plt.subplots(len(top_regions), 1, figsize=(12, 4 * len(top_regions)))
    
    if len(top_regions) == 1:
        axes = [axes]
    
    for idx, region in enumerate(top_regions):
        ax = axes[idx]
        
        # Prepare data for the region
        window_seq = region['sequence']
        positions = list(range(region['start'], region['end'] + 1))
        
        # Calculate per-residue properties
        colors = []
        heights = []
        
        for aa in window_seq:
            if aa in 'DE':
                colors.append('red')
                heights.append(-1)  # Acidic below
            elif aa in 'KRH':
                colors.append('blue')
                heights.append(1)  # Basic above
            elif aa in 'FWY':
                colors.append('purple')
                heights.append(0.5)  # Aromatic
            elif aa in 'AVILM':
                colors.append('green')
                heights.append(0.3)  # Hydrophobic
            else:
                colors.append('gray')
                heights.append(0.1)  # Other
        
        # Create bar plot
        ax.bar(positions, heights, color=colors, width=1.0)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax.set_xlabel('Position')
        ax.set_ylabel('Property')
        ax.set_title(f'Region {region["start"]}-{region["end"]} (SMAP score: {region["smap_score"]})')
        ax.set_xlim(region['start'] - 1, region['end'] + 1)
        
        # Add legend
        if idx == 0:
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor='red', label='Acidic'),
                Patch(facecolor='blue', label='Basic'),
                Patch(facecolor='purple', label='Aromatic'),
                Patch(facecolor='green', label='Hydrophobic'),
                Patch(facecolor='gray', label='Other')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'smap_regions.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Create a summary plot of SMAP scores across the sequence
    fig, ax = plt.subplots(figsize=(12, 4))
    
    all_positions = [r['start'] + r['end'] // 2 for r in regions]
    all_scores = [r['smap_score'] for r in regions]
    
    ax.scatter(all_positions, all_scores, alpha=0.6, s=50)
    ax.set_xlabel('Position (center of window)')
    ax.set_ylabel('SMAP Score')
    ax.set_title('SMAP-like Score Distribution Across Sequence')
    ax.grid(True, alpha=0.3)
    
    # Highlight top regions
    for region in top_regions:
        center = (region['start'] + region['end']) // 2
        ax.scatter(center, region['smap_score'], color='red', s=100, marker='*', 
                  label=f"Top region {region['start']}-{region['end']}")
    
    if top_regions:
        ax.legend()
    
    plt.savefig(output_dir / 'smap_score_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()

def main():
    args = parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Read sequence
    record = next(SeqIO.parse(args.input_fasta, "fasta"))
    sequence = record.seq
    
    print(f"Analyzing potential SMAP regions in: {record.id}")
    print(f"Sequence length: {len(sequence)} amino acids")
    
    # Find potential SMAP regions
    print(f"\n1. Searching for SMAP-like regions (window size: {args.window_size})...")
    regions = find_smap_like_regions(sequence, args.window_size)
    
    # Analyze top candidates
    print("2. Analyzing top candidate regions...")
    top_candidates = regions[:3] if regions else []
    
    detailed_analysis = []
    for i, region in enumerate(top_candidates, 1):
        print(f"   Candidate {i}: positions {region['start']}-{region['end']} (score: {region['smap_score']})")
        
        # Analyze boundaries
        boundaries = analyze_domain_boundaries(sequence, region)
        
        # Compare with known SMAP characteristics
        smap_similarity, features = compare_with_known_smap(region['sequence'])
        
        analysis = {
            'rank': i,
            'region': region,
            'boundaries': boundaries,
            'smap_similarity_score': smap_similarity,
            'smap_features_found': features
        }
        detailed_analysis.append(analysis)
    
    # Generate plots
    print("3. Generating visualizations...")
    plot_smap_analysis(sequence, regions, output_dir)
    
    # Save results
    results = {
        'sequence_id': record.id,
        'sequence_length': len(sequence),
        'window_size': args.window_size,
        'total_regions_analyzed': len(regions),
        'top_candidates': detailed_analysis
    }
    
    # Save JSON results
    with open(output_dir / 'smap_analysis.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    # Generate report
    with open(output_dir / 'smap_report.txt', 'w') as f:
        f.write(f"SMAP Domain Analysis Report\n")
        f.write(f"===========================\n\n")
        f.write(f"Sequence: {record.id}\n")
        f.write(f"Length: {len(sequence)} aa\n")
        f.write(f"Window size: {args.window_size} aa\n\n")
        
        if top_candidates:
            f.write(f"Top Candidate Regions:\n")
            f.write(f"---------------------\n\n")
            
            for analysis in detailed_analysis:
                region = analysis['region']
                f.write(f"Candidate {analysis['rank']}:\n")
                f.write(f"  Position: {region['start']}-{region['end']}\n")
                f.write(f"  SMAP score: {region['smap_score']}/5\n")
                f.write(f"  Net charge: {region['net_charge']}\n")
                f.write(f"  Acidic residues: {region['acidic_residues']} ({region['percent_acidic']:.1f}%)\n")
                f.write(f"  Basic residues: {region['basic_residues']} ({region['percent_basic']:.1f}%)\n")
                f.write(f"  Hydrophobic: {region['hydrophobic_percent']:.1f}%\n")
                
                f.write(f"  SMAP similarity: {analysis['smap_similarity_score']}/5\n")
                if analysis['smap_features_found']:
                    f.write(f"  Features found: {', '.join(analysis['smap_features_found'])}\n")
                
                f.write(f"  Boundary analysis:\n")
                boundaries = analysis['boundaries']
                if boundaries['has_n_terminal_disorder']:
                    f.write(f"    - N-terminal disorder detected\n")
                if boundaries['has_c_terminal_disorder']:
                    f.write(f"    - C-terminal disorder detected\n")
                if boundaries['has_proline_boundary']:
                    f.write(f"    - Proline at domain boundary\n")
                if boundaries['has_charge_transition']:
                    f.write(f"    - Charge transition at boundaries\n")
                
                f.write(f"  Sequence:\n")
                f.write(f"    {region['sequence']}\n\n")
        else:
            f.write("No significant SMAP-like regions detected.\n")
        
        f.write(f"\nNotes:\n")
        f.write(f"------\n")
        f.write(f"- SMAP domains are typically 50-60 aa regions\n")
        f.write(f"- They often contain acidic residues but can have mixed charge\n")
        f.write(f"- The scoring system considers multiple factors including charge, ")
        f.write(f"hydrophobicity, and presence of characteristic residues\n")
        f.write(f"- A higher score suggests more SMAP-like characteristics, but ")
        f.write(f"experimental validation would be needed for confirmation\n")
    
    print(f"\nSMAP analysis complete!")
    print(f"Results saved to {output_dir}/")
    print(f"- smap_analysis.json: Detailed results")
    print(f"- smap_report.txt: Human-readable report")
    print(f"- smap_regions.png: Top candidate visualizations")
    print(f"- smap_score_distribution.png: Score distribution")

if __name__ == "__main__":
    main()