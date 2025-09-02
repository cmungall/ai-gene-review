#!/usr/bin/env python3
"""
Analyze AlphaFold structure for a protein.
"""

import click
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import json

def parse_pdb(pdb_file):
    """Parse PDB file and extract key information."""
    atoms = []
    plddt_scores = []
    residues = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # Parse ATOM line
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                bfactor = float(line[60:66])  # pLDDT score in AlphaFold
                
                if atom_name == 'CA':  # Only keep CA atoms for residue-level analysis
                    atoms.append({'res_num': res_num, 'x': x, 'y': y, 'z': z})
                    plddt_scores.append(bfactor)
                    residues.append(res_name)
    
    return atoms, plddt_scores, residues

def calculate_distances(atoms):
    """Calculate pairwise distances between CA atoms."""
    n = len(atoms)
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            dx = atoms[i]['x'] - atoms[j]['x']
            dy = atoms[i]['y'] - atoms[j]['y']
            dz = atoms[i]['z'] - atoms[j]['z']
            distances[i][j] = np.sqrt(dx**2 + dy**2 + dz**2)
    
    return distances

def analyze_confidence(plddt_scores):
    """Analyze pLDDT confidence scores."""
    plddt_array = np.array(plddt_scores)
    
    confidence_regions = {
        'very_high': int((plddt_array > 90).sum()),
        'confident': int(((plddt_array > 70) & (plddt_array <= 90)).sum()),
        'low': int(((plddt_array > 50) & (plddt_array <= 70)).sum()),
        'very_low': int((plddt_array <= 50).sum())
    }
    
    return {
        'mean_plddt': float(np.mean(plddt_array)),
        'min_plddt': float(np.min(plddt_array)),
        'max_plddt': float(np.max(plddt_array)),
        'confidence_regions': confidence_regions
    }

def identify_structural_features(atoms, plddt_scores, tm_regions=None):
    """Identify structural features from coordinates."""
    features = {}
    
    # Calculate end-to-end distance
    if len(atoms) > 0:
        n_term = atoms[0]
        c_term = atoms[-1]
        end_to_end = np.sqrt(
            (n_term['x'] - c_term['x'])**2 + 
            (n_term['y'] - c_term['y'])**2 + 
            (n_term['z'] - c_term['z'])**2
        )
        features['end_to_end_distance'] = float(end_to_end)
    
    # Analyze TM regions if provided
    if tm_regions:
        tm_confidence = []
        for start, end in tm_regions:
            if start <= len(plddt_scores) and end <= len(plddt_scores):
                tm_plddt = plddt_scores[start-1:end]
                tm_confidence.append({
                    'region': f'{start}-{end}',
                    'mean_plddt': float(np.mean(tm_plddt)),
                    'min_plddt': float(np.min(tm_plddt))
                })
        features['tm_regions_confidence'] = tm_confidence
    
    # Calculate radius of gyration
    if len(atoms) > 0:
        coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
        center = np.mean(coords, axis=0)
        rg = np.sqrt(np.mean(np.sum((coords - center)**2, axis=1)))
        features['radius_of_gyration'] = float(rg)
    
    return features

@click.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.option('--tm-regions', multiple=True, help='TM regions in format start-end')
@click.option('--output-prefix', default='alphafold', help='Prefix for output files')
def main(pdb_file, tm_regions, output_prefix):
    """Analyze AlphaFold structure from PDB file."""
    
    print(f"Analyzing AlphaFold structure: {pdb_file}")
    print("=" * 60)
    
    # Parse PDB file
    atoms, plddt_scores, residues = parse_pdb(pdb_file)
    
    print(f"\nStructure contains {len(atoms)} residues")
    
    # Parse TM regions
    tm_list = []
    if tm_regions:
        for tm in tm_regions:
            start, end = map(int, tm.split('-'))
            tm_list.append((start, end))
    
    # Analyze confidence
    confidence = analyze_confidence(plddt_scores)
    
    print(f"\nConfidence Analysis (pLDDT scores):")
    print("-" * 40)
    print(f"Mean pLDDT: {confidence['mean_plddt']:.1f}")
    print(f"Range: {confidence['min_plddt']:.1f} - {confidence['max_plddt']:.1f}")
    print(f"\nConfidence distribution:")
    for level, count in confidence['confidence_regions'].items():
        print(f"  {level}: {count} residues")
    
    # Analyze structural features
    features = identify_structural_features(atoms, plddt_scores, tm_list)
    
    print(f"\nStructural Features:")
    print("-" * 40)
    print(f"End-to-end distance: {features.get('end_to_end_distance', 0):.1f} Å")
    print(f"Radius of gyration: {features.get('radius_of_gyration', 0):.1f} Å")
    
    if 'tm_regions_confidence' in features:
        print(f"\nTM Regions Confidence:")
        for tm in features['tm_regions_confidence']:
            print(f"  {tm['region']}: mean pLDDT = {tm['mean_plddt']:.1f}")
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: pLDDT scores along sequence
    positions = list(range(1, len(plddt_scores) + 1))
    ax1.plot(positions, plddt_scores, 'b-', linewidth=1.5)
    ax1.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='Very high (>90)')
    ax1.axhline(y=70, color='yellow', linestyle='--', alpha=0.5, label='Confident (>70)')
    ax1.axhline(y=50, color='orange', linestyle='--', alpha=0.5, label='Low (>50)')
    
    # Mark TM regions
    for start, end in tm_list:
        ax1.axvspan(start, end, alpha=0.2, color='gray', label='TM region' if start == tm_list[0][0] else '')
    
    ax1.set_xlabel('Residue Position')
    ax1.set_ylabel('pLDDT Score')
    ax1.set_title('AlphaFold Confidence (pLDDT) Along Sequence')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Distance matrix
    distances = calculate_distances(atoms)
    im = ax2.imshow(distances, cmap='viridis', aspect='auto')
    ax2.set_xlabel('Residue Index')
    ax2.set_ylabel('Residue Index')
    ax2.set_title('CA-CA Distance Matrix')
    plt.colorbar(im, ax=ax2, label='Distance (Å)')
    
    plt.tight_layout()
    plot_file = f'{output_prefix}_analysis.png'
    plt.savefig(plot_file, dpi=150)
    print(f"\nVisualization saved as '{plot_file}'")
    
    # Save results
    results = {
        'pdb_file': pdb_file,
        'num_residues': len(atoms),
        'confidence_analysis': confidence,
        'structural_features': features,
        'tm_regions': tm_list
    }
    
    json_file = f'{output_prefix}_results.json'
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to '{json_file}'")

if __name__ == '__main__':
    main()