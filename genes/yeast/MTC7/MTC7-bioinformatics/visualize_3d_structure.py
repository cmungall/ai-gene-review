#!/usr/bin/env python3
"""
Create 3D visualization of AlphaFold structure with confidence coloring.
"""

import click
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def parse_pdb_for_3d(pdb_file):
    """Parse PDB file for 3D visualization."""
    ca_atoms = []
    all_atoms = []
    plddt_scores = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                res_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                bfactor = float(line[60:66])  # pLDDT score
                
                all_atoms.append({'name': atom_name, 'res': res_num, 'coords': [x, y, z], 'plddt': bfactor})
                
                if atom_name == 'CA':
                    ca_atoms.append({'res': res_num, 'coords': [x, y, z], 'plddt': bfactor})
                    plddt_scores.append(bfactor)
    
    return ca_atoms, all_atoms, plddt_scores

def calculate_secondary_structure_regions(ca_atoms, window=4):
    """Estimate secondary structure regions based on CA distances."""
    regions = []
    
    for i in range(len(ca_atoms) - window):
        # Calculate distances between consecutive CAs
        distances = []
        for j in range(window - 1):
            ca1 = np.array(ca_atoms[i + j]['coords'])
            ca2 = np.array(ca_atoms[i + j + 1]['coords'])
            distances.append(np.linalg.norm(ca2 - ca1))
        
        avg_dist = np.mean(distances)
        
        # Rough estimates: alpha helix ~3.8Å, beta strand ~3.3Å
        if 3.6 <= avg_dist <= 4.0:
            regions.append({'type': 'helix', 'pos': i + window//2})
        elif 3.0 <= avg_dist <= 3.5:
            regions.append({'type': 'strand', 'pos': i + window//2})
    
    return regions

@click.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.option('--tm-regions', multiple=True, help='TM regions in format start-end')
@click.option('--polybasic', multiple=True, help='Polybasic regions in format start-end')
@click.option('--output-prefix', default='structure_3d', help='Prefix for output files')
def main(pdb_file, tm_regions, polybasic, output_prefix):
    """Create 3D visualization of protein structure."""
    
    print(f"Creating 3D visualization of: {pdb_file}")
    print("=" * 60)
    
    # Parse PDB
    ca_atoms, all_atoms, plddt_scores = parse_pdb_for_3d(pdb_file)
    
    # Parse regions
    tm_list = []
    if tm_regions:
        for tm in tm_regions:
            start, end = map(int, tm.split('-'))
            tm_list.append((start, end))
    
    polybasic_list = []
    if polybasic:
        for pb in polybasic:
            start, end = map(int, pb.split('-'))
            polybasic_list.append((start, end))
    
    # Extract coordinates and colors
    coords = np.array([atom['coords'] for atom in ca_atoms])
    colors = [atom['plddt'] for atom in ca_atoms]
    
    # Create figure with multiple views
    fig = plt.figure(figsize=(16, 12))
    
    # 3D structure with confidence coloring
    ax1 = fig.add_subplot(221, projection='3d')
    
    # Plot backbone
    ax1.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'gray', alpha=0.3, linewidth=0.5)
    
    # Plot CA atoms colored by confidence
    scatter = ax1.scatter(coords[:, 0], coords[:, 1], coords[:, 2], 
                         c=colors, cmap='RdYlBu_r', s=50, 
                         vmin=0, vmax=100, alpha=0.8)
    
    # Highlight TM regions
    for start, end in tm_list:
        tm_coords = coords[start-1:end]
        ax1.plot(tm_coords[:, 0], tm_coords[:, 1], tm_coords[:, 2], 
                'red', linewidth=3, alpha=0.7, label=f'TM {start}-{end}')
    
    # Highlight polybasic regions
    for start, end in polybasic_list:
        pb_coords = coords[start-1:end]
        ax1.plot(pb_coords[:, 0], pb_coords[:, 1], pb_coords[:, 2], 
                'blue', linewidth=3, alpha=0.7, label=f'Polybasic {start}-{end}')
    
    ax1.set_xlabel('X (Å)')
    ax1.set_ylabel('Y (Å)')
    ax1.set_zlabel('Z (Å)')
    ax1.set_title('3D Structure with pLDDT Confidence')
    ax1.legend(loc='upper right', fontsize=8)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax1, pad=0.1)
    cbar.set_label('pLDDT Score')
    
    # Top view (XY plane)
    ax2 = fig.add_subplot(222)
    scatter2 = ax2.scatter(coords[:, 0], coords[:, 1], c=colors, 
                          cmap='RdYlBu_r', s=30, vmin=0, vmax=100)
    ax2.plot(coords[:, 0], coords[:, 1], 'gray', alpha=0.3, linewidth=0.5)
    
    # Mark termini
    ax2.scatter(coords[0, 0], coords[0, 1], s=100, c='green', marker='^', label='N-term')
    ax2.scatter(coords[-1, 0], coords[-1, 1], s=100, c='red', marker='v', label='C-term')
    
    ax2.set_xlabel('X (Å)')
    ax2.set_ylabel('Y (Å)')
    ax2.set_title('Top View (XY plane)')
    ax2.legend()
    ax2.axis('equal')
    
    # Side view (XZ plane)
    ax3 = fig.add_subplot(223)
    ax3.scatter(coords[:, 0], coords[:, 2], c=colors, 
               cmap='RdYlBu_r', s=30, vmin=0, vmax=100)
    ax3.plot(coords[:, 0], coords[:, 2], 'gray', alpha=0.3, linewidth=0.5)
    
    # Highlight TM regions in side view
    for start, end in tm_list:
        tm_coords = coords[start-1:end]
        ax3.fill_between(tm_coords[:, 0], tm_coords[:, 2]-2, tm_coords[:, 2]+2, 
                         alpha=0.3, color='red')
    
    ax3.set_xlabel('X (Å)')
    ax3.set_ylabel('Z (Å)')
    ax3.set_title('Side View (XZ plane) - TM regions highlighted')
    ax3.axis('equal')
    
    # pLDDT distribution
    ax4 = fig.add_subplot(224)
    ax4.hist(plddt_scores, bins=20, edgecolor='black', alpha=0.7)
    ax4.axvline(x=90, color='green', linestyle='--', label='Very high (>90)')
    ax4.axvline(x=70, color='yellow', linestyle='--', label='Confident (>70)')
    ax4.axvline(x=50, color='orange', linestyle='--', label='Low (>50)')
    ax4.set_xlabel('pLDDT Score')
    ax4.set_ylabel('Number of Residues')
    ax4.set_title('pLDDT Score Distribution')
    ax4.legend()
    
    plt.suptitle(f'MTC7 AlphaFold Structure Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save figure
    plot_file = f'{output_prefix}.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"3D visualization saved as '{plot_file}'")
    
    # Calculate and print structural metrics
    print("\nStructural Metrics:")
    print("-" * 40)
    
    # End-to-end distance
    end_to_end = np.linalg.norm(coords[-1] - coords[0])
    print(f"End-to-end distance: {end_to_end:.1f} Å")
    
    # Radius of gyration
    center = np.mean(coords, axis=0)
    rg = np.sqrt(np.mean(np.sum((coords - center)**2, axis=1)))
    print(f"Radius of gyration: {rg:.1f} Å")
    
    # Compactness
    max_dist = np.max([np.linalg.norm(coords[i] - coords[j]) 
                       for i in range(len(coords)) 
                       for j in range(i+1, min(i+10, len(coords)))])
    print(f"Maximum dimension: {max_dist:.1f} Å")
    
    # TM region analysis
    if tm_list:
        print("\nTM Region Structural Analysis:")
        for start, end in tm_list:
            tm_coords = coords[start-1:end]
            tm_plddt = plddt_scores[start-1:end]
            
            # Calculate helical parameters
            if len(tm_coords) > 3:
                # Rise per residue (should be ~1.5Å for alpha helix)
                rises = [np.linalg.norm(tm_coords[i+1] - tm_coords[i]) 
                        for i in range(len(tm_coords)-1)]
                avg_rise = np.mean(rises)
                
                # Overall length
                tm_length = np.linalg.norm(tm_coords[-1] - tm_coords[0])
                
                print(f"  TM {start}-{end}:")
                print(f"    Length: {tm_length:.1f} Å")
                print(f"    Avg rise/residue: {avg_rise:.2f} Å")
                print(f"    Mean pLDDT: {np.mean(tm_plddt):.1f}")
    
    # plt.show()  # Comment out for non-interactive use

if __name__ == '__main__':
    main()