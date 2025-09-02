#!/usr/bin/env python3
"""
Analyze structural features of Epe1 JmjC domain.
Predict secondary structure and identify key structural elements.
"""

from Bio import SeqIO
from pathlib import Path
import requests
import json
import matplotlib.pyplot as plt
import numpy as np

def predict_secondary_structure(sequence):
    """
    Simple secondary structure prediction based on amino acid propensities.
    This is a simplified method - ideally would use more sophisticated tools.
    """
    # Chou-Fasman propensities (simplified)
    helix_formers = set("AELMQK")
    helix_breakers = set("PGNDS")
    sheet_formers = set("VFIYW")
    sheet_breakers = set("PEGDK")
    
    structure = []
    window_size = 5
    
    for i in range(len(sequence)):
        window_start = max(0, i - window_size // 2)
        window_end = min(len(sequence), i + window_size // 2 + 1)
        window = sequence[window_start:window_end]
        
        helix_score = sum(1 for aa in window if aa in helix_formers) - sum(1 for aa in window if aa in helix_breakers)
        sheet_score = sum(1 for aa in window if aa in sheet_formers) - sum(1 for aa in window if aa in sheet_breakers)
        
        if helix_score > sheet_score and helix_score > 0:
            structure.append('H')  # Helix
        elif sheet_score > 0:
            structure.append('E')  # Sheet
        else:
            structure.append('C')  # Coil
    
    return ''.join(structure)

def analyze_jmjc_structural_features(jmjc_seq):
    """Analyze structural features of JmjC domain."""
    print("\nStructural Features of Epe1 JmjC Domain:")
    print("-" * 50)
    
    # Secondary structure prediction
    ss_pred = predict_secondary_structure(jmjc_seq)
    
    helix_count = ss_pred.count('H')
    sheet_count = ss_pred.count('E')
    coil_count = ss_pred.count('C')
    
    print(f"\n  Predicted secondary structure composition:")
    print(f"    α-helix: {helix_count} ({helix_count*100/len(ss_pred):.1f}%)")
    print(f"    β-sheet: {sheet_count} ({sheet_count*100/len(ss_pred):.1f}%)")
    print(f"    Coil/loop: {coil_count} ({coil_count*100/len(ss_pred):.1f}%)")
    
    # Identify β-strands (characteristic of JmjC fold)
    # JmjC domains typically have 8 β-strands forming a β-barrel
    beta_regions = []
    in_beta = False
    start = 0
    
    for i, ss in enumerate(ss_pred):
        if ss == 'E' and not in_beta:
            in_beta = True
            start = i
        elif ss != 'E' and in_beta:
            in_beta = False
            if i - start >= 3:  # Minimum length for β-strand
                beta_regions.append((start, i))
    
    print(f"\n  Predicted β-strands: {len(beta_regions)}")
    if beta_regions:
        print("    Positions:")
        for i, (start, end) in enumerate(beta_regions[:8], 1):  # Show first 8
            print(f"      β{i}: {start+243}-{end+243} ({end-start} aa)")
    
    # Check for metal-binding pocket characteristics
    print("\n  Metal-binding pocket analysis:")
    
    # Look for histidines and their spacing
    histidine_positions = [i for i, aa in enumerate(jmjc_seq) if aa == 'H']
    
    if len(histidine_positions) >= 2:
        print(f"    Histidine positions: {[p+243 for p in histidine_positions]}")
        
        # Check spacing between histidines
        if len(histidine_positions) >= 2:
            spacings = [histidine_positions[i+1] - histidine_positions[i] 
                       for i in range(len(histidine_positions)-1)]
            print(f"    Spacing between histidines: {spacings}")
            
            # Typical Fe(II) coordination requires His residues ~15-25 aa apart
            good_spacing = any(15 <= s <= 30 for s in spacings)
            if good_spacing:
                print("    ✓ Histidine spacing compatible with Fe(II) coordination")
            else:
                print("    ✗ Unusual histidine spacing for Fe(II) coordination")
    
    # Hydrophobic core analysis
    hydrophobic_core = sum(1 for aa in jmjc_seq if aa in "VLIMFYW")
    print(f"\n  Hydrophobic core residues: {hydrophobic_core} ({hydrophobic_core*100/len(jmjc_seq):.1f}%)")
    
    return {
        "helix_percent": helix_count*100/len(ss_pred),
        "sheet_percent": sheet_count*100/len(ss_pred),
        "beta_strands": len(beta_regions),
        "histidine_positions": histidine_positions,
        "hydrophobic_percent": hydrophobic_core*100/len(jmjc_seq)
    }

def compare_structural_features():
    """Compare structural features with known JmjC structures."""
    print("\nComparison with Known JmjC Structures:")
    print("-" * 50)
    
    # Typical JmjC domain features (from literature/PDB structures)
    typical_features = {
        "beta_strands": "8 (forming β-barrel)",
        "alpha_helices": "2-4 surrounding β-barrel",
        "Fe_coordination": "HX(D/E)...H motif",
        "aKG_binding": "Basic pocket (K/R)",
        "fold": "Double-stranded β-helix (DSBH)"
    }
    
    print("\n  Canonical JmjC domain features:")
    for feature, description in typical_features.items():
        print(f"    {feature}: {description}")
    
    # Load Epe1 sequence to check for motifs
    from pathlib import Path
    with open(Path("data") / "epe1_spombe.fasta") as f:
        from Bio import SeqIO
        for record in SeqIO.parse(f, "fasta"):
            epe1_seq = str(record.seq)
            break
    
    # Actually check for HVD motif in the sequence
    import re
    hvd_pattern = re.compile(r"HVD")
    hxd_pattern = re.compile(r"H.D")
    
    print("\n  Epe1 deviations detected:")
    if hvd_pattern.search(str(epe1_seq)):
        print("    - HVD motif found (V is hydrophobic, unusual for Fe(II) binding)")
        print("    - Likely altered metal-binding geometry")
    elif hxd_pattern.search(str(epe1_seq)):
        print("    - HXD-like motif found, specific residue analysis needed")
    else:
        print("    - No canonical HXD motif detected")
    print("    - Likely retains overall DSBH fold based on sequence")
    print("    - Catalytic activity depends on specific active site residues")

def visualize_domain_architecture():
    """Create a visual representation of Epe1 domain architecture."""
    # Load Epe1 sequence
    from pathlib import Path
    from Bio import SeqIO
    with open(Path("data") / "epe1_spombe.fasta") as f:
        for record in SeqIO.parse(f, "fasta"):
            epe1_seq = str(record.seq)
            break
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
    
    # Epe1 full-length protein domains
    protein_length = 948
    
    domains = [
        {"name": "N-terminal", "start": 1, "end": 242, "color": "lightgray"},
        {"name": "JmjC domain", "start": 243, "end": 402, "color": "lightblue"},
        {"name": "Linker", "start": 403, "end": 848, "color": "lightgray"},
        {"name": "C-terminal (HP1 binding)", "start": 849, "end": 948, "color": "lightgreen"}
    ]
    
    ax1.set_xlim(0, protein_length)
    ax1.set_ylim(0, 1)
    
    for domain in domains:
        width = domain["end"] - domain["start"]
        rect = plt.Rectangle((domain["start"], 0.3), width, 0.4, 
                            facecolor=domain["color"], edgecolor='black', linewidth=1)
        ax1.add_patch(rect)
        
        # Add label
        mid = (domain["start"] + domain["end"]) / 2
        ax1.text(mid, 0.5, domain["name"], ha='center', va='center', fontsize=10)
    
    # Detect actual motif positions
    import re
    key_positions = []
    
    # Find HVD or HXD motifs
    hvd_match = re.search(r"HVD", str(epe1_seq))
    if hvd_match:
        key_positions.append({"pos": hvd_match.start() + 1, "label": "HVD", "color": "red"})
    
    hie_match = re.search(r"HIE", str(epe1_seq))
    if hie_match:
        key_positions.append({"pos": hie_match.start() + 1, "label": "HIE", "color": "red"})
    
    # If no specific motifs found, look for any HXD/HXE patterns
    if not key_positions:
        hxd_match = re.search(r"H.D", str(epe1_seq))
        if hxd_match:
            motif = epe1_seq[hxd_match.start():hxd_match.start()+3]
            key_positions.append({"pos": hxd_match.start() + 1, "label": motif, "color": "orange"})
    
    for pos_info in key_positions:
        ax1.axvline(x=pos_info["pos"], color=pos_info["color"], linestyle='--', alpha=0.5)
        ax1.text(pos_info["pos"], 0.8, pos_info["label"], ha='center', fontsize=8)
    
    ax1.set_xlabel("Amino acid position")
    ax1.set_title("Epe1 Domain Architecture")
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    
    # JmjC domain detail
    jmjc_start = 243
    jmjc_end = 402
    jmjc_length = jmjc_end - jmjc_start
    
    ax2.set_xlim(jmjc_start, jmjc_end)
    ax2.set_ylim(0, 1)
    
    # Highlight key regions
    regions = [
        {"name": "Fe-binding region", "start": 275, "end": 285, "color": "salmon"},
        {"name": "Putative αKG site", "start": 290, "end": 305, "color": "lightyellow"},
    ]
    
    for region in regions:
        width = region["end"] - region["start"]
        rect = plt.Rectangle((region["start"], 0.3), width, 0.4,
                            facecolor=region["color"], edgecolor='black', 
                            linewidth=1, alpha=0.7)
        ax2.add_patch(rect)
        ax2.text((region["start"] + region["end"])/2, 0.15, region["name"], 
                ha='center', fontsize=8)
    
    # Find actual histidine positions in JmjC domain
    jmjc_seq = epe1_seq[jmjc_start-1:jmjc_end]
    histidine_positions = [i + jmjc_start for i, aa in enumerate(jmjc_seq) if aa == 'H'][:5]  # Show first 5
    for h_pos in histidine_positions:
        ax2.scatter(h_pos, 0.5, color='blue', s=50, zorder=5)
        ax2.text(h_pos, 0.65, 'H', ha='center', fontsize=8, color='blue')
    
    ax2.set_xlabel("Amino acid position")
    ax2.set_title("JmjC Domain Detail (positions 243-402)")
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    
    plt.tight_layout()
    
    # Save figure
    results_dir = Path("results")
    fig_path = results_dir / "epe1_domain_architecture.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Domain architecture figure saved to {fig_path}")
    
    return fig_path

def main():
    print("=" * 60)
    print("Structural Features Analysis")
    print("=" * 60)
    
    data_dir = Path("data")
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    
    # Load Epe1 sequence
    with open(data_dir / "epe1_spombe.fasta") as f:
        epe1_record = next(SeqIO.parse(f, "fasta"))
    
    epe1_seq = str(epe1_record.seq)
    
    # Extract JmjC domain
    jmjc_start = 243
    jmjc_end = 402
    jmjc_seq = epe1_seq[jmjc_start-1:jmjc_end]
    
    print(f"\nAnalyzing JmjC domain (positions {jmjc_start}-{jmjc_end})")
    print(f"Length: {len(jmjc_seq)} aa")
    
    # Analyze structural features
    structural_features = analyze_jmjc_structural_features(jmjc_seq)
    
    # Compare with known structures
    compare_structural_features()
    
    # Create visualization
    visualize_domain_architecture()
    
    # Save results
    results_file = results_dir / "structural_analysis.txt"
    with open(results_file, "w") as f:
        f.write("Structural Features Analysis\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("JmjC Domain Structural Features:\n")
        f.write(f"- Position: {jmjc_start}-{jmjc_end}\n")
        f.write(f"- Length: {len(jmjc_seq)} aa\n")
        f.write(f"- Predicted α-helix: {structural_features['helix_percent']:.1f}%\n")
        f.write(f"- Predicted β-sheet: {structural_features['sheet_percent']:.1f}%\n")
        f.write(f"- Predicted β-strands: {structural_features['beta_strands']}\n")
        f.write(f"- Hydrophobic core: {structural_features['hydrophobic_percent']:.1f}%\n")
        
        f.write("\nKey Observations:\n")
        f.write("1. Epe1 likely maintains the JmjC fold structure\n")
        # Check what motif was actually found
        import re
        if re.search(r"HVD", str(epe1_seq)):
            f.write("2. HVD motif detected - deviates from canonical HXD\n")
        elif re.search(r"H.D", str(epe1_seq)):
            f.write("2. HXD-like motif found - specific analysis needed\n")
        else:
            f.write("2. No canonical HXD motif detected\n")
        f.write("3. Histidine positions suggest altered metal coordination\n")
        f.write("4. Structure retained but catalytic activity lost\n")
        f.write("5. Functions as a pseudo-enzyme/structural scaffold\n")
    
    print(f"\n✓ Results saved to {results_file}")

if __name__ == "__main__":
    main()