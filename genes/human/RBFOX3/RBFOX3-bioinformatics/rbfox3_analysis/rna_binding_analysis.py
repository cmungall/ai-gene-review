#!/usr/bin/env python3
"""
RBFOX3 RNA-Binding Specificity Analysis

This script analyzes RBFOX3's RNA-binding properties including:
- Known target recognition sequences
- RNA-binding domain structure
- Predicted RNA-binding sites
- Analysis of known target genes
"""

import re
import json
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

# RBFOX3 sequence
RBFOX3_SEQUENCE = "MAQPYPPAQYPPPPQNGIPAYAPPPPHPTQDYSGQTPVPTEHGMTLYTPAQTHPEQPGSEASTQPIAGTQTVPQTDEAAQTDSQPLHPSDPTEKQQPKRLHVSNIPFRFRDPDLRQMFGQFGKILDVEIIFNERGSKGFGFVTFETSSDADRAREKLNGTIVEGRKIEVNNATARVMTNKKTGNPYTNGWKLNPVVGAVYGPEFYAVTGFPYPTTGTAVAYRGAHLRGRGRRAVYNTFRAPPPPPIPTYGGAVVYQDGFYGAEIYGGYAAYRYAQPAAAAAYSDSYGRVYAAADPYHHTTIGPAATYSIGTTM"

# Known RBFOX binding motifs from literature
RBFOX_MOTIFS = {
    'canonical': 'UGCAUG',
    'core': 'GCAUG',
    'extended': '(U)GCAUG',
    'degenerate': '[UC]GCA[UC]G'
}

# Known RBFOX3 target genes from literature
KNOWN_TARGETS = {
    'RBFOX2': {
        'description': 'RBFOX3 regulates alternative splicing of RBFOX2',
        'mechanism': 'Promotes exon skipping to enhance nonsense-mediated decay',
        'evidence': 'PMID:21747913'
    },
    'NUMB': {
        'description': 'RBFOX3 regulates Numb alternative splicing',
        'mechanism': 'Controls inclusion of exon 9 during neuronal differentiation',
        'evidence': 'PMID:23420872'
    },
    'GRIA2': {
        'description': 'Glutamate receptor subunit regulation',
        'mechanism': 'Alternative splicing control in neurons',
        'evidence': 'Literature review'
    },
    'CACNA1C': {
        'description': 'Calcium channel alternative splicing',
        'mechanism': 'Neuronal-specific isoform regulation',
        'evidence': 'Literature review'
    }
}

def analyze_rrm_structure():
    """Analyze the RNA Recognition Motif structure in RBFOX3"""
    print("=== RBFOX3 RRM Domain Analysis ===")
    
    # RRM domain boundaries from UniProt
    rrm_start, rrm_end = 99, 175  # 0-indexed
    rrm_sequence = RBFOX3_SEQUENCE[rrm_start:rrm_end]
    
    print(f"RRM domain (aa {rrm_start+1}-{rrm_end}):")
    print(f"  {rrm_sequence}")
    print(f"  Length: {len(rrm_sequence)} amino acids")
    
    # Identify key RRM motifs
    # RNP1 (RNA-binding nucleoprotein 1) - typically K/R-G-F/Y-G-F/Y-V-X-F/Y
    # RNP2 (RNA-binding nucleoprotein 2) - typically L/V-F/Y-V-G-N-L
    
    rnp1_pattern = r'[RK]G[FY]G[FY]'
    rnp2_pattern = r'[LV][FY][VL]G'
    
    rnp1_matches = list(re.finditer(rnp1_pattern, rrm_sequence))
    rnp2_matches = list(re.finditer(rnp2_pattern, rrm_sequence))
    
    print(f"\nRNP motif analysis within RRM:")
    if rnp1_matches:
        for match in rnp1_matches:
            pos = match.start() + rrm_start + 1  # Convert to 1-indexed absolute position
            print(f"  RNP1-like motif: {match.group()} at position {pos}")
    else:
        print("  No canonical RNP1 motif found")
    
    if rnp2_matches:
        for match in rnp2_matches:
            pos = match.start() + rrm_start + 1
            print(f"  RNP2-like motif: {match.group()} at position {pos}")
    else:
        print("  No canonical RNP2 motif found")
    
    # Identify aromatic and basic residues important for RNA binding
    aromatic_residues = 'FYW'
    basic_residues = 'RKH'
    
    aromatic_count = sum(1 for aa in rrm_sequence if aa in aromatic_residues)
    basic_count = sum(1 for aa in rrm_sequence if aa in basic_residues)
    
    print(f"\nRNA-binding relevant residues in RRM:")
    print(f"  Aromatic residues (F,Y,W): {aromatic_count}")
    print(f"  Basic residues (R,K,H): {basic_count}")
    
    return {
        'rrm_sequence': rrm_sequence,
        'rnp1_motifs': [(m.start() + rrm_start + 1, m.group()) for m in rnp1_matches],
        'rnp2_motifs': [(m.start() + rrm_start + 1, m.group()) for m in rnp2_matches],
        'aromatic_count': aromatic_count,
        'basic_count': basic_count
    }

def analyze_binding_specificity():
    """Analyze RBFOX3's RNA-binding specificity"""
    print("\n=== RNA-Binding Specificity Analysis ===")
    
    print("Known RBFOX family binding motifs:")
    for motif_name, sequence in RBFOX_MOTIFS.items():
        print(f"  {motif_name}: {sequence}")
    
    # Analyze the canonical motif UGCAUG
    canonical_motif = 'UGCAUG'
    print(f"\nCanonical RBFOX motif analysis ({canonical_motif}):")
    print("  - High affinity binding site")
    print("  - Purine-rich context enhances binding")
    print("  - Position relative to splice sites affects regulation")
    
    # Context-dependent regulation
    print("\nContext-dependent splicing regulation:")
    print("  - Upstream of exon: Typically promotes inclusion")
    print("  - Downstream of exon: Typically promotes skipping")
    print("  - Distance from splice site affects efficiency")
    
    return {
        'canonical_motif': canonical_motif,
        'binding_mechanism': 'sequence_specific',
        'regulatory_context': 'position_dependent'
    }

def analyze_target_genes():
    """Analyze known RBFOX3 target genes"""
    print("\n=== Known RBFOX3 Target Genes ===")
    
    for gene, info in KNOWN_TARGETS.items():
        print(f"\n{gene}:")
        print(f"  Description: {info['description']}")
        print(f"  Mechanism: {info['mechanism']}")
        print(f"  Evidence: {info['evidence']}")
    
    # Categorize targets by function
    target_categories = {
        'RNA_processing': ['RBFOX2'],
        'neuronal_differentiation': ['NUMB'],
        'synaptic_function': ['GRIA2', 'CACNA1C']
    }
    
    print("\nTarget gene categories:")
    for category, genes in target_categories.items():
        print(f"  {category.replace('_', ' ')}: {', '.join(genes)}")
    
    return {
        'known_targets': KNOWN_TARGETS,
        'target_categories': target_categories
    }

def predict_binding_sites():
    """Predict potential RNA-binding sites in RBFOX3"""
    print("\n=== RNA-Binding Site Prediction ===")
    
    # Based on UniProt experimental evidence
    experimental_sites = [100, 108, 109, 133, 138, 142, 166, 176]
    
    print("Experimentally validated RNA interaction sites:")
    for site in experimental_sites:
        residue = RBFOX3_SEQUENCE[site-1]  # Convert to 0-indexed
        print(f"  Position {site}: {residue}")
    
    # Analyze these sites within domain context
    rrm_sites = [site for site in experimental_sites if 100 <= site <= 175]
    non_rrm_sites = [site for site in experimental_sites if site < 100 or site > 175]
    
    print(f"\nRNA interaction sites in RRM domain: {rrm_sites}")
    print(f"RNA interaction sites outside RRM: {non_rrm_sites}")
    
    # Analyze amino acid properties at binding sites
    binding_residues = [RBFOX3_SEQUENCE[site-1] for site in experimental_sites]
    residue_types = {
        'aromatic': sum(1 for aa in binding_residues if aa in 'FYW'),
        'basic': sum(1 for aa in binding_residues if aa in 'RKH'),
        'polar': sum(1 for aa in binding_residues if aa in 'NQST'),
        'hydrophobic': sum(1 for aa in binding_residues if aa in 'AILV')
    }
    
    print(f"\nAmino acid types at binding sites:")
    for aa_type, count in residue_types.items():
        print(f"  {aa_type}: {count}/{len(binding_residues)}")
    
    return {
        'experimental_sites': experimental_sites,
        'rrm_sites': rrm_sites,
        'binding_residue_types': residue_types
    }

def analyze_splicing_regulation():
    """Analyze RBFOX3's splicing regulation mechanism"""
    print("\n=== Splicing Regulation Mechanism ===")
    
    print("RBFOX3 splicing regulation principles:")
    print("1. Position-dependent regulation:")
    print("   - Binding upstream of exon → inclusion")
    print("   - Binding downstream of exon → skipping")
    print("   - Distance from splice site affects strength")
    
    print("\n2. Cooperative binding:")
    print("   - Multiple RBFOX sites enhance regulation")
    print("   - Can form higher-order complexes")
    print("   - Cooperative with other splicing factors")
    
    print("\n3. Isoform-specific effects:")
    print("   - Nuclear isoforms: Direct splicing regulation")
    print("   - Cytoplasmic isoforms: May affect mRNA stability")
    
    print("\n4. Neuronal specificity:")
    print("   - Tissue-restricted expression")
    print("   - Neuron-specific alternative splicing events")
    print("   - Critical for neuronal differentiation")
    
    return {
        'regulation_type': 'position_dependent',
        'cooperative_binding': True,
        'tissue_specificity': 'neuronal',
        'functional_outcomes': ['exon_inclusion', 'exon_skipping', 'mRNA_stability']
    }

def create_binding_visualization():
    """Create visualization of RNA-binding properties"""
    print("\n=== Creating RNA-Binding Visualization ===")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Domain architecture with RNA-binding sites
    protein_length = len(RBFOX3_SEQUENCE)
    experimental_sites = [100, 108, 109, 133, 138, 142, 166, 176]
    
    # Draw protein backbone
    ax1.plot([1, protein_length], [1, 1], 'k-', linewidth=8, alpha=0.3)
    
    # Draw domains
    domains = [
        ("N-term disordered", 1, 104, "lightgray"),
        ("RRM domain", 100, 175, "red"),
        ("Fox-1 C-domain", 176, 312, "blue")
    ]
    
    for name, start, end, color in domains:
        ax1.plot([start, end], [1, 1], color=color, linewidth=12, alpha=0.8)
        mid_point = (start + end) / 2
        ax1.text(mid_point, 1.3, name, ha='center', va='bottom', fontsize=10)
    
    # Mark RNA interaction sites
    for site in experimental_sites:
        ax1.plot(site, 1, 'go', markersize=8, markeredgecolor='black', markeredgewidth=1)
        ax1.text(site, 0.7, str(site), ha='center', va='top', fontsize=8)
    
    ax1.set_xlim(0, 320)
    ax1.set_ylim(0.5, 1.8)
    ax1.set_xlabel('Amino Acid Position')
    ax1.set_title('RBFOX3 RNA-Binding Sites and Domain Architecture')
    ax1.set_yticks([])
    
    # Plot 2: RNA motif recognition
    motif_data = {
        'UGCAUG': 'Canonical high-affinity motif',
        'GCAUG': 'Core recognition sequence',
        'UCGAUG': 'Alternative variant',
        'UGCAUC': 'Lower affinity variant'
    }
    
    motifs = list(motif_data.keys())
    descriptions = list(motif_data.values())
    
    y_pos = range(len(motifs))
    colors = ['red', 'orange', 'yellow', 'lightblue']
    
    bars = ax2.barh(y_pos, [6]*len(motifs), color=colors, alpha=0.7)
    
    # Add motif sequences on bars
    for i, (motif, bar) in enumerate(zip(motifs, bars)):
        ax2.text(3, i, motif, ha='center', va='center', fontweight='bold', fontsize=12)
    
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(descriptions)
    ax2.set_xlabel('Relative Binding Affinity')
    ax2.set_title('RBFOX3 RNA Recognition Motifs')
    ax2.set_xlim(0, 6)
    
    plt.tight_layout()
    plt.savefig('rbfox3_rna_binding.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("RNA-binding visualization saved as 'rbfox3_rna_binding.png'")

def main():
    """Main analysis function"""
    print("RBFOX3 RNA-Binding Specificity Analysis")
    print("=" * 50)
    
    # RRM structure analysis
    rrm_analysis = analyze_rrm_structure()
    
    # Binding specificity
    specificity_analysis = analyze_binding_specificity()
    
    # Target gene analysis
    target_analysis = analyze_target_genes()
    
    # Binding site prediction
    binding_sites = predict_binding_sites()
    
    # Splicing regulation mechanism
    regulation_analysis = analyze_splicing_regulation()
    
    # Create visualization
    create_binding_visualization()
    
    # Compile results
    results = {
        'rrm_analysis': rrm_analysis,
        'binding_specificity': specificity_analysis,
        'known_targets': target_analysis,
        'binding_sites': binding_sites,
        'regulation_mechanism': regulation_analysis,
        'summary': {
            'primary_function': 'Alternative splicing regulation',
            'binding_motif': 'UGCAUG (canonical)',
            'regulation_type': 'Position-dependent',
            'tissue_specificity': 'Neuronal',
            'key_targets': list(KNOWN_TARGETS.keys())
        }
    }
    
    with open('rna_binding_analysis.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 50)
    print("RNA-binding analysis complete. Results saved to 'rna_binding_analysis.json'")

if __name__ == "__main__":
    main()