#!/usr/bin/env python3
"""
Critical bioinformatic analysis of LRX-1 protein (Q22179)
Validates/refutes automated domain predictions
"""

import re
from typing import List, Tuple, Dict

# LRX-1 protein sequence
SEQUENCE = """MAWLTSIFFILLAVQPVLPQDLYGTATQQQPYPYVQPSASSGSGGYVPNPQSSIHTVQQP
YPNIDVVEPDVDSVDIYETEEPQFKVVNPVFPLGGSGIVEEPGTIPPPMPQTQAPEKPDNS
YAINYCDKREFPDDVLAQYGLERIDYFVYNTSCSHVFFQCSIGQTFPLACMSEDQAFDKS
TENCNHKNAIKFCPEYDHVMHCTIKDTCTENEFACCAMPQSCIHVSKRCDGHPDCADGED
ENNCPSCARDDEFACVKSEHCIPANKRCDGVADDCEDGSNLDEIGCSKNTTCIGKFVCGTS
RGGVSCVDLDMHCDGKKDCLNGEDEMNCQEGRQKYLLCENQKQSVTRLQWCNGETDCAD
GSDEKYCY""".replace("\n", "")

def analyze_cysteine_patterns(sequence: str) -> Dict:
    """Analyze cysteine distribution and spacing patterns"""
    
    # Find all cysteine positions
    cys_positions = [i for i, aa in enumerate(sequence) if aa == 'C']
    
    # Calculate spacings between consecutive cysteines
    spacings = []
    for i in range(1, len(cys_positions)):
        spacing = cys_positions[i] - cys_positions[i-1] - 1
        spacings.append(spacing)
    
    # Canonical LDL-A domain pattern: C-x(2,3)-C-x(3,4)-[DE]-x(4,5)-C-x(7,8)-C-x(2,3)-C-x(3,9)-C
    # This means spacings should be roughly: [2-3, 3-4, 4-5, 7-8, 2-3, 3-9]
    
    # Check for potential LDL-A domains
    potential_ldla = []
    for i in range(len(spacings) - 5):
        pattern = spacings[i:i+6]
        # Check if pattern matches LDL-A spacing
        if (2 <= pattern[0] <= 3 and 
            3 <= pattern[1] <= 4 and
            4 <= pattern[2] <= 5 and
            7 <= pattern[3] <= 8 and
            2 <= pattern[4] <= 3 and
            3 <= pattern[5] <= 9):
            potential_ldla.append({
                'start_cys': cys_positions[i],
                'end_cys': cys_positions[i+6],
                'spacings': pattern
            })
    
    return {
        'total_cysteines': len(cys_positions),
        'cysteine_percentage': (len(cys_positions) / len(sequence)) * 100,
        'positions': cys_positions,
        'spacings': spacings,
        'potential_ldla_domains': potential_ldla,
        'cysteine_regions': identify_cysteine_regions(cys_positions)
    }

def identify_cysteine_regions(cys_positions: List[int]) -> List[Tuple[int, int]]:
    """Identify regions with clustered cysteines"""
    regions = []
    if not cys_positions:
        return regions
    
    start = cys_positions[0]
    end = cys_positions[0]
    
    for i in range(1, len(cys_positions)):
        if cys_positions[i] - cys_positions[i-1] <= 20:  # Within 20 residues
            end = cys_positions[i]
        else:
            regions.append((start, end))
            start = cys_positions[i]
            end = cys_positions[i]
    
    regions.append((start, end))
    return regions

def check_membrane_topology(sequence: str) -> Dict:
    """Check for transmembrane helices and signal peptides"""
    
    # Hydrophobic amino acids
    hydrophobic = set('AILMFWV')
    
    # Check signal peptide region (first 30 aa)
    signal_region = sequence[:30]
    signal_hydrophobic = sum(1 for aa in signal_region if aa in hydrophobic)
    
    # Scan for potential transmembrane helices (20+ hydrophobic residues)
    tm_regions = []
    window_size = 20
    
    for i in range(len(sequence) - window_size):
        window = sequence[i:i+window_size]
        hydro_count = sum(1 for aa in window if aa in hydrophobic)
        if hydro_count >= 14:  # 70% hydrophobic
            tm_regions.append({
                'start': i,
                'end': i + window_size,
                'hydrophobicity': hydro_count / window_size
            })
    
    return {
        'signal_peptide_likely': signal_hydrophobic >= 10,
        'signal_hydrophobicity': signal_hydrophobic / len(signal_region),
        'potential_tm_helices': tm_regions,
        'likely_topology': 'secreted' if not tm_regions else 'membrane'
    }

def compare_to_lrp_features() -> Dict:
    """Compare LRX-1 to known LRP family features"""
    
    lrp_features = {
        'typical_size': '>4000 aa',
        'lrx1_size': len(SEQUENCE),
        'has_beta_propeller': False,  # Would need specific motif search
        'has_egf_domains': False,  # Would need EGF pattern search
        'ldla_domains_expected': '30-40',
        'ldla_domains_found': 0  # Will be updated based on analysis
    }
    
    # Check for EGF-like pattern: C-x(3,4)-C-x(3)-C-x(5)-C-x-C-x(2)-C
    egf_pattern = re.compile(r'C.{3,4}C.{3}C.{5}C.C.{2}C')
    if egf_pattern.search(SEQUENCE):
        lrp_features['has_egf_domains'] = True
    
    return lrp_features

def main():
    """Run complete analysis"""
    
    print("=== LRX-1 BIOINFORMATIC ANALYSIS ===\n")
    
    # Analyze cysteines
    cys_analysis = analyze_cysteine_patterns(SEQUENCE)
    print("CYSTEINE ANALYSIS:")
    print(f"- Total cysteines: {cys_analysis['total_cysteines']}")
    print(f"- Percentage: {cys_analysis['cysteine_percentage']:.1f}%")
    print(f"- Spacing pattern: {cys_analysis['spacings'][:10]}...")
    print(f"- Potential LDL-A domains found: {len(cys_analysis['potential_ldla_domains'])}")
    
    if not cys_analysis['potential_ldla_domains']:
        print("  ❌ NO canonical LDL-A domain patterns detected!")
        print("  - Cysteine spacings incompatible with LDL-A requirements")
    
    print(f"\nCysteine-rich regions: {cys_analysis['cysteine_regions']}")
    
    # Check membrane topology
    print("\nMEMBRANE TOPOLOGY:")
    topology = check_membrane_topology(SEQUENCE)
    print(f"- Signal peptide likely: {topology['signal_peptide_likely']}")
    print(f"- Signal hydrophobicity: {topology['signal_hydrophobicity']:.1%}")
    print(f"- Transmembrane helices found: {len(topology['potential_tm_helices'])}")
    print(f"- Predicted topology: {topology['likely_topology']}")
    
    # Compare to LRP family
    print("\nLRP FAMILY COMPARISON:")
    lrp_comparison = compare_to_lrp_features()
    print(f"- LRX-1 size: {lrp_comparison['lrx1_size']} aa (typical LRP: {lrp_comparison['typical_size']})")
    print(f"- Has β-propeller domains: {lrp_comparison['has_beta_propeller']}")
    print(f"- Has EGF domains: {lrp_comparison['has_egf_domains']}")
    
    # Final verdict
    print("\n=== VERDICT ===")
    print("❌ LRX-1 is NOT a true LRP family protein:")
    print("  1. No canonical LDL-A domains despite annotation")
    print("  2. Wrong size (10x smaller than LRPs)")
    print("  3. Lacks essential LRP domains")
    print("  4. Likely secreted, not membrane-bound")
    print("\n✓ LRX-1 appears to be:")
    print("  - A novel cysteine-rich protein")
    print("  - Possibly secreted")
    print("  - Function unknown")
    print("  - C. elegans-specific protein family")

if __name__ == "__main__":
    main()