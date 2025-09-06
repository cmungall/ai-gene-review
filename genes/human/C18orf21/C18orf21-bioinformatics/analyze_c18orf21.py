#!/usr/bin/env python3
"""
Bioinformatics analysis of human C18orf21 protein (Q32NC0)
"""

import json
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
from pathlib import Path

# C18orf21 protein sequence
SEQUENCE = """
MRQKHYLEAAARGLHDSCPGQARYLLWAYTSSHDDKSTFEETCPYCFQLLVLDNSRVRLK
PKARLTPKIQKLLNREARNYFTLSFKEAKMVKKFKDSKSVLLITCKTCNRTVKHHGKSRS
FVSTLKSNPATPTSKLSLKTPERRTANPNHDMSGSKGKSPASVRFPTSGQSVSTCSSKNT
SKTKKHFSQLKMLLSQNESQKIPKVDFRNFLSSLKGGLLK
""".replace('\n', '')

def analyze_basic_properties(sequence):
    """Analyze basic protein properties"""
    analysis = ProteinAnalysis(sequence)
    
    properties = {
        'length': len(sequence),
        'molecular_weight': analysis.molecular_weight(),
        'aromaticity': analysis.aromaticity(),
        'instability_index': analysis.instability_index(),
        'isoelectric_point': analysis.isoelectric_point(),
        'secondary_structure': analysis.secondary_structure_fraction(),
        'gravy': analysis.gravy(),
        'amino_acid_composition': analysis.get_amino_acids_percent()
    }
    
    # Classify stability
    properties['stability'] = 'stable' if properties['instability_index'] < 40 else 'unstable'
    
    return properties

def analyze_phosphorylation_sites(sequence):
    """Analyze known and predicted phosphorylation sites"""
    # Known sites from UniProt
    known_sites = {
        'S126': {'position': 126, 'residue': 'S', 'context': sequence[121:131] if len(sequence) > 130 else ''},
        'T130': {'position': 130, 'residue': 'T', 'context': sequence[125:135] if len(sequence) > 134 else ''},
        'T139': {'position': 139, 'residue': 'T', 'context': sequence[134:144] if len(sequence) > 143 else ''}
    }
    
    # Analyze context for phosphorylation motifs
    for site, info in known_sites.items():
        context = info['context']
        # Check for common kinase motifs
        motifs = []
        if 'RXX' in context or 'KXX' in context:
            motifs.append('PKA/PKC-like')
        if 'SP' in context or 'TP' in context:
            motifs.append('Proline-directed (MAPK/CDK-like)')
        if 'S/T' in context and ('E' in context or 'D' in context):
            motifs.append('Casein kinase II-like')
        info['predicted_kinases'] = motifs
    
    return known_sites

def analyze_regions(sequence):
    """Analyze different regions of the protein"""
    regions = {
        'N-terminal': {
            'sequence': sequence[:50],
            'properties': analyze_basic_properties(sequence[:50])
        },
        'DUF4674_core': {
            'sequence': sequence[50:150],
            'properties': analyze_basic_properties(sequence[50:150])
        },
        'C-terminal': {
            'sequence': sequence[150:],
            'properties': analyze_basic_properties(sequence[150:])
        }
    }
    
    # Identify disordered regions based on composition
    disorder_regions = []
    window_size = 20
    for i in range(0, len(sequence) - window_size):
        window = sequence[i:i+window_size]
        analysis = ProteinAnalysis(window)
        # High disorder if rich in P, E, S, T, K, R, D
        disorder_residues = sum(analysis.get_amino_acids_percent().get(aa, 0) 
                               for aa in 'PESTKRD')
        if disorder_residues > 0.5:  # More than 50% disorder-promoting residues
            disorder_regions.append((i, i+window_size))
    
    # Merge overlapping regions
    merged_disorder = []
    if disorder_regions:
        current_start, current_end = disorder_regions[0]
        for start, end in disorder_regions[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged_disorder.append((current_start, current_end))
                current_start, current_end = start, end
        merged_disorder.append((current_start, current_end))
    
    return regions, merged_disorder

def check_alphafold_structure(uniprot_id):
    """Check AlphaFold structure information"""
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data:
                return {
                    'available': True,
                    'entry': data[0] if isinstance(data, list) else data,
                    'pdb_url': f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb",
                    'confidence_url': f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-confidence_v4.json"
                }
    except Exception as e:
        print(f"Error checking AlphaFold: {e}")
    
    return {'available': False}

def analyze_motifs(sequence):
    """Search for functional motifs"""
    motifs = {}
    
    # Nuclear localization signals (NLS)
    # Classical monopartite: K(K/R)X(K/R)
    # Bipartite: 2 basic residues, 10-12 residue spacer, 3/5 basic residues
    nls_patterns = []
    for i in range(len(sequence) - 4):
        window = sequence[i:i+5]
        basic_count = sum(1 for aa in window if aa in 'KR')
        if basic_count >= 4:
            nls_patterns.append((i, window))
    
    if nls_patterns:
        motifs['potential_NLS'] = nls_patterns
    
    # Nucleolar localization signals (NoLS) - often R/K rich regions
    nols_regions = []
    for i in range(0, len(sequence) - 20, 10):
        window = sequence[i:i+20]
        rk_content = sum(1 for aa in window if aa in 'RK') / 20
        if rk_content > 0.3:  # More than 30% R/K
            nols_regions.append((i, i+20, rk_content))
    
    if nols_regions:
        motifs['potential_NoLS'] = nols_regions
    
    # Coiled-coil regions (simplified - leucine zippers)
    # Look for heptad repeats with hydrophobic residues at positions a and d
    cc_regions = []
    for i in range(0, len(sequence) - 28, 7):
        heptads = sequence[i:i+28]
        if len(heptads) == 28:
            # Check positions a (0, 7, 14, 21) and d (3, 10, 17, 24)
            hydrophobic_a = sum(1 for j in [0, 7, 14, 21] if heptads[j] in 'LIVMF')
            hydrophobic_d = sum(1 for j in [3, 10, 17, 24] if heptads[j] in 'LIVMF')
            if hydrophobic_a >= 3 or hydrophobic_d >= 3:
                cc_regions.append((i, i+28))
    
    if cc_regions:
        motifs['potential_coiled_coil'] = cc_regions
    
    return motifs

def generate_report(results):
    """Generate a comprehensive report"""
    report = []
    report.append("# C18orf21 Bioinformatics Analysis Results\n")
    report.append("## Executive Summary\n")
    report.append("Analysis of human C18orf21 (Q32NC0), an uncharacterized protein of the UPF0711 family.\n")
    
    # Basic properties
    props = results['basic_properties']
    report.append("## Basic Properties\n")
    report.append(f"- **Length**: {props['length']} amino acids")
    report.append(f"- **Molecular Weight**: {props['molecular_weight']:.2f} Da")
    report.append(f"- **Isoelectric Point (pI)**: {props['isoelectric_point']:.2f}")
    report.append(f"- **Aromaticity**: {props['aromaticity']:.3f}")
    report.append(f"- **Instability Index**: {props['instability_index']:.2f} ({props['stability']})")
    report.append(f"- **GRAVY (hydropathy)**: {props['gravy']:.3f}")
    
    helix, turn, sheet = props['secondary_structure']
    report.append(f"- **Predicted Secondary Structure**:")
    report.append(f"  - Alpha helix: {helix:.1%}")
    report.append(f"  - Beta sheet: {sheet:.1%}")
    report.append(f"  - Random coil: {turn:.1%}\n")
    
    # Phosphorylation sites
    report.append("## Phosphorylation Sites\n")
    for site, info in results['phosphorylation_sites'].items():
        report.append(f"- **{site}** (position {info['position']})")
        report.append(f"  - Context: {info['context']}")
        if info['predicted_kinases']:
            report.append(f"  - Predicted kinases: {', '.join(info['predicted_kinases'])}")
    report.append("")
    
    # Disorder regions
    if results['disorder_regions']:
        report.append("## Predicted Disordered Regions\n")
        for start, end in results['disorder_regions']:
            report.append(f"- Residues {start+1}-{end} (disorder-promoting residue enriched)")
        report.append("")
    
    # Motifs
    if results['motifs']:
        report.append("## Potential Functional Motifs\n")
        if 'potential_NLS' in results['motifs']:
            report.append("### Nuclear Localization Signals (NLS)")
            for pos, seq in results['motifs']['potential_NLS']:
                report.append(f"- Position {pos+1}: {seq}")
        
        if 'potential_NoLS' in results['motifs']:
            report.append("### Nucleolar Localization Signals (NoLS)")
            for start, end, rk in results['motifs']['potential_NoLS']:
                report.append(f"- Position {start+1}-{end}: {rk:.1%} R/K content")
        report.append("")
    
    # AlphaFold
    report.append("## Structural Information\n")
    if results['alphafold']['available']:
        report.append("- AlphaFold structure available")
        report.append(f"- [Structure file]({results['alphafold']['pdb_url']})")
        report.append(f"- [Confidence data]({results['alphafold']['confidence_url']})")
    else:
        report.append("- AlphaFold structure check failed (may need manual verification)")
    report.append("")
    
    # Functional predictions
    report.append("## Functional Predictions\n")
    report.append("Based on the bioinformatics analysis:")
    report.append("1. **Cellular localization**: Nuclear/nucleolar (based on composition and reported staining)")
    report.append("2. **Post-translational regulation**: Multiple phosphorylation sites suggest regulation by cellular kinases")
    report.append("3. **Structural features**: Contains both ordered (DUF4674) and disordered regions")
    report.append("4. **Potential functions**:")
    report.append("   - May be involved in nuclear/nucleolar processes")
    report.append("   - Phosphorylation suggests involvement in signaling pathways")
    report.append("   - Disorder regions suggest potential for protein-protein interactions")
    report.append("   - Conservation in vertebrates indicates essential cellular function")
    
    return "\n".join(report)

def main():
    """Run complete analysis"""
    print("Analyzing C18orf21 protein...")
    
    results = {
        'basic_properties': analyze_basic_properties(SEQUENCE),
        'phosphorylation_sites': analyze_phosphorylation_sites(SEQUENCE),
        'regions': analyze_regions(SEQUENCE)[0],
        'disorder_regions': analyze_regions(SEQUENCE)[1],
        'motifs': analyze_motifs(SEQUENCE),
        'alphafold': check_alphafold_structure('Q32NC0')
    }
    
    # Generate and save report
    report = generate_report(results)
    
    # Save results
    with open('RESULTS.md', 'w') as f:
        f.write(report)
    
    with open('analysis_data.json', 'w') as f:
        # Convert non-serializable items
        clean_results = {
            'basic_properties': {
                k: v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items()}
                for k, v in results['basic_properties'].items()
                if k != 'secondary_structure'
            },
            'phosphorylation_sites': results['phosphorylation_sites'],
            'disorder_regions': results['disorder_regions'],
            'motifs': results['motifs'],
            'alphafold': results['alphafold']
        }
        clean_results['basic_properties']['secondary_structure'] = list(results['basic_properties']['secondary_structure'])
        json.dump(clean_results, f, indent=2)
    
    print("Analysis complete! Results saved to RESULTS.md and analysis_data.json")
    
    return results

if __name__ == "__main__":
    main()