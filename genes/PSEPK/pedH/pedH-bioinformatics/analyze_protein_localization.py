#!/usr/bin/env python3
"""
Generic protein localization analyzer
Analyzes proteins for signal peptides, transmembrane regions, and cellular localization
"""

import requests
from Bio import SeqIO
from io import StringIO
import json
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

def fetch_uniprot_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    # Try new REST API first
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code != 200:
        # Fallback to old API
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        response = requests.get(url)
    
    if response.status_code == 200 and response.text.strip():
        try:
            record = SeqIO.read(StringIO(response.text), "fasta")
            return str(record.seq), record.description
        except:
            return None, None
    return None, None

def read_fasta_file(filepath):
    """Read sequence from FASTA file"""
    with open(filepath, 'r') as f:
        record = SeqIO.read(f, "fasta")
        return str(record.seq), record.description

def analyze_signal_peptide(sequence, verbose=True):
    """Analyze potential signal peptide in N-terminus"""
    results = {
        'has_signal_peptide': False,
        'signal_features': {},
        'cleavage_sites': []
    }
    
    if verbose:
        print("\n=== Signal Peptide Analysis ===")
        print(f"First 60 amino acids:")
        print(sequence[:60])
    
    # Manual analysis based on signal peptide characteristics
    # Typical Sec signal peptide: N-region (positive), H-region (hydrophobic), C-region (cleavage site)
    first_30 = sequence[:30] if len(sequence) >= 30 else sequence
    
    # Check for positive residues in N-terminal
    positive_count = sum(1 for aa in first_30[:10] if aa in 'RK')
    hydrophobic_count = sum(1 for aa in first_30[5:20] if aa in 'AVILMFYW' and 5 < 20 <= len(first_30))
    
    results['signal_features'] = {
        'n_terminal_positive': positive_count,
        'hydrophobic_core': hydrophobic_count,
        'likely_signal': positive_count >= 1 and hydrophobic_count >= 5
    }
    
    if verbose:
        print(f"\nSignal peptide features:")
        print(f"  Positive residues in first 10 aa: {positive_count}")
        print(f"  Hydrophobic residues in positions 5-20: {hydrophobic_count}")
    
    # Look for Ala-X-Ala motif (common cleavage site)
    for i in range(15, min(30, len(sequence)-2)):
        if sequence[i] == 'A' and sequence[i+2] == 'A':
            results['cleavage_sites'].append({
                'position': i,
                'motif': sequence[max(0, i-2):min(len(sequence), i+5)],
                'type': 'Ala-X-Ala'
            })
            if verbose:
                print(f"  Potential cleavage site: {sequence[max(0, i-2):min(len(sequence), i+5)]} at position {i}")
    
    # Also check for other common cleavage motifs
    for i in range(15, min(35, len(sequence)-1)):
        # V-X-A motif
        if i > 0 and sequence[i-1] == 'V' and sequence[i+1] == 'A':
            results['cleavage_sites'].append({
                'position': i,
                'motif': sequence[max(0, i-2):min(len(sequence), i+4)],
                'type': 'V-X-A'
            })
    
    results['has_signal_peptide'] = results['signal_features']['likely_signal']
    
    return results

def check_transmembrane_regions(sequence, verbose=True):
    """Check for transmembrane helices using hydrophobicity analysis"""
    results = {
        'tm_regions': [],
        'is_membrane_protein': False,
        'topology': 'unknown'
    }
    
    if verbose:
        print("\n=== Transmembrane Region Analysis ===")
    
    # Kyte-Doolittle hydrophobicity scale
    hydrophobicity = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    # Calculate hydrophobicity for windows of 19 residues (typical TM helix)
    window_size = 19
    threshold = 1.6  # Threshold for TM helix
    
    tm_regions = []
    for i in range(len(sequence) - window_size):
        window = sequence[i:i+window_size]
        avg_hydrophobicity = sum(hydrophobicity.get(aa, 0) for aa in window) / window_size
        
        if avg_hydrophobicity > threshold:
            tm_regions.append({
                'start': i,
                'end': i + window_size,
                'score': round(avg_hydrophobicity, 2),
                'sequence': window
            })
    
    # Merge overlapping regions
    merged_tm = []
    for region in tm_regions:
        if merged_tm and region['start'] <= merged_tm[-1]['end'] + 5:
            # Merge with previous region
            merged_tm[-1]['end'] = max(merged_tm[-1]['end'], region['end'])
            merged_tm[-1]['score'] = max(merged_tm[-1]['score'], region['score'])
        else:
            merged_tm.append(region)
    
    results['tm_regions'] = merged_tm
    results['is_membrane_protein'] = len(merged_tm) > 0
    
    if verbose:
        if merged_tm:
            print(f"Potential transmembrane regions found: {len(merged_tm)}")
            for i, region in enumerate(merged_tm[:3], 1):  # Show first 3
                print(f"  TM{i}: Position {region['start']}-{region['end']} (score: {region['score']})")
                print(f"       {region['sequence'][:20]}...")
        else:
            print("No strong transmembrane regions detected (threshold > 1.6)")
            print("This suggests the protein is SOLUBLE, not membrane-embedded")
    
    # Classify topology
    if len(merged_tm) == 0:
        results['topology'] = 'soluble'
    elif len(merged_tm) == 1:
        results['topology'] = 'single-pass membrane'
    elif len(merged_tm) <= 3:
        results['topology'] = 'multi-pass membrane'
    else:
        results['topology'] = f'polytopic membrane ({len(merged_tm)} TM)'
    
    return results

def analyze_subcellular_features(sequence, verbose=True):
    """Analyze features that indicate subcellular localization"""
    results = {
        'nuclear_signals': [],
        'mitochondrial_signals': [],
        'peroxisomal_signals': [],
        'er_signals': [],
        'features': {}
    }
    
    if verbose:
        print("\n=== Subcellular Localization Signals ===")
    
    # Nuclear localization signals (NLS)
    # Classical monopartite: K(K/R)X(K/R)
    import re
    nls_mono = re.finditer(r'K[KR].[KR]', sequence)
    for match in nls_mono:
        results['nuclear_signals'].append({
            'type': 'monopartite NLS',
            'position': match.start(),
            'sequence': match.group()
        })
    
    # Bipartite NLS: (K/R)(K/R)X10-12(K/R)3-5
    nls_bi = re.finditer(r'[KR]{2}.{10,12}[KR]{3,5}', sequence)
    for match in nls_bi:
        results['nuclear_signals'].append({
            'type': 'bipartite NLS',
            'position': match.start(),
            'sequence': match.group()[:20] + '...'
        })
    
    # Mitochondrial targeting (typically N-terminal, rich in R, L, S, A)
    if len(sequence) >= 30:
        first_30 = sequence[:30]
        mito_score = sum(1 for aa in first_30 if aa in 'RLSA') / 30
        if mito_score > 0.4:
            results['mitochondrial_signals'].append({
                'type': 'N-terminal targeting',
                'score': round(mito_score, 2)
            })
    
    # Peroxisomal targeting signal (PTS1: SKL or variants at C-terminus)
    if len(sequence) >= 3:
        c_term = sequence[-3:]
        if c_term in ['SKL', 'AKL', 'CKL', 'SKI', 'SRL', 'ARL']:
            results['peroxisomal_signals'].append({
                'type': 'PTS1',
                'sequence': c_term
            })
    
    # ER retention signal (KDEL or variants at C-terminus)
    if len(sequence) >= 4:
        c_term = sequence[-4:]
        if c_term in ['KDEL', 'HDEL', 'RDEL', 'KEEL']:
            results['er_signals'].append({
                'type': 'ER retention',
                'sequence': c_term
            })
    
    # General features
    cys_count = sequence.count('C')
    results['features']['cysteines'] = cys_count
    results['features']['potential_disulfides'] = cys_count // 2
    
    # Charge distribution
    basic = sum(1 for aa in sequence if aa in 'RKH')
    acidic = sum(1 for aa in sequence if aa in 'DE')
    results['features']['basic_residues'] = basic
    results['features']['acidic_residues'] = acidic
    results['features']['net_charge'] = basic - acidic
    results['features']['length'] = len(sequence)
    
    if verbose:
        if results['nuclear_signals']:
            print(f"Nuclear signals: {len(results['nuclear_signals'])}")
        if results['mitochondrial_signals']:
            print(f"Mitochondrial signals: {len(results['mitochondrial_signals'])}")
        if results['peroxisomal_signals']:
            print(f"Peroxisomal signals: {len(results['peroxisomal_signals'])}")
        if results['er_signals']:
            print(f"ER signals: {len(results['er_signals'])}")
        
        print(f"\nProtein features:")
        print(f"  Length: {results['features']['length']} aa")
        print(f"  Cysteines: {results['features']['cysteines']} (potential for {results['features']['potential_disulfides']} disulfide bonds)")
        print(f"  Net charge: {results['features']['net_charge']}")
    
    return results

def predict_localization(signal_results, tm_results, features_results, verbose=True):
    """Predict overall cellular localization based on all features"""
    predictions = {
        'primary_location': 'unknown',
        'confidence': 'low',
        'supporting_evidence': [],
        'go_terms': []
    }
    
    if verbose:
        print("\n=== Localization Prediction ===")
    
    # Decision tree for localization
    has_signal = signal_results['has_signal_peptide']
    has_tm = tm_results['is_membrane_protein']
    topology = tm_results['topology']
    
    # Check for specific organellar signals
    has_nuclear = len(features_results['nuclear_signals']) > 0
    has_mito = len(features_results['mitochondrial_signals']) > 0
    has_perox = len(features_results['peroxisomal_signals']) > 0
    has_er = len(features_results['er_signals']) > 0
    
    # Bacterial proteins
    if has_signal and not has_tm:
        predictions['primary_location'] = 'periplasmic space (bacterial)'
        predictions['confidence'] = 'high'
        predictions['supporting_evidence'].append('Signal peptide present')
        predictions['supporting_evidence'].append('No transmembrane regions')
        predictions['go_terms'].append('GO:0042597 (periplasmic space)')
    elif has_signal and has_tm:
        predictions['primary_location'] = 'membrane (likely inner/plasma membrane)'
        predictions['confidence'] = 'high'
        predictions['supporting_evidence'].append('Signal peptide present')
        predictions['supporting_evidence'].append(f'{len(tm_results["tm_regions"])} TM regions')
        predictions['go_terms'].append('GO:0016020 (membrane)')
    
    # Eukaryotic proteins
    elif has_nuclear:
        predictions['primary_location'] = 'nucleus'
        predictions['confidence'] = 'moderate'
        predictions['supporting_evidence'].append('Nuclear localization signal(s)')
        predictions['go_terms'].append('GO:0005634 (nucleus)')
    elif has_mito:
        predictions['primary_location'] = 'mitochondrion'
        predictions['confidence'] = 'moderate'
        predictions['supporting_evidence'].append('Mitochondrial targeting signal')
        predictions['go_terms'].append('GO:0005739 (mitochondrion)')
    elif has_perox:
        predictions['primary_location'] = 'peroxisome'
        predictions['confidence'] = 'moderate'
        predictions['supporting_evidence'].append('Peroxisomal targeting signal')
        predictions['go_terms'].append('GO:0005777 (peroxisome)')
    elif has_er:
        predictions['primary_location'] = 'endoplasmic reticulum'
        predictions['confidence'] = 'moderate'
        predictions['supporting_evidence'].append('ER retention signal')
        predictions['go_terms'].append('GO:0005783 (endoplasmic reticulum)')
    elif not has_signal and not has_tm:
        predictions['primary_location'] = 'cytoplasm'
        predictions['confidence'] = 'moderate'
        predictions['supporting_evidence'].append('No targeting signals')
        predictions['supporting_evidence'].append('No transmembrane regions')
        predictions['go_terms'].append('GO:0005737 (cytoplasm)')
    
    # Add disulfide bond information
    if features_results['features']['potential_disulfides'] > 0:
        if 'periplasmic' in predictions['primary_location']:
            predictions['supporting_evidence'].append(
                f'{features_results["features"]["potential_disulfides"]} potential disulfide bonds (common in periplasm)')
    
    if verbose:
        print(f"Predicted location: {predictions['primary_location']}")
        print(f"Confidence: {predictions['confidence']}")
        print("\nSupporting evidence:")
        for evidence in predictions['supporting_evidence']:
            print(f"  - {evidence}")
        print("\nRecommended GO terms:")
        for go_term in predictions['go_terms']:
            print(f"  - {go_term}")
    
    return predictions

def main():
    parser = argparse.ArgumentParser(
        description='Analyze protein subcellular localization using sequence features',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze from UniProt ID
  %(prog)s --uniprot Q88JH0
  
  # Analyze from FASTA file
  %(prog)s --fasta protein.fasta
  
  # Save results to JSON
  %(prog)s --uniprot Q88JH0 --output results.json
  
  # Quiet mode (JSON output only)
  %(prog)s --uniprot Q88JH0 --quiet --output results.json
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--uniprot', '-u', help='UniProt ID to analyze')
    input_group.add_argument('--fasta', '-f', help='FASTA file to analyze')
    
    parser.add_argument('--output', '-o', help='Output JSON file for results')
    parser.add_argument('--quiet', '-q', action='store_true', 
                       help='Suppress verbose output')
    
    args = parser.parse_args()
    
    # Get sequence
    if args.uniprot:
        sequence, description = fetch_uniprot_sequence(args.uniprot)
        if not sequence:
            print(f"Error: Failed to fetch sequence for UniProt ID {args.uniprot}", 
                  file=sys.stderr)
            sys.exit(1)
        protein_id = args.uniprot
    else:
        try:
            sequence, description = read_fasta_file(args.fasta)
            protein_id = Path(args.fasta).stem
        except Exception as e:
            print(f"Error reading FASTA file: {e}", file=sys.stderr)
            sys.exit(1)
    
    verbose = not args.quiet
    
    if verbose:
        print("=" * 60)
        print("Protein Localization Analysis")
        print(f"Protein: {protein_id}")
        if description:
            print(f"Description: {description[:80]}...")
        print("=" * 60)
        print(f"\nProtein length: {len(sequence)} amino acids")
    
    # Perform analyses
    signal_results = analyze_signal_peptide(sequence, verbose=verbose)
    tm_results = check_transmembrane_regions(sequence, verbose=verbose)
    features_results = analyze_subcellular_features(sequence, verbose=verbose)
    predictions = predict_localization(signal_results, tm_results, features_results, verbose=verbose)
    
    # Compile all results
    all_results = {
        'protein_id': protein_id,
        'description': description,
        'sequence_length': len(sequence),
        'signal_peptide': signal_results,
        'transmembrane': tm_results,
        'subcellular_features': features_results,
        'localization_prediction': predictions
    }
    
    # Save results if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(all_results, f, indent=2)
        if verbose:
            print(f"\nResults saved to {args.output}")
    
    if verbose:
        print("\n" + "=" * 60)
        print("Analysis complete!")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())