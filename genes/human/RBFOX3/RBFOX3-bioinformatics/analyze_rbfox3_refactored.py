#!/usr/bin/env python3
"""
Generic RNA-binding protein analysis tool
Analyzes RRM domains, disorder regions, and protein interaction motifs
"""

import requests
from Bio import SeqIO
from io import StringIO
import json
import re
import argparse
import sys
from pathlib import Path

def fetch_uniprot_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        record = SeqIO.read(StringIO(response.text), "fasta")
        return str(record.seq), record.description
    return None, None

def read_fasta_file(filepath):
    """Read sequence from FASTA file"""
    with open(filepath, 'r') as f:
        record = SeqIO.read(f, "fasta")
        return str(record.seq), record.description

def analyze_rrm_domain(sequence, verbose=True):
    """Analyze RNA Recognition Motif (RRM) domain"""
    results = {}
    
    # RRM domains are typically 80-90 aa with conserved RNP motifs
    # RNP-1: [RK]-G-[FY]-[GA]-[FY]-[ILV]-X-[FY] (octamer)
    # RNP-2: [ILV]-[FY]-[ILV]-X-N-L (hexamer)
    
    # Simplified RNP-1 pattern
    rnp1_pattern = r'[RK]G[FY][GA][FY][ILV].[FY]'
    rnp1_matches = list(re.finditer(rnp1_pattern, sequence))
    
    results['rnp1_motifs'] = []
    results['predicted_rrm_domains'] = []
    
    if rnp1_matches:
        if verbose:
            print(f"RNP-1 motifs found: {len(rnp1_matches)}")
        for i, match in enumerate(rnp1_matches, 1):
            pos = match.start()
            motif_info = {
                'number': i,
                'start': pos,
                'end': match.end(),
                'sequence': match.group()
            }
            # RRM typically extends ~40 aa before and after RNP-1
            rrm_start = max(0, pos - 40)
            rrm_end = min(len(sequence), pos + 50)
            motif_info['predicted_rrm'] = {'start': rrm_start, 'end': rrm_end}
            results['rnp1_motifs'].append(motif_info)
            
            if verbose:
                print(f"  RNP-1 #{i}: Position {pos}-{match.end()}")
                print(f"    Sequence: {match.group()}")
                print(f"    Predicted RRM domain: {rrm_start}-{rrm_end}")
    
    # Simplified RNP-2 pattern  
    rnp2_pattern = r'[ILV][FY][ILV].NL'
    rnp2_matches = list(re.finditer(rnp2_pattern, sequence))
    
    results['rnp2_motifs'] = []
    if rnp2_matches:
        if verbose:
            print(f"\nRNP-2 motifs found: {len(rnp2_matches)}")
        for i, match in enumerate(rnp2_matches, 1):
            motif_info = {
                'number': i,
                'start': match.start(),
                'end': match.end(),
                'sequence': match.group()
            }
            results['rnp2_motifs'].append(motif_info)
            if verbose:
                print(f"  RNP-2 #{i}: Position {match.start()}-{match.end()}: {match.group()}")
    
    if not (rnp1_matches or rnp2_matches) and verbose:
        print("No canonical RNP motifs detected")
    
    results['has_rrm'] = len(rnp1_matches) > 0
    return results

def analyze_disorder(sequence, window_size=30, disorder_threshold=0.5, verbose=True):
    """Predict intrinsically disordered regions"""
    results = {}
    
    # Simple disorder prediction based on composition
    # Disordered regions are enriched in P, E, S, K, R, Q and depleted in W, C, F, Y, I, V, L
    
    disorder_promoting = set('PESKRQ')
    order_promoting = set('WCFYIVL')
    
    disordered_regions = []
    
    for i in range(0, len(sequence) - window_size):
        window = sequence[i:i+window_size]
        disorder_score = sum(1 for aa in window if aa in disorder_promoting) / window_size
        order_score = sum(1 for aa in window if aa in order_promoting) / window_size
        
        net_disorder = disorder_score - order_score
        
        if net_disorder > disorder_threshold:
            disordered_regions.append((i, i+window_size, net_disorder))
    
    # Merge overlapping regions
    merged = []
    if disordered_regions:
        for start, end, score in disordered_regions:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]), max(score, merged[-1][2]))
            else:
                merged.append((start, end, score))
        
        if verbose:
            print(f"Predicted disordered regions: {len(merged)}")
        
        results['disordered_regions'] = []
        for start, end, score in merged:
            length = end - start
            region_seq = sequence[start:end]
            region_info = {
                'start': start,
                'end': end,
                'length': length,
                'score': round(score, 2),
                'composition': {
                    'S': region_seq.count('S'),
                    'P': region_seq.count('P'),
                    'E': region_seq.count('E'),
                    'K': region_seq.count('K'),
                    'R': region_seq.count('R'),
                    'Q': region_seq.count('Q')
                }
            }
            results['disordered_regions'].append(region_info)
            
            if verbose:
                print(f"  Position {start}-{end} (length: {length} aa, score: {score:.2f})")
                print(f"    Composition: S={region_info['composition']['S']}, "
                      f"P={region_info['composition']['P']}, E={region_info['composition']['E']}")
    else:
        results['disordered_regions'] = []
        if verbose:
            print("No significant disordered regions predicted")
    
    results['has_disorder'] = len(results['disordered_regions']) > 0
    return results

def analyze_rna_binding_features(sequence, verbose=True):
    """Analyze generic RNA binding features"""
    results = {}
    
    if verbose:
        print("\n=== RNA Binding Feature Analysis ===")
        print("Checking for conserved RNA-binding features...")
    
    # Check for basic patches that interact with RNA backbone
    basic_patch_pattern = r'[RK]{2,}'
    basic_patches = list(re.finditer(basic_patch_pattern, sequence))
    
    results['basic_patches'] = []
    if basic_patches:
        if verbose:
            print(f"Basic patches for RNA binding: {len(basic_patches)}")
        for i, match in enumerate(basic_patches, 1):
            patch_info = {
                'number': i,
                'start': match.start(),
                'end': match.end(),
                'sequence': match.group()
            }
            results['basic_patches'].append(patch_info)
            if verbose and i <= 3:
                print(f"  Patch {i}: Position {match.start()}: {match.group()}")
    
    # Check for aromatic residues important for RNA base stacking
    aromatic_clusters = []
    for i in range(0, len(sequence) - 10):
        window = sequence[i:i+10]
        aromatic_count = sum(1 for aa in window if aa in 'FYW')
        if aromatic_count >= 3:
            aromatic_clusters.append((i, i+10, aromatic_count))
    
    results['aromatic_clusters'] = aromatic_clusters
    if aromatic_clusters and verbose:
        print(f"Aromatic-rich regions (potential RNA base stacking): {len(aromatic_clusters)}")
    
    return results

def analyze_protein_interactions(sequence, verbose=True):
    """Analyze potential protein interaction regions"""
    results = {}
    
    if verbose:
        print("\n=== Protein Interaction Region Analysis ===")
    
    # Check for proline-rich regions (often mediate protein interactions)
    proline_rich_pattern = r'P{2,}|P.P.P'
    p_rich = list(re.finditer(proline_rich_pattern, sequence))
    
    results['proline_rich'] = []
    if p_rich:
        if verbose:
            print(f"Proline-rich regions: {len(p_rich)}")
        for match in p_rich:
            region_info = {
                'start': match.start(),
                'end': match.end(),
                'sequence': match.group()
            }
            results['proline_rich'].append(region_info)
            if verbose and len(results['proline_rich']) <= 3:
                print(f"  Position {match.start()}-{match.end()}: {match.group()}")
    
    # Check for SH3 binding motifs (PxxP)
    sh3_pattern = r'P..P'
    sh3_motifs = list(re.finditer(sh3_pattern, sequence))
    
    results['sh3_motifs'] = []
    if sh3_motifs:
        if verbose:
            print(f"\nPotential SH3 binding motifs (PxxP): {len(sh3_motifs)}")
        for match in sh3_motifs:
            results['sh3_motifs'].append({
                'start': match.start(),
                'end': match.end(),
                'sequence': match.group()
            })
        if verbose:
            print(f"  First 3 motifs: {[m.group() for m in sh3_motifs[:3]]}")
    
    # Check for WW domain binding motifs (PPxY)
    ww_pattern = r'PP.Y'
    ww_motifs = list(re.finditer(ww_pattern, sequence))
    
    results['ww_motifs'] = []
    if ww_motifs:
        if verbose:
            print(f"\nPotential WW domain binding motifs (PPxY): {len(ww_motifs)}")
        for match in ww_motifs:
            results['ww_motifs'].append({
                'start': match.start(),
                'end': match.end(),
                'sequence': match.group()
            })
            if verbose:
                print(f"  Position {match.start()}: {match.group()}")
    
    return results

def analyze_post_translational_modifications(sequence, verbose=True):
    """Predict post-translational modification sites"""
    results = {}
    
    if verbose:
        print("\n=== Post-Translational Modification Sites ===")
    
    # Phosphorylation sites
    pka_pattern = r'R.[ST]'  # PKA consensus
    ck2_pattern = r'[ST]..E'  # CK2 consensus
    
    pka_sites = list(re.finditer(pka_pattern, sequence))
    ck2_sites = list(re.finditer(ck2_pattern, sequence))
    
    results['phosphorylation'] = {
        'pka_sites': [{'position': m.start(), 'sequence': m.group()} for m in pka_sites],
        'ck2_sites': [{'position': m.start(), 'sequence': m.group()} for m in ck2_sites]
    }
    
    if verbose:
        print(f"Potential phosphorylation sites:")
        print(f"  PKA sites (R-x-S/T): {len(pka_sites)}")
        if pka_sites:
            print(f"    Positions: {[m.start() for m in pka_sites[:5]]}...")
        
        print(f"  CK2 sites (S/T-x-x-E): {len(ck2_sites)}")
        if ck2_sites:
            print(f"    Positions: {[m.start() for m in ck2_sites[:5]]}...")
    
    # SUMOylation sites (ΨKxE where Ψ is hydrophobic)
    sumo_pattern = r'[VILMF]K.E'
    sumo_sites = list(re.finditer(sumo_pattern, sequence))
    
    results['sumoylation_sites'] = []
    if sumo_sites:
        if verbose:
            print(f"\nPotential SUMOylation sites: {len(sumo_sites)}")
        for match in sumo_sites:
            results['sumoylation_sites'].append({
                'position': match.start(),
                'sequence': match.group()
            })
            if verbose and len(results['sumoylation_sites']) <= 3:
                print(f"  Position {match.start()}: {match.group()}")
    
    return results

def calculate_properties(sequence, verbose=True):
    """Calculate basic protein properties"""
    results = {}
    
    if verbose:
        print("\n=== Protein Properties ===")
    
    length = len(sequence)
    results['length'] = length
    
    if verbose:
        print(f"Length: {length} amino acids")
    
    # Molecular weight
    mw_table = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    mw = sum(mw_table.get(aa, 110) for aa in sequence) / 1000
    results['molecular_weight_kda'] = round(mw, 1)
    
    if verbose:
        print(f"Molecular weight: ~{mw:.1f} kDa")
    
    # Composition
    basic = sum(1 for aa in sequence if aa in 'RKH')
    acidic = sum(1 for aa in sequence if aa in 'DE')
    
    results['basic_residues'] = basic
    results['acidic_residues'] = acidic
    results['net_charge'] = basic - acidic
    
    if verbose:
        print(f"Basic residues: {basic}")
        print(f"Acidic residues: {acidic}")
        print(f"Net charge at pH 7: {basic - acidic}")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Analyze RNA-binding proteins for RRM domains and functional features',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze from UniProt ID
  %(prog)s --uniprot A6NFN3
  
  # Analyze from FASTA file
  %(prog)s --fasta protein.fasta
  
  # Save results to JSON
  %(prog)s --uniprot A6NFN3 --output results.json
  
  # Quiet mode (JSON output only)
  %(prog)s --uniprot A6NFN3 --quiet --output results.json
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--uniprot', '-u', help='UniProt ID to analyze')
    input_group.add_argument('--fasta', '-f', help='FASTA file to analyze')
    
    parser.add_argument('--output', '-o', help='Output JSON file for results')
    parser.add_argument('--quiet', '-q', action='store_true', 
                       help='Suppress verbose output (only show summary)')
    
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
        print("RNA-Binding Protein Analysis")
        print(f"Protein: {protein_id}")
        if description:
            print(f"Description: {description[:80]}...")
        print("=" * 60)
        print(f"\nSequence retrieved: {len(sequence)} amino acids")
    
    # Perform analyses
    all_results = {
        'protein_id': protein_id,
        'description': description,
        'properties': calculate_properties(sequence, verbose=verbose),
        'rrm_analysis': analyze_rrm_domain(sequence, verbose=verbose),
        'disorder_analysis': analyze_disorder(sequence, verbose=verbose),
        'rna_binding_features': analyze_rna_binding_features(sequence, verbose=verbose),
        'protein_interactions': analyze_protein_interactions(sequence, verbose=verbose),
        'ptm_sites': analyze_post_translational_modifications(sequence, verbose=verbose)
    }
    
    # Generate conclusions
    if verbose:
        print("\n" + "=" * 60)
        print("ANALYSIS SUMMARY")
        print("=" * 60)
        
        print("\n1. DOMAIN STRUCTURE:")
        if all_results['rrm_analysis']['has_rrm']:
            print("   ✓ Contains RNA Recognition Motif (RRM) domain")
            print(f"   - {len(all_results['rrm_analysis']['rnp1_motifs'])} RNP-1 motif(s)")
            print(f"   - {len(all_results['rrm_analysis']['rnp2_motifs'])} RNP-2 motif(s)")
        else:
            print("   ✗ No canonical RRM domain detected")
        
        print("\n2. FUNCTIONAL FEATURES:")
        if all_results['disorder_analysis']['has_disorder']:
            print(f"   ✓ Contains {len(all_results['disorder_analysis']['disordered_regions'])} "
                  f"intrinsically disordered region(s)")
        else:
            print("   ✗ No significant disordered regions")
        
        if all_results['rna_binding_features']['basic_patches']:
            print(f"   ✓ {len(all_results['rna_binding_features']['basic_patches'])} "
                  f"basic patches for RNA binding")
        
        print("\n3. PROTEIN INTERACTIONS:")
        interaction_count = (
            len(all_results['protein_interactions']['proline_rich']) +
            len(all_results['protein_interactions']['sh3_motifs']) +
            len(all_results['protein_interactions']['ww_motifs'])
        )
        if interaction_count > 0:
            print(f"   ✓ {interaction_count} potential protein interaction motifs")
        else:
            print("   ✗ No significant interaction motifs detected")
        
        print("\n4. POST-TRANSLATIONAL MODIFICATIONS:")
        ptm_count = (
            len(all_results['ptm_sites']['phosphorylation']['pka_sites']) +
            len(all_results['ptm_sites']['phosphorylation']['ck2_sites']) +
            len(all_results['ptm_sites']['sumoylation_sites'])
        )
        print(f"   ✓ {ptm_count} potential PTM sites predicted")
    
    # Save results if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(all_results, f, indent=2)
        if verbose:
            print(f"\nResults saved to {args.output}")
    
    # Return success
    return 0

if __name__ == "__main__":
    sys.exit(main())