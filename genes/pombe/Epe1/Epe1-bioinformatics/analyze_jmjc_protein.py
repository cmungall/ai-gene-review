#!/usr/bin/env python3
"""
Generic JmjC domain-containing protein analyzer
Analyzes proteins for JmjC domains, demethylase activity features, and chromatin-related motifs
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

def analyze_jmjc_domain(sequence, jmjc_start=None, jmjc_end=None, verbose=True):
    """
    Analyze JmjC (Jumonji C) domain features
    
    Args:
        sequence: Protein sequence
        jmjc_start: Start position of JmjC domain (if known)
        jmjc_end: End position of JmjC domain (if known)
        verbose: Print detailed output
    """
    results = {
        'histidines': [],
        'fe_binding_motifs': [],
        'jmjc_aromatic_content': {},
        'has_jmjc_features': False
    }
    
    if verbose:
        print("\n=== JmjC Domain Analysis ===")
    
    # JmjC domains typically contain Fe(II) and α-ketoglutarate binding residues
    # Key conserved residues: HxD/E...H motif for Fe(II) coordination
    
    # Look for histidine residues that could coordinate Fe(II)
    histidines = [i for i, aa in enumerate(sequence) if aa == 'H']
    results['histidines'] = histidines
    
    if verbose:
        print(f"Total histidine residues: {len(histidines)}")
    
    # Check for HxD or HxE motifs (Fe(II) binding)
    hxd_pattern = r'H.[DE]'
    hxd_matches = list(re.finditer(hxd_pattern, sequence))
    
    for match in hxd_matches:
        pos = match.start()
        motif_info = {
            'position': pos,
            'sequence': match.group(),
            'partner_histidines': []
        }
        
        # Check for second histidine ~100-200 aa downstream (typical spacing)
        downstream_h = [h for h in histidines if pos + 100 <= h <= pos + 200]
        if downstream_h:
            motif_info['partner_histidines'] = downstream_h
        
        results['fe_binding_motifs'].append(motif_info)
    
    if verbose and hxd_matches:
        print(f"\nPotential Fe(II) binding motifs (HxD/E): {len(hxd_matches)}")
        for i, motif in enumerate(results['fe_binding_motifs'][:5], 1):
            print(f"  Motif {i}: Position {motif['position']}: {motif['sequence']}")
            if motif['partner_histidines']:
                print(f"    Potential partner H at: {motif['partner_histidines'][0]}")
    
    # Analyze JmjC region if bounds provided, otherwise scan for aromatic-rich regions
    if jmjc_start is None:
        # Scan for aromatic-rich regions that might be JmjC
        window_size = 200
        best_region = None
        best_score = 0
        
        for start in range(0, len(sequence) - window_size, 50):
            end = min(start + window_size, len(sequence))
            region_seq = sequence[start:end]
            f_count = region_seq.count('F')
            y_count = region_seq.count('Y')
            k_count = region_seq.count('K')
            score = f_count + y_count + (k_count * 0.5)
            
            if score > best_score:
                best_score = score
                best_region = (start, end, f_count, y_count, k_count)
        
        if best_region:
            start, end, f_count, y_count, k_count = best_region
            results['jmjc_aromatic_content'] = {
                'region': f"{start}-{end}",
                'F_count': f_count,
                'Y_count': y_count,
                'K_count': k_count,
                'aromatic_rich': (f_count + y_count) > 5
            }
            
            if verbose:
                print(f"\nPotential JmjC region (positions {start}-{end}):")
                print(f"  Aromatic residues: F={f_count}, Y={y_count}")
                print(f"  Basic residues: K={k_count}")
                if (f_count + y_count) > 5:
                    print("  Aromatic content consistent with JmjC domain")
    else:
        # Use provided JmjC boundaries
        jmjc_end = jmjc_end or min(jmjc_start + 200, len(sequence))
        if jmjc_start < len(sequence):
            jmjc_seq = sequence[jmjc_start:jmjc_end]
            f_count = jmjc_seq.count('F')
            y_count = jmjc_seq.count('Y')
            k_count = jmjc_seq.count('K')
            
            results['jmjc_aromatic_content'] = {
                'region': f"{jmjc_start}-{jmjc_end}",
                'F_count': f_count,
                'Y_count': y_count,
                'K_count': k_count,
                'aromatic_rich': (f_count + y_count) > 5
            }
            
            if verbose:
                print(f"\nJmjC region analysis (positions {jmjc_start}-{jmjc_end}):")
                print(f"  Aromatic residues: F={f_count}, Y={y_count}")
                print(f"  Basic residues: K={k_count}")
    
    results['has_jmjc_features'] = len(results['fe_binding_motifs']) > 0
    return results

def analyze_demethylase_activity(sequence, verbose=True):
    """Analyze potential demethylase activity features"""
    results = {
        'akg_binding_motifs': [],
        'basic_patches': [],
        'catalytic_potential': 'unknown'
    }
    
    if verbose:
        print("\n=== Demethylase Activity Prediction ===")
    
    # Check for α-ketoglutarate binding residues
    # Typically involves R and N residues in specific positions
    akg_pattern = r'R.{5}R'
    akg_matches = list(re.finditer(akg_pattern, sequence))
    
    for match in akg_matches:
        results['akg_binding_motifs'].append({
            'position': match.start(),
            'sequence': match.group()
        })
    
    if verbose and akg_matches:
        print(f"Potential α-ketoglutarate binding motifs: {len(akg_matches)}")
        for motif in results['akg_binding_motifs'][:3]:
            print(f"  Position {motif['position']}: {motif['sequence']}")
    
    # Check for substrate binding groove characteristics
    window = 10
    for i in range(0, len(sequence) - window):
        window_seq = sequence[i:i+window]
        basic_count = sum(1 for aa in window_seq if aa in 'RKH')
        if basic_count >= 4:
            results['basic_patches'].append({
                'position': i,
                'end': i + window,
                'basic_count': basic_count
            })
    
    if verbose and results['basic_patches']:
        print(f"\nBasic patches for histone binding: {len(results['basic_patches'])}")
        print("  Top 3 basic-rich regions:")
        sorted_patches = sorted(results['basic_patches'], 
                              key=lambda x: x['basic_count'], reverse=True)
        for patch in sorted_patches[:3]:
            print(f"    Position {patch['position']}-{patch['end']}: "
                  f"{patch['basic_count']} basic residues")
    
    # Assess catalytic potential
    if len(results['akg_binding_motifs']) > 0 and len(results['basic_patches']) > 5:
        results['catalytic_potential'] = 'high'
    elif len(results['akg_binding_motifs']) > 0 or len(results['basic_patches']) > 3:
        results['catalytic_potential'] = 'moderate'
    else:
        results['catalytic_potential'] = 'low'
    
    if verbose:
        print(f"\nCatalytic potential: {results['catalytic_potential']}")
    
    return results

def analyze_chromatin_features(sequence, verbose=True):
    """Analyze features related to chromatin localization and binding"""
    results = {
        'hp1_binding_motifs': [],
        'aromatic_clusters': [],
        'nuclear_localization_signals': {
            'monopartite': [],
            'bipartite': []
        }
    }
    
    if verbose:
        print("\n=== Chromatin Localization Features ===")
    
    # Check for HP1/Swi6 interaction motifs (PxVxL)
    pxvxl_pattern = r'P.V.L'
    pxvxl_matches = list(re.finditer(pxvxl_pattern, sequence))
    
    for match in pxvxl_matches:
        results['hp1_binding_motifs'].append({
            'position': match.start(),
            'sequence': match.group()
        })
    
    if verbose:
        if pxvxl_matches:
            print(f"HP1/Swi6 binding motifs (PxVxL): {len(pxvxl_matches)}")
            for motif in results['hp1_binding_motifs']:
                print(f"  Position {motif['position']}: {motif['sequence']}")
        else:
            print("No canonical PxVxL motifs found")
    
    # Check for aromatic clusters (methyl-lysine recognition)
    window = 20
    for i in range(0, len(sequence) - window):
        window_seq = sequence[i:i+window]
        aromatic_count = sum(1 for aa in window_seq if aa in 'FYW')
        if aromatic_count >= 4:
            results['aromatic_clusters'].append({
                'position': i,
                'end': i + window,
                'aromatic_count': aromatic_count,
                'sequence': window_seq[:10] + '...'
            })
    
    if verbose and results['aromatic_clusters']:
        print(f"\nAromatic clusters (potential methyl-lysine binding): "
              f"{len(results['aromatic_clusters'])}")
        for cluster in results['aromatic_clusters'][:3]:
            print(f"  Position {cluster['position']}-{cluster['end']}: "
                  f"{cluster['aromatic_count']} aromatic residues")
            print(f"    Sequence: {cluster['sequence']}")
    
    # Nuclear localization signals
    # Monopartite NLS
    mono_nls_pattern = r'K[KR].[KR]'
    mono_matches = list(re.finditer(mono_nls_pattern, sequence))
    
    for match in mono_matches:
        results['nuclear_localization_signals']['monopartite'].append({
            'position': match.start(),
            'sequence': match.group()
        })
    
    # Bipartite NLS
    bipartite_pattern = r'[KR]{2}.{10,12}[KR]{3,5}'
    bipartite_matches = list(re.finditer(bipartite_pattern, sequence))
    
    for match in bipartite_matches:
        results['nuclear_localization_signals']['bipartite'].append({
            'position': match.start(),
            'end': match.end(),
            'sequence': match.group()[:20] + '...' if len(match.group()) > 20 else match.group()
        })
    
    if verbose:
        print("\n=== Nuclear Localization Signals ===")
        if results['nuclear_localization_signals']['monopartite']:
            print(f"Monopartite NLS: {len(results['nuclear_localization_signals']['monopartite'])}")
            for nls in results['nuclear_localization_signals']['monopartite'][:3]:
                print(f"  Position {nls['position']}: {nls['sequence']}")
        
        if results['nuclear_localization_signals']['bipartite']:
            print(f"Bipartite NLS: {len(results['nuclear_localization_signals']['bipartite'])}")
            for nls in results['nuclear_localization_signals']['bipartite'][:2]:
                print(f"  Position {nls['position']}-{nls['end']}")
    
    return results

def calculate_properties(sequence, verbose=True):
    """Calculate basic protein properties"""
    results = {}
    
    if verbose:
        print("\n=== Protein Properties ===")
    
    length = len(sequence)
    results['length'] = length
    
    # Molecular weight
    mw_table = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    mw = sum(mw_table.get(aa, 110) for aa in sequence) / 1000
    results['molecular_weight_kda'] = round(mw, 1)
    
    # Charge
    basic = sum(1 for aa in sequence if aa in 'RKH')
    acidic = sum(1 for aa in sequence if aa in 'DE')
    results['basic_residues'] = basic
    results['acidic_residues'] = acidic
    results['net_charge'] = basic - acidic
    
    # Composition
    results['serine_percent'] = round(sequence.count('S') / length * 100, 1)
    results['proline_percent'] = round(sequence.count('P') / length * 100, 1)
    
    if verbose:
        print(f"Length: {length} amino acids")
        print(f"Molecular weight: ~{mw:.1f} kDa")
        print(f"Basic residues: {basic}")
        print(f"Acidic residues: {acidic}")
        print(f"Net charge at pH 7: {basic - acidic}")
        print(f"\nComposition features:")
        print(f"  Serine content: {results['serine_percent']}% (phosphorylation potential)")
        print(f"  Proline content: {results['proline_percent']}%")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Analyze JmjC domain-containing proteins for demethylase and chromatin features',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze from UniProt ID
  %(prog)s --uniprot O94603
  
  # Analyze from FASTA file
  %(prog)s --fasta protein.fasta
  
  # Specify JmjC domain boundaries if known
  %(prog)s --uniprot O94603 --jmjc-start 400 --jmjc-end 600
  
  # Save results to JSON
  %(prog)s --uniprot O94603 --output results.json
  
  # Quiet mode (JSON output only)
  %(prog)s --uniprot O94603 --quiet --output results.json
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--uniprot', '-u', help='UniProt ID to analyze')
    input_group.add_argument('--fasta', '-f', help='FASTA file to analyze')
    
    parser.add_argument('--jmjc-start', type=int, help='Start position of JmjC domain')
    parser.add_argument('--jmjc-end', type=int, help='End position of JmjC domain')
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
        print("JmjC Domain Protein Analysis")
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
        'jmjc_analysis': analyze_jmjc_domain(sequence, 
                                            jmjc_start=args.jmjc_start,
                                            jmjc_end=args.jmjc_end,
                                            verbose=verbose),
        'demethylase_activity': analyze_demethylase_activity(sequence, verbose=verbose),
        'chromatin_features': analyze_chromatin_features(sequence, verbose=verbose)
    }
    
    # Generate conclusions
    if verbose:
        print("\n" + "=" * 60)
        print("ANALYSIS SUMMARY")
        print("=" * 60)
        
        print("\n1. DOMAIN STRUCTURE:")
        if all_results['jmjc_analysis']['has_jmjc_features']:
            print("   ✓ Contains JmjC domain features")
            print(f"   - {len(all_results['jmjc_analysis']['fe_binding_motifs'])} "
                  f"Fe(II) binding motifs")
            if all_results['jmjc_analysis']['jmjc_aromatic_content']:
                if all_results['jmjc_analysis']['jmjc_aromatic_content']['aromatic_rich']:
                    print("   - Aromatic content consistent with JmjC domain")
        else:
            print("   ✗ No clear JmjC domain features detected")
        
        print("\n2. CATALYTIC ACTIVITY:")
        catalytic = all_results['demethylase_activity']['catalytic_potential']
        if catalytic == 'high':
            print("   ✓ High potential for demethylase activity")
        elif catalytic == 'moderate':
            print("   ⚠ Moderate potential for demethylase activity")
        else:
            print("   ✗ Low potential for demethylase activity")
        
        print(f"   - {len(all_results['demethylase_activity']['akg_binding_motifs'])} "
              f"α-ketoglutarate binding motifs")
        print(f"   - {len(all_results['demethylase_activity']['basic_patches'])} "
              f"basic patches for substrate binding")
        
        print("\n3. CHROMATIN FEATURES:")
        chrom_features = all_results['chromatin_features']
        if chrom_features['hp1_binding_motifs']:
            print(f"   ✓ {len(chrom_features['hp1_binding_motifs'])} HP1 binding motifs")
        if chrom_features['aromatic_clusters']:
            print(f"   ✓ {len(chrom_features['aromatic_clusters'])} "
                  f"aromatic clusters for methyl-lysine recognition")
        
        nls_count = (len(chrom_features['nuclear_localization_signals']['monopartite']) +
                    len(chrom_features['nuclear_localization_signals']['bipartite']))
        if nls_count > 0:
            print(f"   ✓ {nls_count} nuclear localization signals")
        
        print("\n4. PROTEIN CLASS:")
        if all_results['jmjc_analysis']['has_jmjc_features']:
            if catalytic in ['high', 'moderate']:
                print("   Likely active histone demethylase")
            else:
                print("   Likely JmjC reader protein (non-catalytic)")
        else:
            print("   Not a typical JmjC domain protein")
    
    # Save results if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(all_results, f, indent=2)
        if verbose:
            print(f"\nResults saved to {args.output}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())