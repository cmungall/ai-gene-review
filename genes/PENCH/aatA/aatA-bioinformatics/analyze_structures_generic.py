#!/usr/bin/env python3
"""
Generic structural analysis tool for PDB structures
Analyzes protein structures for specified residues and cleavage sites
"""

import requests
import pandas as pd
from Bio import PDB
from pathlib import Path
import json
from datetime import datetime
import warnings
import argparse
import sys
import yaml
warnings.filterwarnings('ignore')

def download_pdb(pdb_id, output_dir="."):
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_file = Path(output_dir) / f"{pdb_id}.pdb"
        with open(pdb_file, 'w') as f:
            f.write(response.text)
        return pdb_file
    else:
        raise Exception(f"Failed to download {pdb_id}")

def analyze_structure(pdb_id, pdb_file, residues_of_interest=None):
    """
    Analyze PDB structure for specified residues.
    
    Args:
        pdb_id: PDB identifier
        pdb_file: Path to PDB file
        residues_of_interest: Dict of {position: expected_residue} to check
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    
    results = {
        "pdb_id": pdb_id,
        "chains": [],
        "residues_analyzed": {},
        "num_chains": 0,
        "total_residues": 0
    }
    
    # Analyze each chain
    for model in structure:
        for chain in model:
            chain_residues = list(chain.get_residues())
            chain_info = {
                "chain_id": chain.get_id(),
                "num_residues": len(chain_residues),
                "start_residue": chain_residues[0].get_id()[1] if chain_residues else None,
                "end_residue": chain_residues[-1].get_id()[1] if chain_residues else None
            }
            results["chains"].append(chain_info)
            results["total_residues"] += len(chain_residues)
            
            # Look for specified residues
            if residues_of_interest:
                for residue in chain:
                    res_id = residue.get_id()[1]
                    res_name = residue.get_resname()
                    
                    if res_id in residues_of_interest:
                        key = f"position_{res_id}"
                        if key not in results["residues_analyzed"]:
                            results["residues_analyzed"][key] = []
                        
                        results["residues_analyzed"][key].append({
                            "chain": chain.get_id(),
                            "residue": res_name,
                            "expected": residues_of_interest[res_id],
                            "matches": res_name == residues_of_interest[res_id] or 
                                      (isinstance(residues_of_interest[res_id], list) and 
                                       res_name in residues_of_interest[res_id])
                        })
    
    results["num_chains"] = len(results["chains"])
    
    # Determine if this might be a cleaved form
    if results["num_chains"] > 1:
        results["form"] = "multi-chain (possibly cleaved)"
    else:
        results["form"] = "single-chain"
    
    return results

def analyze_modified_residues(pdb_file):
    """Check for modified residues in PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    pdb_id = Path(pdb_file).stem
    structure = parser.get_structure(pdb_id, pdb_file)
    
    modified_residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                # Check for non-standard amino acids
                standard_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                              'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                              'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
                
                if res_name not in standard_aa and not res_name.startswith('HOH'):
                    modified_residues.append({
                        'chain': chain.get_id(),
                        'position': residue.get_id()[1],
                        'residue': res_name
                    })
    
    return modified_residues

def generate_report(all_results, analysis_config):
    """Generate markdown report of findings."""
    report = []
    report.append(f"# Structural Analysis Report")
    report.append(f"\n**Analysis Date:** {datetime.now().strftime('%Y-%m-%d')}")
    
    if 'title' in analysis_config:
        report.append(f"\n## {analysis_config['title']}")
    
    if 'description' in analysis_config:
        report.append(f"\n{analysis_config['description']}")
    
    report.append("\n## PDB Structures Analyzed\n")
    
    for result in all_results:
        pdb_info = analysis_config.get('pdb_structures', {}).get(result['pdb_id'], {})
        description = pdb_info.get('description', 'No description provided')
        
        report.append(f"### {result['pdb_id']}: {description}")
        report.append(f"- **Form:** {result['form']}")
        report.append(f"- **Number of chains:** {result['num_chains']}")
        report.append(f"- **Total residues:** {result['total_residues']}")
        
        if result.get('modified_residues'):
            report.append(f"- **Modified residues:** {len(result['modified_residues'])}")
            for mod in result['modified_residues'][:3]:  # Show first 3
                report.append(f"  - Chain {mod['chain']}, Position {mod['position']}: {mod['residue']}")
        
        if result['residues_analyzed']:
            report.append("\n**Key Residues:**")
            for pos, residues in result['residues_analyzed'].items():
                for res_info in residues:
                    status = "✓" if res_info['matches'] else "✗"
                    report.append(f"- {pos}: {res_info['residue']} in chain {res_info['chain']} "
                                f"(expected: {res_info['expected']}) {status}")
    
    report.append("\n## Analysis Summary\n")
    
    # Check for consistency across structures
    if 'expected_findings' in analysis_config:
        report.append("### Expected Findings:")
        for finding in analysis_config['expected_findings']:
            report.append(f"- {finding}")
    
    # Summarize key residue findings
    residue_summary = {}
    for result in all_results:
        for pos, residues in result['residues_analyzed'].items():
            if pos not in residue_summary:
                residue_summary[pos] = []
            for res_info in residues:
                residue_summary[pos].append({
                    'pdb': result['pdb_id'],
                    'residue': res_info['residue'],
                    'matches': res_info['matches']
                })
    
    if residue_summary:
        report.append("\n### Residue Conservation:")
        for pos, observations in residue_summary.items():
            conserved = all(obs['matches'] for obs in observations)
            status = "✓ Conserved" if conserved else "⚠ Variable"
            report.append(f"- {pos}: {status}")
            for obs in observations:
                report.append(f"  - {obs['pdb']}: {obs['residue']}")
    
    # Check for cleavage patterns
    single_chain = sum(1 for r in all_results if r['form'] == 'single-chain')
    multi_chain = sum(1 for r in all_results if 'multi-chain' in r['form'])
    
    if single_chain > 0 and multi_chain > 0:
        report.append("\n### Structural Forms:")
        report.append(f"- Single-chain structures: {single_chain}")
        report.append(f"- Multi-chain structures: {multi_chain}")
        report.append("- This suggests possible proteolytic processing")
    
    if 'conclusions' in analysis_config:
        report.append("\n## Conclusions\n")
        for conclusion in analysis_config['conclusions']:
            report.append(f"- {conclusion}")
    
    if 'references' in analysis_config:
        report.append("\n## References\n")
        for ref in analysis_config['references']:
            report.append(f"- {ref}")
    
    return "\n".join(report)

def load_config(config_file):
    """Load analysis configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(
        description='Analyze PDB structures for specified residues and features',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Configuration file format (YAML):
---
title: "Analysis Title"
description: "Brief description"

pdb_structures:
  2X1C:
    description: "Mutant form"
  2X1D:
    description: "Wild-type"

residues_of_interest:
  103: CYS  # or [CYS, ALA] for multiple possibilities
  102: GLY
  309: SER

expected_findings:
  - "Finding 1"
  - "Finding 2"

conclusions:
  - "Conclusion 1"

references:
  - "Reference 1"

Examples:
  # Analyze with configuration file
  %(prog)s --config analysis.yaml
  
  # Analyze specific PDB IDs with residue checks
  %(prog)s --pdb 2X1C 2X1D --residues 103:CYS 102:GLY
  
  # Save results to specific directory
  %(prog)s --config analysis.yaml --output results/
        """
    )
    
    parser.add_argument('--config', '-c', help='Configuration file (YAML)')
    parser.add_argument('--pdb', nargs='+', help='PDB IDs to analyze')
    parser.add_argument('--residues', nargs='+', 
                       help='Residues to check (format: position:expected, e.g., 103:CYS)')
    parser.add_argument('--output', '-o', default='.', 
                       help='Output directory (default: current directory)')
    parser.add_argument('--report', '-r', default='structural_analysis_report.md',
                       help='Report filename (default: structural_analysis_report.md)')
    
    args = parser.parse_args()
    
    # Load configuration
    if args.config:
        config = load_config(args.config)
        pdb_ids = list(config.get('pdb_structures', {}).keys())
        residues_of_interest = config.get('residues_of_interest', {})
    else:
        if not args.pdb:
            print("Error: Either --config or --pdb must be specified", file=sys.stderr)
            sys.exit(1)
        
        pdb_ids = args.pdb
        residues_of_interest = {}
        
        # Parse residue specifications
        if args.residues:
            for res_spec in args.residues:
                if ':' in res_spec:
                    pos, expected = res_spec.split(':', 1)
                    try:
                        residues_of_interest[int(pos)] = expected
                    except ValueError:
                        print(f"Warning: Invalid residue specification: {res_spec}", 
                              file=sys.stderr)
        
        # Create minimal config
        config = {
            'pdb_structures': {pdb_id: {'description': f'Structure {pdb_id}'} 
                             for pdb_id in pdb_ids}
        }
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Starting structural analysis of {len(pdb_ids)} structures...")
    
    all_results = []
    
    for pdb_id in pdb_ids:
        print(f"\nAnalyzing {pdb_id}...")
        
        # Download PDB file
        print(f"  Downloading PDB file...")
        try:
            pdb_file = download_pdb(pdb_id, output_dir)
        except Exception as e:
            print(f"  Error downloading {pdb_id}: {e}", file=sys.stderr)
            continue
        
        # Analyze structure
        print(f"  Analyzing structure...")
        results = analyze_structure(pdb_id, pdb_file, residues_of_interest)
        
        # Check for modified residues
        print(f"  Checking for modified residues...")
        results['modified_residues'] = analyze_modified_residues(pdb_file)
        
        all_results.append(results)
        
        # Save intermediate results
        with open(output_dir / f"{pdb_id}_analysis.json", 'w') as f:
            json.dump(results, f, indent=2)
    
    # Generate report
    print("\nGenerating report...")
    report = generate_report(all_results, config)
    
    # Save report
    report_file = output_dir / args.report
    with open(report_file, 'w') as f:
        f.write(report)
    
    print(f"\nAnalysis complete! Results saved to {report_file}")
    
    # Print summary to console
    print("\n" + "="*50)
    print("ANALYSIS SUMMARY:")
    print("="*50)
    for result in all_results:
        print(f"\n{result['pdb_id']}: {result['form']}")
        print(f"  Chains: {result['num_chains']}, Residues: {result['total_residues']}")
        if result['residues_analyzed']:
            matches = sum(1 for pos_data in result['residues_analyzed'].values() 
                         for res in pos_data if res['matches'])
            total = sum(len(pos_data) for pos_data in result['residues_analyzed'].values())
            print(f"  Key residues: {matches}/{total} match expected")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())