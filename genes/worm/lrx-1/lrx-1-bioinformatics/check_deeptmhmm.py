#!/usr/bin/env python3
"""
DeepTMHMM data preparation and result parsing for LRX-1 transmembrane topology analysis

This script prepares data for DeepTMHMM analysis and can parse results if provided.
DeepTMHMM is the most accurate method for predicting transmembrane helices and topology.
"""

import json
from typing import Dict, List, Optional
from pathlib import Path
import sys

def prepare_fasta() -> str:
    """Prepare FASTA format for DeepTMHMM submission"""
    
    # Load from FASTA file
    fasta_path = Path(__file__).parent.parent / "lrx-1.fasta"
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found at {fasta_path}")
    
    with open(fasta_path, 'r') as f:
        fasta = f.read()
    
    return fasta

def parse_deeptmhmm_output(output_file: str) -> Dict:
    """
    Parse DeepTMHMM output file if available
    
    DeepTMHMM output format includes:
    - Topology prediction (SP, TM, Globular, etc.)
    - Per-residue predictions
    - Probability scores
    """
    
    results = {
        "file": output_file,
        "topology": None,
        "signal_peptide": None,
        "tm_helices": [],
        "probabilities": {}
    }
    
    try:
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        # Parse the output based on DeepTMHMM format
        # This is a template for when actual output is available
        for line in lines:
            line = line.strip()
            if line.startswith('#'):
                continue
            
            # Parse topology line
            if 'Topology:' in line:
                results["topology"] = line.split('Topology:')[1].strip()
            
            # Parse signal peptide
            elif 'Signal peptide:' in line:
                results["signal_peptide"] = line.split('Signal peptide:')[1].strip()
            
            # Parse TM regions
            elif 'TM helix:' in line:
                tm_info = line.split('TM helix:')[1].strip()
                results["tm_helices"].append(tm_info)
        
    except FileNotFoundError:
        results["error"] = f"Output file not found: {output_file}"
    except Exception as e:
        results["error"] = f"Error parsing file: {str(e)}"
    
    return results

def analyze_topology_implications(topology: Optional[str]) -> Dict:
    """
    Analyze the biological implications of the topology prediction
    
    Returns objective analysis without predetermined conclusions
    """
    
    implications = {
        "SP+Globular": {
            "localization": "Likely secreted or extracellular",
            "characteristics": "Signal peptide followed by globular domain",
            "functional_implication": "Protein is processed through ER and likely secreted"
        },
        "SP+TM": {
            "localization": "Membrane-anchored with extracellular domain",
            "characteristics": "Signal peptide followed by transmembrane helix",
            "functional_implication": "Type I membrane protein topology"
        },
        "TM": {
            "localization": "Membrane-embedded",
            "characteristics": "Contains transmembrane helices without signal peptide",
            "functional_implication": "Integral membrane protein"
        },
        "Globular": {
            "localization": "Cytoplasmic or nuclear",
            "characteristics": "No signal peptide or TM helices",
            "functional_implication": "Soluble intracellular protein"
        }
    }
    
    if topology and topology in implications:
        return implications[topology]
    else:
        return {"note": "Topology prediction needed for functional analysis"}

def main():
    """Main function for DeepTMHMM analysis"""
    
    print("=== DeepTMHMM Data Preparation for LRX-1 ===\n")
    
    # Prepare FASTA
    try:
        fasta = prepare_fasta()
        fasta_output = "lrx1_for_deeptmhmm.fasta"
        with open(fasta_output, "w") as f:
            f.write(fasta)
        print(f"✓ FASTA file prepared: {fasta_output}")
        print(f"  Sequence length: {len([c for c in fasta if c.isalpha()])} aa")
    except Exception as e:
        print(f"✗ Error preparing FASTA: {e}")
        return
    
    # Check if there's an output file to parse
    if len(sys.argv) > 1:
        output_file = sys.argv[1]
        print(f"\n=== Parsing DeepTMHMM Output ===")
        results = parse_deeptmhmm_output(output_file)
        
        if "error" in results:
            print(f"✗ {results['error']}")
        else:
            print(f"✓ Parsed results from: {output_file}")
            print(f"  Topology: {results.get('topology', 'Not found')}")
            print(f"  Signal peptide: {results.get('signal_peptide', 'Not found')}")
            print(f"  TM helices: {len(results.get('tm_helices', []))}")
            
            # Analyze implications
            if results.get('topology'):
                implications = analyze_topology_implications(results['topology'])
                print("\n=== Biological Implications ===")
                for key, value in implications.items():
                    print(f"  {key}: {value}")
            
            # Save parsed results
            with open("deeptmhmm_parsed_results.json", "w") as f:
                json.dump(results, f, indent=2)
            print("\n✓ Results saved to: deeptmhmm_parsed_results.json")
    else:
        print("\n=== Next Steps ===")
        print("1. Submit the FASTA file to DeepTMHMM for analysis")
        print("2. Save the output to a file")
        print("3. Run: python check_deeptmhmm.py <output_file>")
        print("   to parse and analyze the results")
    
    print("\n=== Analysis Information ===")
    print("DeepTMHMM predicts:")
    print("- Signal peptides (SP)")
    print("- Transmembrane helices (TM)")
    print("- Overall topology classification")
    print("- Position-specific probabilities")
    print("\nThe tool will objectively analyze whatever topology is predicted.")

if __name__ == "__main__":
    main()