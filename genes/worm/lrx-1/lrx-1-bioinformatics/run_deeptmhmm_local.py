#!/usr/bin/env python3
"""
Run DeepTMHMM locally using BioLib
This actually runs the tool, not just instructions!
"""

import json
from typing import Dict, Optional

# Check if biolib is installed
try:
    import biolib
    BIOLIB_AVAILABLE = True
except ImportError:
    BIOLIB_AVAILABLE = False
    print("BioLib not installed. Install with: pip install pybiolib")

def run_deeptmhmm_analysis(sequence: str, sequence_id: str = "LRX1") -> Optional[Dict]:
    """
    Run DeepTMHMM on a protein sequence using BioLib
    
    Args:
        sequence: Protein sequence (amino acids)
        sequence_id: Identifier for the sequence
    
    Returns:
        Dictionary with prediction results or None if failed
    """
    
    if not BIOLIB_AVAILABLE:
        return {"error": "BioLib not installed. Run: pip install pybiolib"}
    
    try:
        # Load DeepTMHMM from BioLib
        print("Loading DeepTMHMM from BioLib...")
        deeptmhmm = biolib.load('DTU/DeepTMHMM')
        
        # Prepare input in FASTA format
        fasta_input = f">{sequence_id}\n{sequence}"
        
        # Run DeepTMHMM
        print("Running DeepTMHMM prediction...")
        # Write FASTA to temporary file for biolib
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_input)
            fasta_file = f.name
        
        try:
            # Run with file input
            result = deeptmhmm.cli(args=['--fasta', fasta_file])
            
            # Get output (use get_stdout() as per warning)
            if hasattr(result, 'get_stdout'):
                output = result.get_stdout()
            else:
                output = result.stdout
            
            # Parse the results
            parsed = parse_deeptmhmm_output(output)
            
        finally:
            # Clean up temp file
            import os
            if os.path.exists(fasta_file):
                os.unlink(fasta_file)
        
        return parsed
        
    except Exception as e:
        return {"error": f"Failed to run DeepTMHMM: {str(e)}"}

def parse_deeptmhmm_output(output: str) -> Dict:
    """
    Parse DeepTMHMM output
    
    Expected format varies but typically includes:
    - Topology prediction (GLOB, TM, SP+GLOB, SP+TM)
    - Per-residue predictions
    """
    
    lines = output.strip().split('\n')
    
    result = {
        "raw_output": output,
        "topology": None,
        "signal_peptide": None,
        "tm_helices": [],
        "regions": []
    }
    
    for line in lines:
        line = line.strip()
        
        # Look for topology prediction
        if 'GLOB' in line or 'TM' in line or 'SP' in line:
            if 'SP+GLOB' in line:
                result["topology"] = "SP+GLOBULAR (secreted)"
            elif 'SP+TM' in line:
                result["topology"] = "SP+TM (membrane)"
            elif 'GLOB' in line:
                result["topology"] = "GLOBULAR"
            elif 'TM' in line:
                result["topology"] = "TRANSMEMBRANE"
        
        # Look for signal peptide region
        if 'Signal peptide' in line or 'SP' in line:
            # Extract position info if available
            import re
            pos_match = re.search(r'(\d+)-(\d+)', line)
            if pos_match:
                result["signal_peptide"] = {
                    "start": int(pos_match.group(1)),
                    "end": int(pos_match.group(2))
                }
        
        # Look for TM regions
        if 'Transmembrane' in line or 'TM' in line:
            import re
            pos_match = re.search(r'(\d+)-(\d+)', line)
            if pos_match:
                result["tm_helices"].append({
                    "start": int(pos_match.group(1)),
                    "end": int(pos_match.group(2))
                })
    
    return result

def main():
    """Run DeepTMHMM on LRX-1 sequence"""
    
    # LRX-1 sequence
    sequence = """MAWLTSIFFILLAVQPVLPQDLYGTATQQQPYPYVQPSASSGSGGYVPNPQSSIHTVQQP
YPNIDVVEPDVDSVDIYETEEPQFKVVNPVFPLGGSGIVEEPGTIPPPMPQTQAPEKPDNS
YAINYCDKREFPDDVLAQYGLERIDYFVYNTSCSHVFFQCSIGQTFPLACMSEDQAFDKS
TENCNHKNAIKFCPEYDHVMHCTIKDTCTENEFACCAMPQSCIHVSKRCDGHPDCADGED
ENNCPSCARDDEFACVKSEHCIPANKRCDGVADDCEDGSNLDEIGCSKNTTCIGKFVCGTS
RGGVSCVDLDMHCDGKKDCLNGEDEMNCQEGRQKYLLCENQKQSVTRLQWCNGETDCAD
GSDEKYCY""".replace("\n", "").replace(" ", "")
    
    print("=== DeepTMHMM Analysis for LRX-1 (Q22179) ===\n")
    
    # Run analysis
    result = run_deeptmhmm_analysis(sequence, "Q22179_LRX1")
    
    if result and "error" in result:
        print(f"Error: {result['error']}")
        print("\nTo fix:")
        print("1. Install BioLib: pip install pybiolib")
        print("2. You may need to authenticate: biolib login")
        return
    
    # Display results
    print("\n=== Results ===")
    print(f"Topology Prediction: {result.get('topology', 'Unknown')}")
    
    if result.get("signal_peptide"):
        sp = result["signal_peptide"]
        print(f"Signal Peptide: YES (positions {sp['start']}-{sp['end']})")
    else:
        print("Signal Peptide: Not detected or not reported")
    
    if result.get("tm_helices"):
        print(f"Transmembrane Helices: {len(result['tm_helices'])}")
        for i, tm in enumerate(result["tm_helices"], 1):
            print(f"  TM{i}: positions {tm['start']}-{tm['end']}")
    else:
        print("Transmembrane Helices: NONE")
    
    # Save results
    with open("deeptmhmm_results.json", "w") as f:
        json.dump(result, f, indent=2)
    print("\n✓ Results saved to deeptmhmm_results.json")
    
    # Interpretation
    print("\n=== Interpretation ===")
    if result.get("topology") == "SP+GLOBULAR (secreted)":
        print("✓ LRX-1 is a SECRETED protein")
        print("✓ Has signal peptide but NO transmembrane helices")
        print("✓ UniProt's membrane annotation is INCORRECT")
    elif "TM" in str(result.get("topology", "")):
        print("✗ LRX-1 appears to be a MEMBRANE protein")
        print("✗ This contradicts our sequence analysis")
        print("✗ Check the raw output for details")

if __name__ == "__main__":
    main()