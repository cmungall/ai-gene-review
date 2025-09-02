#!/usr/bin/env python3
"""
Parse DeepTMHMM results for LRX-1
Save DeepTMHMM output to 'deeptmhmm_results.txt' and run this script
"""

import json
from pathlib import Path
from typing import Dict

def parse_deeptmhmm_results(filepath: str = "deeptmhmm_results.txt") -> Dict:
    """
    Parse DeepTMHMM output file
    
    Expected format example:
    >Q22179
    M	INSIDE	0.00001	0.00002	0.99997
    A	INSIDE	0.00001	0.00001	0.99998
    W	SIGNAL	0.00001	0.99997	0.00002
    ...
    
    Or summary format:
    >Q22179 SP+GLOBULAR
    """
    
    if not Path(filepath).exists():
        return {"error": f"File {filepath} not found. Please save DeepTMHMM results to this file."}
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    result = {
        "protein_id": None,
        "topology_prediction": None,
        "signal_peptide": {"present": False, "start": None, "end": None},
        "tm_helices": [],
        "regions": [],
        "per_residue": []
    }
    
    current_region = None
    position = 0
    
    for line in lines:
        line = line.strip()
        
        # Header line
        if line.startswith('>'):
            parts = line.split()
            result["protein_id"] = parts[0][1:]
            if len(parts) > 1:
                result["topology_prediction"] = parts[1]
            continue
        
        # Skip empty lines
        if not line:
            continue
        
        # Parse per-residue predictions
        parts = line.split('\t')
        if len(parts) >= 3:
            position += 1
            parts[0]
            prediction = parts[1]
            
            # Track regions
            if prediction != current_region:
                if current_region:
                    # Close previous region
                    result["regions"][-1]["end"] = position - 1
                
                # Start new region
                current_region = prediction
                result["regions"].append({
                    "type": prediction,
                    "start": position,
                    "end": None
                })
                
                # Track specific features
                if prediction == "SIGNAL":
                    result["signal_peptide"]["present"] = True
                    result["signal_peptide"]["start"] = position
                elif prediction in ["TRANSMEM", "TM", "MEMBRANE"]:
                    result["tm_helices"].append({"start": position, "end": None})
            
            # Update signal peptide end
            if prediction == "SIGNAL":
                result["signal_peptide"]["end"] = position
            
            # Update TM helix end
            if prediction in ["TRANSMEM", "TM", "MEMBRANE"] and result["tm_helices"]:
                result["tm_helices"][-1]["end"] = position
    
    # Close last region
    if result["regions"] and result["regions"][-1]["end"] is None:
        result["regions"][-1]["end"] = position
    
    # Determine topology if not explicitly given
    if not result["topology_prediction"]:
        if result["signal_peptide"]["present"] and not result["tm_helices"]:
            result["topology_prediction"] = "SP+GLOBULAR"
        elif result["signal_peptide"]["present"] and result["tm_helices"]:
            result["topology_prediction"] = "SP+TM"
        elif result["tm_helices"]:
            result["topology_prediction"] = "TM"
        else:
            result["topology_prediction"] = "GLOBULAR"
    
    return result

def interpret_results(parsed: Dict) -> str:
    """Interpret parsed DeepTMHMM results"""
    
    interpretation = []
    
    interpretation.append(f"Protein: {parsed.get('protein_id', 'Unknown')}")
    interpretation.append(f"Overall Topology: {parsed.get('topology_prediction', 'Unknown')}")
    
    if parsed["signal_peptide"]["present"]:
        sp = parsed["signal_peptide"]
        interpretation.append(f"Signal Peptide: YES (positions {sp['start']}-{sp['end']})")
    else:
        interpretation.append("Signal Peptide: NO")
    
    if parsed["tm_helices"]:
        interpretation.append(f"Transmembrane Helices: {len(parsed['tm_helices'])}")
        for i, tm in enumerate(parsed["tm_helices"], 1):
            interpretation.append(f"  TM{i}: positions {tm['start']}-{tm['end']}")
    else:
        interpretation.append("Transmembrane Helices: NONE")
    
    # Objective Analysis
    interpretation.append("\nTopology Analysis:")
    if parsed["topology_prediction"] == "SP+GLOBULAR":
        interpretation.append("• Predicted topology: Signal peptide followed by globular domain")
        interpretation.append("• Characteristic of: Secreted or extracellular proteins")
        interpretation.append("• Processing: Likely through ER/Golgi secretory pathway")
    elif parsed["topology_prediction"] == "SP+TM":
        interpretation.append("• Predicted topology: Signal peptide followed by transmembrane helix")
        interpretation.append("• Characteristic of: Type I membrane proteins")
        interpretation.append("• Localization: Membrane-anchored with extracellular domain")
    elif parsed["topology_prediction"] == "TM":
        interpretation.append("• Predicted topology: Transmembrane protein without signal peptide")
        interpretation.append("• Characteristic of: Multi-pass membrane proteins")
        interpretation.append("• Localization: Integral membrane protein")
    elif parsed["topology_prediction"] == "GLOBULAR":
        interpretation.append("• Predicted topology: Globular protein without signal peptide")
        interpretation.append("• Characteristic of: Cytoplasmic or nuclear proteins")
        interpretation.append("• Localization: Intracellular")
    else:
        interpretation.append(f"• Topology prediction: {parsed['topology_prediction']}")
        interpretation.append("• Further analysis needed to determine localization")
    
    # Note about comparing with other annotations
    interpretation.append("\nNote: Compare this prediction with:")
    interpretation.append("• Experimental evidence from literature")
    interpretation.append("• Other prediction tools (SignalP, TMHMM, Phobius)")
    interpretation.append("• Functional studies and cellular localization data")
    
    return "\n".join(interpretation)

def main():
    """Parse and interpret DeepTMHMM results"""
    
    print("=== DeepTMHMM Results Parser for LRX-1 ===\n")
    
    # Try to parse results file
    parsed = parse_deeptmhmm_results()
    
    if "error" in parsed:
        print(parsed["error"])
        print("\nTo use this script:")
        print("1. Go to https://dtu.biolib.com/DeepTMHMM")
        print("2. Submit the sequence from lrx1_for_deeptmhmm.fasta")
        print("3. Save the results to 'deeptmhmm_results.txt'")
        print("4. Run this script again")
        return
    
    # Output parsed results
    print("Parsed Results:")
    print(json.dumps(parsed, indent=2))
    
    print("\n" + "="*50)
    print("Interpretation:")
    print("="*50)
    print(interpret_results(parsed))
    
    # Save interpreted results
    with open("deeptmhmm_interpretation.json", "w") as f:
        json.dump(parsed, f, indent=2)
    
    print("\n✓ Results saved to deeptmhmm_interpretation.json")

if __name__ == "__main__":
    main()