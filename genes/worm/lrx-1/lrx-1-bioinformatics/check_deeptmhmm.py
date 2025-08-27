#!/usr/bin/env python3
"""
Check DeepTMHMM predictions for LRX-1 transmembrane topology
DeepTMHMM is the most accurate method for predicting transmembrane helices and topology
"""

import json
from typing import Dict

def prepare_fasta() -> str:
    """Prepare FASTA format for DeepTMHMM submission"""
    
    sequence = """MAWLTSIFFILLAVQPVLPQDLYGTATQQQPYPYVQPSASSGSGGYVPNPQSSIHTVQQP
YPNIDVVEPDVDSVDIYETEEPQFKVVNPVFPLGGSGIVEEPGTIPPPMPQTQAPEKPDNS
YAINYCDKREFPDDVLAQYGLERIDYFVYNTSCSHVFFQCSIGQTFPLACMSEDQAFDKS
TENCNHKNAIKFCPEYDHVMHCTIKDTCTENEFACCAMPQSCIHVSKRCDGHPDCADGED
ENNCPSCARDDEFACVKSEHCIPANKRCDGVADDCEDGSNLDEIGCSKNTTCIGKFVCGTS
RGGVSCVDLDMHCDGKKDCLNGEDEMNCQEGRQKYLLCENQKQSVTRLQWCNGETDCAD
GSDEKYCY""".replace("\n", "").replace(" ", "")
    
    fasta = ">Q22179|LRX-1|C.elegans\n"
    # Add sequence in 60-character lines
    for i in range(0, len(sequence), 60):
        fasta += sequence[i:i+60] + "\n"
    
    return fasta

def interpret_deeptmhmm_results() -> Dict:
    """
    Interpret DeepTMHMM results for LRX-1
    
    Note: This function documents what to look for in DeepTMHMM output.
    Actual results need to be obtained from: https://dtu.biolib.com/DeepTMHMM
    """
    
    interpretation_guide = {
        "what_to_check": {
            "signal_peptide": "Look for 'SP' (signal peptide) prediction in first ~20-30 residues",
            "tm_helices": "Check for 'TM' (transmembrane) regions - especially after signal peptide",
            "topology": "Check overall prediction: 'Globular', 'SP+Globular', 'TM', or 'SP+TM'",
            "confidence": "DeepTMHMM provides probability scores for each prediction"
        },
        "expected_for_membrane_protein": {
            "pattern": "SP (signal) followed by TM (transmembrane) region",
            "topology": "Type Ia: N-term outside, C-term inside",
            "example": "Positions 1-20: SP, Positions 21-45: TM"
        },
        "expected_for_secreted_protein": {
            "pattern": "SP (signal) followed by non-TM regions",
            "topology": "SP+Globular",
            "example": "Positions 1-20: SP, Rest: Outside/Globular"
        },
        "lrx1_hypothesis": {
            "prediction": "SP+Globular (secreted)",
            "reasoning": [
                "Has signal peptide (1-19)",
                "No hydrophobic region after signal",
                "Cysteine-rich suggesting extracellular",
                "AlphaFold shows no TM helix"
            ]
        }
    }
    
    return interpretation_guide

def generate_submission_instructions() -> str:
    """Generate instructions for DeepTMHMM submission"""
    
    fasta = prepare_fasta()
    
    instructions = f"""
=== DeepTMHMM Analysis Instructions for LRX-1 ===

1. Go to: https://dtu.biolib.com/DeepTMHMM

2. Paste this FASTA sequence:
{fasta}

3. Click "Run" and wait for results (~30 seconds)

4. Check the results for:
   - Signal peptide (SP) prediction
   - Transmembrane helix (TM) predictions
   - Overall topology classification
   - Probability plot showing confidence

5. Key questions to answer:
   - Is there a TM helix after the signal peptide?
   - What is the overall topology prediction?
   - Does it predict "SP+Globular" (secreted) or "SP+TM" (membrane)?

6. Expected result based on our analysis:
   - Signal peptide: YES (positions ~1-19)
   - TM helices: NONE
   - Topology: SP+Globular (secreted protein)
   
7. If DeepTMHMM predicts membrane protein:
   - Check the probability scores
   - Look at the specific region predicted as TM
   - Compare hydrophobicity of that region
"""
    
    return instructions

def main():
    """Generate DeepTMHMM analysis setup"""
    
    print("=== DeepTMHMM Analysis Setup for LRX-1 ===\n")
    
    # Generate FASTA
    fasta = prepare_fasta()
    with open("lrx1_for_deeptmhmm.fasta", "w") as f:
        f.write(fasta)
    print("✓ FASTA file created: lrx1_for_deeptmhmm.fasta")
    
    # Show interpretation guide
    guide = interpret_deeptmhmm_results()
    print("\n=== What to Look For ===")
    print(json.dumps(guide, indent=2))
    
    # Print instructions
    instructions = generate_submission_instructions()
    print(instructions)
    
    print("\n=== Summary ===")
    print("DeepTMHMM is the gold standard for membrane topology prediction.")
    print("It uses deep learning trained on experimentally validated proteins.")
    print("This analysis will definitively determine if LRX-1 is membrane-bound or secreted.")

if __name__ == "__main__":
    main()