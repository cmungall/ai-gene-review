#!/usr/bin/env python3
"""Extract protein sequence from UniProt file."""

import sys
from pathlib import Path

def extract_sequence_from_uniprot(uniprot_file):
    """Extract the protein sequence from a UniProt text file."""
    with open(uniprot_file, 'r') as f:
        lines = f.readlines()
    
    sequence_lines = []
    in_sequence = False
    
    for line in lines:
        if line.startswith('SQ'):
            in_sequence = True
            continue
        if in_sequence:
            if line.startswith('//'):
                break
            # Remove line numbers and spaces from sequence lines
            seq_part = line.strip()
            if seq_part:
                # Split by spaces and join to get clean sequence
                parts = seq_part.split()
                # Filter out line numbers (first element is usually the line number)
                seq_parts = [p for p in parts if p.replace(' ', '').isalpha()]
                sequence_lines.extend(seq_parts)
    
    return ''.join(sequence_lines)

def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else "../CG6051-uniprot.txt"
    input_path = Path(input_file)
    
    if not input_path.exists():
        print(f"Error: Input file {input_path} not found")
        sys.exit(1)
    
    sequence = extract_sequence_from_uniprot(input_path)
    
    # Save to FASTA format
    output_file = "CG6051.fasta"
    with open(output_file, 'w') as f:
        f.write(">Q9VB70|LST2_DROME Lateral signaling target protein 2 homolog\n")
        # Write sequence in 60-character lines
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
    
    print(f"Extracted sequence ({len(sequence)} aa) saved to {output_file}")
    print(f"First 50 aa: {sequence[:50]}")
    print(f"Last 50 aa: {sequence[-50:]}")

if __name__ == "__main__":
    main()