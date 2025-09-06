#!/usr/bin/env python3
"""
Fetch and analyze human LST2/ZFYVE28 protein for comparison with Drosophila CG6051.
"""

import sys
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
import json

def fetch_uniprot_sequence(uniprot_id: str) -> tuple:
    """Fetch protein sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        fasta_io = StringIO(response.text)
        record = next(SeqIO.parse(fasta_io, "fasta"))
        return str(record.id), str(record.seq), str(record.description)
    else:
        print(f"Error fetching {uniprot_id}: {response.status_code}")
        return None, None, None

def save_fasta(seq_id, sequence, description, filename):
    """Save sequence to FASTA file."""
    with open(filename, 'w') as f:
        f.write(f">{seq_id} {description}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

def main():
    # Fetch human LST2/ZFYVE28
    human_id = "Q9HCC9"
    print(f"Fetching human LST2/ZFYVE28 ({human_id})...")
    
    seq_id, sequence, description = fetch_uniprot_sequence(human_id)
    
    if sequence:
        # Save to file
        save_fasta(seq_id, sequence, description, "human_LST2.fasta")
        
        print(f"Retrieved human LST2:")
        print(f"  ID: {seq_id}")
        print(f"  Length: {len(sequence)} aa")
        print(f"  Description: {description}")
        
        # Extract key regions based on UniProt annotations
        # TOS motif should be around position 401-405
        tos_region = sequence[395:410] if len(sequence) > 410 else "N/A"
        print(f"\n  TOS motif region (396-410): {tos_region}")
        
        # Check for FDIDI motif
        if "FDIDI" in sequence:
            pos = sequence.find("FDIDI") + 1
            print(f"  FDIDI motif found at position {pos}")
            context = sequence[max(0, pos-6):min(len(sequence), pos+9)]
            print(f"  Context: {context}")
        
        # Check FYVE domain region (should be around 815-875)
        if len(sequence) > 875:
            fyve_region = sequence[814:875]
            cys_count = fyve_region.count('C')
            print(f"\n  FYVE domain region (815-875): {len(fyve_region)} aa")
            print(f"  Cysteines in FYVE: {cys_count}")
            print(f"  First 30 aa of FYVE: {fyve_region[:30]}")
        
        print(f"\nSequence saved to human_LST2.fasta")
    else:
        print("Failed to fetch human LST2 sequence")

if __name__ == "__main__":
    main()