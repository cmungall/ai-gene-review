#!/usr/bin/env python3
"""
Fetch RBFOX family protein sequences from UniProt
"""

import requests
from pathlib import Path

def fetch_uniprot_fasta(uniprot_id):
    """Fetch FASTA sequence from UniProt API"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Failed to fetch {uniprot_id}: {response.status_code}")

def main():
    # RBFOX family UniProt IDs
    rbfox_proteins = {
        "RBFOX1": "Q9NWB1",  # Human RBFOX1
        "RBFOX2": "O43251",  # Human RBFOX2  
        "RBFOX3": "A6NFN3"   # Human RBFOX3 (already have this)
    }
    
    output_dir = Path(__file__).parent
    
    for protein_name, uniprot_id in rbfox_proteins.items():
        fasta_file = output_dir / f"{protein_name}.fasta"
        
        # Skip RBFOX3 if it already exists in parent directory
        if protein_name == "RBFOX3":
            parent_fasta = output_dir.parent.parent / "RBFOX3.fasta"
            if parent_fasta.exists():
                print(f"✓ {protein_name} already exists at {parent_fasta}")
                continue
        
        print(f"Fetching {protein_name} ({uniprot_id})...")
        try:
            fasta_content = fetch_uniprot_fasta(uniprot_id)
            
            with open(fasta_file, 'w') as f:
                f.write(fasta_content)
            
            print(f"✓ Saved {protein_name} to {fasta_file}")
            
            # Print first line and sequence length for verification
            lines = fasta_content.strip().split('\n')
            seq_length = len(''.join(lines[1:]))
            print(f"  Header: {lines[0]}")
            print(f"  Length: {seq_length} aa")
            
        except Exception as e:
            print(f"✗ Error fetching {protein_name}: {e}")
    
    print("\nDone! FASTA files ready for analysis.")

if __name__ == "__main__":
    main()