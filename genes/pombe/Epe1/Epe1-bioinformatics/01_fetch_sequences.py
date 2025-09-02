#!/usr/bin/env python3
"""
Fetch protein sequences for Epe1 and related proteins for analysis.
This script fetches:
1. S. pombe Epe1 (O94603)
2. Known active JmjC demethylases for comparison
3. Homologs from other species
"""

import requests
import json
from pathlib import Path
import sys

def fetch_uniprot_sequence(uniprot_id):
    """Fetch protein sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to fetch {uniprot_id}: {response.status_code}")
        return None

def fetch_uniprot_json(uniprot_id):
    """Fetch full UniProt entry as JSON."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to fetch JSON for {uniprot_id}: {response.status_code}")
        return None

def main():
    # Create data directory
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)
    
    # Primary target: S. pombe Epe1
    print("Fetching S. pombe Epe1 sequence...")
    epe1_fasta = fetch_uniprot_sequence("O94603")
    if epe1_fasta:
        with open(data_dir / "epe1_spombe.fasta", "w") as f:
            f.write(epe1_fasta)
        print("✓ Saved Epe1 sequence")
    
    # Fetch full UniProt entry for detailed analysis
    epe1_json = fetch_uniprot_json("O94603")
    if epe1_json:
        with open(data_dir / "epe1_uniprot.json", "w") as f:
            json.dump(epe1_json, f, indent=2)
        print("✓ Saved Epe1 UniProt data")
    
    # Known active JmjC demethylases for comparison
    active_demethylases = {
        "P84027": "KDM4A_HUMAN",  # JMJD2A - H3K9me3/H3K36me3 demethylase
        "Q9Y2K7": "KDM2A_HUMAN",  # FBXL11 - H3K36me2 demethylase
        "Q6ZMT4": "KDM5C_HUMAN",  # JARID1C - H3K4me3/me2 demethylase
        "Q92833": "KDM3A_HUMAN",  # JMJD1A - H3K9me2/me1 demethylase
        "P41229": "KDM5B_HUMAN",  # JARID1B - H3K4me3/me2 demethylase
    }
    
    print("\nFetching active JmjC demethylases for comparison...")
    for uniprot_id, name in active_demethylases.items():
        fasta = fetch_uniprot_sequence(uniprot_id)
        if fasta:
            with open(data_dir / f"{name.lower()}.fasta", "w") as f:
                f.write(fasta)
            print(f"✓ Saved {name}")
    
    # Try to find Epe1 homologs using UniProt search
    print("\nSearching for Epe1 homologs...")
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": "gene:epe1 AND taxonomy:fungi",
        "format": "json",
        "size": 10
    }
    
    response = requests.get(search_url, params=params)
    if response.status_code == 200:
        results = response.json()
        if results.get("results"):
            homologs_file = data_dir / "epe1_homologs.fasta"
            with open(homologs_file, "w") as f:
                for entry in results["results"]:
                    accession = entry["primaryAccession"]
                    organism = entry.get("organism", {}).get("scientificName", "Unknown")
                    if accession != "O94603":  # Skip our primary target
                        fasta = fetch_uniprot_sequence(accession)
                        if fasta:
                            f.write(fasta)
                            print(f"✓ Added homolog from {organism} ({accession})")
    
    print("\nSequence fetching complete!")
    print(f"Data saved in: {data_dir.absolute()}")

if __name__ == "__main__":
    main()