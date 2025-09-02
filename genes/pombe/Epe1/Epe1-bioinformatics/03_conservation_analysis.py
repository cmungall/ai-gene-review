#!/usr/bin/env python3
"""
Perform conservation analysis of Epe1 across species.
Focus on JmjC domain and critical residues.
"""

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import requests
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def search_homologs():
    """Search for Epe1 homologs across different species."""
    print("Searching for Epe1 homologs in different species...")
    
    # Use BLAST-like search via UniProt
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # Search for proteins similar to Epe1
    params = {
        "query": "(gene:epe1 OR name:\"Jmjc domain-containing histone demethylation protein\") AND reviewed:true",
        "format": "json",
        "size": 25
    }
    
    response = requests.get(search_url, params=params)
    homologs = []
    
    if response.status_code == 200:
        results = response.json()
        for entry in results.get("results", []):
            organism = entry.get("organism", {}).get("scientificName", "")
            accession = entry["primaryAccession"]
            gene = entry.get("genes", [{}])[0].get("geneName", {}).get("value", "")
            
            homologs.append({
                "accession": accession,
                "organism": organism,
                "gene": gene
            })
    
    # Also search for JmjC proteins in fission yeast relatives
    yeast_search = {
        "query": "jmjc AND taxonomy:\"schizosaccharomyces\"",
        "format": "json",
        "size": 10
    }
    
    response = requests.get(search_url, params=yeast_search)
    if response.status_code == 200:
        results = response.json()
        for entry in results.get("results", []):
            organism = entry.get("organism", {}).get("scientificName", "")
            accession = entry["primaryAccession"]
            gene = entry.get("genes", [{}])[0].get("geneName", {}).get("value", "")
            
            if accession not in [h["accession"] for h in homologs]:
                homologs.append({
                    "accession": accession,
                    "organism": organism,
                    "gene": gene
                })
    
    return homologs

def fetch_and_save_homologs(homologs):
    """Fetch sequences for homologs."""
    sequences = []
    data_dir = Path("data")
    
    for homolog in homologs:
        url = f"https://rest.uniprot.org/uniprotkb/{homolog['accession']}.fasta"
        response = requests.get(url)
        
        if response.status_code == 200:
            # Parse the FASTA
            lines = response.text.strip().split('\n')
            if len(lines) > 1:
                header = lines[0]
                sequence = ''.join(lines[1:])
                
                # Create simplified ID
                org_short = homolog['organism'].replace(' ', '_')[:20]
                gene = homolog['gene'] if homolog['gene'] else homolog['accession']
                seq_id = f"{gene}_{org_short}"
                
                record = SeqRecord(
                    Seq(sequence),
                    id=seq_id,
                    description=f"{homolog['accession']} {homolog['organism']}"
                )
                sequences.append(record)
                print(f"  Fetched: {gene} from {homolog['organism']}")
    
    # Save all sequences
    if sequences:
        output_file = data_dir / "epe1_homologs_extended.fasta"
        with open(output_file, "w") as f:
            SeqIO.write(sequences, f, "fasta")
        print(f"\nSaved {len(sequences)} homolog sequences to {output_file}")
    
    return sequences

def extract_jmjc_domains(sequences):
    """Extract JmjC domain regions from sequences."""
    jmjc_domains = []
    
    # For Epe1, we know the JmjC domain is approximately 243-402
    # We'll use pattern matching to find similar regions in homologs
    
    for seq_record in sequences:
        sequence = str(seq_record.seq)
        
        # Look for JmjC domain characteristics
        # Simplified approach: look for region with histidines and acidic residues
        best_score = 0
        best_start = 0
        window_size = 160  # Approximate JmjC domain size
        
        for i in range(len(sequence) - window_size):
            window = sequence[i:i+window_size]
            # Score based on JmjC characteristics
            h_count = window.count('H')
            de_count = window.count('D') + window.count('E')
            
            # Simple scoring
            score = h_count + de_count
            
            if score > best_score:
                best_score = score
                best_start = i
        
        # Extract the best matching region
        jmjc_seq = sequence[best_start:best_start+window_size]
        
        jmjc_record = SeqRecord(
            Seq(jmjc_seq),
            id=seq_record.id,
            description=f"JmjC domain {best_start}-{best_start+window_size}"
        )
        jmjc_domains.append(jmjc_record)
    
    return jmjc_domains

def analyze_conservation(alignment_file):
    """Analyze conservation from alignment."""
    # This would typically use an alignment file
    # For now, we'll do a simple analysis
    
    sequences = []
    with open(alignment_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append(str(record.seq))
    
    if not sequences:
        return None
    
    # Calculate conservation scores
    seq_len = len(sequences[0])
    conservation = []
    
    for pos in range(seq_len):
        column = [seq[pos] if pos < len(seq) else '-' for seq in sequences]
        # Simple conservation: fraction of most common residue
        if column:
            most_common = max(set(column), key=column.count)
            if most_common != '-':
                score = column.count(most_common) / len(column)
            else:
                score = 0
        else:
            score = 0
        conservation.append(score)
    
    return conservation

def identify_critical_positions():
    """Identify critical positions based on literature."""
    # Based on active JmjC demethylases structure
    critical_positions = {
        "Fe_binding_1": "HXD/HXE motif - First histidine",
        "Fe_binding_2": "HXD/HXE motif - Aspartate/Glutamate",
        "Fe_binding_3": "Additional histidine for Fe coordination",
        "aKG_binding": "Lysine or Serine for α-ketoglutarate"
    }
    
    return critical_positions

def main():
    print("=" * 60)
    print("Conservation Analysis of Epe1")
    print("=" * 60)
    
    data_dir = Path("data")
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    
    # Load Epe1 sequence
    with open(data_dir / "epe1_spombe.fasta") as f:
        epe1_record = next(SeqIO.parse(f, "fasta"))
    
    epe1_seq = str(epe1_record.seq)
    
    # Extract JmjC domain from Epe1 (positions 243-402 from UniProt)
    epe1_jmjc = epe1_seq[242:402]  # 0-indexed
    
    print("\n1. Epe1 JmjC Domain:")
    print("-" * 40)
    print(f"  Length: {len(epe1_jmjc)} aa")
    print(f"  Position: 243-402")
    
    # Search for homologs
    print("\n2. Searching for Homologs:")
    print("-" * 40)
    homologs = search_homologs()
    
    if homologs:
        print(f"  Found {len(homologs)} potential homologs")
        sequences = fetch_and_save_homologs(homologs[:10])  # Limit to 10 for analysis
    else:
        print("  No additional homologs found via API search")
        sequences = []
    
    # Load active demethylases for comparison
    print("\n3. Loading Active Demethylases for Comparison:")
    print("-" * 40)
    
    active_seqs = []
    for fasta_file in data_dir.glob("kdm*.fasta"):
        with open(fasta_file) as f:
            record = next(SeqIO.parse(f, "fasta"))
            active_seqs.append(record)
            print(f"  Loaded: {fasta_file.stem.upper()}")
    
    # Combine all sequences for analysis
    all_sequences = [epe1_record] + sequences + active_seqs
    
    # Extract JmjC domains
    print("\n4. Extracting JmjC Domains:")
    print("-" * 40)
    jmjc_domains = extract_jmjc_domains(all_sequences)
    
    # Save JmjC domains
    jmjc_file = results_dir / "jmjc_domains.fasta"
    with open(jmjc_file, "w") as f:
        SeqIO.write(jmjc_domains, f, "fasta")
    print(f"  Saved {len(jmjc_domains)} JmjC domains to {jmjc_file}")
    
    # Analyze specific positions in Epe1
    print("\n5. Critical Residue Analysis in Epe1 JmjC:")
    print("-" * 40)
    
    # Check for canonical Fe(II) binding residues
    print("\n  Checking for canonical catalytic residues:")
    
    # Position ~280 should have HXD motif
    pos_280_region = epe1_jmjc[35:45]  # Around position 280 in full sequence
    print(f"    Region around expected HXD (275-285): {pos_280_region}")
    
    # Check if we have HXD or HXE
    import re
    hxd_pattern = re.compile(r"H.[DE]")
    matches = list(hxd_pattern.finditer(epe1_jmjc))
    
    if matches:
        print(f"    Found {len(matches)} HX[DE] motifs:")
        for match in matches:
            abs_pos = match.start() + 243  # Convert to absolute position
            print(f"      - Position {abs_pos}: {match.group()}")
    else:
        print("    ⚠️  No canonical HXD/HXE motifs found!")
    
    # Count key residues
    h_count = epe1_jmjc.count('H')
    d_count = epe1_jmjc.count('D')
    e_count = epe1_jmjc.count('E')
    k_count = epe1_jmjc.count('K')
    
    print(f"\n  Residue composition in JmjC domain:")
    print(f"    Histidines (H): {h_count}")
    print(f"    Aspartates (D): {d_count}")
    print(f"    Glutamates (E): {e_count}")
    print(f"    Lysines (K): {k_count}")
    
    # Compare with typical active demethylase
    if active_seqs:
        # Use KDM4A as reference (well-characterized)
        kdm4a = next((s for s in active_seqs if "KDM4A" in s.id), None)
        if kdm4a:
            kdm4a_seq = str(kdm4a.seq)
            # Extract approximate JmjC domain (varies by protein)
            kdm4a_jmjc_approx = kdm4a_seq[150:310]  # Approximate
            
            print(f"\n  Comparison with KDM4A (active demethylase):")
            print(f"    KDM4A JmjC - H: {kdm4a_jmjc_approx.count('H')}")
            print(f"    KDM4A JmjC - D: {kdm4a_jmjc_approx.count('D')}")
            print(f"    KDM4A JmjC - E: {kdm4a_jmjc_approx.count('E')}")
    
    # Analyze conservation if we have homologs
    if len(jmjc_domains) > 1:
        print("\n6. Conservation Analysis:")
        print("-" * 40)
        
        conservation = analyze_conservation(jmjc_file)
        if conservation:
            # Find highly conserved positions
            high_conservation = [(i, score) for i, score in enumerate(conservation) if score > 0.8]
            
            print(f"  Highly conserved positions (>80% identity): {len(high_conservation)}")
            
            # Check conservation at critical positions
            # Position 38 in JmjC domain (approximate HXD location)
            if len(conservation) > 38:
                print(f"  Conservation at HXD region (pos ~38): {conservation[38]:.2f}")
    
    # Generate summary
    summary_file = results_dir / "conservation_analysis.txt"
    with open(summary_file, "w") as f:
        f.write("Conservation Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Epe1 JmjC Domain Features:\n")
        f.write(f"- Domain position: 243-402\n")
        f.write(f"- Length: {len(epe1_jmjc)} aa\n")
        f.write(f"- Histidines: {h_count}\n")
        f.write(f"- Aspartates: {d_count}\n")
        f.write(f"- Glutamates: {e_count}\n")
        f.write(f"- HX[DE] motifs: {len(matches)}\n")
        
        if matches:
            f.write("\nHX[DE] motifs found:\n")
            for match in matches:
                abs_pos = match.start() + 243
                f.write(f"  - Position {abs_pos}: {match.group()}\n")
        else:
            f.write("\n⚠️ No canonical HXD/HXE motifs found\n")
        
        # Count HXD/HXE motifs from earlier analysis
        hxd_count = len(hxd_matches) if 'hxd_matches' in locals() else 0
        hxe_count = len(hxe_matches) if 'hxe_matches' in locals() else 0
        
        f.write("\nConclusion based on analysis:\n")
        if hxd_count > 0:
            f.write(f"Found {hxd_count} HXD-like motif(s) in Epe1 JmjC domain.\n")
            f.write("Detailed analysis of the specific residues (e.g., HVD vs HXD)\n")
            f.write("is needed to determine if catalytic activity is retained.\n")
        elif hxe_count > 0:
            f.write(f"Found {hxe_count} HXE motif(s) but no HXD motifs.\n")
            f.write("This suggests altered or absent catalytic function.\n")
        else:
            f.write("The Epe1 JmjC domain lacks canonical HXD/HXE motifs,\n")
            f.write("indicating loss of demethylase catalytic activity.\n")
    
    print(f"\n✓ Analysis complete. Results saved to {summary_file}")

if __name__ == "__main__":
    main()