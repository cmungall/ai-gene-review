#!/usr/bin/env python3
"""
Analyze JmjC domain architecture and catalytic residues in Epe1.
Compare with active JmjC demethylases to identify missing residues.
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from pathlib import Path
import json
import re

def load_sequences():
    """Load all sequences from data directory."""
    sequences = {}
    data_dir = Path("data")
    
    # Load Epe1
    with open(data_dir / "epe1_spombe.fasta") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences["Epe1_spombe"] = str(record.seq)
    
    # Load active demethylases
    for fasta_file in data_dir.glob("kdm*.fasta"):
        with open(fasta_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                name = fasta_file.stem.upper()
                sequences[name] = str(record.seq)
    
    return sequences

def extract_jmjc_region(sequence, start_motif="HX[DE]", end_length=200):
    """
    Extract approximate JmjC domain region based on conserved motifs.
    JmjC domains typically contain HXD/HXE motifs for Fe(II) coordination.
    """
    # Look for HXD or HXE motif
    pattern = re.compile(r"H.[DE]")
    matches = list(pattern.finditer(sequence))
    
    if matches:
        # Take region around first match
        start = max(0, matches[0].start() - 50)
        end = min(len(sequence), start + end_length)
        return sequence[start:end], start
    return None, None

def identify_catalytic_residues(sequence):
    """
    Identify key catalytic residues in JmjC domain.
    Active JmjC demethylases require:
    1. HXD/HXE motif for Fe(II) binding
    2. Additional histidine for Fe(II) coordination
    3. Lysine or serine for α-ketoglutarate binding
    """
    catalytic_features = {
        "HXD_motif": [],
        "HXE_motif": [],
        "histidines": [],
        "lysines": [],
        "serines": [],
        "aspartates": [],
        "glutamates": []
    }
    
    # Find HXD/HXE motifs
    hxd_pattern = re.compile(r"H.D")
    hxe_pattern = re.compile(r"H.E")
    
    for match in hxd_pattern.finditer(sequence):
        catalytic_features["HXD_motif"].append({
            "position": match.start() + 1,  # 1-based
            "motif": match.group()
        })
    
    for match in hxe_pattern.finditer(sequence):
        catalytic_features["HXE_motif"].append({
            "position": match.start() + 1,
            "motif": match.group()
        })
    
    # Find all potentially important residues
    for i, aa in enumerate(sequence):
        pos = i + 1  # 1-based position
        if aa == 'H':
            catalytic_features["histidines"].append(pos)
        elif aa == 'K':
            catalytic_features["lysines"].append(pos)
        elif aa == 'S':
            catalytic_features["serines"].append(pos)
        elif aa == 'D':
            catalytic_features["aspartates"].append(pos)
        elif aa == 'E':
            catalytic_features["glutamates"].append(pos)
    
    return catalytic_features

def analyze_epe1_domains():
    """Analyze Epe1 domain structure from UniProt data."""
    with open("data/epe1_uniprot.json") as f:
        uniprot_data = json.load(f)
    
    # Extract domain information
    domains = {}
    if "features" in uniprot_data:
        for feature in uniprot_data["features"]:
            if feature.get("type") == "Domain":
                location = feature.get("location", {})
                start = location.get("start", {}).get("value")
                end = location.get("end", {}).get("value")
                description = feature.get("description", "")
                domains[description] = {"start": start, "end": end}
    
    return domains

def compare_jmjc_domains(sequences):
    """Compare JmjC domains across proteins."""
    results = {}
    
    for protein_name, sequence in sequences.items():
        print(f"\nAnalyzing {protein_name}...")
        
        # Extract JmjC region (approximate)
        jmjc_region, start_pos = extract_jmjc_region(sequence)
        
        if jmjc_region:
            # Identify catalytic residues
            catalytic = identify_catalytic_residues(jmjc_region)
            
            results[protein_name] = {
                "jmjc_start": start_pos,
                "jmjc_length": len(jmjc_region),
                "catalytic_residues": catalytic,
                "has_HXD": len(catalytic["HXD_motif"]) > 0,
                "has_HXE": len(catalytic["HXE_motif"]) > 0,
                "total_histidines": len(catalytic["histidines"]),
                "jmjc_sequence": jmjc_region[:50] + "..." if len(jmjc_region) > 50 else jmjc_region
            }
        else:
            results[protein_name] = {
                "jmjc_start": None,
                "error": "Could not identify JmjC domain"
            }
    
    return results

def main():
    print("=" * 60)
    print("JmjC Domain Architecture Analysis")
    print("=" * 60)
    
    # Load sequences
    sequences = load_sequences()
    print(f"Loaded {len(sequences)} sequences")
    
    # Analyze Epe1 domains from UniProt
    print("\n1. Epe1 Domain Structure (from UniProt):")
    print("-" * 40)
    domains = analyze_epe1_domains()
    for domain_name, info in domains.items():
        print(f"  {domain_name}: {info['start']}-{info['end']}")
    
    # Get Epe1 sequence for detailed analysis
    epe1_seq = sequences.get("Epe1_spombe", "")
    
    # Check for JmjC domain in Epe1 specifically
    print("\n2. Epe1 JmjC Domain Analysis:")
    print("-" * 40)
    
    # The JmjC domain in Epe1 is approximately at positions 200-400
    # Based on literature, it should be around residues 200-350
    jmjc_start = 200
    jmjc_end = 350
    epe1_jmjc = epe1_seq[jmjc_start-1:jmjc_end]
    
    print(f"  JmjC domain region (approx): {jmjc_start}-{jmjc_end}")
    print(f"  Sequence length: {len(epe1_jmjc)} aa")
    
    # Analyze catalytic residues in Epe1 JmjC domain
    epe1_catalytic = identify_catalytic_residues(epe1_jmjc)
    
    print("\n  Catalytic residue analysis:")
    print(f"    HXD motifs: {len(epe1_catalytic['HXD_motif'])}")
    if epe1_catalytic['HXD_motif']:
        for motif in epe1_catalytic['HXD_motif']:
            abs_pos = motif['position'] + jmjc_start - 1
            print(f"      - Position {abs_pos}: {motif['motif']}")
    
    print(f"    HXE motifs: {len(epe1_catalytic['HXE_motif'])}")
    if epe1_catalytic['HXE_motif']:
        for motif in epe1_catalytic['HXE_motif']:
            abs_pos = motif['position'] + jmjc_start - 1
            print(f"      - Position {abs_pos}: {motif['motif']}")
    
    print(f"    Total Histidines: {len(epe1_catalytic['histidines'])}")
    print(f"    Total Aspartates: {len(epe1_catalytic['aspartates'])}")
    print(f"    Total Glutamates: {len(epe1_catalytic['glutamates'])}")
    
    # Compare with active demethylases
    print("\n3. Comparison with Active JmjC Demethylases:")
    print("-" * 40)
    
    comparison_results = compare_jmjc_domains(sequences)
    
    # Summary table
    print("\n  Summary of Fe(II)-binding motifs:")
    print("  " + "-" * 50)
    print(f"  {'Protein':<20} {'HXD':<6} {'HXE':<6} {'Histidines':<10}")
    print("  " + "-" * 50)
    
    for protein, data in comparison_results.items():
        if "error" not in data:
            has_hxd = "Yes" if data["has_HXD"] else "No"
            has_hxe = "Yes" if data["has_HXE"] else "No"
            hist_count = data["total_histidines"]
            print(f"  {protein:<20} {has_hxd:<6} {has_hxe:<6} {hist_count:<10}")
    
    # Detailed comparison of critical positions
    print("\n4. Critical Residue Conservation:")
    print("-" * 40)
    
    # Known critical positions in active demethylases (approximate)
    # These typically include:
    # - Fe(II) binding: H-X-D/E, additional H
    # - α-KG binding: K or S
    
    active_proteins = [p for p in comparison_results.keys() if p.startswith("KDM")]
    
    if active_proteins:
        reference = active_proteins[0]
        ref_data = comparison_results[reference]
        
        print(f"\n  Using {reference} as reference for active demethylase")
        print(f"  Reference has:")
        if ref_data.get("catalytic_residues"):
            cat_res = ref_data["catalytic_residues"]
            if cat_res["HXD_motif"]:
                print(f"    - HXD motif at position(s): {[m['position'] for m in cat_res['HXD_motif']]}")
            if cat_res["HXE_motif"]:
                print(f"    - HXE motif at position(s): {[m['position'] for m in cat_res['HXE_motif']]}")
    
    # Check for C-terminal region
    print("\n5. C-terminal Region Analysis (HP1/Swi6 binding):")
    print("-" * 40)
    
    c_term_start = len(epe1_seq) - 100
    c_term = epe1_seq[c_term_start:]
    
    print(f"  C-terminal region ({c_term_start}-{len(epe1_seq)}):")
    print(f"  Length: {len(c_term)} aa")
    
    # Look for potential protein-protein interaction motifs
    # HP1 binding often involves hydrophobic patches
    hydrophobic_count = sum(1 for aa in c_term if aa in "FWYLMIV")
    basic_count = sum(1 for aa in c_term if aa in "KRH")
    
    print(f"  Hydrophobic residues: {hydrophobic_count} ({hydrophobic_count*100/len(c_term):.1f}%)")
    print(f"  Basic residues: {basic_count} ({basic_count*100/len(c_term):.1f}%)")
    
    # Save results
    results_file = Path("results") / "jmjc_domain_analysis.txt"
    results_file.parent.mkdir(exist_ok=True)
    
    with open(results_file, "w") as f:
        f.write("JmjC Domain Analysis Results\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Epe1 JmjC Domain Features:\n")
        f.write(f"- Domain location: {jmjc_start}-{jmjc_end}\n")
        f.write(f"- HXD motifs found: {len(epe1_catalytic['HXD_motif'])}\n")
        if epe1_catalytic['HXD_motif']:
            for motif in epe1_catalytic['HXD_motif']:
                abs_pos = motif['position'] + jmjc_start - 1
                f.write(f"  - Position {abs_pos}: {motif['motif']}\n")
        f.write(f"- HXE motifs found: {len(epe1_catalytic['HXE_motif'])}\n")
        if epe1_catalytic['HXE_motif']:
            for motif in epe1_catalytic['HXE_motif']:
                abs_pos = motif['position'] + jmjc_start - 1
                f.write(f"  - Position {abs_pos}: {motif['motif']}\n")
        f.write(f"- Total histidines in domain: {len(epe1_catalytic['histidines'])}\n")
        
        f.write("\nAnalysis Results:\n")
        if epe1_catalytic['HXD_motif']:
            # Check if it's actually HVD (valine instead of small/polar)
            motif_seq = epe1_catalytic['HXD_motif'][0]['motif']
            if 'V' in motif_seq:
                f.write(f"- Found atypical HVD motif at position {epe1_catalytic['HXD_motif'][0]['position'] + jmjc_start - 1}\n")
                f.write("- Valine (V) is hydrophobic and bulky, incompatible with Fe(II) coordination\n")
                f.write("- This substitution likely abolishes catalytic activity\n")
            else:
                f.write(f"- Found canonical HXD motif: {motif_seq}\n")
        elif epe1_catalytic['HXE_motif']:
            f.write(f"- Found HXE motif but no HXD motif\n")
            f.write("- Missing critical Fe(II)-binding residues\n")
        else:
            f.write("- No canonical HXD/HXE motifs found\n")
            f.write("- Lacks required Fe(II) binding residues for catalytic activity\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Comparison with Active Demethylases:\n")
        for protein, data in comparison_results.items():
            if "error" not in data:
                f.write(f"\n{protein}:\n")
                f.write(f"  Has HXD: {data['has_HXD']}\n")
                f.write(f"  Has HXE: {data['has_HXE']}\n")
                f.write(f"  Histidines: {data['total_histidines']}\n")
    
    print(f"\n✓ Results saved to: {results_file}")

if __name__ == "__main__":
    main()