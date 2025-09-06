#!/usr/bin/env python3
"""
Basic sequence analysis for C18orf21 protein
Analyzes composition, hydrophobicity, and basic properties
"""

import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
import json

def analyze_sequence(fasta_file):
    """Perform basic sequence analysis"""
    
    # Read the sequence
    record = next(SeqIO.parse(fasta_file, "fasta"))
    sequence = str(record.seq)
    
    # Create ProtParam object for analysis
    analyzer = ProtParam.ProteinAnalysis(sequence)
    
    results = {
        "id": record.id,
        "description": record.description,
        "length": len(sequence),
        "sequence": sequence,
        "molecular_weight": analyzer.molecular_weight(),
        "isoelectric_point": analyzer.isoelectric_point(),
        "instability_index": analyzer.instability_index(),
        "gravy": analyzer.gravy(),
        "aromaticity": analyzer.aromaticity(),
        "amino_acid_composition": analyzer.get_amino_acids_percent(),
        "secondary_structure_fraction": analyzer.secondary_structure_fraction(),
        "classification": "Unstable" if analyzer.instability_index() > 40 else "Stable"
    }
    
    # Analyze charged residues
    positive_charged = sequence.count('K') + sequence.count('R') + sequence.count('H')
    negative_charged = sequence.count('D') + sequence.count('E')
    results["positive_charged"] = positive_charged
    results["negative_charged"] = negative_charged
    results["net_charge"] = positive_charged - negative_charged
    
    # Identify potential features
    features = []
    
    # Look for lysine-rich regions (potential nuclear localization)
    for i in range(len(sequence) - 10):
        window = sequence[i:i+10]
        if window.count('K') + window.count('R') >= 4:
            features.append({
                "type": "lysine_rich_region",
                "start": i + 1,
                "end": i + 10,
                "sequence": window
            })
    
    # Look for phosphorylation sites
    phospho_sites = []
    for i, aa in enumerate(sequence):
        if aa in 'STY':  # Serine, Threonine, Tyrosine
            context = sequence[max(0, i-3):min(len(sequence), i+4)]
            phospho_sites.append({
                "position": i + 1,
                "residue": aa,
                "context": context
            })
    
    results["features"] = features
    results["potential_phosphorylation_sites"] = phospho_sites
    
    return results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 01_sequence_analysis.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = Path(sys.argv[1])
    if not fasta_file.exists():
        print(f"Error: File {fasta_file} not found")
        sys.exit(1)
    
    results = analyze_sequence(fasta_file)
    
    # Save results
    output_file = Path("results/01_sequence_analysis.json")
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Print summary
    print(f"Sequence Analysis for {results['id']}")
    print(f"Length: {results['length']} amino acids")
    print(f"Molecular Weight: {results['molecular_weight']:.2f} Da")
    print(f"Isoelectric Point: {results['isoelectric_point']:.2f}")
    print(f"Instability Index: {results['instability_index']:.2f} ({results['classification']})")
    print(f"GRAVY: {results['gravy']:.2f}")
    print(f"Net Charge: {results['net_charge']}")
    print(f"Found {len(results['features'])} lysine-rich regions")
    print(f"Found {len(results['potential_phosphorylation_sites'])} potential phosphorylation sites")
    print(f"\nResults saved to {output_file}")