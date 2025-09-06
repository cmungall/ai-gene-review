#!/usr/bin/env python3
"""
Motif and domain analysis for C18orf21
Searches for conserved motifs and functional patterns
"""

import sys
import re
from pathlib import Path
from Bio import SeqIO
import json

def search_motifs(sequence):
    """Search for known functional motifs and patterns"""
    
    motifs = []
    
    # Nuclear localization signals (NLS)
    nls_patterns = [
        (r'[KR]{2,}.{0,2}[KR]{2,}', 'Bipartite NLS'),
        (r'[KR]{4,}', 'Classical NLS'),
        (r'P[KR]{2,4}', 'Pat7 NLS'),
    ]
    
    for pattern, name in nls_patterns:
        for match in re.finditer(pattern, sequence):
            motifs.append({
                "type": "NLS",
                "name": name,
                "start": match.start() + 1,
                "end": match.end(),
                "sequence": match.group(),
                "pattern": pattern
            })
    
    # Phosphorylation motifs
    phospho_patterns = [
        (r'R..[ST]', 'PKA/PKC motif'),
        (r'[ST]P', 'Proline-directed kinase motif'),
        (r'[DE].{0,2}[ST]', 'CK2 motif'),
        (r'[ST]..[DE]', 'CK1 motif'),
        (r'[RK].{2}[ST]', 'CAMK2 motif'),
        (r'[ST]Q', 'ATM/ATR motif'),
    ]
    
    for pattern, name in phospho_patterns:
        for match in re.finditer(pattern, sequence):
            motifs.append({
                "type": "Phosphorylation",
                "name": name,
                "start": match.start() + 1,
                "end": match.end(),
                "sequence": match.group(),
                "pattern": pattern
            })
    
    # Protein-protein interaction motifs
    ppi_patterns = [
        (r'P..P', 'SH3 binding motif'),
        (r'[ILV]..L[ILV]', 'Leucine-rich motif'),
        (r'[FWY]..[FWY]', 'Aromatic interaction motif'),
        (r'L.{0,4}L.{0,4}L', 'Leucine zipper-like'),
    ]
    
    for pattern, name in ppi_patterns:
        for match in re.finditer(pattern, sequence):
            motifs.append({
                "type": "PPI",
                "name": name,
                "start": match.start() + 1,
                "end": match.end(),
                "sequence": match.group(),
                "pattern": pattern
            })
    
    # DNA/RNA binding motifs
    binding_patterns = [
        (r'[KR]{2,3}.{1,3}[KR]{2,3}', 'Basic DNA binding motif'),
        (r'C.{2}C.{3,20}C.{2}C', 'Zinc finger-like motif'),
        (r'[KR]G[KR]', 'RGG-like RNA binding'),
    ]
    
    for pattern, name in binding_patterns:
        for match in re.finditer(pattern, sequence):
            motifs.append({
                "type": "Nucleic_acid_binding",
                "name": name,
                "start": match.start() + 1,
                "end": match.end(),
                "sequence": match.group(),
                "pattern": pattern
            })
    
    # Degradation signals
    degradation_patterns = [
        (r'R..L.{4,5}[ILMV].[ILMV]', 'Destruction box (D-box)'),
        (r'[RK].L.{2,3}[ILMV].N', 'KEN box'),
        (r'[DE]SG.{2,3}[ST]', 'PEST motif'),
    ]
    
    for pattern, name in degradation_patterns:
        for match in re.finditer(pattern, sequence):
            motifs.append({
                "type": "Degradation",
                "name": name,
                "start": match.start() + 1,
                "end": match.end(),
                "sequence": match.group(),
                "pattern": pattern
            })
    
    return motifs

def analyze_composition_bias(sequence):
    """Analyze amino acid composition biases"""
    
    aa_counts = {}
    for aa in sequence:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1
    
    total = len(sequence)
    aa_freq = {aa: count/total for aa, count in aa_counts.items()}
    
    # Expected frequencies (from UniProt human proteome)
    expected_freq = {
        'A': 0.070, 'R': 0.056, 'N': 0.036, 'D': 0.047, 'C': 0.023,
        'Q': 0.048, 'E': 0.071, 'G': 0.066, 'H': 0.026, 'I': 0.044,
        'L': 0.100, 'K': 0.058, 'M': 0.021, 'F': 0.037, 'P': 0.063,
        'S': 0.083, 'T': 0.054, 'W': 0.011, 'Y': 0.027, 'V': 0.060
    }
    
    enriched = []
    depleted = []
    
    for aa, freq in aa_freq.items():
        expected = expected_freq.get(aa, 0.05)
        ratio = freq / expected
        
        if ratio > 1.5:
            enriched.append({
                "amino_acid": aa,
                "frequency": freq,
                "expected": expected,
                "enrichment": ratio
            })
        elif ratio < 0.67:
            depleted.append({
                "amino_acid": aa,
                "frequency": freq,
                "expected": expected,
                "depletion": ratio
            })
    
    # Check for specific biases
    charged_freq = sum(aa_freq.get(aa, 0) for aa in 'DEKR')
    hydrophobic_freq = sum(aa_freq.get(aa, 0) for aa in 'AILMFWV')
    polar_freq = sum(aa_freq.get(aa, 0) for aa in 'STYNQ')
    
    biases = {
        "enriched_aa": sorted(enriched, key=lambda x: x["enrichment"], reverse=True),
        "depleted_aa": sorted(depleted, key=lambda x: x["depletion"]),
        "charged_fraction": charged_freq,
        "hydrophobic_fraction": hydrophobic_freq,
        "polar_fraction": polar_freq,
        "proline_fraction": aa_freq.get('P', 0),
        "lysine_fraction": aa_freq.get('K', 0),
        "serine_fraction": aa_freq.get('S', 0)
    }
    
    return biases

def find_repeats(sequence, min_length=3, max_length=10):
    """Find tandem repeats and repetitive sequences"""
    
    repeats = []
    
    for length in range(min_length, min(max_length + 1, len(sequence) // 2)):
        for i in range(len(sequence) - 2 * length + 1):
            pattern = sequence[i:i+length]
            
            # Check for tandem repeat
            count = 1
            j = i + length
            while j <= len(sequence) - length:
                if sequence[j:j+length] == pattern:
                    count += 1
                    j += length
                else:
                    break
            
            if count >= 2:
                repeats.append({
                    "type": "tandem",
                    "pattern": pattern,
                    "start": i + 1,
                    "end": i + count * length,
                    "count": count,
                    "length": length
                })
    
    # Remove overlapping repeats, keeping the longest
    filtered = []
    for r in sorted(repeats, key=lambda x: x["end"] - x["start"], reverse=True):
        if not any(f["start"] <= r["start"] <= f["end"] or 
                  f["start"] <= r["end"] <= f["end"] 
                  for f in filtered):
            filtered.append(r)
    
    return sorted(filtered, key=lambda x: x["start"])

def analyze_motifs(fasta_file):
    """Main motif analysis function"""
    
    # Read sequence
    record = next(SeqIO.parse(fasta_file, "fasta"))
    sequence = str(record.seq)
    
    # Search for motifs
    motifs = search_motifs(sequence)
    
    # Analyze composition bias
    composition_bias = analyze_composition_bias(sequence)
    
    # Find repeats
    repeats = find_repeats(sequence)
    
    # Organize motifs by type
    motifs_by_type = {}
    for motif in motifs:
        motif_type = motif["type"]
        if motif_type not in motifs_by_type:
            motifs_by_type[motif_type] = []
        motifs_by_type[motif_type].append(motif)
    
    # Check for specific features relevant to known phosphorylation sites
    known_phospho_sites = [126, 130, 139]
    phospho_context = []
    
    for pos in known_phospho_sites:
        if pos <= len(sequence):
            context_start = max(0, pos - 6)
            context_end = min(len(sequence), pos + 5)
            context = sequence[context_start:context_end]
            
            phospho_context.append({
                "position": pos,
                "residue": sequence[pos-1],
                "context": context,
                "context_positions": f"{context_start+1}-{context_end}"
            })
    
    results = {
        "id": record.id,
        "sequence_length": len(sequence),
        "total_motifs": len(motifs),
        "motifs_by_type": motifs_by_type,
        "composition_bias": composition_bias,
        "repeats": repeats,
        "known_phosphorylation_context": phospho_context,
        "summary": {
            "nls_found": len([m for m in motifs if m["type"] == "NLS"]),
            "phospho_motifs": len([m for m in motifs if m["type"] == "Phosphorylation"]),
            "ppi_motifs": len([m for m in motifs if m["type"] == "PPI"]),
            "binding_motifs": len([m for m in motifs if m["type"] == "Nucleic_acid_binding"]),
            "degradation_signals": len([m for m in motifs if m["type"] == "Degradation"]),
            "tandem_repeats": len(repeats)
        }
    }
    
    return results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 03_motif_analysis.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = Path(sys.argv[1])
    if not fasta_file.exists():
        print(f"Error: File {fasta_file} not found")
        sys.exit(1)
    
    results = analyze_motifs(fasta_file)
    
    # Save results
    output_file = Path("results/03_motif_analysis.json")
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Print summary
    print(f"Motif Analysis for {results['id']}")
    print(f"Total motifs found: {results['total_motifs']}")
    print(f"\nMotif summary:")
    for key, count in results['summary'].items():
        if count > 0:
            print(f"  {key}: {count}")
    
    print("\nComposition bias:")
    if results['composition_bias']['enriched_aa']:
        print("  Enriched amino acids:")
        for aa in results['composition_bias']['enriched_aa'][:3]:
            print(f"    {aa['amino_acid']}: {aa['enrichment']:.2f}x expected")
    
    print(f"\nKnown phosphorylation sites context:")
    for site in results['known_phosphorylation_context']:
        print(f"  {site['residue']}{site['position']}: {site['context']}")
    
    if results['motifs_by_type'].get('NLS'):
        print(f"\nNuclear localization signals found:")
        for nls in results['motifs_by_type']['NLS']:
            print(f"  {nls['name']} at {nls['start']}-{nls['end']}: {nls['sequence']}")
    
    print(f"\nResults saved to {output_file}")