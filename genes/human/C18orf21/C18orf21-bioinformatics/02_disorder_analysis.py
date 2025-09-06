#!/usr/bin/env python3
"""
Disorder and low-complexity region analysis for C18orf21
Predicts intrinsically disordered regions and binding sites
"""

import sys
from pathlib import Path
from Bio import SeqIO
import numpy as np
import json
import matplotlib.pyplot as plt

def calculate_disorder_propensity(sequence):
    """Calculate disorder propensity using a simplified method"""
    
    # Disorder-promoting amino acids
    disorder_aa = {'A': 0.06, 'R': 0.18, 'N': 0.007, 'D': 0.192, 
                   'C': -0.02, 'Q': 0.318, 'E': 0.736, 'G': 0.166,
                   'H': 0.303, 'I': -0.486, 'L': -0.326, 'K': 0.586,
                   'M': -0.397, 'F': -0.697, 'P': 0.987, 'S': 0.341,
                   'T': 0.059, 'W': -0.884, 'Y': -0.510, 'V': -0.121}
    
    window_size = 15
    scores = []
    
    for i in range(len(sequence)):
        start = max(0, i - window_size // 2)
        end = min(len(sequence), i + window_size // 2 + 1)
        window = sequence[start:end]
        
        score = sum(disorder_aa.get(aa, 0) for aa in window) / len(window)
        scores.append(score)
    
    return scores

def find_low_complexity_regions(sequence, window_size=20, threshold=2.5):
    """Find low-complexity regions using Shannon entropy"""
    
    regions = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        
        # Calculate Shannon entropy
        aa_counts = {}
        for aa in window:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        entropy = 0
        for count in aa_counts.values():
            p = count / window_size
            if p > 0:
                entropy -= p * np.log2(p)
        
        if entropy < threshold:
            regions.append({
                "start": i + 1,
                "end": i + window_size,
                "entropy": entropy,
                "sequence": window
            })
    
    # Merge overlapping regions
    merged = []
    for region in sorted(regions, key=lambda x: x["start"]):
        if merged and region["start"] <= merged[-1]["end"]:
            merged[-1]["end"] = max(merged[-1]["end"], region["end"])
            merged[-1]["entropy"] = min(merged[-1]["entropy"], region["entropy"])
        else:
            merged.append(region)
    
    return merged

def predict_binding_regions(sequence, disorder_scores):
    """Predict potential protein binding regions (MoRFs)"""
    
    # MoRFs are typically short disordered regions that undergo disorder-to-order transition
    morfs = []
    
    # Look for moderately disordered regions flanked by more ordered regions
    for i in range(5, len(sequence) - 5):
        if 0.2 < disorder_scores[i] < 0.6:  # Moderate disorder
            # Check if flanking regions are more ordered
            left_flank = np.mean(disorder_scores[max(0, i-5):i])
            right_flank = np.mean(disorder_scores[i+1:min(len(sequence), i+6)])
            
            if left_flank < disorder_scores[i] and right_flank < disorder_scores[i]:
                # Potential MoRF region - extend to find boundaries
                start = i
                while start > 0 and disorder_scores[start-1] > 0.1:
                    start -= 1
                end = i
                while end < len(sequence) - 1 and disorder_scores[end+1] > 0.1:
                    end += 1
                
                if 5 <= (end - start + 1) <= 25:  # Typical MoRF length
                    morfs.append({
                        "start": start + 1,
                        "end": end + 1,
                        "sequence": sequence[start:end+1],
                        "avg_disorder": np.mean(disorder_scores[start:end+1])
                    })
    
    # Remove duplicates
    unique_morfs = []
    for morf in morfs:
        if not any(m["start"] == morf["start"] and m["end"] == morf["end"] for m in unique_morfs):
            unique_morfs.append(morf)
    
    return unique_morfs

def analyze_disorder(fasta_file):
    """Main disorder analysis function"""
    
    # Read sequence
    record = next(SeqIO.parse(fasta_file, "fasta"))
    sequence = str(record.seq)
    
    # Calculate disorder propensity
    disorder_scores = calculate_disorder_propensity(sequence)
    
    # Find disordered regions (score > 0.5)
    disordered_regions = []
    in_disorder = False
    start = 0
    
    for i, score in enumerate(disorder_scores):
        if score > 0.5 and not in_disorder:
            start = i
            in_disorder = True
        elif score <= 0.5 and in_disorder:
            disordered_regions.append({
                "start": start + 1,
                "end": i,
                "length": i - start,
                "avg_score": np.mean(disorder_scores[start:i]),
                "sequence": sequence[start:i]
            })
            in_disorder = False
    
    # Handle if ends in disorder
    if in_disorder:
        disordered_regions.append({
            "start": start + 1,
            "end": len(sequence),
            "length": len(sequence) - start,
            "avg_score": np.mean(disorder_scores[start:]),
            "sequence": sequence[start:]
        })
    
    # Find low-complexity regions
    low_complexity = find_low_complexity_regions(sequence)
    
    # Predict binding regions
    binding_regions = predict_binding_regions(sequence, disorder_scores)
    
    # Calculate overall statistics
    disordered_residues = sum(r["length"] for r in disordered_regions)
    
    results = {
        "id": record.id,
        "sequence_length": len(sequence),
        "disorder_scores": disorder_scores,
        "disordered_regions": disordered_regions,
        "disordered_fraction": disordered_residues / len(sequence),
        "low_complexity_regions": low_complexity,
        "predicted_binding_regions": binding_regions,
        "avg_disorder_score": np.mean(disorder_scores),
        "max_disorder_score": max(disorder_scores),
        "min_disorder_score": min(disorder_scores)
    }
    
    # Plot disorder profile
    plt.figure(figsize=(12, 6))
    
    positions = list(range(1, len(sequence) + 1))
    plt.plot(positions, disorder_scores, 'b-', linewidth=1, label='Disorder propensity')
    plt.axhline(y=0.5, color='r', linestyle='--', alpha=0.5, label='Disorder threshold')
    
    # Highlight disordered regions
    for region in disordered_regions:
        plt.axvspan(region["start"], region["end"], alpha=0.2, color='red')
    
    # Highlight binding regions
    for region in binding_regions:
        plt.axvspan(region["start"], region["end"], alpha=0.3, color='green')
    
    plt.xlabel('Residue Position')
    plt.ylabel('Disorder Propensity')
    plt.title(f'Disorder Profile - {record.id}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plot_file = Path("results/02_disorder_profile.png")
    plot_file.parent.mkdir(exist_ok=True)
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Remove disorder_scores from JSON (too large)
    json_results = {k: v for k, v in results.items() if k != 'disorder_scores'}
    
    return json_results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 02_disorder_analysis.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = Path(sys.argv[1])
    if not fasta_file.exists():
        print(f"Error: File {fasta_file} not found")
        sys.exit(1)
    
    results = analyze_disorder(fasta_file)
    
    # Save results
    output_file = Path("results/02_disorder_analysis.json")
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Print summary
    print(f"Disorder Analysis for {results['id']}")
    print(f"Sequence length: {results['sequence_length']} aa")
    print(f"Disordered fraction: {results['disordered_fraction']:.2%}")
    print(f"Number of disordered regions: {len(results['disordered_regions'])}")
    print(f"Number of low-complexity regions: {len(results['low_complexity_regions'])}")
    print(f"Predicted binding regions (MoRFs): {len(results['predicted_binding_regions'])}")
    
    if results['disordered_regions']:
        print("\nDisordered regions:")
        for r in results['disordered_regions']:
            print(f"  {r['start']}-{r['end']} ({r['length']} aa): {r['sequence'][:20]}...")
    
    if results['predicted_binding_regions']:
        print("\nPredicted binding regions:")
        for r in results['predicted_binding_regions']:
            print(f"  {r['start']}-{r['end']}: {r['sequence']}")
    
    print(f"\nResults saved to {output_file}")
    print(f"Plot saved to results/02_disorder_profile.png")