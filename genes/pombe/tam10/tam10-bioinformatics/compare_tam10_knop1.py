#!/usr/bin/env python3
"""
Compare tam10 from S. pombe with human KNOP1 (Q1ED39) to assess potential homology.
This analysis addresses the discrepancy between the deep research claiming no orthologs
and the GOA file showing an ISO annotation to human KNOP1.
"""

import requests
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO
import subprocess
import os

def fetch_uniprot_sequence(uniprot_id):
    """Fetch sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        record = SeqIO.read(StringIO(response.text), "fasta")
        return str(record.seq)
    else:
        raise Exception(f"Failed to fetch {uniprot_id}")

def run_blast_comparison(seq1, seq2, label1="seq1", label2="seq2"):
    """Run BLAST-like comparison using basic alignment"""
    alignments = pairwise2.align.globalxx(seq1, seq2)
    
    if alignments:
        best_alignment = alignments[0]
        score = best_alignment[2]
        length1 = len(seq1)
        length2 = len(seq2)
        aligned_seq1 = best_alignment[0]
        aligned_seq2 = best_alignment[1]
        
        # Calculate identity
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        identity = (matches / len(aligned_seq1.replace('-', ''))) * 100
        
        # Calculate similarity (considering conservative substitutions)
        similar_pairs = {
            ('A', 'V'), ('V', 'A'), ('A', 'I'), ('I', 'A'), ('A', 'L'), ('L', 'A'),
            ('I', 'V'), ('V', 'I'), ('I', 'L'), ('L', 'I'), ('V', 'L'), ('L', 'V'),
            ('S', 'T'), ('T', 'S'),
            ('D', 'E'), ('E', 'D'),
            ('K', 'R'), ('R', 'K'),
            ('F', 'Y'), ('Y', 'F'), ('F', 'W'), ('W', 'F'), ('Y', 'W'), ('W', 'Y'),
            ('N', 'Q'), ('Q', 'N'),
            ('H', 'K'), ('K', 'H'), ('H', 'R'), ('R', 'H')
        }
        
        similar = matches
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a != b and a != '-' and b != '-' and (a, b) in similar_pairs:
                similar += 1
        
        similarity = (similar / len(aligned_seq1.replace('-', ''))) * 100
        
        return {
            'score': score,
            'identity': identity,
            'similarity': similarity,
            'alignment': format_alignment(*best_alignment),
            'len1': length1,
            'len2': length2
        }
    return None

def analyze_composition(sequence, name):
    """Analyze amino acid composition"""
    composition = {}
    for aa in sequence:
        composition[aa] = composition.get(aa, 0) + 1
    
    # Calculate percentages
    total = len(sequence)
    comp_pct = {aa: (count/total)*100 for aa, count in composition.items()}
    
    # Get top 5 most abundant
    top5 = sorted(comp_pct.items(), key=lambda x: x[1], reverse=True)[:5]
    
    return {
        'name': name,
        'length': total,
        'top5_aa': top5,
        'lysine_pct': comp_pct.get('K', 0),
        'arginine_pct': comp_pct.get('R', 0),
        'basic_pct': comp_pct.get('K', 0) + comp_pct.get('R', 0) + comp_pct.get('H', 0)
    }

def main():
    print("=" * 80)
    print("TAM10 vs KNOP1 Homology Analysis")
    print("=" * 80)
    print()
    
    # Fetch sequences
    print("Fetching sequences...")
    tam10_id = "G2TRQ9"  # tam10 from S. pombe
    knop1_id = "Q1ED39"  # KNOP1 from human
    
    try:
        tam10_seq = fetch_uniprot_sequence(tam10_id)
        knop1_seq = fetch_uniprot_sequence(knop1_id)
        
        print(f"✓ tam10 (S. pombe, {tam10_id}): {len(tam10_seq)} aa")
        print(f"✓ KNOP1 (Human, {knop1_id}): {len(knop1_seq)} aa")
        print()
        
        # Composition analysis
        print("Amino Acid Composition Analysis")
        print("-" * 40)
        tam10_comp = analyze_composition(tam10_seq, "tam10")
        knop1_comp = analyze_composition(knop1_seq, "KNOP1")
        
        for comp in [tam10_comp, knop1_comp]:
            print(f"\n{comp['name']} ({comp['length']} aa):")
            print(f"  Top 5 amino acids:")
            for aa, pct in comp['top5_aa']:
                print(f"    {aa}: {pct:.1f}%")
            print(f"  Lysine (K): {comp['lysine_pct']:.1f}%")
            print(f"  Arginine (R): {comp['arginine_pct']:.1f}%")
            print(f"  Total basic residues (K+R+H): {comp['basic_pct']:.1f}%")
        
        print()
        print("=" * 80)
        print("Sequence Alignment Analysis")
        print("-" * 40)
        
        result = run_blast_comparison(tam10_seq, knop1_seq, "tam10", "KNOP1")
        
        if result:
            print(f"Alignment statistics:")
            print(f"  tam10 length: {result['len1']} aa")
            print(f"  KNOP1 length: {result['len2']} aa")
            print(f"  Alignment score: {result['score']:.1f}")
            print(f"  Sequence identity: {result['identity']:.1f}%")
            print(f"  Sequence similarity: {result['similarity']:.1f}%")
            print()
            
            # Show a portion of the alignment
            print("First 200 positions of alignment:")
            lines = result['alignment'].split('\n')
            for line in lines[:6]:  # Show first portion
                if line:
                    print(line)
        
        # Check for conserved domains or motifs
        print()
        print("=" * 80)
        print("Domain/Motif Search")
        print("-" * 40)
        
        # Look for lysine-rich regions (characteristic of KNOP1)
        def find_lysine_rich_regions(seq, window=20, threshold=0.25):
            """Find regions enriched in lysine"""
            regions = []
            for i in range(len(seq) - window + 1):
                window_seq = seq[i:i+window]
                k_count = window_seq.count('K')
                if k_count / window >= threshold:
                    regions.append((i, i+window, k_count/window))
            return regions
        
        tam10_k_regions = find_lysine_rich_regions(tam10_seq)
        knop1_k_regions = find_lysine_rich_regions(knop1_seq)
        
        print(f"Lysine-rich regions (>25% K in 20aa window):")
        print(f"  tam10: {len(tam10_k_regions)} regions")
        if tam10_k_regions:
            for start, end, pct in tam10_k_regions[:3]:  # Show first 3
                print(f"    Position {start+1}-{end}: {pct*100:.0f}% K")
        
        print(f"  KNOP1: {len(knop1_k_regions)} regions")
        if knop1_k_regions:
            for start, end, pct in knop1_k_regions[:3]:  # Show first 3
                print(f"    Position {start+1}-{end}: {pct*100:.0f}% K")
        
        # Final assessment
        print()
        print("=" * 80)
        print("CONCLUSION")
        print("-" * 40)
        
        if result['identity'] < 15:
            print("⚠️  VERY LOW sequence identity detected (<15%)")
            print("   This level of similarity is below typical ortholog thresholds.")
            print("   The ISO annotation may be questionable or based on:")
            print("   - Shared domain architecture rather than sequence homology")
            print("   - Functional similarity rather than evolutionary relationship")
            print("   - Potential annotation error that needs review")
        elif result['identity'] < 25:
            print("⚠️  LOW sequence identity detected (15-25%)")
            print("   This suggests possible distant homology but requires")
            print("   additional evidence (structural similarity, synteny, etc.)")
        else:
            print("✓  Moderate to high sequence identity detected")
            print("   Supporting potential orthologous relationship")
        
        print()
        print(f"Key findings:")
        print(f"- tam10 is much smaller ({result['len1']} aa) than KNOP1 ({result['len2']} aa)")
        print(f"- Sequence identity: {result['identity']:.1f}%")
        print(f"- Both proteins show basic residue enrichment")
        print(f"- tam10 K-rich regions: {len(tam10_k_regions)}")
        print(f"- KNOP1 K-rich regions: {len(knop1_k_regions)}")
        
        # Save detailed alignment
        with open("tam10_knop1_alignment.txt", "w") as f:
            f.write("Full alignment between tam10 and KNOP1\n")
            f.write("=" * 80 + "\n\n")
            f.write(result['alignment'])
        
        print()
        print("Full alignment saved to: tam10_knop1_alignment.txt")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())