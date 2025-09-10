#!/usr/bin/env python3
"""
Correct calculation of sequence identity between tam10 and KNOP1
"""

import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from io import StringIO

def fetch_uniprot_sequence(uniprot_id):
    """Fetch sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        record = SeqIO.read(StringIO(response.text), "fasta")
        return str(record.seq)
    else:
        raise Exception(f"Failed to fetch {uniprot_id}")

def calculate_identity_properly(seq1, seq2):
    """Calculate sequence identity using proper global alignment"""
    
    # Create aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.1
    
    # Perform alignment
    alignments = aligner.align(seq1, seq2)
    
    if alignments:
        alignment = alignments[0]
        
        # Parse alignment to count matches
        aln_lines = str(alignment).split('\n')
        seq1_aln = aln_lines[0]  # First sequence with gaps
        seq2_aln = aln_lines[2]  # Second sequence with gaps
        
        matches = 0
        total_aligned = 0
        
        for c1, c2 in zip(seq1_aln, seq2_aln):
            if c1 != '-' or c2 != '-':  # Count all positions that aren't double gaps
                total_aligned += 1
                if c1 == c2 and c1 != '-':  # Count matches
                    matches += 1
        
        # Calculate identity over full alignment length
        identity_over_alignment = (matches / len(seq1_aln)) * 100 if len(seq1_aln) > 0 else 0
        
        # Calculate identity over shorter sequence
        shorter_len = min(len(seq1), len(seq2))
        identity_over_shorter = (matches / shorter_len) * 100 if shorter_len > 0 else 0
        
        # Calculate identity over aligned regions only (excluding gaps)
        non_gap_positions = sum(1 for c1, c2 in zip(seq1_aln, seq2_aln) 
                               if c1 != '-' and c2 != '-')
        identity_aligned_only = (matches / non_gap_positions) * 100 if non_gap_positions > 0 else 0
        
        return {
            'matches': matches,
            'alignment_length': len(seq1_aln),
            'non_gap_positions': non_gap_positions,
            'identity_over_alignment': identity_over_alignment,
            'identity_over_shorter': identity_over_shorter,
            'identity_aligned_only': identity_aligned_only,
            'seq1_len': len(seq1),
            'seq2_len': len(seq2),
            'alignment': alignment
        }
    
    return None

def main():
    print("CORRECT SEQUENCE IDENTITY CALCULATION")
    print("=" * 60)
    
    # Fetch sequences
    tam10_seq = fetch_uniprot_sequence("G2TRQ9")
    knop1_seq = fetch_uniprot_sequence("Q1ED39")
    
    print(f"tam10: {len(tam10_seq)} aa")
    print(f"KNOP1: {len(knop1_seq)} aa")
    print()
    
    result = calculate_identity_properly(tam10_seq, knop1_seq)
    
    if result:
        print("ALIGNMENT STATISTICS:")
        print("-" * 40)
        print(f"Total alignment length: {result['alignment_length']} positions")
        print(f"Non-gap aligned positions: {result['non_gap_positions']} positions")
        print(f"Identical residues: {result['matches']} positions")
        print()
        
        print("IDENTITY CALCULATIONS:")
        print("-" * 40)
        print(f"1. Over full alignment (incl. gaps): {result['identity_over_alignment']:.1f}%")
        print(f"   ({result['matches']}/{result['alignment_length']})")
        print()
        print(f"2. Over shorter sequence: {result['identity_over_shorter']:.1f}%")
        print(f"   ({result['matches']}/{min(result['seq1_len'], result['seq2_len'])})")
        print()
        print(f"3. Over aligned regions only: {result['identity_aligned_only']:.1f}%")
        print(f"   ({result['matches']}/{result['non_gap_positions']})")
        print()
        
        # Show alignment excerpt
        aln_str = str(result['alignment'])
        lines = aln_str.split('\n')
        
        print("ALIGNMENT EXCERPT (first 300 chars):")
        print("-" * 40)
        for i, line in enumerate(lines):
            if i < 3:  # Show first 3 lines
                if len(line) > 100:
                    print(line[:100])
                else:
                    print(line)
        
        print()
        print("INTERPRETATION:")
        print("-" * 40)
        
        # Use the most conservative measure (over shorter sequence)
        identity = result['identity_over_shorter']
        
        if identity < 15:
            print(f"❌ NO SIGNIFICANT HOMOLOGY (identity: {identity:.1f}%)")
            print()
            print("This level of sequence identity is essentially random.")
            print("The proteins are NOT homologous.")
            print()
            print("The ISO annotation linking tam10 to KNOP1 appears INCORRECT.")
            print("The deep research stating 'no orthologs' is CORRECT.")
        elif identity < 25:
            print(f"⚠️  MARGINAL SIMILARITY (identity: {identity:.1f}%)")
            print()
            print("This is in the twilight zone - possibly very distant homology")
            print("or convergent evolution. Additional evidence needed.")
        else:
            print(f"✓ POTENTIAL HOMOLOGY (identity: {identity:.1f}%)")
            print()
            print("This level of identity suggests possible evolutionary relationship.")
        
        # Save detailed results
        with open("final_verdict.txt", "w") as f:
            f.write("TAM10 vs KNOP1 HOMOLOGY ANALYSIS - FINAL VERDICT\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Sequence lengths:\n")
            f.write(f"  tam10: {result['seq1_len']} aa\n")
            f.write(f"  KNOP1: {result['seq2_len']} aa\n\n")
            f.write(f"Alignment statistics:\n")
            f.write(f"  Matches: {result['matches']}\n")
            f.write(f"  Total alignment: {result['alignment_length']} positions\n")
            f.write(f"  Non-gap positions: {result['non_gap_positions']}\n\n")
            f.write(f"Identity measurements:\n")
            f.write(f"  Over full alignment: {result['identity_over_alignment']:.1f}%\n")
            f.write(f"  Over shorter sequence: {result['identity_over_shorter']:.1f}%\n")
            f.write(f"  Over aligned regions: {result['identity_aligned_only']:.1f}%\n\n")
            f.write(f"VERDICT: ")
            if identity < 15:
                f.write("NO HOMOLOGY - ISO annotation appears incorrect\n")
            elif identity < 25:
                f.write("QUESTIONABLE - Requires additional evidence\n")
            else:
                f.write("POSSIBLE HOMOLOGY - ISO annotation may be valid\n")
        
        print()
        print("Results saved to: final_verdict.txt")

if __name__ == "__main__":
    main()