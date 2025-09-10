#!/usr/bin/env python3
"""
Advanced homology analysis between tam10 and KNOP1 using proper BLAST alignment.
"""

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import tempfile
import subprocess
import os

def fetch_uniprot_sequence(uniprot_id):
    """Fetch sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        record = SeqIO.read(StringIO(response.text), "fasta")
        return str(record.seq), record.description
    else:
        raise Exception(f"Failed to fetch {uniprot_id}")

def run_emboss_needle(seq1, seq2, id1="seq1", id2="seq2"):
    """Run EMBOSS needle for global alignment if available"""
    try:
        # Check if needle is available
        result = subprocess.run(['which', 'needle'], capture_output=True, text=True)
        if result.returncode != 0:
            return None
        
        # Create temp files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f1:
            f1.write(f">{id1}\n{seq1}\n")
            file1 = f1.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f2:
            f2.write(f">{id2}\n{seq2}\n")
            file2 = f2.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.needle', delete=False) as out:
            outfile = out.name
        
        # Run needle
        cmd = ['needle', file1, file2, outfile, '-gapopen', '10', '-gapextend', '0.5']
        subprocess.run(cmd, capture_output=True)
        
        # Parse output
        with open(outfile, 'r') as f:
            content = f.read()
        
        # Clean up
        os.unlink(file1)
        os.unlink(file2)
        os.unlink(outfile)
        
        return content
    except:
        return None

def calculate_real_identity(seq1, seq2):
    """Calculate identity without alignment (for verification)"""
    # This gives us a baseline - how many identical residues if we just compare position by position
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    return (matches / min_len) * 100

def run_blast_api(query_seq, subject_seq):
    """Use NCBI BLAST API for alignment"""
    # Note: This would require actual BLAST setup or API access
    # For now, we'll use local alignment
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    from Bio.SubsMat import MatrixInfo as matlist
    
    # Use BLOSUM62 matrix for more accurate protein alignment
    matrix = matlist.blosum62
    
    # Perform alignment with gap penalties
    alignments = pairwise2.align.globalds(query_seq, subject_seq, matrix, -10, -0.5)
    
    if alignments:
        best = alignments[0]
        aln_seq1 = best[0]
        aln_seq2 = best[1]
        score = best[2]
        
        # Calculate proper identity
        matches = 0
        aln_len = 0
        for a, b in zip(aln_seq1, aln_seq2):
            if a != '-' and b != '-':
                aln_len += 1
                if a == b:
                    matches += 1
        
        if aln_len > 0:
            identity = (matches / aln_len) * 100
        else:
            identity = 0
        
        # Calculate coverage
        coverage_q = (len(aln_seq1.replace('-', '')) / len(query_seq)) * 100
        coverage_s = (len(aln_seq2.replace('-', '')) / len(subject_seq)) * 100
        
        return {
            'identity': identity,
            'matches': matches,
            'alignment_length': aln_len,
            'query_coverage': coverage_q,
            'subject_coverage': coverage_s,
            'score': score,
            'alignment': format_alignment(*best)
        }
    
    return None

def main():
    print("=" * 80)
    print("ADVANCED HOMOLOGY ANALYSIS: tam10 vs KNOP1")
    print("=" * 80)
    print()
    
    # Fetch sequences
    tam10_id = "G2TRQ9"
    knop1_id = "Q1ED39"
    
    tam10_seq, tam10_desc = fetch_uniprot_sequence(tam10_id)
    knop1_seq, knop1_desc = fetch_uniprot_sequence(knop1_id)
    
    print(f"tam10: {len(tam10_seq)} aa")
    print(f"KNOP1: {len(knop1_seq)} aa")
    print()
    
    # First, check direct positional identity (no alignment)
    print("Direct Positional Comparison (no alignment):")
    print("-" * 40)
    direct_identity = calculate_real_identity(tam10_seq, knop1_seq)
    print(f"Position-by-position identity (first {min(len(tam10_seq), len(knop1_seq))} aa): {direct_identity:.1f}%")
    print()
    
    # Run proper alignment
    print("Global Sequence Alignment (BLOSUM62):")
    print("-" * 40)
    result = run_blast_api(tam10_seq, knop1_seq)
    
    if result:
        print(f"Identity: {result['identity']:.1f}% ({result['matches']}/{result['alignment_length']} residues)")
        print(f"Query coverage: {result['query_coverage']:.1f}%")
        print(f"Subject coverage: {result['subject_coverage']:.1f}%")
        print(f"Alignment score: {result['score']:.1f}")
        print()
        
        # Show first part of alignment
        print("Alignment excerpt (first 500 characters):")
        print(result['alignment'][:500])
        print()
        
        # Save full alignment
        with open("detailed_alignment.txt", "w") as f:
            f.write("DETAILED ALIGNMENT: tam10 vs KNOP1\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"tam10 ({tam10_id}): {tam10_desc}\n")
            f.write(f"KNOP1 ({knop1_id}): {knop1_desc}\n\n")
            f.write(f"Statistics:\n")
            f.write(f"  Identity: {result['identity']:.1f}%\n")
            f.write(f"  Matches: {result['matches']}/{result['alignment_length']}\n")
            f.write(f"  Query coverage: {result['query_coverage']:.1f}%\n")
            f.write(f"  Subject coverage: {result['subject_coverage']:.1f}%\n\n")
            f.write("Full Alignment:\n")
            f.write(result['alignment'])
    
    # Check for shared motifs
    print("=" * 80)
    print("SHARED FEATURES ANALYSIS")
    print("-" * 40)
    
    # Check N-terminal region (often more conserved)
    n_term_len = 50
    if len(tam10_seq) >= n_term_len and len(knop1_seq) >= n_term_len:
        n_term_result = run_blast_api(tam10_seq[:n_term_len], knop1_seq[:n_term_len])
        if n_term_result:
            print(f"N-terminal (first {n_term_len} aa) identity: {n_term_result['identity']:.1f}%")
    
    # Check composition similarity
    def get_composition(seq):
        comp = {}
        for aa in seq:
            comp[aa] = comp.get(aa, 0) + 1
        total = len(seq)
        return {aa: (count/total)*100 for aa, count in comp.items()}
    
    tam10_comp = get_composition(tam10_seq)
    knop1_comp = get_composition(knop1_seq)
    
    # Compare lysine content (key feature of both proteins)
    print(f"\nLysine content:")
    print(f"  tam10: {tam10_comp.get('K', 0):.1f}%")
    print(f"  KNOP1: {knop1_comp.get('K', 0):.1f}%")
    
    # Calculate composition similarity
    shared_aa = set(tam10_comp.keys()) & set(knop1_comp.keys())
    comp_diff = sum(abs(tam10_comp.get(aa, 0) - knop1_comp.get(aa, 0)) for aa in shared_aa)
    comp_similarity = 100 - (comp_diff / 2)  # Normalize
    print(f"\nComposition similarity: {comp_similarity:.1f}%")
    
    print()
    print("=" * 80)
    print("ASSESSMENT")
    print("-" * 40)
    
    if result['identity'] < 20:
        print("⚠️  VERY WEAK HOMOLOGY")
        print()
        print("The sequence identity is well below typical ortholog thresholds (usually >25-30%).")
        print("This suggests:")
        print("1. These proteins are likely NOT orthologs in the classical sense")
        print("2. The ISO annotation may be based on:")
        print("   - Functional similarity (both are nucleolar, lysine-rich)")
        print("   - Domain architecture similarity")
        print("   - Computational prediction that needs experimental validation")
        print("3. The deep research statement about 'no orthologs' appears more accurate")
        print()
        print("RECOMMENDATION: The ISO annotation should be reviewed and possibly removed")
        print("unless there is strong experimental evidence for functional conservation.")
    elif result['identity'] < 30:
        print("⚠️  WEAK HOMOLOGY")
        print("Possible distant relationship but requires additional evidence.")
    else:
        print("✓  SIGNIFICANT HOMOLOGY")
        print("Supporting potential orthologous relationship.")
    
    print()
    print("Files created:")
    print("  - detailed_alignment.txt")

if __name__ == "__main__":
    main()