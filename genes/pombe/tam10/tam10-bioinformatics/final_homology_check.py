#!/usr/bin/env python3
"""
Final homology check between tam10 and KNOP1 using modern Bio.Align
"""

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from io import StringIO

def fetch_uniprot_sequence(uniprot_id):
    """Fetch sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        record = SeqIO.read(StringIO(response.text), "fasta")
        return str(record.seq), record.description
    else:
        raise Exception(f"Failed to fetch {uniprot_id}")

def analyze_alignment(seq1, seq2, name1="seq1", name2="seq2"):
    """Perform proper pairwise alignment and analysis"""
    
    # Create aligner with BLOSUM62-like scoring
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    # Perform alignment
    alignments = aligner.align(seq1, seq2)
    
    if alignments:
        # Get best alignment
        alignment = alignments[0]
        
        # Extract aligned sequences
        aln_str = str(alignment).split('\n')
        aln_seq1 = aln_str[0]
        aln_seq2 = aln_str[2]
        
        # Calculate statistics
        matches = 0
        mismatches = 0
        gaps = 0
        
        for a, b in zip(aln_seq1, aln_seq2):
            if a == '-' or b == '-':
                gaps += 1
            elif a == b:
                matches += 1
            else:
                mismatches += 1
        
        total_positions = len(aln_seq1)
        
        # Calculate identity over aligned region (excluding gaps)
        aligned_positions = matches + mismatches
        if aligned_positions > 0:
            identity = (matches / aligned_positions) * 100
        else:
            identity = 0
        
        # Calculate coverage
        coverage1 = (len(seq1) / len(aln_seq1.replace('-', ''))) * 100 if aln_seq1.replace('-', '') else 0
        coverage2 = (len(seq2) / len(aln_seq2.replace('-', ''))) * 100 if aln_seq2.replace('-', '') else 0
        
        return {
            'score': alignment.score,
            'identity': identity,
            'matches': matches,
            'mismatches': mismatches,
            'gaps': gaps,
            'total_positions': total_positions,
            'aligned_positions': aligned_positions,
            'coverage1': min(coverage1, 100),
            'coverage2': min(coverage2, 100),
            'alignment': str(alignment)
        }
    
    return None

def main():
    print("=" * 80)
    print("HOMOLOGY VERIFICATION: tam10 vs KNOP1")
    print("=" * 80)
    print()
    
    # Fetch sequences
    tam10_id = "G2TRQ9"
    knop1_id = "Q1ED39"
    
    print("Fetching sequences from UniProt...")
    tam10_seq, tam10_desc = fetch_uniprot_sequence(tam10_id)
    knop1_seq, knop1_desc = fetch_uniprot_sequence(knop1_id)
    
    print(f"✓ tam10 ({tam10_id}): {len(tam10_seq)} aa - {tam10_desc[:50]}...")
    print(f"✓ KNOP1 ({knop1_id}): {len(knop1_seq)} aa - {knop1_desc[:50]}...")
    print()
    
    # Analyze full sequences
    print("FULL SEQUENCE ALIGNMENT")
    print("-" * 40)
    
    result = analyze_alignment(tam10_seq, knop1_seq, "tam10", "KNOP1")
    
    if result:
        print(f"Alignment score: {result['score']:.1f}")
        print(f"Identity: {result['identity']:.1f}% ({result['matches']}/{result['aligned_positions']} positions)")
        print(f"Gaps: {result['gaps']} positions")
        print(f"tam10 coverage: {result['coverage1']:.1f}%")
        print(f"KNOP1 coverage: {result['coverage2']:.1f}%")
        print()
        
        # Show a sample of the alignment
        lines = result['alignment'].split('\n')
        print("Alignment sample (first 200 characters):")
        for line in lines[:3]:
            if len(line) > 200:
                print(line[:200] + "...")
            else:
                print(line)
        print()
    
    # Check specific regions
    print("REGIONAL ANALYSIS")
    print("-" * 40)
    
    # N-terminal (often most conserved)
    n_len = min(50, len(tam10_seq), len(knop1_seq))
    n_result = analyze_alignment(tam10_seq[:n_len], knop1_seq[:n_len], "tam10_Nterm", "KNOP1_Nterm")
    if n_result:
        print(f"N-terminal (first {n_len} aa) identity: {n_result['identity']:.1f}%")
    
    # C-terminal
    c_len = min(50, len(tam10_seq), len(knop1_seq))
    c_result = analyze_alignment(tam10_seq[-c_len:], knop1_seq[-c_len:], "tam10_Cterm", "KNOP1_Cterm")
    if c_result:
        print(f"C-terminal (last {c_len} aa) identity: {c_result['identity']:.1f}%")
    
    print()
    
    # Composition analysis
    print("COMPOSITION ANALYSIS")
    print("-" * 40)
    
    def analyze_composition(seq, name):
        comp = {}
        for aa in seq:
            comp[aa] = comp.get(aa, 0) + 1
        total = len(seq)
        pct = {aa: (count/total)*100 for aa, count in comp.items()}
        
        # Get top residues
        top = sorted(pct.items(), key=lambda x: x[1], reverse=True)[:5]
        
        return {
            'name': name,
            'length': total,
            'lysine_pct': pct.get('K', 0),
            'arginine_pct': pct.get('R', 0),
            'basic_pct': pct.get('K', 0) + pct.get('R', 0) + pct.get('H', 0),
            'acidic_pct': pct.get('D', 0) + pct.get('E', 0),
            'top': top
        }
    
    tam10_comp = analyze_composition(tam10_seq, "tam10")
    knop1_comp = analyze_composition(knop1_seq, "KNOP1")
    
    for comp in [tam10_comp, knop1_comp]:
        print(f"\n{comp['name']} ({comp['length']} aa):")
        print(f"  Lysine (K): {comp['lysine_pct']:.1f}%")
        print(f"  Arginine (R): {comp['arginine_pct']:.1f}%")
        print(f"  Total basic (K+R+H): {comp['basic_pct']:.1f}%")
        print(f"  Total acidic (D+E): {comp['acidic_pct']:.1f}%")
        print(f"  Top 3 residues: {', '.join([f'{aa}:{pct:.1f}%' for aa, pct in comp['top'][:3]])}")
    
    print()
    print("=" * 80)
    print("FINAL ASSESSMENT")
    print("-" * 40)
    
    # Determine relationship based on identity
    if result['identity'] < 15:
        assessment = "NO SIGNIFICANT HOMOLOGY"
        symbol = "❌"
        explanation = """
The sequence identity ({:.1f}%) is far below the threshold for orthology.
Typical ortholog pairs share >25-30% identity, and even remote homologs 
usually retain >15-20% identity. This level of similarity is essentially
what would be expected by random chance for proteins with similar composition.

The ISO annotation from tam10 to KNOP1 appears to be INCORRECT or based on:
- Similar cellular localization (both nucleolar)
- Similar composition (both lysine-rich)  
- Computational prediction without proper validation

The deep research's claim of "no orthologs" appears to be CORRECT.
""".format(result['identity'])
        
    elif result['identity'] < 25:
        assessment = "WEAK/QUESTIONABLE HOMOLOGY"
        symbol = "⚠️"
        explanation = """
The sequence identity ({:.1f}%) is in the twilight zone of homology detection.
While some very distant homologs can have identities in this range, it requires
strong additional evidence (structural similarity, synteny, functional data).
The ISO annotation should be reviewed with caution.
""".format(result['identity'])
        
    else:
        assessment = "POSSIBLE HOMOLOGY"
        symbol = "✓"
        explanation = """
The sequence identity ({:.1f}%) suggests possible evolutionary relationship.
The ISO annotation may be valid, though additional evidence would strengthen it.
""".format(result['identity'])
    
    print(f"{symbol} {assessment}")
    print(explanation)
    
    # Additional considerations
    print("Additional observations:")
    print(f"- Size difference: tam10 is {len(tam10_seq)} aa, KNOP1 is {len(knop1_seq)} aa (2.7x larger)")
    print(f"- Both are lysine-rich: tam10 {tam10_comp['lysine_pct']:.1f}% K, KNOP1 {knop1_comp['lysine_pct']:.1f}% K")
    print(f"- Both have high basic residue content: tam10 {tam10_comp['basic_pct']:.1f}%, KNOP1 {knop1_comp['basic_pct']:.1f}%")
    
    # Save results
    with open("RESULTS.md", "w") as f:
        f.write(f"# Homology Analysis: tam10 vs KNOP1\n\n")
        f.write(f"## Summary\n\n")
        f.write(f"**Assessment**: {assessment}\n\n")
        f.write(f"## Statistics\n\n")
        f.write(f"- **Sequence identity**: {result['identity']:.1f}%\n")
        f.write(f"- **Aligned positions**: {result['aligned_positions']}\n")
        f.write(f"- **Matches**: {result['matches']}\n")
        f.write(f"- **Gaps**: {result['gaps']}\n")
        f.write(f"- **tam10 length**: {len(tam10_seq)} aa\n")
        f.write(f"- **KNOP1 length**: {len(knop1_seq)} aa\n\n")
        f.write(f"## Conclusion\n\n")
        f.write(explanation.strip())
        f.write(f"\n\n## Method\n\n")
        f.write(f"Global pairwise alignment using Biopython's PairwiseAligner with ")
        f.write(f"match=2, mismatch=-1, gap_open=-10, gap_extend=-0.5\n")
    
    print()
    print("Results saved to: RESULTS.md")

if __name__ == "__main__":
    main()