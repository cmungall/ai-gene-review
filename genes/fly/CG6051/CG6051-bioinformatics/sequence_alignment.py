#!/usr/bin/env python3
"""
Perform sequence alignment between Drosophila CG6051 and human LST2/ZFYVE28.
Focuses on domain conservation and functional motifs.
"""

import sys
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import json

def load_sequences():
    """Load both protein sequences."""
    sequences = {}
    
    # Load fly sequence
    fly_file = "CG6051.fasta"
    if Path(fly_file).exists():
        record = next(SeqIO.parse(fly_file, "fasta"))
        sequences["fly"] = {
            "id": str(record.id),
            "seq": str(record.seq),
            "length": len(record.seq)
        }
    
    # Load human sequence
    human_file = "human_LST2.fasta"
    if Path(human_file).exists():
        record = next(SeqIO.parse(human_file, "fasta"))
        sequences["human"] = {
            "id": str(record.id), 
            "seq": str(record.seq),
            "length": len(record.seq)
        }
    
    return sequences

def align_domains(seq1, seq2, name1="Seq1", name2="Seq2"):
    """Perform pairwise alignment of domain sequences."""
    # Use BLOSUM62 matrix for protein alignment
    matrix = substitution_matrices.load("BLOSUM62")
    
    # Create aligner
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'
    
    # Perform alignment
    alignments = aligner.align(seq1, seq2)
    
    if alignments:
        best = alignments[0]
        aligned_seq1 = str(best.target)
        aligned_seq2 = str(best.query)
        score = best.score
        
        # Calculate identity and similarity
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        identity = (matches / len(aligned_seq1.replace('-', ''))) * 100
        
        # Calculate similarity (including conservative substitutions)
        similar = 0
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a != '-' and b != '-':
                if a == b:
                    similar += 1
                elif matrix.get((a, b), 0) > 0:
                    similar += 1
        
        similarity = (similar / len(aligned_seq1.replace('-', ''))) * 100
        
        return {
            "aligned_seq1": aligned_seq1,
            "aligned_seq2": aligned_seq2,
            "score": score,
            "identity": identity,
            "similarity": similarity,
            "length": len(aligned_seq1)
        }
    
    return None

def compare_fyve_domains(fly_seq, human_seq):
    """Compare FYVE domains between species."""
    # Extract FYVE domains based on UniProt annotations
    fly_fyve = fly_seq[908:969]  # 909-969 in 1-based
    human_fyve = human_seq[814:875] if len(human_seq) > 875 else None
    
    if not human_fyve:
        return None
    
    print("\nFYVE Domain Comparison:")
    print("-" * 50)
    print(f"Fly FYVE (909-969): {len(fly_fyve)} aa")
    print(f"Human FYVE (815-875): {len(human_fyve)} aa")
    
    # Align FYVE domains
    alignment = align_domains(fly_fyve, human_fyve, "Fly_FYVE", "Human_FYVE")
    
    if alignment:
        print(f"Identity: {alignment['identity']:.1f}%")
        print(f"Similarity: {alignment['similarity']:.1f}%")
        
        # Show alignment
        print("\nAlignment:")
        for i in range(0, len(alignment['aligned_seq1']), 60):
            print(f"Fly:   {alignment['aligned_seq1'][i:i+60]}")
            print(f"Human: {alignment['aligned_seq2'][i:i+60]}")
            print()
        
        # Check conserved cysteines
        fly_cys = [i+1 for i, aa in enumerate(fly_fyve) if aa == 'C']
        human_cys = [i+1 for i, aa in enumerate(human_fyve) if aa == 'C']
        print(f"Fly cysteines at positions: {fly_cys}")
        print(f"Human cysteines at positions: {human_cys}")
    
    return alignment

def search_tos_in_alignment(fly_seq, human_seq):
    """Search for TOS motifs in both sequences."""
    print("\nTOS Motif Analysis:")
    print("-" * 50)
    
    # Canonical TOS: F[DE][LIVM][DE][LIVM]
    import re
    tos_pattern = r"F[DE][LIVM][DE][LIVM]"
    
    # Search in human
    human_tos = []
    for match in re.finditer(tos_pattern, human_seq):
        human_tos.append({
            "position": match.start() + 1,
            "sequence": match.group(),
            "context": human_seq[max(0, match.start()-5):min(len(human_seq), match.end()+5)]
        })
    
    # Search in fly
    fly_tos = []
    for match in re.finditer(tos_pattern, fly_seq):
        fly_tos.append({
            "position": match.start() + 1,
            "sequence": match.group(),
            "context": fly_seq[max(0, match.start()-5):min(len(fly_seq), match.end()+5)]
        })
    
    print(f"Human TOS motifs found: {len(human_tos)}")
    for tos in human_tos:
        print(f"  Position {tos['position']}: {tos['sequence']} in context '{tos['context']}'")
    
    print(f"\nFly TOS motifs found: {len(fly_tos)}")
    for tos in fly_tos:
        print(f"  Position {tos['position']}: {tos['sequence']} in context '{tos['context']}'")
    
    # Also check for the specific human FDIDI sequence
    if "FDIDI" in human_seq:
        pos = human_seq.find("FDIDI") + 1
        context = human_seq[max(0, pos-6):min(len(human_seq), pos+14)]
        print(f"\nHuman FDIDI motif at position {pos}: context '{context}'")
        
        # Look for similar sequences in fly
        # FDIDI, xDIDI, FDIxx patterns
        patterns = ["FDIDI", r"[FY]DI[DE][LIVM]", r"[FY][DE]I[DE][LIVM]"]
        print("\nSearching for FDIDI-like patterns in fly:")
        for pattern in patterns:
            if pattern == "FDIDI":
                if pattern in fly_seq:
                    pos = fly_seq.find(pattern) + 1
                    print(f"  Exact match '{pattern}' at position {pos}")
            else:
                matches = re.finditer(pattern, fly_seq)
                for match in matches:
                    print(f"  Pattern '{pattern}' match at position {match.start()+1}: {match.group()}")
    
    return {"human": human_tos, "fly": fly_tos}

def overall_alignment(fly_seq, human_seq):
    """Perform and analyze overall sequence alignment."""
    print("\nOverall Sequence Alignment:")
    print("-" * 50)
    
    # For large sequences, we'll calculate overall stats without full alignment display
    matrix = substitution_matrices.load("BLOSUM62")
    
    # Create aligner for local alignment
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'local'
    
    # Use local alignment for better detection of conserved regions
    alignments = aligner.align(fly_seq, human_seq)
    
    if alignments:
        best = alignments[0]
        aligned_seq1 = str(best.target)
        aligned_seq2 = str(best.query)
        score = best.score
        
        # Calculate statistics
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        gaps1 = aligned_seq1.count('-')
        gaps2 = aligned_seq2.count('-')
        
        print(f"Alignment score: {score:.1f}")
        print(f"Aligned region length: {len(aligned_seq1)}")
        print(f"Identical positions: {matches}")
        print(f"Gaps in fly: {gaps1}")
        print(f"Gaps in human: {gaps2}")
        
        # Find highly conserved regions (stretches of >10 identical residues)
        conserved_regions = []
        in_conserved = False
        start_pos = 0
        
        for i, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2)):
            if a == b and a != '-':
                if not in_conserved:
                    in_conserved = True
                    start_pos = i
            else:
                if in_conserved and i - start_pos >= 10:
                    conserved_regions.append((start_pos, i))
                in_conserved = False
        
        if conserved_regions:
            print(f"\nHighly conserved regions (≥10 continuous matches):")
            for start, end in conserved_regions[:5]:  # Show first 5
                print(f"  Position {start}-{end}: {aligned_seq1[start:end]}")
        
        return {
            "score": score,
            "matches": matches,
            "length": len(aligned_seq1),
            "identity": (matches / min(len(fly_seq), len(human_seq))) * 100
        }
    
    return None

def main():
    print("CG6051 vs Human LST2/ZFYVE28 Comparison")
    print("=" * 50)
    
    # Load sequences
    sequences = load_sequences()
    
    if "fly" not in sequences or "human" not in sequences:
        print("Error: Could not load both sequences")
        sys.exit(1)
    
    fly = sequences["fly"]
    human = sequences["human"]
    
    print(f"Fly CG6051: {fly['length']} aa")
    print(f"Human LST2: {human['length']} aa")
    print(f"Length difference: {fly['length'] - human['length']} aa")
    
    results = {
        "fly_length": fly['length'],
        "human_length": human['length'],
        "analyses": {}
    }
    
    # 1. FYVE domain comparison
    fyve_result = compare_fyve_domains(fly['seq'], human['seq'])
    if fyve_result:
        results["analyses"]["fyve_alignment"] = {
            "identity": fyve_result['identity'],
            "similarity": fyve_result['similarity']
        }
    
    # 2. TOS motif search
    tos_results = search_tos_in_alignment(fly['seq'], human['seq'])
    results["analyses"]["tos_motifs"] = tos_results
    
    # 3. Overall alignment
    overall_result = overall_alignment(fly['seq'], human['seq'])
    if overall_result:
        results["analyses"]["overall_alignment"] = overall_result
        print(f"\nOverall identity: {overall_result['identity']:.1f}%")
    
    # Save results
    with open("alignment_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 50)
    print("Analysis complete. Results saved to alignment_results.json")

if __name__ == "__main__":
    main()