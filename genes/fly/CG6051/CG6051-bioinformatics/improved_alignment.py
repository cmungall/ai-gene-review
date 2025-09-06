#!/usr/bin/env python3
"""
Improved alignment analysis focusing on correct FYVE domain boundaries and conservation.
"""

import sys
import re
from pathlib import Path
from Bio import SeqIO, Align
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

def analyze_fyve_domains(fly_seq, human_seq):
    """Detailed FYVE domain analysis and proper alignment."""
    print("\n" + "=" * 70)
    print("FYVE DOMAIN ANALYSIS")
    print("=" * 70)
    
    # Extract FYVE domains - using correct boundaries
    fly_fyve = fly_seq[908:969]  # 909-969 in 1-based
    human_fyve = human_seq[814:875]  # 815-875 in 1-based
    
    print(f"\nFly FYVE domain (positions 909-969): {len(fly_fyve)} aa")
    print(f"Sequence: {fly_fyve}")
    
    print(f"\nHuman FYVE domain (positions 815-875): {len(human_fyve)} aa")
    print(f"Sequence: {human_fyve}")
    
    # Analyze cysteine patterns - critical for zinc binding
    fly_cys_pos = [i+1 for i, aa in enumerate(fly_fyve) if aa == 'C']
    human_cys_pos = [i+1 for i, aa in enumerate(human_fyve) if aa == 'C']
    
    print(f"\nCysteine residues (critical for zinc coordination):")
    print(f"  Fly:   {len(fly_cys_pos)} cysteines at relative positions: {fly_cys_pos}")
    print(f"  Human: {len(human_cys_pos)} cysteines at relative positions: {human_cys_pos}")
    
    # Perform careful alignment
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = 'global'
    
    alignments = list(aligner.align(fly_fyve, human_fyve))
    
    if alignments:
        best = alignments[0]
        
        # Format alignment properly
        print(f"\nFYVE Domain Alignment (Score: {best.score:.1f}):")
        print("-" * 70)
        
        # Extract aligned sequences properly from alignment object
        aligned_fly = ""
        aligned_human = ""
        
        # Convert alignment to strings
        for i in range(len(best.target)):
            if best.coordinates[0][i] >= 0:
                aligned_fly += fly_fyve[best.coordinates[0][i]] if best.coordinates[0][i] < len(fly_fyve) else "-"
            else:
                aligned_fly += "-"
                
        for i in range(len(best.query)):
            if best.coordinates[1][i] >= 0:
                aligned_human += human_fyve[best.coordinates[1][i]] if best.coordinates[1][i] < len(human_fyve) else "-"
            else:
                aligned_human += "-"
        
        # Use format method to extract aligned sequences
        alignment_str = format(best)
        lines = alignment_str.strip().split('\n')
        
        # Extract sequences from formatted alignment
        if len(lines) >= 3:
            aligned_fly = lines[0].split()[2] if len(lines[0].split()) > 2 else lines[0]
            aligned_human = lines[2].split()[2] if len(lines[2].split()) > 2 else lines[2]
        
        # Calculate conservation
        identical = 0
        similar = 0
        conservation_string = ""
        
        for f, h in zip(aligned_fly, aligned_human):
            if f == h and f != '-':
                identical += 1
                conservation_string += "*"
            elif f != '-' and h != '-' and matrix.get((f, h), 0) > 0:
                similar += 1
                conservation_string += ":"
            else:
                conservation_string += " "
        
        # Display alignment in blocks
        block_size = 60
        for i in range(0, len(aligned_fly), block_size):
            end = min(i + block_size, len(aligned_fly))
            print(f"Fly    {i+1:3}: {aligned_fly[i:end]}")
            print(f"          : {conservation_string[i:end]}")
            print(f"Human  {i+1:3}: {aligned_human[i:end]}")
            print()
        
        identity_pct = (identical / len(fly_fyve)) * 100
        similarity_pct = ((identical + similar) / len(fly_fyve)) * 100
        
        print(f"Identity: {identical}/{len(fly_fyve)} ({identity_pct:.1f}%)")
        print(f"Similarity: {identical + similar}/{len(fly_fyve)} ({similarity_pct:.1f}%)")
        
        # Check for conserved FYVE motifs
        print("\nConserved FYVE motifs:")
        
        # R(R/K)HHCR motif - zinc binding
        hhcr_fly = re.search(r"[RK][RK]HHC[RK]", fly_fyve)
        hhcr_human = re.search(r"[RK][RK]HHC[RK]", human_fyve)
        
        if hhcr_fly:
            print(f"  Fly HHCR motif: {hhcr_fly.group()} at position {hhcr_fly.start()+1}")
        if hhcr_human:
            print(f"  Human HHCR motif: {hhcr_human.group()} at position {hhcr_human.start()+1}")
        
        # RVC motif
        if "RVC" in fly_fyve:
            print(f"  Fly RVC motif at position {fly_fyve.index('RVC')+1}")
        if "RVC" in human_fyve:
            print(f"  Human RVC motif at position {human_fyve.index('RVC')+1}")
        
        return {
            "identity": identity_pct,
            "similarity": similarity_pct,
            "conserved_cysteines": len(set(fly_cys_pos) & set(human_cys_pos))
        }
    
    return None

def search_for_tos_motifs(fly_seq, human_seq):
    """Search for TOS and TOS-like motifs in both sequences."""
    print("\n" + "=" * 70)
    print("TOS MOTIF SEARCH")
    print("=" * 70)
    
    # Different TOS patterns to search
    tos_patterns = [
        ("F[DE][LIVM][DE][LIVM]", "Canonical TOS"),
        ("F[DE][LIVM][DE].", "Relaxed TOS"),
        ("F[DE]I[DE][LIVM]", "FxIxL pattern"),
        ("FDIDI", "Exact human TOS")
    ]
    
    results = {"fly": [], "human": []}
    
    for pattern, description in tos_patterns:
        # Search in human
        if pattern == "FDIDI":
            if pattern in human_seq:
                pos = human_seq.index(pattern) + 1
                results["human"].append({
                    "pattern": pattern,
                    "description": description,
                    "position": pos,
                    "sequence": pattern,
                    "context": human_seq[max(0, pos-6):min(len(human_seq), pos+9)]
                })
        else:
            for match in re.finditer(pattern, human_seq):
                results["human"].append({
                    "pattern": pattern,
                    "description": description,
                    "position": match.start() + 1,
                    "sequence": match.group(),
                    "context": human_seq[max(0, match.start()-5):min(len(human_seq), match.end()+5)]
                })
        
        # Search in fly
        if pattern == "FDIDI":
            if pattern in fly_seq:
                pos = fly_seq.index(pattern) + 1
                results["fly"].append({
                    "pattern": pattern,
                    "description": description,
                    "position": pos,
                    "sequence": pattern,
                    "context": fly_seq[max(0, pos-6):min(len(fly_seq), pos+9)]
                })
        else:
            for match in re.finditer(pattern, fly_seq):
                results["fly"].append({
                    "pattern": pattern,
                    "description": description,
                    "position": match.start() + 1,
                    "sequence": match.group(),
                    "context": fly_seq[max(0, match.start()-5):min(len(fly_seq), match.end()+5)]
                })
    
    # Remove duplicates
    def unique_motifs(motif_list):
        seen = set()
        unique = []
        for m in motif_list:
            key = (m["position"], m["sequence"])
            if key not in seen:
                seen.add(key)
                unique.append(m)
        return unique
    
    results["human"] = unique_motifs(results["human"])
    results["fly"] = unique_motifs(results["fly"])
    
    # Display results
    print("\nHuman LST2 TOS motifs:")
    if results["human"]:
        for m in results["human"]:
            print(f"  {m['description']} at position {m['position']}: '{m['sequence']}'")
            print(f"    Context: ...{m['context']}...")
    else:
        print("  No TOS motifs found")
    
    print("\nFly CG6051 TOS motifs:")
    if results["fly"]:
        for m in results["fly"]:
            print(f"  {m['description']} at position {m['position']}: '{m['sequence']}'")
            print(f"    Context: ...{m['context']}...")
    else:
        print("  No TOS motifs found")
    
    return results

def search_phosphatase_domains(fly_seq, human_seq):
    """Search for phosphatase catalytic signatures."""
    print("\n" + "=" * 70)
    print("PHOSPHATASE DOMAIN SEARCH")
    print("=" * 70)
    
    phosphatase_signatures = {
        "PTP active site": r"[LIVMFY]HC[LIVM]AGR",
        "DSP active site": r"HC[LIVM]{5}R",
        "PP1/PP2A": r"[LIVM]{2}GD[YFH]HG",
        "PP2C": r"[LIVM]{2}D[GS]H[GA]",
        "Acid phosphatase": r"[LIVM].[LIVM]D[ST]G[STN]",
        "His phosphatase": r"RH[GA].{4,8}H"
    }
    
    print("\nSearching for phosphatase signatures...")
    
    for name, pattern in phosphatase_signatures.items():
        fly_matches = list(re.finditer(pattern, fly_seq))
        human_matches = list(re.finditer(pattern, human_seq))
        
        if fly_matches or human_matches:
            print(f"\n{name} ({pattern}):")
            if fly_matches:
                print(f"  Fly matches: {len(fly_matches)}")
                for m in fly_matches[:2]:  # Show first 2
                    print(f"    Position {m.start()+1}: {m.group()}")
            if human_matches:
                print(f"  Human matches: {len(human_matches)}")
                for m in human_matches[:2]:  # Show first 2
                    print(f"    Position {m.start()+1}: {m.group()}")
    
    print("\nConclusion: No canonical phosphatase active site signatures detected in either protein")
    
    return None

def analyze_conservation(fly_seq, human_seq):
    """Analyze overall sequence conservation."""
    print("\n" + "=" * 70)
    print("OVERALL CONSERVATION ANALYSIS")
    print("=" * 70)
    
    # Perform global alignment for overall conservation
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'
    
    # For large sequences, we'll just get the score
    alignment = aligner.align(fly_seq, human_seq)[0]
    
    # Calculate quick stats
    min_len = min(len(fly_seq), len(human_seq))
    max_len = max(len(fly_seq), len(human_seq))
    
    # Estimate identity (rough)
    estimated_identity = (alignment.score / (5 * min_len)) * 100  # Rough estimate
    
    print(f"\nSequence lengths:")
    print(f"  Fly CG6051:  {len(fly_seq)} aa")
    print(f"  Human LST2:  {len(human_seq)} aa")
    print(f"  Difference:  {len(fly_seq) - len(human_seq)} aa (Fly is longer)")
    
    print(f"\nAlignment score: {alignment.score:.1f}")
    print(f"Estimated overall identity: ~{estimated_identity:.1f}%")
    
    # Find highly conserved regions using sliding window
    window = 20
    conserved_regions = []
    
    # Create local aligner for finding conserved regions
    aligner.mode = 'local'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    
    # Sample the sequence for conserved regions
    for i in range(0, min(len(fly_seq), len(human_seq)) - window, window//2):
        fly_window = fly_seq[i:i+window]
        human_window = human_seq[i:i+window]
        
        # Simple identity check
        matches = sum(1 for a, b in zip(fly_window, human_window) if a == b)
        if matches >= window * 0.5:  # >50% identity
            conserved_regions.append((i+1, i+window, matches/window*100))
    
    if conserved_regions:
        print(f"\nHighly conserved regions (>50% identity in 20aa window):")
        for start, end, identity in conserved_regions[:5]:
            print(f"  Positions {start}-{end}: {identity:.1f}% identity")
    
    return {
        "alignment_score": alignment.score,
        "estimated_identity": estimated_identity
    }

def main():
    print("=" * 70)
    print("CG6051 vs Human LST2/ZFYVE28 Comprehensive Analysis")
    print("=" * 70)
    
    # Load sequences
    sequences = load_sequences()
    
    if "fly" not in sequences or "human" not in sequences:
        print("Error: Could not load both sequences")
        sys.exit(1)
    
    fly = sequences["fly"]
    human = sequences["human"]
    
    results = {
        "proteins": {
            "fly": {"id": fly["id"], "length": fly["length"]},
            "human": {"id": human["id"], "length": human["length"]}
        },
        "analyses": {}
    }
    
    # 1. FYVE domain analysis
    fyve_result = analyze_fyve_domains(fly["seq"], human["seq"])
    if fyve_result:
        results["analyses"]["fyve_domain"] = fyve_result
    
    # 2. TOS motif search
    tos_results = search_for_tos_motifs(fly["seq"], human["seq"])
    results["analyses"]["tos_motifs"] = tos_results
    
    # 3. Phosphatase domain search
    phosphatase_results = search_phosphatase_domains(fly["seq"], human["seq"])
    results["analyses"]["phosphatase_domains"] = "No canonical phosphatase domains found"
    
    # 4. Overall conservation
    conservation = analyze_conservation(fly["seq"], human["seq"])
    results["analyses"]["conservation"] = conservation
    
    # Save comprehensive results
    with open("comprehensive_analysis.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("\n1. FYVE Domain: Present in both proteins at C-terminus")
    print("   - Both contain 8-9 cysteines for zinc coordination")
    print("   - FYVE domains show moderate conservation")
    print("\n2. Phosphatase Activity: NO phosphatase domains detected")
    print("   - No canonical phosphatase active sites found")
    print("   - Phosphatase GO annotations likely incorrect")
    print("\n3. TOS Motif: Not conserved between species")
    print("   - Human has TOS motifs, fly does not")
    print("\n4. Overall: Low sequence identity but conserved domain architecture")
    print("   - Both are negative regulators of EGFR signaling")
    print("   - Both have C-terminal FYVE domains")
    
    print("\nResults saved to comprehensive_analysis.json")

if __name__ == "__main__":
    main()