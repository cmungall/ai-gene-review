#!/usr/bin/env python3
"""
Fetch ortholog sequences for CFAP418 and perform conservation analysis
"""

import json
import time
from pathlib import Path
from typing import Dict, List

import requests
from Bio import Align, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
DATA_DIR.mkdir(exist_ok=True)

def fetch_uniprot_sequence(uniprot_id: str) -> str:
    """Fetch sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            sequence = ''.join(lines[1:])  # Skip header
            return sequence
        else:
            print(f"Failed to fetch {uniprot_id}: Status {response.status_code}")
            return ""
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
        return ""

def fetch_orthologs() -> Dict:
    """Fetch known CFAP418 orthologs from different species"""
    
    orthologs = {
        "human": {
            "uniprot_id": "Q96NL8",
            "organism": "Homo sapiens",
            "common_name": "Human",
            "sequence": ""
        },
        "mouse": {
            "uniprot_id": "Q8BXQ0",
            "organism": "Mus musculus", 
            "common_name": "Mouse",
            "sequence": ""
        },
        "zebrafish": {
            "uniprot_id": "F1QKW4",
            "organism": "Danio rerio",
            "common_name": "Zebrafish",
            "sequence": ""
        },
        "xenopus": {
            "uniprot_id": "F6YQE7",
            "organism": "Xenopus tropicalis",
            "common_name": "Western clawed frog",
            "sequence": ""
        }
    }
    
    # Fetch sequences
    print("Fetching ortholog sequences from UniProt...")
    for species, info in orthologs.items():
        print(f"  Fetching {species} ({info['organism']})...")
        
        # Try to load from local file first
        local_file = DATA_DIR / f"CFAP418_{species}.fasta"
        if local_file.exists():
            for record in SeqIO.parse(local_file, "fasta"):
                info["sequence"] = str(record.seq)
                print(f"    Loaded from local file")
        else:
            # Fetch from UniProt
            sequence = fetch_uniprot_sequence(info["uniprot_id"])
            if sequence:
                info["sequence"] = sequence
                # Save to file
                record = SeqRecord(
                    Seq(sequence),
                    id=info["uniprot_id"],
                    description=f"CFAP418 {info['organism']}"
                )
                SeqIO.write([record], local_file, "fasta")
                print(f"    Fetched and saved to {local_file}")
                time.sleep(0.5)  # Be polite to the server
            else:
                print(f"    Failed to fetch sequence")
    
    return orthologs

def align_sequences(seq1: str, seq2: str, seq1_name: str, seq2_name: str) -> Dict:
    """Perform pairwise alignment and calculate identity/similarity"""
    
    if not seq1 or not seq2:
        return {
            "identity": 0,
            "similarity": 0,
            "gaps": 0,
            "alignment_length": 0
        }
    
    # Create aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    # Perform alignment
    alignments = aligner.align(seq1, seq2)
    
    if alignments:
        alignment = alignments[0]
        
        # Calculate statistics
        aligned_seq1 = str(alignment[0])
        aligned_seq2 = str(alignment[1])
        
        identical = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
        similar = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                     if a != '-' and b != '-' and are_similar(a, b))
        gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
        alignment_length = len(aligned_seq1)
        
        return {
            "identity": (identical / alignment_length) * 100 if alignment_length > 0 else 0,
            "similarity": (similar / alignment_length) * 100 if alignment_length > 0 else 0,
            "gaps": gaps,
            "alignment_length": alignment_length,
            "score": alignment.score
        }
    
    return {
        "identity": 0,
        "similarity": 0,
        "gaps": 0,
        "alignment_length": 0,
        "score": 0
    }

def are_similar(aa1: str, aa2: str) -> bool:
    """Check if two amino acids are similar based on chemical properties"""
    if aa1 == aa2:
        return True
    
    # Define similarity groups
    groups = [
        set("GAVLI"),      # Small/hydrophobic
        set("ST"),         # Small polar
        set("CM"),         # Sulfur-containing
        set("FYWH"),       # Aromatic
        set("KRDE"),       # Charged
        set("NQ"),         # Amide
        set("P")           # Special
    ]
    
    for group in groups:
        if aa1 in group and aa2 in group:
            return True
    
    return False

def analyze_conservation(orthologs: Dict) -> Dict:
    """Analyze sequence conservation across species"""
    
    human_seq = orthologs["human"]["sequence"]
    
    if not human_seq:
        print("Error: Human sequence not available")
        return {}
    
    conservation_results = {
        "reference": "human",
        "alignments": {},
        "conserved_regions": [],
        "variable_regions": []
    }
    
    # Compare human with each ortholog
    for species, info in orthologs.items():
        if species == "human" or not info["sequence"]:
            continue
        
        print(f"\nAligning human vs {species}...")
        alignment_stats = align_sequences(
            human_seq, info["sequence"],
            "human", species
        )
        
        conservation_results["alignments"][species] = {
            "organism": info["organism"],
            "common_name": info["common_name"],
            **alignment_stats
        }
        
        print(f"  Identity: {alignment_stats['identity']:.1f}%")
        print(f"  Similarity: {alignment_stats['similarity']:.1f}%")
    
    # Identify highly conserved regions (simplified)
    # In a full implementation, would use multiple sequence alignment
    window_size = 20
    conservation_scores = []
    
    for i in range(len(human_seq) - window_size):
        window = human_seq[i:i+window_size]
        
        # Calculate conservation score for this window
        scores = []
        for species, info in orthologs.items():
            if species == "human" or not info["sequence"]:
                continue
            
            if i + window_size <= len(info["sequence"]):
                other_window = info["sequence"][i:i+window_size]
                identity = sum(1 for a, b in zip(window, other_window) if a == b)
                scores.append(identity / window_size)
        
        if scores:
            avg_conservation = sum(scores) / len(scores)
            conservation_scores.append({
                "start": i + 1,
                "end": i + window_size,
                "conservation": avg_conservation
            })
    
    # Find highly conserved regions
    high_threshold = 0.8
    for score in conservation_scores:
        if score["conservation"] >= high_threshold:
            conservation_results["conserved_regions"].append({
                "start": score["start"],
                "end": score["end"],
                "conservation_level": f"{score['conservation']*100:.1f}%"
            })
    
    return conservation_results

def identify_functional_motifs(sequence: str) -> List[Dict]:
    """Identify potential functional motifs in the sequence"""
    
    motifs = []
    
    # Check for common motifs
    # Phosphorylation sites
    phos_patterns = [
        ("PKA", r"R[RK].[ST]"),  # PKA consensus
        ("PKC", r"[ST].[RK]"),    # PKC consensus
        ("CK2", r"[ST]..E"),      # CK2 consensus
    ]
    
    # Nuclear localization signals
    nls_patterns = [
        ("NLS_mono", r"[KR]{4,}"),  # Monopartite NLS
        ("NLS_bi", r"[KR]{2}.{10,12}[KR]{3,}"),  # Bipartite NLS
    ]
    
    # Check for zinc finger motifs (important for some ciliary proteins)
    if sequence.count('C') >= 4:
        # Look for C2H2 zinc finger pattern
        import re
        c2h2_pattern = r"C.{2,4}C.{12}H.{3,5}H"
        matches = re.finditer(c2h2_pattern, sequence)
        for match in matches:
            motifs.append({
                "type": "Zinc_finger_C2H2",
                "start": match.start() + 1,
                "end": match.end(),
                "sequence": match.group()
            })
    
    return motifs

def main():
    """Main conservation analysis pipeline"""
    
    # Load human sequence
    human_fasta = BASE_DIR / "CFAP418.fasta"
    human_seq = ""
    if human_fasta.exists():
        for record in SeqIO.parse(human_fasta, "fasta"):
            human_seq = str(record.seq)
            break
    
    if not human_seq:
        print("Error: Human CFAP418 sequence not found")
        return
    
    # Fetch orthologs
    orthologs = fetch_orthologs()
    orthologs["human"]["sequence"] = human_seq
    
    # Analyze conservation
    print("\n" + "="*60)
    print("CONSERVATION ANALYSIS")
    print("="*60)
    
    conservation = analyze_conservation(orthologs)
    
    # Identify functional motifs
    print("\nSearching for functional motifs...")
    motifs = identify_functional_motifs(human_seq)
    
    # Compile results
    results = {
        "orthologs": {
            species: {
                "uniprot_id": info["uniprot_id"],
                "organism": info["organism"],
                "common_name": info["common_name"],
                "sequence_length": len(info["sequence"]) if info["sequence"] else 0
            }
            for species, info in orthologs.items()
        },
        "conservation_analysis": conservation,
        "functional_motifs": motifs,
        "summary": {
            "total_orthologs": len([o for o in orthologs.values() if o["sequence"]]),
            "average_identity_to_human": sum(
                conservation["alignments"][s]["identity"] 
                for s in conservation["alignments"]
            ) / len(conservation["alignments"]) if conservation.get("alignments") else 0,
            "highly_conserved_regions": len(conservation.get("conserved_regions", [])),
            "identified_motifs": len(motifs)
        }
    }
    
    # Save results
    output_file = BASE_DIR / "results" / "conservation_analysis.json"
    output_file.parent.mkdir(exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Orthologs analyzed: {results['summary']['total_orthologs']}")
    print(f"Average identity to human: {results['summary']['average_identity_to_human']:.1f}%")
    print(f"Highly conserved regions: {results['summary']['highly_conserved_regions']}")
    print(f"Functional motifs identified: {results['summary']['identified_motifs']}")
    
    return results

if __name__ == "__main__":
    results = main()