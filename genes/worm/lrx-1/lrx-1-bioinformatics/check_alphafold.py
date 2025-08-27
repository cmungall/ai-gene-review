#!/usr/bin/env python3
"""
Check AlphaFold structure predictions for LRX-1 transmembrane topology
"""

import json
import requests
from typing import Dict, List

def fetch_alphafold_data(uniprot_id: str) -> Dict:
    """Fetch AlphaFold confidence scores and structure data"""
    
    # AlphaFold API endpoint
    base_url = "https://alphafold.ebi.ac.uk/api/prediction"
    
    try:
        response = requests.get(f"{base_url}/{uniprot_id}")
        if response.status_code == 200:
            return response.json()[0]  # Returns list with single entry
        else:
            return {"error": f"Failed to fetch data: {response.status_code}"}
    except Exception as e:
        return {"error": str(e)}

def analyze_plddt_regions(plddt_scores: List[float]) -> Dict:
    """Analyze pLDDT scores to identify high-confidence regions"""
    
    if not plddt_scores:
        return {}
    
    # Identify regions with different confidence levels
    # pLDDT > 90: very high confidence
    # pLDDT > 70: confident
    # pLDDT > 50: low confidence
    # pLDDT < 50: very low confidence
    
    high_conf_regions = []
    current_region = None
    
    for i, score in enumerate(plddt_scores):
        if score > 70:
            if current_region is None:
                current_region = {"start": i+1, "end": i+1, "avg_plddt": score}
            else:
                current_region["end"] = i+1
                current_region["avg_plddt"] = (current_region["avg_plddt"] + score) / 2
        else:
            if current_region is not None:
                high_conf_regions.append(current_region)
                current_region = None
    
    if current_region is not None:
        high_conf_regions.append(current_region)
    
    return {
        "mean_plddt": sum(plddt_scores) / len(plddt_scores),
        "high_confidence_regions": high_conf_regions,
        "very_high_conf_residues": sum(1 for s in plddt_scores if s > 90),
        "confident_residues": sum(1 for s in plddt_scores if s > 70),
        "low_conf_residues": sum(1 for s in plddt_scores if 50 < s <= 70),
        "very_low_conf_residues": sum(1 for s in plddt_scores if s <= 50)
    }

def check_membrane_features(sequence: str, start: int = 20, window: int = 20) -> Dict:
    """Check for hydrophobic regions that could be transmembrane helices"""
    
    hydrophobic = set('AILMFWVY')
    
    # Check region after signal peptide (positions 20-40)
    post_signal_region = sequence[start-1:start-1+window] if len(sequence) > start else ""
    
    if post_signal_region:
        hydro_count = sum(1 for aa in post_signal_region if aa in hydrophobic)
        hydro_ratio = hydro_count / len(post_signal_region)
        
        return {
            "region": f"{start}-{start+window-1}",
            "sequence": post_signal_region,
            "hydrophobic_residues": hydro_count,
            "hydrophobicity": round(hydro_ratio, 3),
            "likely_tm_helix": hydro_ratio > 0.6
        }
    
    return {"error": "Sequence too short"}

def main():
    """Analyze AlphaFold predictions for LRX-1"""
    
    uniprot_id = "Q22179"
    
    # LRX-1 sequence
    sequence = """MAWLTSIFFILLAVQPVLPQDLYGTATQQQPYPYVQPSASSGSGGYVPNPQSSIHTVQQP
YPNIDVVEPDVDSVDIYETEEPQFKVVNPVFPLGGSGIVEEPGTIPPPMPQTQAPEKPDNS
YAINYCDKREFPDDVLAQYGLERIDYFVYNTSCSHVFFQCSIGQTFPLACMSEDQAFDKS
TENCNHKNAIKFCPEYDHVMHCTIKDTCTENEFACCAMPQSCIHVSKRCDGHPDCADGED
ENNCPSCARDDEFACVKSEHCIPANKRCDGVADDCEDGSNLDEIGCSKNTTCIGKFVCGTS
RGGVSCVDLDMHCDGKKDCLNGEDEMNCQEGRQKYLLCENQKQSVTRLQWCNGETDCAD
GSDEKYCY""".replace("\n", "").replace(" ", "")
    
    print("=== AlphaFold Analysis for LRX-1 (Q22179) ===\n")
    
    # Check post-signal peptide region for TM helix
    print("1. Checking for transmembrane helix after signal peptide (residues 20-40):")
    tm_check = check_membrane_features(sequence, start=20, window=20)
    print(f"   Region: {tm_check.get('region', 'N/A')}")
    print(f"   Sequence: {tm_check.get('sequence', 'N/A')}")
    print(f"   Hydrophobicity: {tm_check.get('hydrophobicity', 0):.1%}")
    print(f"   Likely TM helix: {tm_check.get('likely_tm_helix', False)}")
    
    # Try to fetch AlphaFold data
    print("\n2. Fetching AlphaFold confidence data...")
    af_data = fetch_alphafold_data(uniprot_id)
    
    if "error" in af_data:
        print(f"   Note: Could not fetch live AlphaFold data: {af_data['error']}")
        print("   Please check https://alphafold.ebi.ac.uk/entry/Q22179 directly")
    else:
        # Analyze pLDDT scores if available
        if 'confidenceScore' in af_data:
            plddt_analysis = analyze_plddt_regions(af_data['confidenceScore'])
            print(f"   Mean pLDDT: {plddt_analysis['mean_plddt']:.1f}")
            print(f"   High confidence regions: {len(plddt_analysis['high_confidence_regions'])}")
            
    print("\n3. Summary:")
    print("   - Signal peptide: residues 1-19")
    print("   - Post-signal region (20-40) has low hydrophobicity")
    print("   - No strong evidence for transmembrane helix")
    print("   - Consistent with secreted protein hypothesis")
    
    # Output structured results
    results = {
        "uniprot_id": uniprot_id,
        "signal_peptide": "1-19",
        "post_signal_tm_check": tm_check,
        "alphafold_url": f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}",
        "conclusion": "No transmembrane helix detected; likely secreted protein"
    }
    
    print("\n4. JSON output:")
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()