#!/usr/bin/env python3
"""
Run DeepTMHMM analysis on LRX-1 protein sequence
"""

import requests
import json
from typing import Dict, Optional

def run_deeptmhmm(sequence: str) -> Optional[Dict]:
    """
    Attempt to run DeepTMHMM via BioLib API
    
    Note: BioLib requires authentication for API access.
    This script documents the API structure but may require API keys.
    """
    
    # BioLib API endpoint for DeepTMHMM
    base_url = "https://api.biolib.com/v1"
    app_id = "DTU/DeepTMHMM"
    
    # Prepare the request
    headers = {
        "Content-Type": "application/json",
        # API key would go here if available
        # "Authorization": "Bearer YOUR_API_KEY"
    }
    
    payload = {
        "input_fasta": f">LRX1\n{sequence}",
        "format": "json"
    }
    
    try:
        # Attempt to submit job
        response = requests.post(
            f"{base_url}/apps/{app_id}/jobs",
            headers=headers,
            json=payload,
            timeout=30
        )
        
        if response.status_code == 401:
            return {
                "error": "API authentication required",
                "message": "BioLib API requires authentication. Cannot run automatically.",
                "manual_url": "https://dtu.biolib.com/DeepTMHMM"
            }
        elif response.status_code == 200:
            return response.json()
        else:
            return {
                "error": f"API request failed with status {response.status_code}",
                "response": response.text
            }
            
    except requests.exceptions.RequestException as e:
        return {
            "error": "Network request failed",
            "message": str(e)
        }

def parse_deeptmhmm_output(result: str) -> Dict:
    """
    Parse DeepTMHMM output format
    
    Expected format:
    # ID  Type  Positions
    LRX1  SP    1-19
    LRX1  GLOB  20-369
    """
    
    parsed = {
        "signal_peptide": None,
        "tm_helices": [],
        "topology": "unknown",
        "regions": []
    }
    
    if not result:
        return parsed
    
    lines = result.strip().split('\n')
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
            
        parts = line.split('\t')
        if len(parts) >= 3:
            region_type = parts[1]
            positions = parts[2]
            
            if region_type == 'SP':
                parsed["signal_peptide"] = positions
            elif region_type == 'TM':
                parsed["tm_helices"].append(positions)
            elif region_type in ['GLOB', 'GLOBULAR']:
                parsed["topology"] = "globular"
                
            parsed["regions"].append({
                "type": region_type,
                "positions": positions
            })
    
    # Determine overall topology
    if parsed["signal_peptide"] and not parsed["tm_helices"]:
        parsed["topology"] = "SP+Globular (secreted)"
    elif parsed["signal_peptide"] and parsed["tm_helices"]:
        parsed["topology"] = "SP+TM (membrane)"
    elif parsed["tm_helices"]:
        parsed["topology"] = "TM (membrane)"
    else:
        parsed["topology"] = "Globular"
    
    return parsed

def main():
    """Run DeepTMHMM analysis"""
    
    # LRX-1 sequence
    sequence = """MAWLTSIFFILLAVQPVLPQDLYGTATQQQPYPYVQPSASSGSGGYVPNPQSSIHTVQQP
YPNIDVVEPDVDSVDIYETEEPQFKVVNPVFPLGGSGIVEEPGTIPPPMPQTQAPEKPDNS
YAINYCDKREFPDDVLAQYGLERIDYFVYNTSCSHVFFQCSIGQTFPLACMSEDQAFDKS
TENCNHKNAIKFCPEYDHVMHCTIKDTCTENEFACCAMPQSCIHVSKRCDGHPDCADGED
ENNCPSCARDDEFACVKSEHCIPANKRCDGVADDCEDGSNLDEIGCSKNTTCIGKFVCGTS
RGGVSCVDLDMHCDGKKDCLNGEDEMNCQEGRQKYLLCENQKQSVTRLQWCNGETDCAD
GSDEKYCY""".replace("\n", "").replace(" ", "")
    
    print("=== DeepTMHMM Analysis for LRX-1 ===\n")
    
    # Try to run via API
    result = run_deeptmhmm(sequence)
    
    if result and "error" in result:
        print(f"Cannot run automatically: {result['error']}")
        print(f"Message: {result.get('message', '')}")
        if "manual_url" in result:
            print(f"\nManual submission required at: {result['manual_url']}")
        return
    
    # If we somehow got results, parse them
    if result:
        parsed = parse_deeptmhmm_output(result.get("output", ""))
        print("Results:")
        print(json.dumps(parsed, indent=2))
    else:
        print("No results obtained")

if __name__ == "__main__":
    main()