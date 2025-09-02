#!/usr/bin/env python3
"""
Check AlphaFold structure predictions for protein transmembrane topology.

This script fetches and analyzes AlphaFold predictions to assess:
- pLDDT confidence scores
- High confidence regions
- Membrane topology predictions
"""

import json
import requests
import click
from typing import Dict, List, Tuple
from pathlib import Path
from Bio import SeqIO


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


@click.command()
@click.argument('fasta_file', type=click.Path(exists=True))
@click.option('--signal-end', default=19, help='End position of signal peptide (default: 19)')
@click.option('--check-region-start', default=20, help='Start position to check for TM helix (default: 20)')
@click.option('--check-region-length', default=20, help='Length of region to check (default: 20)')
@click.option('--output', '-o', help='Output JSON file (optional)')
def main(fasta_file, signal_end, check_region_start, check_region_length, output):
    """
    Analyze AlphaFold predictions for a protein.
    
    FASTA_FILE: Path to the input FASTA file containing the protein sequence
    """
    
    # Load sequence from FASTA file
    with open(fasta_file, 'r') as f:
        record = next(SeqIO.parse(f, "fasta"))
    
    sequence = str(record.seq)
    
    # Extract UniProt ID from header
    header = record.description
    if '|' in header:
        parts = header.split('|')
        uniprot_id = parts[0]
    else:
        uniprot_id = record.id
    
    protein_name = record.description.split('|')[1] if '|' in record.description else record.id
    
    click.echo(f"=== AlphaFold Analysis for {protein_name} ({uniprot_id}) ===\n")
    
    # Check post-signal peptide region for TM helix
    click.echo(f"1. Checking for transmembrane helix after signal peptide (residues {check_region_start}-{check_region_start+check_region_length-1}):")
    tm_check = check_membrane_features(sequence, start=check_region_start, window=check_region_length)
    click.echo(f"   Region: {tm_check.get('region', 'N/A')}")
    click.echo(f"   Sequence: {tm_check.get('sequence', 'N/A')}")
    click.echo(f"   Hydrophobicity: {tm_check.get('hydrophobicity', 0):.1%}")
    click.echo(f"   Likely TM helix: {tm_check.get('likely_tm_helix', False)}")
    
    # Try to fetch AlphaFold data
    click.echo("\n2. Fetching AlphaFold confidence data...")
    af_data = fetch_alphafold_data(uniprot_id)
    
    plddt_analysis = None
    if "error" in af_data:
        click.echo(f"   Note: Could not fetch live AlphaFold data: {af_data['error']}")
        click.echo(f"   Please check https://alphafold.ebi.ac.uk/entry/{uniprot_id} directly")
    else:
        # Analyze pLDDT scores if available
        if 'confidenceScore' in af_data:
            plddt_analysis = analyze_plddt_regions(af_data['confidenceScore'])
            click.echo(f"   Mean pLDDT: {plddt_analysis['mean_plddt']:.1f}")
            click.echo(f"   High confidence regions: {len(plddt_analysis['high_confidence_regions'])}")
            click.echo(f"   Very high confidence residues: {plddt_analysis['very_high_conf_residues']}")
            click.echo(f"   Confident residues: {plddt_analysis['confident_residues']}")
            click.echo(f"   Low confidence residues: {plddt_analysis['low_conf_residues']}")
            click.echo(f"   Very low confidence residues: {plddt_analysis['very_low_conf_residues']}")
    
    click.echo("\n3. Analysis Summary:")
    click.echo(f"   - Predicted signal peptide region: residues 1-{signal_end}")
    click.echo(f"   - Post-signal region ({check_region_start}-{check_region_start+check_region_length-1}) hydrophobicity: {tm_check.get('hydrophobicity', 0):.1%}")
    click.echo(f"   - Transmembrane helix likelihood: {'Low' if not tm_check.get('likely_tm_helix', False) else 'High'}")
    
    if plddt_analysis:
        click.echo(f"   - Overall structure confidence: {plddt_analysis['mean_plddt']:.1f} (mean pLDDT)")
    
    click.echo("   - Topology prediction requires additional evidence")
    
    # Output structured results
    results = {
        "uniprot_id": uniprot_id,
        "protein_name": protein_name,
        "sequence_length": len(sequence),
        "signal_peptide": f"1-{signal_end}",
        "post_signal_tm_check": tm_check,
        "alphafold_url": f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}",
        "plddt_analysis": plddt_analysis if plddt_analysis else None,
        "analysis": "Membrane topology prediction based on hydrophobicity and AlphaFold confidence"
    }
    
    if output:
        with open(output, 'w') as f:
            json.dump(results, f, indent=2)
        click.echo(f"\n4. Results saved to '{output}'")
    else:
        click.echo("\n4. JSON output:")
        click.echo(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()