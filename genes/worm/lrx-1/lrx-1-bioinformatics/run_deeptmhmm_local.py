#!/usr/bin/env python3
"""
Run DeepTMHMM locally using BioLib.

This script runs DeepTMHMM membrane topology prediction on a protein sequence
from a FASTA file. Requires BioLib Python package to be installed.
"""

import json
import click
from typing import Dict, Optional
from Bio import SeqIO
from pathlib import Path

# Check if biolib is installed
try:
    import biolib
    BIOLIB_AVAILABLE = True
except ImportError:
    BIOLIB_AVAILABLE = False


def run_deeptmhmm_analysis(sequence: str, sequence_id: str = "protein") -> Optional[Dict]:
    """
    Run DeepTMHMM on a protein sequence using BioLib
    
    Args:
        sequence: Protein sequence (amino acids)
        sequence_id: Identifier for the sequence
    
    Returns:
        Dictionary with prediction results or None if failed
    """
    
    if not BIOLIB_AVAILABLE:
        return {"error": "BioLib not installed. Run: pip install pybiolib"}
    
    try:
        # Load DeepTMHMM from BioLib
        print("Loading DeepTMHMM from BioLib...")
        deeptmhmm = biolib.load('DTU/DeepTMHMM')
        
        # Prepare input in FASTA format
        fasta_input = f">{sequence_id}\n{sequence}"
        
        # Run DeepTMHMM
        print("Running DeepTMHMM prediction...")
        # Write FASTA to temporary file for biolib
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(fasta_input)
            fasta_file = f.name
        
        try:
            # Run with file input
            result = deeptmhmm.cli(args=['--fasta', fasta_file])
            
            # Get output (use get_stdout() as per warning)
            if hasattr(result, 'get_stdout'):
                output = result.get_stdout()
            else:
                output = result.stdout
            
            # Parse the results
            parsed = parse_deeptmhmm_output(output)
            
        finally:
            # Clean up temp file
            import os
            if os.path.exists(fasta_file):
                os.unlink(fasta_file)
        
        return parsed
        
    except Exception as e:
        return {"error": f"Failed to run DeepTMHMM: {str(e)}"}


def parse_deeptmhmm_output(output: str) -> Dict:
    """
    Parse DeepTMHMM output
    
    Expected format varies but typically includes:
    - Topology prediction (GLOB, TM, SP+GLOB, SP+TM)
    - Per-residue predictions
    """
    
    lines = output.strip().split('\n')
    
    result = {
        "raw_output": output,
        "topology": None,
        "signal_peptide": None,
        "tm_helices": [],
        "regions": []
    }
    
    for line in lines:
        line = line.strip()
        
        # Look for topology prediction
        if 'GLOB' in line or 'TM' in line or 'SP' in line:
            if 'SP+GLOB' in line:
                result["topology"] = "SP+GLOBULAR (secreted)"
            elif 'SP+TM' in line:
                result["topology"] = "SP+TM (membrane)"
            elif 'GLOB' in line:
                result["topology"] = "GLOBULAR"
            elif 'TM' in line:
                result["topology"] = "TRANSMEMBRANE"
        
        # Look for signal peptide region
        if 'Signal peptide' in line or 'SP' in line:
            # Extract position info if available
            import re
            pos_match = re.search(r'(\d+)-(\d+)', line)
            if pos_match:
                result["signal_peptide"] = {
                    "start": int(pos_match.group(1)),
                    "end": int(pos_match.group(2))
                }
        
        # Look for TM regions
        if 'Transmembrane' in line or 'TM' in line:
            import re
            pos_match = re.search(r'(\d+)-(\d+)', line)
            if pos_match:
                result["tm_helices"].append({
                    "start": int(pos_match.group(1)),
                    "end": int(pos_match.group(2))
                })
    
    return result


@click.command()
@click.argument('fasta_file', type=click.Path(exists=True))
@click.option('--output', '-o', default='deeptmhmm_results.json',
              help='Output JSON file (default: deeptmhmm_results.json)')
@click.option('--raw-output', default='deeptmhmm_raw.txt',
              help='Raw output file (default: deeptmhmm_raw.txt)')
def main(fasta_file, output, raw_output):
    """
    Run DeepTMHMM membrane topology prediction on a protein sequence.
    
    FASTA_FILE: Path to the input FASTA file containing the protein sequence
    
    Requires BioLib to be installed:
        pip install pybiolib
    
    May require authentication:
        biolib login
    """
    
    # Check BioLib availability
    if not BIOLIB_AVAILABLE:
        click.echo("Error: BioLib not installed", err=True)
        click.echo("\nTo install BioLib:", err=True)
        click.echo("  pip install pybiolib", err=True)
        click.echo("\nYou may also need to authenticate:", err=True)
        click.echo("  biolib login", err=True)
        return 1
    
    # Read sequence from FASTA file
    with open(fasta_file, 'r') as f:
        record = next(SeqIO.parse(f, "fasta"))
    
    sequence = str(record.seq)
    
    # Extract protein info
    protein_id = record.id
    protein_name = record.description.split('|')[1] if '|' in record.description else record.id
    
    click.echo(f"=== DeepTMHMM Analysis for {protein_name} ({protein_id}) ===\n")
    click.echo(f"Sequence length: {len(sequence)} aa")
    
    # Run analysis
    result = run_deeptmhmm_analysis(sequence, protein_id)
    
    if result and "error" in result:
        click.echo(f"Error: {result['error']}", err=True)
        click.echo("\nTroubleshooting:", err=True)
        click.echo("1. Install BioLib: pip install pybiolib", err=True)
        click.echo("2. Authenticate if needed: biolib login", err=True)
        click.echo("3. Check internet connection", err=True)
        return 1
    
    # Display results
    click.echo("\n=== Results ===")
    click.echo(f"Topology Prediction: {result.get('topology', 'Unknown')}")
    
    if result.get("signal_peptide"):
        sp = result["signal_peptide"]
        click.echo(f"Signal Peptide: YES (positions {sp['start']}-{sp['end']})")
    else:
        click.echo("Signal Peptide: Not detected or not reported")
    
    if result.get("tm_helices"):
        click.echo(f"Transmembrane Helices: {len(result['tm_helices'])}")
        for i, tm in enumerate(result["tm_helices"], 1):
            click.echo(f"  TM{i}: positions {tm['start']}-{tm['end']}")
    else:
        click.echo("Transmembrane Helices: NONE")
    
    # Save results
    with open(output, 'w') as f:
        json.dump(result, f, indent=2)
    click.echo(f"\n✓ Results saved to {output}")
    
    # Save raw output if available
    if result.get("raw_output"):
        with open(raw_output, 'w') as f:
            f.write(result["raw_output"])
        click.echo(f"✓ Raw output saved to {raw_output}")
    
    # Interpretation
    click.echo("\n=== Topology Analysis ===")
    topology = result.get("topology", "Unknown")
    if topology == "SP+GLOBULAR (secreted)":
        click.echo("• Predicted topology: Signal peptide + Globular domain")
        click.echo("• Characteristic of: Secreted/extracellular proteins")
        click.echo("• Processing: Through ER/Golgi secretory pathway")
    elif "SP+TM" in str(topology):
        click.echo("• Predicted topology: Signal peptide + Transmembrane helix")
        click.echo("• Characteristic of: Type I membrane proteins")
        click.echo("• Localization: Membrane-anchored with extracellular domain")
    elif "TM" in str(topology):
        click.echo("• Predicted topology: Transmembrane protein")
        click.echo("• Contains transmembrane helix(es)")
        click.echo("• Review raw output for specific helix positions")
    else:
        click.echo(f"• Topology prediction: {topology}")
        click.echo("• Further analysis needed")
    
    return 0


if __name__ == "__main__":
    exit(main())