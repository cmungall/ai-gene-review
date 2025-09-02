#!/usr/bin/env python3
"""
Prepare FASTA file for DeepTMHMM analysis.

This script reads a protein sequence from a FASTA file and prepares it
for submission to the DeepTMHMM web service or local installation.
"""

import click
from Bio import SeqIO
from pathlib import Path


@click.command()
@click.argument('input_fasta', type=click.Path(exists=True))
@click.option('--output', '-o', default='deeptmhmm_input.fasta', 
              help='Output FASTA file for DeepTMHMM (default: deeptmhmm_input.fasta)')
@click.option('--max-length', default=5000, 
              help='Maximum sequence length for DeepTMHMM (default: 5000)')
def main(input_fasta, output, max_length):
    """
    Prepare protein sequence for DeepTMHMM analysis.
    
    INPUT_FASTA: Path to the input FASTA file containing the protein sequence
    
    DeepTMHMM predicts transmembrane topology including:
    - Signal peptides
    - Transmembrane helices
    - Inside/outside orientation
    
    After running this script, submit the output file to:
    https://dtu.biolib.com/DeepTMHMM
    """
    
    # Read the sequence from input FASTA
    with open(input_fasta, 'r') as f:
        record = next(SeqIO.parse(f, "fasta"))
    
    sequence = str(record.seq)
    
    # Extract protein info from header
    header = record.description
    if '|' in header:
        parts = header.split('|')
        protein_id = parts[0]
        protein_name = parts[1] if len(parts) > 1 else "Unknown"
        organism = parts[2] if len(parts) > 2 else "Unknown"
    else:
        protein_id = record.id
        protein_name = "Unknown"
        organism = "Unknown"
    
    click.echo(f"Processing: {protein_name} ({protein_id})")
    click.echo(f"Organism: {organism}")
    click.echo(f"Sequence length: {len(sequence)} aa")
    
    # Check sequence length
    if len(sequence) > max_length:
        click.echo(f"WARNING: Sequence length ({len(sequence)}) exceeds DeepTMHMM limit ({max_length})")
        click.echo(f"Truncating to first {max_length} amino acids")
        sequence = sequence[:max_length]
    
    # Create formatted header for DeepTMHMM
    # Format: >ID|Name|Organism
    formatted_header = f"{protein_id}|{protein_name}|{organism}"
    
    # Write output FASTA
    with open(output, 'w') as f:
        f.write(f">{formatted_header}\n")
        # Write sequence in chunks of 60 characters
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + '\n')
    
    click.echo(f"\nOutput written to: {output}")
    click.echo("\nNext steps:")
    click.echo("1. Go to https://dtu.biolib.com/DeepTMHMM")
    click.echo(f"2. Upload the file: {output}")
    click.echo("3. Run the analysis")
    click.echo("4. Download results to 'biolib_results/' directory")
    click.echo("5. Run 'python check_deeptmhmm.py' to parse the results")
    
    # Display first 30 amino acids for signal peptide preview
    if len(sequence) >= 30:
        click.echo(f"\nFirst 30 amino acids (potential signal peptide region):")
        click.echo(f"  {sequence[:30]}")
        
        # Quick hydrophobicity check
        hydrophobic = set('AILMFWVY')
        hydro_count = sum(1 for aa in sequence[:30] if aa in hydrophobic)
        click.echo(f"  Hydrophobicity: {hydro_count}/30 ({hydro_count/30:.1%})")


if __name__ == "__main__":
    main()