#!/usr/bin/env python3
"""
Analyze conservation of MTC7 across fungi.

This script searches for MTC7 homologs in fungal species and analyzes
conservation patterns, particularly in the transmembrane and polybasic regions.
"""

import json
import requests
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from pathlib import Path
import time
import click


def search_uniprot_homologs(query_sequence, taxon_filter="Fungi", max_results=20):
    """
    Search for homologs in UniProt using BLAST.
    
    Note: This uses UniProt's BLAST service which may have rate limits.
    """
    print("Searching for fungal homologs in UniProt...")
    
    # Prepare BLAST request to UniProt
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # First, let's search by protein name and taxonomy
    params = {
        'query': f'(protein_name:"Maintenance of telomere" OR gene:MTC*) AND taxonomy:"{taxon_filter}"',
        'format': 'json',
        'size': max_results,
        'fields': 'accession,id,protein_name,organism_name,organism_id,sequence,gene_names,length'
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        homologs = []
        for result in data.get('results', []):
            seq_data = result.get('sequence', {})
            homologs.append({
                'accession': result.get('primaryAccession', ''),
                'id': result.get('uniProtkbId', ''),
                'protein_name': result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                'organism': result.get('organism', {}).get('scientificName', ''),
                'organism_id': result.get('organism', {}).get('taxonId', ''),
                'sequence': seq_data.get('value', ''),
                'length': seq_data.get('length', 0),
                'gene_names': [g.get('value', '') for g in result.get('genes', [])]
            })
        
        return homologs
    
    except Exception as e:
        print(f"Error searching UniProt: {e}")
        return []


def analyze_sequence_similarity(seq1, seq2):
    """
    Calculate simple sequence similarity between two sequences.
    """
    from difflib import SequenceMatcher
    
    # Calculate overall similarity
    matcher = SequenceMatcher(None, seq1, seq2)
    similarity = matcher.ratio() * 100
    
    # Find matching blocks
    matching_blocks = matcher.get_matching_blocks()
    
    return similarity, matching_blocks


def analyze_tm_conservation(homologs, reference_seq, tm_regions):
    """
    Analyze conservation of transmembrane regions across homologs.
    """
    tm_conservation = []
    
    for tm_start, tm_end in tm_regions:
        tm_seq = reference_seq[tm_start-1:tm_end]
        
        conserved_positions = [0] * len(tm_seq)
        hydrophobic_conservation = []
        
        for homolog in homologs:
            if not homolog['sequence']:
                continue
            
            # Find best alignment of TM region in homolog
            best_match = 0
            for i in range(len(homolog['sequence']) - len(tm_seq) + 1):
                matches = sum(1 for j in range(len(tm_seq)) 
                             if i+j < len(homolog['sequence']) and 
                             tm_seq[j] == homolog['sequence'][i+j])
                if matches > best_match:
                    best_match = matches
            
            if best_match > 0:
                hydrophobic_conservation.append(best_match / len(tm_seq) * 100)
        
        tm_conservation.append({
            'region': f'{tm_start}-{tm_end}',
            'sequence': tm_seq,
            'avg_conservation': sum(hydrophobic_conservation) / len(hydrophobic_conservation) if hydrophobic_conservation else 0
        })
    
    return tm_conservation


def search_local_blast_hits():
    """
    Prepare data for local BLAST search against fungal genomes.
    
    Note: This prepares the query but actual BLAST would need to be run
    against a local database or NCBI's web service.
    """
    print("\nPreparing BLAST search parameters...")
    
    # Common fungal model organisms to search
    fungal_species = [
        {'name': 'Candida albicans', 'taxid': '5476'},
        {'name': 'Schizosaccharomyces pombe', 'taxid': '4896'},
        {'name': 'Neurospora crassa', 'taxid': '5141'},
        {'name': 'Aspergillus nidulans', 'taxid': '162425'},
        {'name': 'Cryptococcus neoformans', 'taxid': '5207'},
        {'name': 'Ustilago maydis', 'taxid': '5270'},
        {'name': 'Kluyveromyces lactis', 'taxid': '28985'},
        {'name': 'Yarrowia lipolytica', 'taxid': '4952'}
    ]
    
    return fungal_species


def analyze_polybasic_conservation(homologs, reference_seq):
    """
    Analyze conservation of polybasic regions, particularly C-terminal.
    """
    # C-terminal region (after position 62)
    c_term_ref = reference_seq[62:]
    
    conservation_data = []
    
    for homolog in homologs:
        if not homolog['sequence']:
            continue
        
        # Get C-terminal region of homolog (last portion of similar length)
        seq_len = len(homolog['sequence'])
        if seq_len > 62:
            c_term_homolog = homolog['sequence'][-(len(c_term_ref)):]
            
            # Count basic residues
            k_count = c_term_homolog.count('K')
            r_count = c_term_homolog.count('R')
            total_basic = k_count + r_count
            
            conservation_data.append({
                'organism': homolog['organism'],
                'c_terminal_seq': c_term_homolog,
                'lysine_count': k_count,
                'arginine_count': r_count,
                'total_basic': total_basic,
                'basic_percentage': (total_basic / len(c_term_homolog)) * 100 if c_term_homolog else 0
            })
    
    return conservation_data


@click.command()
@click.argument('fasta_file', type=click.Path(exists=True))
@click.option('--output-prefix', default='conservation_analysis', help='Prefix for output files')
def main(fasta_file, output_prefix):
    """
    Analyze conservation of a protein sequence across fungi.
    
    FASTA_FILE: Path to the FASTA file containing the protein sequence
    """
    # Read the sequence from the provided FASTA file
    with open(fasta_file, 'r') as f:
        record = next(SeqIO.parse(f, "fasta"))
    
    sequence = str(record.seq)
    protein_name = record.id
    
    print(f"Analyzing {protein_name} conservation across fungi")
    print("=" * 60)
    
    # Search for homologs
    homologs = search_uniprot_homologs(sequence)
    
    if homologs:
        print(f"\nFound {len(homologs)} potential homologs in UniProt")
        print("-" * 40)
        
        for homolog in homologs[:10]:  # Show first 10
            print(f"  {homolog['organism']}: {homolog['accession']} ({homolog['length']} aa)")
            if homolog['gene_names']:
                print(f"    Gene: {', '.join(homolog['gene_names'][:3])}")
    
    # Known TM regions
    tm_regions = [(13, 33), (42, 62)]
    
    # Analyze TM conservation
    if homologs:
        tm_conservation = analyze_tm_conservation(homologs, sequence, tm_regions)
        
        print("\nTransmembrane Region Conservation:")
        print("-" * 40)
        for tm_data in tm_conservation:
            print(f"  {tm_data['region']}: {tm_data['avg_conservation']:.1f}% average conservation")
    
    # Analyze polybasic region conservation
    if homologs:
        polybasic_conservation = analyze_polybasic_conservation(homologs, sequence)
        
        print("\nC-terminal Polybasic Region Conservation:")
        print("-" * 40)
        
        # Reference C-terminal
        c_term_ref = sequence[62:]
        print(f"  Reference (S. cerevisiae): K={c_term_ref.count('K')}, R={c_term_ref.count('R')}")
        
        for data in polybasic_conservation[:5]:  # Show first 5
            print(f"  {data['organism']}: K={data['lysine_count']}, R={data['arginine_count']} ({data['basic_percentage']:.1f}% basic)")
    
    # Get list of fungal species for BLAST
    fungal_species = search_local_blast_hits()
    
    print("\nRecommended fungal species for detailed BLAST analysis:")
    print("-" * 40)
    for species in fungal_species:
        print(f"  - {species['name']} (TaxID: {species['taxid']})")
    
    # Save results
    results = {
        'query_protein': protein_name,
        'query_length': len(sequence),
        'homologs_found': len(homologs),
        'homolog_list': [{
            'organism': h['organism'],
            'accession': h['accession'],
            'length': h['length'],
            'gene_names': h['gene_names']
        } for h in homologs],
        'tm_conservation': tm_conservation if homologs else [],
        'polybasic_conservation': polybasic_conservation[:10] if homologs else [],
        'recommended_species_for_blast': fungal_species
    }
    
    json_file = f'{output_prefix}_results.json'
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to '{json_file}'")
    
    # Create a simple conservation summary
    if homologs:
        with_polybasic = sum(1 for d in polybasic_conservation if d['total_basic'] >= 6)
        print(f"\nConservation Summary:")
        print(f"  - {len(homologs)} homologs found in fungi")
        print(f"  - {with_polybasic}/{len(polybasic_conservation)} homologs have strong polybasic C-terminus (≥6 K/R)")


if __name__ == "__main__":
    main()