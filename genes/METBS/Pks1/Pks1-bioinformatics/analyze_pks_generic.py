#!/usr/bin/env python3
"""
Generic polyketide synthase (PKS) analyzer
Analyzes PKS proteins for domain architecture, active sites, and biosynthetic potential
"""

import re
import json
import requests
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from io import StringIO
from collections import Counter
import numpy as np
import argparse
import sys
from pathlib import Path

def fetch_uniprot_sequence(uniprot_id):
    """Fetch protein sequence from UniProt"""
    # Try new REST API first
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code != 200:
        # Fallback to old API
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        response = requests.get(url)
    
    if response.status_code == 200 and response.text.strip():
        try:
            record = SeqIO.read(StringIO(response.text), "fasta")
            return str(record.seq), record.description
        except:
            return None, None
    return None, None

def read_fasta_file(filepath):
    """Read sequence from FASTA file"""
    with open(filepath, 'r') as f:
        record = SeqIO.read(f, "fasta")
        return str(record.seq), record.description

class PKSAnalyzer:
    def __init__(self, sequence, protein_id="PKS"):
        self.sequence = sequence.replace('\n', '').replace(' ', '')
        self.length = len(self.sequence)
        self.protein_id = protein_id
        self.domains = []
        self.active_sites = []
        
    def scan_for_domains(self) -> List[Dict]:
        """Scan sequence for PKS domains using conserved motifs."""
        domains = []
        
        # KS (Ketosynthase) domain - look for DTACSS motif
        ks_pattern = r'D[TS]ACS[SG]'
        for match in re.finditer(ks_pattern, self.sequence):
            pos = match.start()
            # KS domains are typically ~400-450 aa
            domains.append({
                "name": "KS",
                "motif_start": pos,
                "estimated_start": max(0, pos - 200),
                "estimated_end": min(len(self.sequence), pos + 250),
                "confidence": "high" if match.group() == "DTACSS" else "moderate",
                "function": "Ketosynthase - chain elongation",
                "color": "#4ECDC4"
            })
        
        # AT/MAT (Acyltransferase) domain - look for GHSLG motif
        at_pattern = r'GHS[LI]G'
        for match in re.finditer(at_pattern, self.sequence):
            pos = match.start()
            # AT domains are typically ~300 aa
            domains.append({
                "name": "AT/MAT",
                "motif_start": pos,
                "estimated_start": max(0, pos - 100),
                "estimated_end": min(len(self.sequence), pos + 200),
                "confidence": "high",
                "function": "Acyltransferase",
                "color": "#45B7D1"
            })
        
        # ACP (Acyl Carrier Protein) - look for DSL motif around phosphopantetheine site
        acp_pattern = r'D[SG]L'
        for match in re.finditer(acp_pattern, self.sequence):
            pos = match.start()
            # Check for surrounding context - ACPs are small (~80-100 aa)
            window = self.sequence[max(0, pos-40):min(len(self.sequence), pos+40)]
            if window.count('S') >= 3:  # ACPs are serine-rich
                domains.append({
                    "name": "ACP",
                    "motif_start": pos,
                    "estimated_start": max(0, pos - 40),
                    "estimated_end": min(len(self.sequence), pos + 40),
                    "confidence": "moderate",
                    "function": "Acyl carrier protein",
                    "color": "#FFEAA7"
                })
        
        # DH (Dehydratase) - look for HxxxGxxxxP motif
        dh_pattern = r'H...G....P'
        for match in re.finditer(dh_pattern, self.sequence):
            pos = match.start()
            domains.append({
                "name": "DH",
                "motif_start": pos,
                "estimated_start": max(0, pos - 100),
                "estimated_end": min(len(self.sequence), pos + 150),
                "confidence": "moderate",
                "function": "Dehydratase",
                "color": "#96CEB4"
            })
        
        # TE (Thioesterase) - look for GxSxG motif
        te_pattern = r'G.S.G'
        for match in re.finditer(te_pattern, self.sequence):
            pos = match.start()
            # Check if this is likely in C-terminal region
            if pos > len(self.sequence) * 0.7:  # TE usually at C-terminus
                domains.append({
                    "name": "TE",
                    "motif_start": pos,
                    "estimated_start": max(0, pos - 100),
                    "estimated_end": min(len(self.sequence), pos + 150),
                    "confidence": "moderate",
                    "function": "Thioesterase - chain release",
                    "color": "#DDA0DD"
                })
        
        # KR (Ketoreductase) - look for NADPH binding motif
        kr_pattern = r'G.GG.[GA]'
        for match in re.finditer(kr_pattern, self.sequence):
            pos = match.start()
            domains.append({
                "name": "KR",
                "motif_start": pos,
                "estimated_start": max(0, pos - 100),
                "estimated_end": min(len(self.sequence), pos + 200),
                "confidence": "low",
                "function": "Ketoreductase",
                "color": "#FFB6C1"
            })
        
        # Sort domains by position
        domains.sort(key=lambda x: x['motif_start'])
        self.domains = domains
        return domains
    
    def predict_active_sites(self) -> List[Dict]:
        """Predict active sites based on domain positions and conserved residues."""
        active_sites = []
        
        for domain in self.domains:
            if domain['name'] == 'KS':
                # KS active site cysteine is typically ~25 aa before DTACSS
                cys_pos = domain['motif_start'] - 25
                if cys_pos > 0 and self.sequence[cys_pos] == 'C':
                    active_sites.append({
                        "position": cys_pos + 1,  # 1-indexed
                        "residue": "C",
                        "domain": "KS",
                        "function": "Catalytic cysteine",
                        "confidence": "high"
                    })
                
                # Histidines in KS domain
                for i in range(domain['estimated_start'], min(domain['estimated_end'], len(self.sequence))):
                    if self.sequence[i] == 'H':
                        active_sites.append({
                            "position": i + 1,
                            "residue": "H",
                            "domain": "KS",
                            "function": "Catalytic histidine",
                            "confidence": "moderate"
                        })
                        break
            
            elif domain['name'] == 'AT/MAT':
                # AT active site serine in GHSLG motif
                ser_pos = domain['motif_start'] + 2  # S in GHSLG
                if ser_pos < len(self.sequence):
                    active_sites.append({
                        "position": ser_pos + 1,
                        "residue": "S",
                        "domain": "AT/MAT",
                        "function": "Catalytic serine",
                        "confidence": "high"
                    })
            
            elif domain['name'] == 'ACP':
                # ACP phosphopantetheine attachment site
                # Look for serine near DSL motif
                for i in range(max(0, domain['motif_start'] - 10), 
                              min(len(self.sequence), domain['motif_start'] + 10)):
                    if self.sequence[i] == 'S':
                        active_sites.append({
                            "position": i + 1,
                            "residue": "S",
                            "domain": "ACP",
                            "function": "Phosphopantetheine attachment",
                            "confidence": "moderate"
                        })
                        break
            
            elif domain['name'] == 'TE':
                # TE active site serine in GxSxG
                ser_pos = domain['motif_start'] + 2
                if ser_pos < len(self.sequence):
                    active_sites.append({
                        "position": ser_pos + 1,
                        "residue": "S",
                        "domain": "TE",
                        "function": "Catalytic serine",
                        "confidence": "moderate"
                    })
        
        self.active_sites = active_sites
        return active_sites
    
    def analyze_composition(self) -> Dict:
        """Analyze amino acid composition and physicochemical properties."""
        aa_count = Counter(self.sequence)
        total = sum(aa_count.values())
        
        # Calculate properties
        hydrophobic = sum(aa_count[aa] for aa in 'AILMFWV')
        charged = sum(aa_count[aa] for aa in 'DEKR')
        polar = sum(aa_count[aa] for aa in 'STNQCY')
        aromatic = sum(aa_count[aa] for aa in 'FWY')
        
        # Molecular weight calculation (approximate)
        mw_dict = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'E': 147, 'Q': 146,
                   'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165,
                   'P': 115, 'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117}
        
        mw = sum(aa_count[aa] * mw_dict.get(aa, 0) for aa in aa_count) - 18 * (total - 1)
        
        return {
            "length": total,
            "molecular_weight_da": mw,
            "molecular_weight_kda": round(mw / 1000, 1),
            "hydrophobic_percent": round(100 * hydrophobic / total, 2),
            "charged_percent": round(100 * charged / total, 2),
            "polar_percent": round(100 * polar / total, 2),
            "aromatic_percent": round(100 * aromatic / total, 2),
            "most_common": aa_count.most_common(5)
        }
    
    def classify_pks_type(self) -> Dict:
        """Classify the PKS type based on domain architecture."""
        domain_names = [d['name'] for d in self.domains]
        unique_domains = list(dict.fromkeys(domain_names))  # Preserve order
        
        classification = {
            "domain_count": len(self.domains),
            "unique_domains": unique_domains,
            "domain_order": domain_names
        }
        
        # Type I PKS characteristics
        if len(self.sequence) > 1500 and 'KS' in domain_names and 'AT/MAT' in domain_names:
            if 'KR' in domain_names or 'DH' in domain_names:
                classification["type"] = "Type I reducing PKS"
                classification["subtype"] = "Modular"
            else:
                classification["type"] = "Type I non-reducing PKS"
                classification["subtype"] = "Iterative"
            
            # Check for ACP domains
            acp_count = domain_names.count('ACP')
            if acp_count > 1:
                classification["acp_count"] = acp_count
                classification["note"] = f"Multiple ACP domains ({acp_count}) suggest iterative processing"
        
        # Type II PKS (smaller, separate proteins)
        elif len(self.sequence) < 500 and len(unique_domains) <= 2:
            classification["type"] = "Possible Type II PKS component"
            classification["subtype"] = "Dissociated"
        
        # Type III PKS (no ACP)
        elif 'ACP' not in domain_names and 'KS' in domain_names:
            classification["type"] = "Possible Type III PKS"
            classification["subtype"] = "ACP-independent"
        
        else:
            classification["type"] = "Unclassified PKS"
            classification["subtype"] = "Unknown"
        
        return classification
    
    def predict_product_class(self) -> Dict:
        """Predict the class of polyketide product based on domains."""
        predictions = {
            "reducing_domains": [],
            "cyclization_potential": False,
            "product_predictions": []
        }
        
        domain_names = [d['name'] for d in self.domains]
        
        # Check for reducing domains
        if 'KR' in domain_names:
            predictions["reducing_domains"].append("KR (ketoreductase)")
        if 'DH' in domain_names:
            predictions["reducing_domains"].append("DH (dehydratase)")
        if 'ER' in domain_names:
            predictions["reducing_domains"].append("ER (enoylreductase)")
        
        # Predict product type
        if len(predictions["reducing_domains"]) == 0:
            predictions["product_predictions"].append("Aromatic polyketide (no reduction)")
            predictions["cyclization_potential"] = True
        elif len(predictions["reducing_domains"]) >= 2:
            predictions["product_predictions"].append("Highly reduced polyketide")
            predictions["product_predictions"].append("Possible fatty acid-like product")
        else:
            predictions["product_predictions"].append("Partially reduced polyketide")
        
        # Check for specific product indicators
        if 'TE' in domain_names:
            predictions["product_predictions"].append("Macrocyclic product possible (TE-mediated)")
        
        if domain_names.count('ACP') > 1:
            predictions["product_predictions"].append("Complex iterative product")
        
        return predictions
    
    def visualize_architecture(self, output_file=None):
        """Create a visual representation of domain architecture."""
        if not self.domains:
            print("No domains found to visualize")
            return
        
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Draw protein backbone
        ax.add_patch(patches.Rectangle((0, 0.45), self.length, 0.1, 
                                       facecolor='lightgray', edgecolor='black'))
        
        # Draw domains
        for domain in self.domains:
            start = domain['estimated_start']
            width = domain['estimated_end'] - start
            color = domain.get('color', '#808080')
            
            # Draw domain box
            rect = patches.Rectangle((start, 0.3), width, 0.4,
                                    facecolor=color, edgecolor='black', alpha=0.7)
            ax.add_patch(rect)
            
            # Add domain label
            ax.text(start + width/2, 0.5, domain['name'],
                   ha='center', va='center', fontsize=8, fontweight='bold')
            
            # Mark motif position
            ax.plot(domain['motif_start'], 0.5, 'r*', markersize=8)
        
        # Mark active sites
        for site in self.active_sites:
            if site['confidence'] == 'high':
                ax.plot(site['position'], 0.2, 'ro', markersize=6)
                ax.text(site['position'], 0.1, f"{site['residue']}{site['position']}", 
                       ha='center', fontsize=6, rotation=45)
        
        # Formatting
        ax.set_xlim(-50, self.length + 50)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Position (aa)', fontsize=12)
        ax.set_title(f'PKS Domain Architecture - {self.protein_id}', fontsize=14, fontweight='bold')
        ax.set_yticks([])
        
        # Add legend
        legend_elements = []
        seen_domains = set()
        for domain in self.domains:
            if domain['name'] not in seen_domains:
                seen_domains.add(domain['name'])
                color = domain.get('color', '#808080')
                legend_elements.append(patches.Patch(facecolor=color, alpha=0.7,
                                                    label=f"{domain['name']}: {domain['function']}"))
        
        ax.legend(handles=legend_elements, loc='upper center', 
                 bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=8)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Domain architecture saved to {output_file}")
        else:
            plt.show()
        
        plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Analyze polyketide synthase proteins for domain architecture and biosynthetic potential',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze from UniProt ID
  %(prog)s --uniprot Q4WF66
  
  # Analyze from FASTA file
  %(prog)s --fasta pks_protein.fasta
  
  # Save visualization
  %(prog)s --uniprot Q4WF66 --plot pks_domains.png
  
  # Save results to JSON
  %(prog)s --uniprot Q4WF66 --output results.json
  
  # Quiet mode
  %(prog)s --uniprot Q4WF66 --quiet --output results.json
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--uniprot', '-u', help='UniProt ID to analyze')
    input_group.add_argument('--fasta', '-f', help='FASTA file to analyze')
    
    parser.add_argument('--output', '-o', help='Output JSON file for results')
    parser.add_argument('--plot', '-p', help='Output file for domain architecture plot')
    parser.add_argument('--quiet', '-q', action='store_true', 
                       help='Suppress verbose output')
    
    args = parser.parse_args()
    
    # Get sequence
    if args.uniprot:
        sequence, description = fetch_uniprot_sequence(args.uniprot)
        if not sequence:
            print(f"Error: Failed to fetch sequence for UniProt ID {args.uniprot}", 
                  file=sys.stderr)
            sys.exit(1)
        protein_id = args.uniprot
    else:
        try:
            sequence, description = read_fasta_file(args.fasta)
            protein_id = Path(args.fasta).stem
        except Exception as e:
            print(f"Error reading FASTA file: {e}", file=sys.stderr)
            sys.exit(1)
    
    verbose = not args.quiet
    
    if verbose:
        print("=" * 60)
        print("Polyketide Synthase Analysis")
        print(f"Protein: {protein_id}")
        if description:
            print(f"Description: {description[:80]}...")
        print("=" * 60)
    
    # Initialize analyzer
    analyzer = PKSAnalyzer(sequence, protein_id)
    
    # Perform analyses
    if verbose:
        print(f"\nSequence length: {analyzer.length} amino acids")
        print("\nScanning for PKS domains...")
    
    domains = analyzer.scan_for_domains()
    
    if verbose:
        print(f"Found {len(domains)} potential domain regions")
        for domain in domains:
            print(f"  {domain['name']} at {domain['estimated_start']}-{domain['estimated_end']} "
                  f"(confidence: {domain['confidence']})")
    
    # Predict active sites
    if verbose:
        print("\nPredicting active sites...")
    
    active_sites = analyzer.predict_active_sites()
    
    if verbose and active_sites:
        print(f"Predicted {len(active_sites)} active sites:")
        for site in active_sites:
            if site['confidence'] == 'high':
                print(f"  {site['residue']}{site['position']} in {site['domain']} - {site['function']}")
    
    # Analyze composition
    composition = analyzer.analyze_composition()
    
    if verbose:
        print(f"\nProtein composition:")
        print(f"  Molecular weight: {composition['molecular_weight_kda']} kDa")
        print(f"  Hydrophobic: {composition['hydrophobic_percent']}%")
        print(f"  Charged: {composition['charged_percent']}%")
        print(f"  Aromatic: {composition['aromatic_percent']}%")
    
    # Classify PKS type
    classification = analyzer.classify_pks_type()
    
    if verbose:
        print(f"\nPKS Classification:")
        print(f"  Type: {classification['type']}")
        print(f"  Subtype: {classification['subtype']}")
        print(f"  Domain order: {' -> '.join(classification['domain_order'][:5])}...")
    
    # Predict product
    product_pred = analyzer.predict_product_class()
    
    if verbose:
        print(f"\nProduct predictions:")
        for pred in product_pred['product_predictions']:
            print(f"  - {pred}")
    
    # Create visualization if requested
    if args.plot:
        if verbose:
            print(f"\nGenerating domain architecture plot...")
        analyzer.visualize_architecture(args.plot)
    
    # Compile all results
    all_results = {
        'protein_id': protein_id,
        'description': description,
        'sequence_length': analyzer.length,
        'domains': domains,
        'active_sites': active_sites,
        'composition': composition,
        'classification': classification,
        'product_predictions': product_pred
    }
    
    # Save results if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(all_results, f, indent=2)
        if verbose:
            print(f"\nResults saved to {args.output}")
    
    if verbose:
        print("\n" + "=" * 60)
        print("Analysis complete!")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())