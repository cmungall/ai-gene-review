# lrx-1 Gene Research Notes

## Gene Overview
- **Gene Symbol**: lrx-1
- **Alternative Names**: T04H1.6, egf-5, CELE_T04H1.6
- **Organism**: *Caenorhabditis elegans*
- **UniProt ID**: Q22179 (TrEMBL)
- **Chromosome Location**: V:12252623-12255450
- **WormBase ID**: WBGene00003075

## Protein Structure
- **Protein Name**: LRP X(Cross)-hybridizing
- **Length**: 368 amino acids
- **Key Domains**: Contains four Class A low-density lipoprotein receptor (LDLR) domains
- **Membrane Topology**: Single-pass type I membrane protein
- **Subcellular Localization**: Endomembrane system, membrane-bound

## Molecular Function
The protein structure with four LDLR Class A domains suggests receptor or receptor-like function, likely involved in:
- Ligand binding
- Signal transduction
- Protein-protein interactions

## Protein Interactions
High-throughput yeast two-hybrid screens have identified interactions with:

1. **RSP-4**: A protein involved in mRNA splicing, embryo development, locomotion, and cell proliferation regulation [BioGRID interaction data]

2. **PFD-3**: A putative prefoldin subunit required for normal alpha-tubulin synthesis, microtubule growth, mitotic spindle formation and positioning, and early embryonic cell division [BioGRID interaction data]

## Related LRP Proteins in C. elegans
While specific functional data for lrx-1 is limited, related lipoprotein receptors provide context:

### LRP-1 and LRP-2
- Major lipoprotein receptors in C. elegans with well-characterized functions
- LRP-1 essential for molting and cuticle shedding [PMID:9876188 "mutations confer a striking defect, an inability to shed and degrade all of the old cuticle at each of the larval molts"]
- LRP-2 regulates glutamatergic signaling and GLR-1 receptor trafficking [PMID:37179968 "lrp-2 loss-of-function mutants have defects in glutamatergic mechanosensory nose-touch behavior"]
- Both involved in vulval development and Wnt signaling [PMID:32550404]

### Functional Context from LRP Family
- Lipoprotein receptors regulate EGL-17/FGF export during vulval development [PMID:14630941]
- Involved in endocytosis and protein trafficking pathways
- Play roles in developmental signaling and cell fate specification

## Current GO Annotations
Based on electronic annotations:
- **Cellular Component**: Endomembrane system (GO:0012505), Membrane (GO:0016020)
- **Biological Process**: Vesicle-mediated transport (GO:0016192)

## Potential Functions (Inferred)
Based on protein structure and interactions:
1. **Protein trafficking**: Interaction with PFD-3 (prefoldin) suggests role in protein folding/quality control
2. **RNA processing**: Interaction with RSP-4 suggests potential involvement in post-transcriptional regulation
3. **Signal transduction**: LDLR domains typically involved in ligand binding and signaling
4. **Membrane transport**: Localization and GO annotations support role in vesicular trafficking

## Research Gaps
- No direct experimental characterization of lrx-1 function
- No mutant phenotype data available
- Expression pattern not characterized
- Specific ligands or binding partners unknown beyond Y2H interactions
- Developmental or physiological role undefined

## Notes on Gene Discovery
- Originally identified in C. elegans genome sequencing project [PMID:9851916 "Genome sequence of the nematode C. elegans: a platform for investigating biology"]
- Named "LRP X(Cross)-hybridizing" suggesting it was identified through cross-hybridization with LRP probes
- Currently only TrEMBL (unreviewed) UniProt entry, not Swiss-Prot reviewed

## CRITICAL BIOINFORMATIC ANALYSIS (2024)

### Major Findings - Annotations are INCORRECT
1. **NOT an LRP family protein** - lacks β-propeller domains, EGF-like domains, wrong size (369 aa vs >4000 aa for true LRPs)
2. **NO functional LDL-A domains** - cysteine spacing incompatible with canonical LDL-A pattern
3. **NOT a membrane protein** - no strong transmembrane helix, likely secreted
4. **ARBA predictions unreliable** - based on false domain annotations

### Actual Protein Features
- 30 cysteines (8.1% of sequence) in non-canonical patterns
- Likely contains novel C. elegans-specific cysteine-rich domains
- Probable secreted protein of unknown function
- UniProt correctly notes "lacks conserved residues" for LDL domains