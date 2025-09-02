# DeepTMHMM Analysis for LRX-1

## Overview

DeepTMHMM is a state-of-the-art deep learning tool for predicting:
- Signal peptides (SP)
- Transmembrane helices (TM)
- Protein topology (inside/outside/membrane)
- Overall architecture (Globular, SP+Globular, TM, SP+TM)

## Sequence Analysis Observations

### Signal Peptide Region (positions 1-19)
- **Sequence**: MAWLTSIFFILLAVQPVLP
- **Characteristics**: Hydrophobic residues typical of signal peptides
- **Analysis**: Consistent with ER-targeting signal sequence

### Post-Signal Region (positions 20-40)
- **Sequence**: QDLYGTATQQQPYPYVQPSA
- **Hydrophobicity**: 35% (7/20 hydrophobic residues)
- **Characteristics**: Contains multiple polar/charged residues (Q, D, Y, T)
- **Analysis**: Low hydrophobicity, not typical of transmembrane helices

## Topology Predictions and Interpretations

### If predicted as "SP+Globular":
- Indicates signal peptide followed by globular domain
- Characteristic of secreted or extracellular proteins
- Protein would be processed through ER/Golgi pathway
- Mature protein would be extracellular

### If predicted as "SP+TM":
- Indicates signal peptide followed by transmembrane helix
- Characteristic of Type I membrane proteins
- Would suggest membrane anchoring after signal peptide
- Specific TM region location and confidence scores should be examined

### If predicted as "TM" only:
- Would indicate transmembrane protein without signal peptide
- Characteristic of multi-pass membrane proteins
- Would require review of specific helix positions

### If predicted as "Globular":
- Would indicate cytoplasmic/nuclear protein
- No signal peptide or TM helices
- Would contradict N-terminal hydrophobic region analysis

## Additional Evidence

### Sequence-based observations:
1. **Cysteine distribution**: 30 cysteines (8.2% of sequence)
   - High cysteine content typical of extracellular proteins
   - Suggests disulfide bond formation

2. **Protein size**: 368 amino acids
   - Much smaller than typical LRP family proteins (4000-6000 aa)
   - Size consistent with secreted proteins

3. **AlphaFold predictions**:
   - No clear transmembrane helix structure
   - Low overall confidence (mean pLDDT: 54.17)

## Analysis Approach

When interpreting DeepTMHMM results:
1. Check the overall topology prediction
2. Examine confidence scores for each region
3. Review specific positions of predicted features
4. Compare with other prediction tools and experimental evidence
5. Consider biological context and known protein functions

## Data Files

- Input sequence: `lrx1_for_deeptmhmm.fasta`
- DeepTMHMM results: `biolib_results/predicted_topologies.3line`
- Parsed results: `deeptmhmm_parsed_results.json`

## Note on Interpretation

Topology predictions should be considered alongside:
- Experimental evidence from literature
- Other computational predictions
- Functional studies
- Cellular localization data

No single prediction tool should be considered definitive without supporting evidence.