# DeepTMHMM Analysis for LRX-1

## Submission Instructions

1. Visit: https://dtu.biolib.com/DeepTMHMM
2. Paste the sequence from `lrx1_for_deeptmhmm.fasta`
3. Click "Run" to get predictions

## What DeepTMHMM Predicts

DeepTMHMM is the current state-of-the-art for predicting:
- Signal peptides (SP)
- Transmembrane helices (TM)
- Protein topology (inside/outside/membrane)
- Overall architecture (Globular, SP+Globular, TM, SP+TM)

## Expected Results for LRX-1

Based on our sequence analysis, we expect DeepTMHMM to predict:

### Signal Peptide
- **Expected**: YES (positions 1-19)
- **Sequence**: MAWLTSIFFILLAVQPVLP
- **Confidence**: Should be high (this is a typical signal peptide)

### Transmembrane Helices
- **Expected**: NONE
- **Reasoning**: 
  - Post-signal region (20-40) has low hydrophobicity
  - Sequence: QDLYGTATQQQPYPYVQPSA (contains Q, D, Y, T - hydrophilic)
  - No 20+ residue hydrophobic stretch

### Overall Topology
- **Expected**: SP+Globular (secreted protein)
- **Not expected**: SP+TM (single-pass membrane protein)

## Interpretation Guide

### If DeepTMHMM confirms "SP+Globular":
- LRX-1 is definitively a secreted protein
- Signal peptide directs to ER
- Protein is released to extracellular space
- All membrane-related GO annotations are incorrect

### If DeepTMHMM predicts "SP+TM":
- Check the specific TM region predicted
- Examine probability scores (low confidence?)
- Review the sequence of the predicted TM region
- Consider that even DeepTMHMM can be wrong for edge cases

## Why This Matters

DeepTMHMM uses deep learning trained on thousands of experimentally validated proteins. It's particularly good at:
- Distinguishing signal peptides from TM helices
- Predicting topology of multi-pass membrane proteins
- Identifying secreted vs membrane proteins

For LRX-1, this analysis will provide definitive evidence about whether UniProt's "single-pass type I membrane protein" annotation is correct or not.

## Current Evidence Summary

All analyses point to LRX-1 being secreted:
1. **Sequence analysis**: No hydrophobic region after signal peptide
2. **AlphaFold**: Low confidence, no TM helix visible
3. **Domain analysis**: No LRP domains, wrong size
4. **Cysteine pattern**: Consistent with extracellular protein

DeepTMHMM will provide the final confirmation.