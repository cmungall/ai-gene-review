# PedH Cellular Localization Analysis Results

## Summary
Bioinformatics analysis confirms that PedH (Q88JH0) from Pseudomonas putida KT2440 is a **soluble periplasmic enzyme**, not membrane-associated.

## Key Findings

### 1. Signal Peptide Analysis
- **Present**: First ~25 amino acids form a classic Sec signal peptide
- **Cleavage site**: Around position 24-25 (AHA-AV)
- **Function**: Directs export to periplasm, then cleaved off

### 2. Transmembrane Region Analysis
- **No TM regions in mature protein**
- Hydrophobic regions detected at positions 6-27 are part of the signal peptide
- After signal peptide cleavage, no membrane-spanning helices remain
- Conclusion: PedH is SOLUBLE, not membrane-embedded

### 3. Protein Architecture
- **Mature protein**: ~570 amino acids (after signal peptide cleavage)
- **Domain**: Eight-bladed β-propeller structure typical of PQQ-ADHs
- **Cofactors**: PQQ binding site and lanthanide metal binding site
- **No membrane anchors**: Lacks transmembrane helices or lipid anchors

### 4. Comparison with Related Enzymes
All characterized PQQ-dependent alcohol dehydrogenases are soluble periplasmic:
- ExaA (P. aeruginosa): Soluble periplasmic
- PedE (P. putida): Soluble periplasmic  
- MxaF (methylotrophs): Soluble periplasmic
- XoxF (methylotrophs): Soluble periplasmic

### 5. Functional Implications
- PedH functions throughout the periplasmic space
- Freely diffusible, encounters substrates anywhere in periplasm
- Transfers electrons to cytochrome c, also in periplasm
- Not restricted to membrane boundaries

## GO Term Recommendation

**CORRECT**: GO:0042597 (periplasmic space)
- Accurately describes soluble periplasmic localization
- Consistent with functional data

**INCORRECT**: GO:0030288 (outer membrane-bounded periplasmic space)
- Implies membrane association that doesn't exist
- Too specific for a freely diffusible enzyme

## Methods
- Signal peptide prediction based on sequence analysis
- Hydrophobicity analysis (Kyte-Doolittle scale)
- Domain architecture from UniProt
- Comparative analysis with characterized PQQ-ADHs

## References
- UniProt Q88JH0
- Wehrmann et al. (2017) mBio - PMID:28655819
- Mückschel et al. (2012) Appl Environ Microbiol - PMID:23023748