# Rgn (Regucalcin/SMP30) GO Annotation Curation Summary

## Overview
Systematic review of 47 existing GO annotations for rat Rgn gene based on literature evidence and functional understanding.

## Key Curation Decisions

### Core Functions Identified

1. **Gluconolactonase Activity & Vitamin C Biosynthesis**
   - GO:0004341 (gluconolactonase activity) - ACCEPTED (multiple evidence codes)
   - GO:0019853 (L-ascorbic acid biosynthetic process) - ACCEPTED
   - Strong experimental evidence from PMID:16585534 showing this is the penultimate enzyme in vitamin C synthesis
   - Knockout mice develop scurvy, confirming essential role

2. **Calcium Binding & Homeostasis**
   - GO:0005509 (calcium ion binding) - ACCEPTED (multiple evidence codes)
   - GO:0006874 (intracellular calcium ion homeostasis) - ACCEPTED
   - GO:0050848 (regulation of calcium-mediated signaling) - ACCEPTED
   - Well-established since original 1978 discovery (PMID:699201)
   - Regulates Ca2+-ATPase activity (PMID:16786169)

3. **Macromolecule Biosynthesis Suppression**
   - GO:0010558 (negative regulation of macromolecule biosynthetic process) - ACCEPTED
   - GO:2000279 (negative regulation of DNA biosynthetic process) - ACCEPTED
   - GO:1902679 (negative regulation of RNA biosynthetic process) - ACCEPTED
   - Inhibits protein synthesis via aminoacyl-tRNA synthetase (PMID:2280766)
   - Suppresses DNA (PMID:11500948) and RNA (PMID:12397604) synthesis

4. **Cytoprotective/Anti-apoptotic Function**
   - GO:0043066 (negative regulation of apoptotic process) - ACCEPTED
   - GO:1903625 (negative regulation of DNA catabolic process) - ACCEPTED
   - Multiple studies show protection from apoptosis (PMID:15806309, PMID:16167335)

5. **Additional Core Functions**
   - GO:0008270 (zinc ion binding) - ACCEPTED (cofactor for enzymatic activity)
   - GO:0032781 (positive regulation of ATP-dependent activity) - MODIFIED to GO:0060590 (ATPase regulator activity)
   - GO:0050680 (negative regulation of epithelial cell proliferation) - ACCEPTED
   - GO:1903052 (positive regulation of proteolysis) - ACCEPTED
   - GO:0045019 (negative regulation of nitric oxide biosynthetic process) - ACCEPTED

### Non-Core/Secondary Functions (KEEP_AS_NON_CORE)

1. **Developmental Processes**
   - GO:0001822 (kidney development) - expression pattern, not core function
   - GO:0001889 (liver development) - expression pattern, not core function
   - GO:0097421 (liver regeneration) - context-specific response
   - GO:1903011 (negative regulation of bone development) - secondary effect of calcium dysregulation

2. **Metabolic Regulation**
   - GO:0010867 (positive regulation of triglyceride biosynthetic process)
   - GO:0045723 (positive regulation of fatty acid biosynthetic process)
   - GO:0010907 (positive regulation of glucose metabolic process)
   - These are secondary metabolic effects from overexpression studies

3. **Reproductive Functions**
   - GO:0007283 (spermatogenesis) - tissue-specific function
   - GO:1901318 (negative regulation of flagellated sperm motility) - tissue-specific effect

### Removed Annotations

1. **Overly Broad Terms**
   - GO:0016787 (hydrolase activity) - REMOVED (too general, gluconolactonase is specific)
   - GO:0046872 (metal ion binding) - REMOVED (too general, calcium and zinc binding are specific)

### Modified Annotations

1. **Enzyme Regulator Activity**
   - GO:0030234 (enzyme regulator activity) - MODIFIED to GO:0060590 (ATPase regulator activity)
   - More specific term better describes Ca2+-ATPase regulation

2. **ATP-dependent Activity**
   - GO:0032781 (positive regulation of ATP-dependent activity) - MODIFIED to GO:0060590 (ATPase regulator activity)
   - More specific to actual function

### Cellular Localization

- GO:0005737 (cytoplasm) - ACCEPTED (predominant location)
- GO:0005634 (nucleus) - ACCEPTED (regulates nuclear processes)
- GO:0005739 (mitochondrion) - Added for Ca2+-ATPase regulation function

## Summary Statistics

- Total annotations reviewed: 47
- ACCEPTED: 35
- KEEP_AS_NON_CORE: 8
- REMOVED: 2
- MODIFIED: 2
- MARK_AS_OVER_ANNOTATED: 0
- UNDECIDED: 0

## Key Literature Support

- PMID:16585534 - Gluconolactonase activity and vitamin C biosynthesis
- PMID:16786169 - Ca2+-ATPase regulation in mitochondria
- PMID:699201 - Original calcium-binding protein discovery
- PMID:2280766 - Protein synthesis inhibition
- PMID:11500948 - DNA synthesis regulation
- PMID:12397604 - RNA synthesis regulation
- PMID:15806309, PMID:16167335 - Anti-apoptotic effects
- PMID:23615721 - Sperm function effects

## Conclusion

The curation confirms Rgn/SMP30 as a multifunctional protein with four major core functions:
1. Essential enzyme in vitamin C biosynthesis (gluconolactonase)
2. Calcium-binding protein regulating calcium homeostasis
3. Suppressor of macromolecule biosynthesis (protein, DNA, RNA)
4. Cytoprotective factor with anti-apoptotic effects

Secondary functions in development, metabolism, and reproduction were appropriately marked as non-core. Overly broad terms were removed or replaced with more specific annotations.