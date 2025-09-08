# Surf1 (mouse) Gene Review Notes

## Key Facts
- SURF1 is a nuclear-encoded assembly factor for cytochrome c oxidase (Complex IV)
- Located in mitochondrial inner membrane with two transmembrane domains
- Component of MITRAC complex (mitochondrial translation regulation assembly intermediate)
- Essential for proper COX assembly, particularly COX1 subunit maturation and heme a3 incorporation

## Experimental Evidence from Literature

### PMID:12566387 - Constitutive knockout study
- Surf1 KO mice show:
  - High embryonic lethality (~90%)
  - Early-onset mortality in surviving mice
  - Profound COX deficiency in skeletal muscle and liver
  - Mitochondrial proliferation in skeletal muscle
  - Motor performance deficits
  - NO neurological symptoms despite COX deficiency
- Confirms COX assembly function [PMID:12566387, "profound and isolated defect of COX activity in skeletal muscle and liver"]

### PMID:14651853 - Mitochondrial proteomics
- Confirms mitochondrial localization via high-throughput proteomics [PMID:14651853, "proteomic survey of mitochondria from mouse brain, heart, kidney, and liver"]
- Part of comprehensive mitochondrial protein inventory

### PMID:18614015 - MitoCarta compendium
- Identified as mitochondrial protein in comprehensive MS study [PMID:18614015, "mitochondrial compendium of 1098 genes"]
- Confirmed mitochondrial localization across multiple tissues

### PMID:23838831 - Paradoxical phenotype study  
- Surf1 KO mice show paradoxical phenotype:
  - 49% reduction in COX activity
  - 24% reduction in state 3 respiration  
  - 24% increase in H2O2 production
  - BUT: enhanced memory, increased lifespan
  - Increased cerebral blood flow and glucose metabolism
  - Elevated HIF-1α and pCREB levels
- Direct evidence for COX activity enabling function [PMID:23838831, "49% reduction in COX activity"]

## Mouse vs Human Differences
- Mouse Surf1 KO: viable with metabolic adaptation, extended lifespan
- Human SURF1 mutations: severe Leigh syndrome, early mortality
- Mouse has better compensatory mechanisms than humans

## GO Annotation Review Summary

### Core Functions (well-supported):
1. **GO:0033617 - mitochondrial respiratory chain complex IV assembly** - Strong experimental support from multiple sources
2. **GO:0004129 - cytochrome-c oxidase activity** - Enables COX activity through assembly function (IMP evidence)
3. **GO:0005739 - mitochondrion** - Well-supported localization (HDA evidence)

### Questionable/Over-annotations:
1. **GO:1902600 - proton transmembrane transport** - Indirect inference from COX activity, not a direct function of SURF1
2. **GO:0016020 - membrane** - Too general, should be more specific (mitochondrial inner membrane)
3. **GO:0031966 - mitochondrial membrane** - Correct but could be more specific (inner membrane)

## Core Functions Summary
Based on all evidence, SURF1's core functions are:
1. Assembly factor for cytochrome c oxidase (Complex IV)
2. Component of MITRAC complex facilitating COX1 maturation
3. Facilitates heme a3 incorporation into COX active site
4. Stabilizes COX assembly intermediates

NOT core functions (avoid over-annotation):
- Direct proton transport (this is COX's function, not SURF1's)
- Direct electron transport (again, COX's function)
- ATP synthesis (downstream effect)