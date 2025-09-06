# LPL1 GO Annotation Curation Summary

## Overview
Systematic review of 17 existing GO annotations for yeast LPL1 (Q08448), a lipid droplet phospholipase B enzyme.

## Key Findings

### Accepted Annotations (8)
- **GO:0006629** (lipid metabolic process) - IBA, IEA - Correctly captures broad function
- **GO:0005811** (lipid droplet) - Multiple evidence codes (IBA, IEA, IDA x3) - Well-supported localization
- **GO:0016042** (lipid catabolic process) - IEA - Accurate for phospholipid breakdown
- **GO:0055088** (lipid homeostasis) - IMP - Supported by mutant phenotypes

### Annotations Requiring Modification (4)
- **GO:0004622** (phosphatidylcholine lysophospholipase activity) - IBA, IDA, IMP
  - **Issue**: Too narrow; captures only second step of phospholipase B activity
  - **Proposed replacement**: GO:0102545 (phospholipase B activity)
  
- **GO:0004623** (phospholipase A2 activity) - IEA
  - **Issue**: Captures only sn-2 cleavage, missing sn-1 activity
  - **Proposed replacement**: GO:0102545 (phospholipase B activity)

### Annotations to Remove (1)
- **GO:0047372** (monoacylglycerol lipase activity) - IBA
  - **Rationale**: Incorrect substrate specificity; LPL1 acts on phospholipids, not monoacylglycerols
  - Likely mis-transferred from paralog ROG1

### Over-annotated Terms (2)
- **GO:0016020** (membrane) - IEA - Too general; lipid droplet is more specific
- **GO:0016787** (hydrolase activity) - IEA - Uninformatively broad

### Non-core Annotations (1)
- **GO:0005783** (endoplasmic reticulum) - HDA - Transient/minor localization

## Critical Missing Annotations

Based on literature evidence, particularly Weisshaar et al. 2017 (PMID:28218240/PMC5349779), the following important annotations are missing:

### 1. Protein Quality Control Function
- **GO:0043161** (proteasome-mediated ubiquitin-dependent protein catabolic process) or related regulation term
- **Evidence**: LPL1 is part of Rpn4 regulon, induced under proteotoxic stress, and lpl1Δ hac1Δ double mutants accumulate ubiquitinated proteins and show defective proteasomal degradation
- **Suggested annotation**: "positive regulation of proteasome-mediated ubiquitin-dependent protein catabolic process"

### 2. Correct Molecular Function
- **GO:0102545** (phospholipase B activity)
- **Evidence**: Direct biochemical assays show sequential cleavage at both sn-1 and sn-2 positions
- Should replace the narrower lysophospholipase and phospholipase A2 annotations

### 3. Lipid Droplet Organization
- **GO:0034389** (lipid droplet organization)
- **Evidence**: lpl1Δ mutants show aberrant lipid droplet morphology with enlarged/fused droplets

### 4. Response to Stress
- **GO:0006950** (response to stress) or more specific stress response terms
- **Evidence**: LPL1 expression induced by Rpn4 under proteotoxic stress conditions

## Recommendations

1. **Consolidate molecular function annotations**: Replace multiple partial activity terms with comprehensive phospholipase B activity (GO:0102545)

2. **Add protein quality control annotations**: The literature strongly supports a role in proteostasis that is completely missing from current annotations

3. **Remove incorrect annotations**: The monoacylglycerol lipase activity has no experimental support

4. **Reduce redundancy**: Multiple duplicate annotations with different evidence codes could be streamlined while maintaining evidence trails

5. **Add missing biological processes**: Lipid droplet organization and stress response functions are well-supported but not annotated

## Literature Gaps
- The key Weisshaar et al. 2017 paper (showing proteostasis role) appears to not be fully incorporated into GO annotations
- The phospholipase B activity characterization from Selvaraju et al. 2014 is incompletely represented

## Conclusion
Current annotations capture LPL1's lipid metabolism and droplet localization well but miss its important role in protein quality control and use overly narrow enzymatic activity terms. The annotation set would benefit from updating to reflect the full breadth of LPL1 function based on available experimental evidence.