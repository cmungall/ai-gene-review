# LPL1 (YOR059C) Gene Review Notes

## Summary
LPL1 encodes a lipid droplet phospholipase B with dual roles in lipid metabolism and protein quality control in S. cerevisiae.

## Key Findings

### Molecular Function
- **Phospholipase B activity** - cleaves both sn-1 and sn-2 positions [PMID:25014274 "The purified Lpl1p showed phospholipase activity with broader substrate specificity, acting on all glycerophospholipids primarily at sn-2 position and later at sn-1 position"]
- Broad substrate specificity: PE, PC, PS, PA, PG
- Contains conserved GXSXG lipase motif essential for activity

### Cellular Localization
- **Lipid droplet** - primary localization confirmed by multiple studies [PMID:25014274, PMID:10515935, PMID:24868093]
- Particularly abundant in stationary phase cells
- Minor ER association detected by high-throughput study [PMID:26928762]

### Biological Roles

#### 1. Lipid Droplet Homeostasis
- Deletion causes enlarged/aberrant lipid droplets [PMID:25014274 "deletion of LPL1 resulted in altered morphology of LDs"]
- Regulates droplet phospholipid composition through hydrolysis
- Overexpression decreases glycerophospholipids and increases free fatty acids

#### 2. Protein Quality Control
- Part of Rpn4-mediated proteotoxic stress response [Weisshaar et al. 2017]
- Required for efficient proteasomal degradation under stress
- Double mutant hac1Δ lpl1Δ shows severe stress sensitivity and accumulates ubiquitinated proteins
- Links lipid droplet function with proteostasis

### Regulation
- Induced by Rpn4 transcription factor under proteotoxic stress
- Contains PACE element in promoter for Rpn4 binding
- Expression increases in stationary phase

## Annotation Review Decisions

### Key Changes Made:
1. **MODIFY**: Phosphatidylcholine lysophospholipase activity → Phospholipase B activity (GO:0102545)
2. **REMOVE**: Monoacylglycerol lipase activity (no evidence, likely from paralog ROG1)
3. **ACCEPT**: Lipid droplet localization (multiple supporting studies)
4. **KEEP_AS_NON_CORE**: ER localization (minor/transient)

### Missing Annotations Identified:
- Cellular response to misfolded protein (GO:0071218)
- Phospholipase B activity (GO:0102545) - more accurate than current narrow terms
- Lipid droplet organization (GO:0034389)

## Evolutionary Context
- Paralog of ROG1 (monoacylglycerol lipase) in yeast
- Distant homologs in mammals (FAM135A/B) but function not characterized
- Member of conserved ROG1 lipase family

## References
- Selvaraju et al. 2014 (PMID:25014274) - Biochemical characterization
- Weisshaar et al. 2017 (PMID:28100635) - Proteostasis role
- Athenstaedt et al. 1999 (PMID:10515935) - Early lipid droplet localization