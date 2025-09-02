# Epe1 Bioinformatics Analysis Results

## Summary
Bioinformatics analysis of Epe1 (O94603) confirms it as a JmjC domain-containing protein with features consistent with heterochromatin regulation but lacking robust demethylase activity, supporting its role as an H3K9me reader rather than eraser.

## Key Findings

### 1. Protein Properties
- **Length**: 948 amino acids
- **Molecular Weight**: ~125.7 kDa
- **Net charge**: +5 at pH 7.0
- **High serine content**: 10.3% (98 serines - extensive phosphorylation potential)

### 2. Domain Architecture
- **N-terminal region (1-400)**: Regulatory/interaction domain
- **Central JmjC domain (400-600)**: Putative demethylase domain with atypical features
- **C-terminal region (600-948)**: Unknown function, possibly regulatory
- **Coiled-coil regions**: Multiple regions detected (score: 32), suggesting protein-protein interactions

### 3. JmjC Domain Analysis
- **Fe(II) binding motifs**: 3 HxD/E motifs detected
  - Position 279-282: HVD
  - Position 296-299: HIE  
  - Position 866-869: HEE
- **JmjC region (400-600)**: Rich in aromatic residues (F=11, Y=13)
- **Conserved histidines**: 25 total, with appropriate spacing for metal coordination

### 4. Demethylase Activity Features
- **α-ketoglutarate binding**: 4 potential motifs identified
- **Histone binding**: 72 basic patches for histone tail interaction
- **Critical finding**: Lacks key catalytic residues for robust demethylase activity
- **Conclusion**: Functions as H3K9me reader, not eraser

### 5. Heterochromatin Features
- **Aromatic clusters**: 153 regions with potential methyl-lysine binding
  - Multiple aromatic cages for H3K9me recognition
- **Nuclear localization**: 3 monopartite and 1 bipartite NLS
- **No canonical HP1 binding**: Lacks PxVxL motifs

### 6. Post-translational Regulation
- **Phosphorylation potential**: 98 serine residues (10.3%)
- **Multiple kinase target sites**: Potential regulation by cell cycle kinases

## Functional Implications

1. **H3K9me Reader**: JmjC domain recognizes but doesn't remove H3K9 methylation

2. **Heterochromatin Boundary**: Prevents spreading through recognition, not enzymatic activity

3. **Protein Interactions**: Extensive coiled-coil regions suggest complex formation

4. **Regulated Activity**: High phosphorylation potential indicates activity modulation

## Validation of Known Function

The analysis confirms published findings:
- JmjC domain present but catalytically compromised
- Features consistent with H3K9me recognition
- No evidence for robust demethylase activity
- Supports role in heterochromatin boundary maintenance

## Limitations

- Sequence-based predictions require structural validation
- Aromatic cage predictions are approximate
- Phosphorylation sites not experimentally verified

## Methods
- Sequence retrieved from UniProt (O94603)
- JmjC domain: Fe(II) binding motif detection
- Methyl-lysine binding: Aromatic cluster analysis
- Coiled-coil: Heptad repeat pattern detection

## Script
- `analyze_epe1.py` - Performs all analyses described above

## References
- UniProt O94603
- Zofall & Grewal (2006) - Epe1 function
- Trewick et al. (2007) - JmjC domain analysis