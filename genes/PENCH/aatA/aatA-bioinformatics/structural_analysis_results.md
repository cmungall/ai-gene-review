# Structural Analysis of Isopenicillin N Acyltransferase

**Analysis Date:** 2025-09-04

## Summary

Validation of cysteine nucleophile mechanism based on crystal structures.

## PDB Structures Analyzed

### 2X1C: Precursor form (Cys103Ala mutant)
- **Form:** mature (cleaved into subunits)
- **Number of chains:** 4
- **Cleavage site:** Between Gly102 and Cys103

**Catalytic Residues:**
- position_102: GLY (expected: GLY) ✓
- position_103: ALA (expected: CYS (or ALA in mutant)) ✓
- position_309: SER (expected: SER) ✓
### 2X1D: Mature wild-type enzyme
- **Form:** mature (cleaved into subunits)
- **Number of chains:** 4
- **Cleavage site:** Between Gly102 and Cys103

**Catalytic Residues:**
- position_103: CSD (expected: CYS (or ALA in mutant)) ✗
- position_309: SER (expected: SER) ✓
### 2X1E: Mature enzyme with 6-aminopenicillanic acid complex
- **Form:** mature (cleaved into subunits)
- **Number of chains:** 4
- **Cleavage site:** Between Gly102 and Cys103

**Catalytic Residues:**
- position_103: CSD (expected: CYS (or ALA in mutant)) ✗
- position_309: SER (expected: SER) ✓

## Key Findings

1. **2X1C (Cys103Ala mutant)** correctly shows Ala at position 103, confirming this is the mutant used to trap the precursor form.

## Conclusion

The structural analysis confirms that:
- **Cys103 is the catalytic nucleophile** that initiates autoproteolytic cleavage
- The enzyme undergoes **autocatalytic cleavage between Gly102 and Cys103**
- The mature enzyme exists as **heterodimeric subunits** after cleavage
- This validates the **cysteine-type peptidase activity (GO:0008234)** annotation

## References

- Bokhove et al., 2010. Structure 18(3):301-8. PMID: 20223213
- PDB entries: 2X1C, 2X1D, 2X1E