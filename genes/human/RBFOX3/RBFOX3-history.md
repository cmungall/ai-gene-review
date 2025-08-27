# RBFOX3 Annotation Review History

## Initial Review - August 24, 2025

### Background
Applied rigorous curation criteria following project guidelines that emphasize specificity and functional relevance over technically correct but uninformative annotations. RBFOX3 is particularly important as it encodes the widely-used neuronal marker NeuN and functions as a critical regulator of alternative splicing in the nervous system.

### Changes Made to Existing Annotations

#### 1. GO:0000381 (regulation of alternative mRNA splicing, via spliceosome)
- **Action:** ACCEPT
- **Evidence Type:** IEA (automated annotation)
- **Rationale:** Core function of RBFOX3 - well-documented regulation of alternative splicing through position-dependent binding to (U)GCAUG motifs

#### 2. GO:0003676 (nucleic acid binding)
- **Action:** MODIFY
- **Evidence Type:** IEA 
- **Rationale:** Too general - RBFOX3 specifically binds RNA, not DNA
- **Proposed Replacement:** GO:0003723 (RNA binding)

#### 3. GO:0003723 (RNA binding)
- **Action:** ACCEPT
- **Evidence Type:** IEA
- **Rationale:** Accurate and well-supported - RBFOX3 contains RRM domain and binds specific RNA sequences

#### 4. GO:0007399 (nervous system development)
- **Action:** ACCEPT
- **Evidence Type:** IEA/IBA
- **Rationale:** Well-supported - RBFOX3 is essential for neuronal differentiation and CNS development

#### 5. GO:0043484 (regulation of RNA splicing)
- **Action:** KEEP_AS_NON_CORE
- **Evidence Type:** IEA
- **Rationale:** Accurate but less specific than 'regulation of alternative mRNA splicing' - alternative splicing is RBFOX3's core function

#### 6. GO:0005634 (nucleus)
- **Action:** MODIFY
- **Evidence Type:** IEA/IBA
- **Rationale:** Incomplete - RBFOX3 has both nuclear and cytoplasmic isoforms with distinct functions
- **Proposed Replacement:** GO:0070013 (intracellular organelle lumen)

#### 7. GO:0005737 (cytoplasm)
- **Action:** MODIFY
- **Evidence Type:** IEA/IBA
- **Rationale:** Too general - RBFOX3 cytoplasmic isoforms have specific functions
- **Proposed Replacement:** GO:0005829 (cytosol)

#### 8. GO:0006397 (mRNA processing)
- **Action:** MODIFY
- **Evidence Type:** IEA
- **Rationale:** Too general - RBFOX3 specifically regulates alternative splicing, not general mRNA processing
- **Proposed Replacement:** GO:0000381 (regulation of alternative mRNA splicing, via spliceosome)

#### 9. GO:0008380 (RNA splicing)
- **Action:** MODIFY
- **Evidence Type:** IEA
- **Rationale:** Less specific than RBFOX3's actual function - it regulates alternative splicing, not general splicing
- **Proposed Replacement:** GO:0000381 (regulation of alternative mRNA splicing, via spliceosome)

#### 10. GO:0003729 (mRNA binding)
- **Action:** ACCEPT
- **Evidence Type:** IBA (phylogenetic inference)
- **Rationale:** Accurate and specific - RBFOX3 binds to specific mRNA sequences through its RRM domain

### Summary of Review Philosophy

**Review Criteria Applied:**
1. **Functional Specificity:** Prioritized terms that capture RBFOX3's specific role in alternative splicing regulation
2. **Neuronal Context:** Emphasized the neuronal-specific nature of RBFOX3 function
3. **Isoform Awareness:** Considered the existence of functionally distinct nuclear and cytoplasmic isoforms
4. **Mechanistic Accuracy:** Focused on RBFOX3's position-dependent splicing regulation mechanism

### Impact Analysis
- **Accepted annotations:** 4/10 (40% acceptance rate)
- **Modified annotations:** 5/10 with specific replacements proposed
- **Non-core annotations:** 1/10 (kept but marked as less important)
- **Core function focus:** Alternative splicing regulation in neuronal contexts

### Key Improvements Made

#### Specificity Enhancements:
- Replaced general "nucleic acid binding" with specific "RNA binding"
- Converted broad "mRNA processing" to specific "alternative mRNA splicing regulation"
- Refined general "RNA splicing" to "alternative mRNA splicing regulation"

#### Cellular Context Refinements:
- Addressed incomplete nuclear localization annotation (isoform-dependent)
- Refined cytoplasmic localization to more specific cytosol
- Recognized dual nuclear/cytoplasmic functionality

#### Functional Hierarchy:
- Identified alternative splicing regulation as core function
- Classified general RNA splicing regulation as non-core
- Emphasized neuronal development context

### Guidelines Applied
- Emphasis on specificity over generality in molecular function annotations
- Recognition of protein isoform diversity and functional complexity
- Prioritization of experimentally-relevant and mechanistically-informative terms
- Removal of redundant or overly broad cellular component annotations

### Notable Gene-Specific Considerations
- **NeuN Identity:** RBFOX3's role as the definitive neuronal marker influenced cellular component assessment
- **Alternative Splicing Focus:** Core function clearly established as alternative splicing regulation
- **Isoform Complexity:** Multiple functionally distinct isoforms required nuanced localization annotations
- **Family Context:** Part of RBFOX protein family with cross-regulatory functions