# CFAP300 Annotation Review History

## Major Review Update - August 24, 2025

### Background
Applied more rigorous curation criteria following project guidelines that emphasize specificity and functional relevance over technically correct but uninformative annotations.

### Changes Made to Existing Annotations

#### 1. GO:0005515 (protein binding)
- **Previous Action:** ACCEPT
- **New Action:** MODIFY
- **Rationale:** Following explicit project guideline that "protein binding doesn't tell us anything about actual function"
- **Proposed Replacement:** GO:0030674 (protein-macromolecule adaptor activity)
- **Evidence:** Strong experimental evidence for CFAP300-DNAAF2 interaction supports adapter function

#### 2. GO:0005737 (cytoplasm)  
- **Previous Action:** ACCEPT
- **New Action:** MARK_AS_OVER_ANNOTATED
- **Rationale:** Too broad and non-specific for a protein whose primary functional location is in motile cilia
- **Evidence Type:** ISS (inferred from sequence similarity)

#### 3. GO:0005856 (cytoskeleton)
- **Previous Action:** ACCEPT  
- **New Action:** MODIFY
- **Rationale:** Too general - lacks specificity for CFAP300's actual functional location
- **Proposed Replacement:** GO:0035085 (cilium axoneme)
- **Evidence Type:** IEA (automated annotation)

#### 4. GO:0031514 (motile cilium)
- **Action:** ACCEPT (unchanged)
- **Rationale:** Most accurate and experimentally supported annotation
- **Evidence Type:** ISS

#### 5. GO:0042995 (cell projection)
- **Previous Action:** ACCEPT
- **New Action:** REMOVE
- **Rationale:** Uninformative, redundant with motile cilium annotation, provides no functional insight
- **Evidence Type:** IEA (automated annotation)

### Summary of Review Philosophy Changes

**Previous approach:** Accepted technically correct annotations even if they were overly broad or uninformative

**New approach:** Applied "ruthless" criteria prioritizing:
1. **Functional specificity** over technical correctness
2. **Informative annotations** that provide biological insight
3. **Core function relevance** aligned with experimental evidence
4. **Removal of redundant** or overly general terms

### Impact
- Reduced accepted annotations from 5/5 to 1/5 (20% acceptance rate)
- Increased specificity through 2 MODIFY actions with proposed replacements
- Eliminated 1 over-annotation and 1 uninformative term
- Maintained focus on CFAP300's core ciliary dynein arm assembly function

### Guidelines Applied
- "Avoid the term `protein binding`, this doesn't tell us anything about the actual function"
- Emphasis on specificity over generality
- Removal of annotations that don't contribute to functional understanding
- Prioritization of experimentally supported, functionally relevant terms