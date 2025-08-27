---
gene_symbol: JAK1
species: Homo sapiens
initial_search_date: 2025-08-25
last_search_date: 2025-08-25
---

# Human JAK1 (Janus Kinase 1) - Research Notes

## Gene Overview

**Gene Symbol:** JAK1  
**Official Name:** Tyrosine-protein kinase JAK1  
**Alternative Names:** Janus kinase 1, JAK-1, JAK1A, JAK1B  
**Location:** Chromosome 1p31.3  
**UniProt ID:** P23458  
**Protein Size:** 1154 amino acids  
**EC Number:** EC 2.7.10.2

## Molecular Function

JAK1 is a **non-receptor tyrosine kinase** that plays a pivotal role in cytokine signaling [file:JAK1-deep-research.md, "JAK1 is a non-receptor tyrosine kinase that plays a pivotal role in cytokine signaling"]. It is one of four Janus kinases in humans (JAK1, JAK2, JAK3, TYK2) and serves as a critical component of the JAK-STAT signaling cascade.

### Key Functions:
- **Protein Tyrosine Kinase Activity**: Catalyzes phosphorylation of tyrosine residues on target proteins [file:JAK1-deep-research.md, "JAK1's kinase activity (protein tyrosine kinase activity, GO:0004715) is essential for transmitting signals from extracellular cytokines to changes in gene expression"]
- **Cytokine Signal Transduction**: Essential mediator of cytokine-mediated signaling pathways 
- **STAT Activation**: Phosphorylates STAT transcription factors to induce dimerization and nuclear translocation
- **Receptor Complex Assembly**: Required for assembly and function of multiple cytokine receptor complexes

### Protein Structure and Domains:
JAK1 contains seven JAK homology (JH1–JH7) regions forming four major domains [file:JAK1-deep-research.md, "The JAK1 protein is composed of 1154 amino acids and features a multi-domain architecture characteristic of the Janus kinase family"]:

- **N-terminal FERM domain (JH6–JH7)**: Mediates binding to cytokine receptor cytoplasmic tails [file:JAK1-deep-research.md, "A band 4.1, ezrin, radixin, moesin (FERM) domain at the N-terminus mediates binding to cytokine receptor cytoplasmic tails"]
- **SH2-like domain (JH3–JH5)**: Further aids in receptor interaction and positioning
- **Pseudokinase domain (JH2)**: Regulatory domain that maintains JAK1 in an off state until activation [file:JAK1-deep-research.md, "The pseudokinase exerts an autoinhibitory effect – it maintains JAK1 in an off state until the correct activation signal is received"]
- **C-terminal kinase domain (JH1)**: Active tyrosine kinase domain responsible for enzymatic activity [file:JAK1-deep-research.md, "The C-terminus of JAK1 is a tyrosine kinase domain responsible for its enzymatic activity"]

### Molecular Interactions:
- **Cytokine Receptors**: Associates with Type I and Type II cytokine receptors including IFN-α/β receptor (IFNAR), IFN-γ receptor, IL-6 family receptors (gp130), and γc-dependent receptors
- **JAK Partners**: Forms heterodimers with TYK2 (Type I IFN signaling) and JAK2 (IFN-γ signaling)
- **STAT Proteins**: Phosphorylates multiple STAT family members (STAT1, STAT3, STAT5, etc.)
- **PI3K Pathway**: Links to PI3K-Akt signaling through p85 subunit phosphorylation [file:JAK1-deep-research.md, "it associates with the interleukin-2 receptor β chain and is required for efficient recruitment and phosphorylation of the PI3K p85 subunit"]

## Cellular Localization and Subcellular Distribution

JAK1 is predominantly localized to the **cytoplasm** and **cytoplasmic side of plasma membrane** [file:JAK1-deep-research.md, "JAK1 is predominantly an intracellular cytosolic protein localized at the cytoplasmic side of the plasma membrane"]:

- **Primary Location**: Cytosol (GO:0005829) associated with membrane-proximal cytokine receptor complexes
- **Membrane Association**: Extrinsic component of cytoplasmic face of plasma membrane through receptor binding
- **Nuclear Translocation**: Can translocate to nucleus in certain contexts, though primary function is cytoplasmic [file:JAK1-deep-research.md, "there is evidence that JAK1 (and JAK2) can also translocate to the nucleus in certain contexts"]

## Expression Pattern

JAK1 shows **ubiquitous expression** across human tissues [file:JAK1-deep-research.md, "JAK1 is ubiquitously expressed in human tissues, consistent with its fundamental role in multiple cytokine systems"]:

### Tissue Distribution:
- **Low tissue specificity**: Expressed in all examined tissues
- **High expression**: Immune organs (spleen, lymph nodes, thymus) and immune cells
- **Broad expression**: Non-immune tissues (lung, liver, CNS) enabling cytokine responses
- **Developmental**: Expressed from early embryonic stages

### Regulation:
- **Constitutive expression**: Generally stable baseline expression
- **Post-transcriptional control**: Regulated by SOCS proteins, microRNAs (miR-15/16), and proteasomal degradation
- **Feedback regulation**: SOCS1/3 proteins provide negative feedback by targeting JAK1 for ubiquitination

## Biological Processes and Signaling Pathways

### Core Signaling Pathways:
JAK1 is essential for multiple cytokine signaling pathways [file:JAK1-deep-research.md, "JAK1 is involved in a broad array of biological processes, predominantly those related to cytokine signaling and immune function"]:

#### Type I Interferon Signaling:
- **IFN-α/β Responses**: Absolutely required - JAK1-deficient cells are completely unresponsive [file:JAK1-deep-research.md, "cells lacking JAK1 are completely unresponsive to both type I interferons (IFN-α/β) and type II interferon (IFN-γ)"]
- **Antiviral Defense**: Critical for innate immune responses to viral infections
- **Partner JAK**: Works with TYK2 in IFNAR complex

#### Type II Interferon Signaling:
- **IFN-γ Responses**: Essential component of IFN-γ receptor signaling
- **Immune Regulation**: Controls inflammatory and immune surveillance functions
- **Partner JAK**: Heterodimerizes with JAK2 in IFN-γ receptor complex

#### IL-6 Family Cytokine Signaling:
- **gp130-dependent cytokines**: Major kinase for IL-6, IL-11, LIF signaling [file:JAK1-deep-research.md, "JAK1 is the major kinase that phosphorylates the gp130 receptor subunit and activates STAT1/STAT3"]
- **STAT1/3 Activation**: Primary kinase for STAT phosphorylation in IL-6 pathway
- **Inflammatory Responses**: Mediates acute phase and inflammatory gene expression

#### γc-dependent Cytokine Signaling:
- **IL-2 family cytokines**: Required for IL-2, IL-4, IL-7, IL-9, IL-15, IL-21 responses
- **Lymphocyte Development**: Essential for T-cell and B-cell development through IL-7
- **Immune Cell Function**: Controls lymphocyte proliferation, survival, and differentiation

### Downstream Biological Processes:
- **Immune system development** (GO:0002520): Via IL-7 and other cytokines for lymphocyte development
- **Inflammatory response** (GO:0006954): Through IL-6, IL-27, IFN-γ pathways
- **Antiviral defense**: Type I interferon-mediated responses
- **Cell differentiation** (GO:0030154): Multiple cytokine-driven differentiation programs
- **Hematopoiesis**: Control of blood cell development and function

## Disease Associations and Clinical Significance

### Loss-of-Function Effects:
- **Mouse Studies**: JAK1⁻/⁻ mice are perinatally lethal, fail to nurse, and show severe runting [file:JAK1-deep-research.md, "Jak1⁻/⁻ mice are runted, fail to nurse, and die perinatally"]
- **Immunodeficiency**: Complete loss would cause severe immunodeficiency affecting multiple cytokine pathways
- **Human Mutations**: No viable complete JAK1 deficiency reported in humans (likely lethal)

### Gain-of-Function and Cancer:
#### Oncogenic Mutations:
- **Acute Leukemias**: Somatic activating mutations in T-cell ALL and other hematological cancers
- **JAK1 V658F**: Analogous to JAK2 V617F, causes constitutive signaling [file:JAK1-deep-research.md, "A well-known example is the JAK1 V658F mutation in the pseudokinase domain"]
- **Immune Evasion**: Tumor cells acquire inactivating JAK1 mutations to resist IFN-γ and immune surveillance

#### Inflammatory Disorders:
- **Germline GOF mutations**: Cause severe allergic inflammation and eosinophilia [file:JAK1-deep-research.md, "rare germline GOF mutations in JAK1 have recently been shown to cause an immunological disorder"]
- **Atopic Diseases**: JAK1 hyperactivation leads to extreme atopic dermatitis and asthma

### Therapeutic Targeting:
- **JAK Inhibitors**: Tofacitinib, upadacitinib, filgotinib for autoimmune/inflammatory diseases
- **Clinical Applications**: Rheumatoid arthritis, ulcerative colitis, other cytokine-driven disorders
- **Cancer Therapy**: Under investigation for tumors with aberrant JAK-STAT activation

## Key Experimental Evidence

### Landmark Studies:
- **Müller et al. (1993)**: Demonstrated absolute requirement for IFN-α/β and IFN-γ signaling [file:JAK1-deep-research.md, "created a mutant human cell line lacking JAK1 and found it was completely unresponsive to interferon-α/β and interferon-γ"]
- **Rodig et al. (1998)**: JAK1 knockout mouse phenotype showing perinatal lethality and selective cytokine unresponsiveness [file:JAK1-deep-research.md, "JAK1⁻/⁻ mice were born alive but failed to nurse and died perinatally"]
- **Hornakova et al. (2008-2010)**: Discovery of oncogenic JAK1 mutations in leukemia [file:JAK1-deep-research.md, "somatic JAK1 mutations in T-ALL patients that led to cytokine-independent JAK-STAT activation"]

## Evolutionary Conservation

JAK1 is highly conserved across vertebrates and invertebrates [file:JAK1-deep-research.md, "JAK1 and the JAK family are well conserved through evolution among multicellular organisms"]:

- **Vertebrate Conservation**: ~86% amino acid identity between human and mouse JAK1
- **Domain Architecture**: FERM, SH2, pseudokinase, kinase domains conserved from insects to mammals
- **Functional Conservation**: Key catalytic residues and regulatory motifs are invariant across species
- **Invertebrate Homologs**: Drosophila Hopscotch gene shows similar domain organization and function

## Research Implications

### Current Therapeutic Applications:
- **Precision Medicine**: JAK1 mutation status guides treatment decisions in hematological malignancies
- **Drug Development**: Multiple JAK1-selective inhibitors in clinical development
- **Biomarker Potential**: JAK1 activity as readout for cytokine pathway function

### Future Research Directions:
1. **Tissue-Specific Functions**: Understanding JAK1 roles beyond immune cells
2. **Nuclear Functions**: Investigating non-canonical nuclear JAK1 activities
3. **Combination Therapies**: JAK1 inhibitors with other targeted agents
4. **Resistance Mechanisms**: Understanding how tumors evade JAK1-dependent immune responses

## References and Data Sources

### Key Research Articles:
- **PMID:8232552** - Müller et al. - JAK1 requirement for interferon signaling (1993)
- **PMID:9590172** - Rodig et al. - JAK1 knockout mouse phenotype (1998) 
- **PMID:7537214** - IL-6/gp130 pathway JAK1 dominance study
- **PMID:20868368** - Hornakova et al. - Oncogenic JAK1 mutations in leukemia
- **PMID:21840487** - JAK1 in tumor cell invasion and metastasis

### Deep Research Sources:
- **file:JAK1-deep-research.md** - Comprehensive AI-generated research report using OpenAI Deep Research API, incorporating current literature and web sources (2025)
- **file:JAK1-citations.md** - Citation list from deep research analysis

### Database Sources:
- **UniProt** (P23458) - Tyrosine-protein kinase JAK1
- **Gene Ontology Annotation** - JAK1-goa.tsv with current GO annotations

*Last updated: August 2025*