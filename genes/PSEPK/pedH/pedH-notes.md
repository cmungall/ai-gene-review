# pedH Gene Research Notes

## Gene Overview
- **Gene Symbol**: pedH (PP_2679)
- **Organism**: Pseudomonas putida KT2440
- **UniProt ID**: Q88JH0
- **Protein Name**: Quinoprotein alcohol dehydrogenase PedH
- **Alternative Names**: Lanthanide-dependent PQQ-ADH, qedH-II

## Core Function
PedH encodes a lanthanide-dependent pyrroloquinoline quinone (PQQ)-dependent alcohol dehydrogenase that catalyzes the oxidation of alcohols and aldehydes in the periplasm. This enzyme is crucial for bacterial growth on volatile organic compounds (VOCs) as carbon and energy sources.

## Enzymatic Properties

### Cofactor Requirements
1. **Pyrroloquinoline Quinone (PQQ)**: Essential cofactor bound 1:1 per subunit
2. **Lanthanides**: Absolute requirement for enzymatic activity
   - Active lanthanides: La³⁺, Ce³⁺, Pr³⁺, Nd³⁺, Sm³⁺, Gd³⁺, Tb³⁺
   - Highest activity with Pr³⁺ and Nd³⁺
   - Activity range: 10 nM to 100 μM lanthanides, peak at 1 μM
   - Responds to ecologically relevant concentrations (1-10 nM La³⁺)

### Substrate Specificity
[PMID:28655819 "Functional Role of Lanthanides in Enzymatic Activity and Transcriptional Regulation of Pyrroloquinoline Quinone-Dependent Alcohol Dehydrogenases in Pseudomonas putida KT2440", "PedH exhibits enzyme activity on a range of substrates including linear and aromatic primary and secondary alcohols, as well as aldehydes"]

**Primary Alcohols**:
- Ethanol → Acetaldehyde (KM = 177 μM, Vmax = 10.6 μmol/min/mg)
- Butan-1-ol → Butanal
- Hexan-1-ol → Hexanal  
- Octan-1-ol → Octanal
- 2-Phenylethanol → 2-Phenylacetaldehyde (KM = 329 μM, Vmax = 11.8 μmol/min/mg)
- Cinnamyl alcohol → Cinnamaldehyde
- Farnesol → Farnesal

**Secondary Alcohols**:
- Butan-2-ol → Butan-2-one

**Aldehydes**:
- Acetaldehyde → Acetate (KM = 2261 μM, Vmax = 8.4 μmol/min/mg)
- Butanal → Butanoate
- Hexanal → Hexanoate
- Octanal → Octanoate

### Electron Transfer
PedH uses cytochrome c550 as its natural electron acceptor. The enzyme follows the typical PQQ-dependent mechanism:
1. Two-electron oxidation of substrate by PQQ
2. Two sequential one-electron transfers to cytochrome c via PQQ semiquinone radical
3. Reoxidation of PQQ for next catalytic cycle

## Subcellular Localization
- **Primary Location**: Periplasm [PMID:28655819 "Functional Role of Lanthanides", "periplasmic oxidation system"]
- **Signal Peptide**: Present (residues 1-27)
- **Mature Protein**: Residues 28-595

## Regulation and Expression
[PMID:28655819 "Functional Role of Lanthanides", "The underlying regulatory network is responsive to as little as 1 to 10 nM lanthanum, a concentration assumed to be of ecological relevance"]

### Transcriptional Control
- **Induction**: Presence of lanthanides
- **Repression**: Absence of lanthanides (pedE is induced instead)  
- **REE-switch**: Inverse regulation with pedE (calcium-dependent counterpart)
- **Sensory Function**: PedH itself acts as lanthanide sensory module for transcriptional regulation

### Two-Component System
- **Regulatory System**: PedS2/PedR2 two-component system
- **Response Range**: 5 nM to 10 μM La³⁺ depending on medium

## Physiological Role
[PMID:28655819 "Functional Role of Lanthanides", "These enzymes are crucial for efficient bacterial growth with a variety of volatile alcohols"]

1. **VOC Metabolism**: Enables growth on volatile organic compounds as sole carbon/energy source
2. **Environmental Adaptation**: Responds to lanthanide availability in soil environments
3. **Metabolic Flexibility**: Provides functional redundancy with pedE for alcohol oxidation
4. **Ecological Strategy**: Adaptive response to variable rare earth element availability

## Structural Features
- **Crystal Structures**: PDB entries 6ZCW (wild-type) and 6ZCV (F412V/W561A mutant)
- **Resolution**: 1.65-1.70 Å
- **Domain Architecture**: Eight-bladed β-propeller fold typical of PQQ dehydrogenases
- **Disulfide Bond**: Cys-131/Cys-132 essential for electron transfer to cytochrome c550
- **PQQ Binding Sites**: Multiple residues coordinate PQQ cofactor

## Comparison with PedE
- **PedE**: Calcium-dependent PQQ-ADH (PP_2674)
- **PedH**: Lanthanide-dependent PQQ-ADH (PP_2679)
- **Substrate Overlap**: Similar range of alcohols and aldehydes
- **Regulation**: Inversely regulated (REE-switch mechanism)
- **Function**: Provides metabolic redundancy for VOC utilization

## Biotechnology Applications  
[From structural studies: "attractive biocatalyst for oxidation of 5-(hydroxymethyl)furoic acid into 5-formylfuroic acid"]
- **Biocatalysis**: Engineered variants for bio-based chemical production
- **FDCA Production**: Key enzyme in furan-2,5-dicarboxylic acid synthesis pathway
- **Green Chemistry**: Alternative to terephthalic acid in polymer production

## Key Research Papers
- [PMID:28655819] Wehrmann et al. 2017 - Initial functional characterization
- Crystal structure studies (2020) - Protein engineering and structural biology
- [PMID:31736923] 2019 - REE transport and iron availability effects

## Implications for GO Annotation
Based on this research, key functional aspects for GO annotation include:
- Specific alcohol dehydrogenase activity (cytochrome c dependent)
- PQQ cofactor binding
- Lanthanide metal ion binding (NOT calcium!)
- Periplasmic localization
- Role in volatile organic compound catabolism
- Transcriptional regulatory function