# CFAP418 (C8orf37) Gene Review Notes

## Gene Overview
- **Gene Symbol**: CFAP418 (also known as C8orf37, smalltalk)
- **UniProt ID**: Q96NL8
- **Protein**: Cilia- and flagella-associated protein 418
- **Size**: 207 amino acids
- **Domain**: Contains RMP (Retinal Maintenance Protein) domain (pfam14996, aa 63-175)

## Disease Associations
1. **Cone-rod dystrophy 16 (CORD16)** [MIM:614500]
   - Autosomal recessive
   - Early macular involvement
   - Cone loss precedes rod degeneration
   
2. **Retinitis pigmentosa 64 (RP64)** [MIM:614500]
   - Autosomal recessive
   - Rod loss precedes cone degeneration
   - Progressive peripheral vision loss

3. **Bardet-Biedl syndrome 21 (BBS21)** [MIM:617406]
   - Syndromic ciliopathy
   - Features: retinal degeneration, obesity, polydactyly, renal malformations, intellectual disability
   - First functional evidence from zebrafish studies [PMID:27008867 "C8orf37 knockdown reproduced cardinal BBS phenotypes"]

## Key Pathogenic Variants
- **R177W**: Associated with CORD16 and BBS21 [PMID:22177090 "c.529C>T [p.Arg177Trp]"; PMID:36233334 "does not affect interaction with FAM161A"]
- **Q182R**: Associated with RP64 [PMID:22177090 "c.545A>G [p.Gln182Arg]"; PMID:36233334 "does not affect interaction with FAM161A"]
- **L166***: Nonsense mutation in RP patient [PMID:22177090 "c.497T>A [p.Leu166(∗)]"]
- **c.156-2A>G**: Splice site mutation, associated with postaxial polydactyly [PMID:22177090 "two CRD siblings with the c.156−2A>G mutation also showed unilateral postaxial polydactyly"]

## Protein Localization
- **Primary cilium base**: Localized at basal body/transition zone in RPE1 cells [PMID:22177090 "C8orf37 localization at the base of the primary cilium of human retinal pigment epithelium cells"]
- **Photoreceptor connecting cilium**: Enriched at junction between inner and outer segments [PMID:22177090 "at the base of connecting cilia of mouse photoreceptors"]
- **Photoreceptor inner segment**: Present throughout inner segment [PMID:36233334 "C8orf37 immunoreactivity was enriched at the inner segment, including the ciliary base"]
- **Cytoplasm**: Diffuse cytoplasmic localization also observed
- **NOT in outer segment**: Absent from photoreceptor outer segment itself

## Protein Interactions
### FAM161A Interaction
- **Direct interaction confirmed**: Y2H, co-IP, proximity ligation assays [PMID:36233334 "C8orf37 interacted with FAM161A"]
- **Interaction domains**:
  - CFAP418 N-terminus (aa 1-75) required for binding [PMID:36233334 "N-terminal aa 1–75 region on C8orf37 was sufficient and required for its interaction with FAM161A"]
  - FAM161A UPF0564 domain (aa 341-517) [PMID:36233334 "aa 341–517 of FAM161A were sufficient for interaction with C8orf37"]
- **Pathogenic mutations do not disrupt interaction**: R177W and Q182R maintain FAM161A binding [PMID:36233334 "these mutations did not affect interactions between C8orf37 and FAM161A"]

### Other Interactions
- **CAPNS1**: Calpain small subunit 1 (IntAct database)

## Functional Evidence

### Mouse Knockout Studies
- Progressive photoreceptor degeneration (rods and cones) [Deep research: "C8orf37 knockout mouse, the absence of CFAP418 caused disorganized photoreceptor outer segment discs"]
- **Disorganized outer segment discs**: Key phenotype, suggests role in disc morphogenesis [Deep research: "severely disorganized outer segment discs in photoreceptors lacking CFAP418"]
- Normal connecting cilium structure
- No systemic BBS features in mice (species difference)

### Zebrafish Knockdown
- Visual impairment
- Kupffer's vesicle defects (ciliary organ)
- Delayed retrograde intraflagellar transport [PMID:27008867 "defects in retrograde melanosome transport"]
- Left-right asymmetry defects

## Molecular Function
- No enzymatic domains identified
- Likely scaffolding/adaptor protein at ciliary base
- May regulate:
  - Photoreceptor outer segment disc morphogenesis
  - Protein trafficking through connecting cilium
  - Intraflagellar transport (IFT)

## Expression Pattern
- **Ubiquitous expression** with enrichment in:
  - Retina (photoreceptors)
  - Brain
  - Heart
- Consistent with ciliary protein expression pattern

## Evolutionary Conservation
- Highly conserved across ciliated eukaryotes
- Absent in non-ciliated organisms (plants, fungi)
- C-terminal two-thirds most conserved
- Mouse ortholog 82% identical to human

## GO Annotation Assessment

### Cellular Component Annotations
1. **GO:0001917 (photoreceptor inner segment)**: Well-supported by experimental evidence
2. **GO:0005737 (cytoplasm)**: Supported, though broad term
3. **GO:0097546 (ciliary base)**: Strong experimental support from PMID:22177090

### Molecular Function Annotations
1. **GO:0005515 (protein binding)**: Too general, should specify FAM161A interaction

### Biological Process Annotations
1. **GO:0008594 (photoreceptor cell morphogenesis)**: Supported by mouse knockout data

## Core Functions Summary
Based on the evidence, CFAP418 functions as:
1. **Ciliary base scaffold protein** essential for photoreceptor survival
2. **Regulator of photoreceptor outer segment disc morphogenesis**
3. **Component of ciliary protein trafficking machinery** (via FAM161A interaction)
4. **Contributor to ciliary transport processes** (retrograde IFT)

## Key Supporting Literature
- PMID:22177090 - Initial disease gene identification, localization studies
- PMID:27008867 - BBS link, zebrafish functional studies
- PMID:36233334 - FAM161A interaction, domain mapping
- PMC5884456 - Mouse knockout, outer segment disc phenotype