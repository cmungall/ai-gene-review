
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q9BRQ4|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|CFAP300|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|C11orf70|
|description[?](https://w3id.org/ai4curation/gene_review/description)|CFAP300 (Cilia And Flagella Associated Protein 300) is a small (~267 amino acids) cytoplasmic protein essential for the assembly and transport of axonemal dynein arms. Functions as both a chaperone/assembly factor for dynein complexes and as a cargo adaptor linking dynein components to intraflagellar transport (IFT) machinery. Contains a conserved DUF4498 domain critical for protein-protein interactions. Mutations in CFAP300  cause Primary Ciliary Dyskinesia 38 (CILD38), characterized by complete loss of outer and inner dynein arms, resulting in immotile cilia, chronic respiratory infections, male infertility, and laterality defects.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|findings|
|---|---|---|
|PMID:29727692|C11orf70 mutations disrupting the intraflagellar transport-dependent assembly of multiple axonemal dyneins cause primary ciliary dyskinesia||
|PMID:29727693|Mutations in C11orf70 cause primary ciliary dyskinesia with randomization of left/right body asymmetry due to defects of outer and inner dynein arms||
|file:human/CFAP300/CFAP300-deep-research.md|Deep Research Report: CFAP300 comprehensive analysis||
|GO_REF:0000024|Manual transfer of experimentally-verified manual GO annotation data to orthologs by curator judgment of sequence similarity.||
|GO_REF:0000043|Gene Ontology annotation based on UniProtKB/Swiss-Prot keyword mapping||
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:29727692](PMID:29727692)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|DNAAF2 interaction supported by both publications but experimental details need verification|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:29727693 states "C11orf70 shows an interaction with cytoplasmic ODA/IDA assembly factor DNAAF2" but publications don't detail experimental methods in abstracts. However, 'protein binding' is uninformative per curation guidelines. CFAP300 functions as an adapter/scaffolding protein for dynein arm assembly.|

#### proposed_replacement_terms


|Slot|Value|
|---|---|
|additional_reference_ids[?](https://w3id.org/ai4curation/gene_review/additional_reference_ids)|[PMID:29727693](PMID:29727693)|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cytoplasmic localization directly supported by experimental evidence in publications|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:29727692 demonstrates that "Tagged C11orf70 in Paramecium and Chlamydomonas localizes mainly in the cytoplasm with a small amount in the ciliary component." PMID:29727693 confirms involvement in "cytoplasmic assembly of dynein arms." This validates the sequence similarity inference.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005856|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoskeleton|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on keywords|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Too general - while technically correct that cilia are part of cytoskeleton, this annotation lacks specificity for CFAP300's actual functional location|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0031514|
|label[?](https://w3id.org/ai4curation/gene_review/label)|motile cilium|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Motile cilium localization confirmed by experimental evidence showing ciliary transport|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:29727692 shows "Tagged C11orf70...localizes mainly in the cytoplasm with a small amount in the ciliary component" and demonstrates "its transport within cilia is IFT dependent." Both papers show CFAP300 mutations cause "immotile respiratory cilia" confirming functional relevance to motile cilia.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042995|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell projection|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on keywords|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Uninformative and redundant with motile cilium annotation. Provides no functional insight and is too general to be useful for understanding CFAP300 function|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CFAP300 is involved in cytoplasmic assembly of axonemal dynein arms before their intraflagellar transport to cilia. However, its primary functional location is in the cytoplasm as part of the dynein assembly machinery, making this cellular component annotation appropriate.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|CFAP300 functions in cytoplasmic dynein arm assembly before transport to cilia. The cytoplasmic localization is functionally relevant and well-supported.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0070286|
|label[?](https://w3id.org/ai4curation/gene_review/label)|axonemal dynein complex assembly|
|ISS|[file:human/CFAP300/CFAP300-deep-research.md](file:human/CFAP300/CFAP300-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CFAP300 is essential for assembly of both outer and inner dynein arm complexes in the cytoplasm before their transport to cilia.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This specific biological process is CFAP300's primary function - it acts as a dynein arm assembly factor that is required for both ODA and IDA formation. Mutations cause complete loss of both dynein arm types.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003351|
|label[?](https://w3id.org/ai4curation/gene_review/label)|epithelial cilium movement involved in extracellular fluid movement|
|ISS|[file:human/CFAP300/CFAP300-deep-research.md](file:human/CFAP300/CFAP300-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CFAP300 enables ciliary beating in respiratory epithelium for mucociliary clearance and in ependymal cells for CSF circulation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Loss of CFAP300 causes immotile respiratory cilia leading to chronic infections due to failed mucociliary clearance. This specific ciliary function is distinct from general cilium movement.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030317|
|label[?](https://w3id.org/ai4curation/gene_review/label)|flagellated sperm motility|
|ISS|[file:human/CFAP300/CFAP300-deep-research.md](file:human/CFAP300/CFAP300-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CFAP300 is essential for sperm flagellar motility through assembly of flagellar dynein arms. Mutations cause complete sperm immotility and male infertility.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Direct evidence that CFAP300 mutations cause asthenozoospermia with completely immotile sperm due to absence of flagellar dynein arms.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007368|
|label[?](https://w3id.org/ai4curation/gene_review/label)|determination of left/right symmetry|
|ISS|[file:human/CFAP300/CFAP300-deep-research.md](file:human/CFAP300/CFAP300-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CFAP300 is required for nodal ciliary beating that establishes left-right body asymmetry. Mutations frequently cause situs inversus and laterality defects.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|CFAP300 mutations result in randomized organ laterality due to nodal cilia dysfunction during embryonic development.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060285|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cilium-dependent cell motility|
|ISS|[PMID:29727692](PMID:29727692)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CFAP300 is required for cell motility in Paramecium and Chlamydomonas through its role in ciliary dynein arm assembly.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Model organism studies show CFAP300 knockdown causes reduced ciliary beating and impaired swimming, demonstrating its role in cilium-dependent cell locomotion.|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|anatomical_locations|substrates|
|---|---|---|---|---|---|---|
|Cytoplasmic assembly and intraflagellar transport of axonemal dynein arms||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030674|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein-macromolecule adaptor activity|
|||||

## suggested_questions


|question|
|---|
|How does CFAP300 coordinate with other cilia and flagella associated proteins to ensure proper axoneme assembly and stability?|
|What are the molecular mechanisms by which CFAP300 contributes to ciliary motility and sperm flagellar function?|
|How do mutations in CFAP300 lead to primary ciliary dyskinesia and male infertility at the cellular and molecular level?|
|What role does CFAP300 play in the radial spoke complex and how does it interact with dynein arms?|

## suggested_experiments


|description|
|---|
|Cryo-electron tomography of cilia and sperm flagella from CFAP300-deficient cells to visualize ultrastructural defects in the axoneme|
|High-resolution live-cell imaging to track CFAP300 dynamics during ciliogenesis and correlate with ciliary beating patterns|
|Proteomics analysis of the radial spoke complex to map CFAP300 protein-protein interactions and assembly hierarchy|
|CRISPR-mediated knock-in of disease-associated CFAP300 variants to study structure-function relationships in ciliary motility|
