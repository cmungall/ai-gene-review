
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q8N441|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|FGFRL1|
|description[?](https://w3id.org/ai4curation/gene_review/description)|FGFRL1 is the fifth member of the FGFR family that acts as a non-signaling decoy receptor. It has three extracellular Ig-like domains for binding FGF ligands and heparin, but lacks the intracellular tyrosine kinase domain. FGFRL1 functions primarily in ligand sequestration, cell adhesion via constitutive homodimers, and cell-cell fusion. It is essential for kidney development, diaphragm muscle formation, and skeletal development.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|findings|
|---|---|---|
|GO_REF:0000033|Annotation inferences using phylogenetic trees||
|GO_REF:0000044|Gene Ontology annotation based on UniProtKB/Swiss-Prot Subcellular Location vocabulary mapping, accompanied by conservative changes to GO terms applied by UniProt.||
|GO_REF:0000107|Automatic transfer of experimentally verified manual GO annotation data to orthologs using Ensembl Compara.||
|GO_REF:0000108|Automatic assignment of GO terms using logical inference, based on on inter-ontology links.||
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.||
|PMID:12813049|Characterization of FGFRL1, a novel fibroblast growth factor (FGF) receptor preferentially expressed in skeletal tissues.||
|PMID:18061161|The cell surface receptor FGFRL1 forms constitutive dimers that promote cell adhesion.||
|Reactome:R-HSA-5654510|FGFRL1 binds SPRED1/2||
|Reactome:R-HSA-5654511|FGFRL1 dimer binds FGFs||
|FGFRL1-deep-research.md|FGFRL1 Structure, Primary Function, and Evolutionary Perspective||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005007|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor receptor activity|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 binds FGF ligands through its extracellular Ig-like domains but lacks the intracellular kinase domain for signal transduction. It acts as a decoy receptor that sequesters FGFs without transducing signals.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While FGFRL1 binds FGF ligands, calling it a "receptor activity" is misleading because it does not transduce signals. A more accurate term would focus on its binding activity.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 is a single-pass transmembrane protein that localizes to the plasma membrane where it forms constitutive homodimers and mediates cell adhesion.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This cellular component annotation is well-supported by experimental evidence showing FGFRL1 localizes to the plasma membrane.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0017134|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor binding|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 directly binds multiple FGF ligands including FGF2, FGF3, FGF4, FGF8, FGF10, FGF18, and FGF22 through its extracellular Ig-like domains.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This molecular function is well-established and represents a core function of FGFRL1 as demonstrated by multiple experimental studies.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008543|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor receptor signaling pathway|
|IEA|[GO_REF:0000108](GO_REF:0000108)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 does not participate in canonical FGF receptor signaling as it lacks the intracellular kinase domain. Instead, it acts as a negative regulator by sequestering FGF ligands.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|FGFRL1 does not transduce FGF signals but rather negatively regulates the pathway by acting as a decoy receptor. A more accurate term would be negative regulation of the pathway.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016020|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane|
|IEA|[GO_REF:0000044](GO_REF:0000044)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 is a transmembrane protein that localizes to plasma membrane and can also be found in Golgi and transport vesicles during trafficking.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This is a valid but general cellular component term. More specific terms like plasma membrane are preferable but this broad term is not incorrect.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001501|
|label[?](https://w3id.org/ai4curation/gene_review/label)|skeletal system development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 is essential for proper skeletal development, including cartilage and bone formation. Knockout mice show cranial abnormalities and humans with FGFRL1 mutations develop craniosynostosis.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This biological process annotation is well-supported by developmental phenotypes in knockout mice and human patients with FGFRL1 mutations.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003179|
|label[?](https://w3id.org/ai4curation/gene_review/label)|heart valve morphogenesis|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|While FGFRL1 is expressed in heart tissue, there is no specific experimental evidence for its role in heart valve morphogenesis. This appears to be an over-annotation based on computational inference.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This annotation lacks specific experimental support for heart valve morphogenesis. FGFRL1 expression in heart tissue alone does not justify this specific developmental process annotation.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005007|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor receptor activity|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a duplicate of the IBA annotation above. FGFRL1 binds FGFs but does not have receptor activity in the classical sense due to lack of kinase domain.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Same reasoning as the IBA annotation - FGFRL1 binds FGFs but lacks signal transduction activity. Should be changed to binding activity.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a duplicate of the IBA annotation above. FGFRL1 is correctly localized to the plasma membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate cellular component annotation supported by experimental evidence.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008285|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of cell population proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 overexpression reduces cell proliferation and DNA synthesis in cultured cells by acting as a decoy receptor that sequesters growth-promoting FGF signals.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This biological process annotation is well-supported by experimental evidence showing FGFRL1 overexpression inhibits cell proliferation.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0017134|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor binding|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a duplicate of the IBA annotation above. FGFRL1 binds multiple FGF ligands through its extracellular domains.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate molecular function annotation representing a core function of FGFRL1.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060412|
|label[?](https://w3id.org/ai4curation/gene_review/label)|ventricular septum morphogenesis|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a very specific cardiac developmental process annotation with no experimental support. While FGFRL1 is expressed in heart, there is no evidence for its specific role in ventricular septum formation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This annotation is overly specific and lacks experimental support. FGFRL1 expression in heart tissue does not justify this precise developmental process annotation.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060539|
|label[?](https://w3id.org/ai4curation/gene_review/label)|diaphragm development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 is essential for diaphragm development. Knockout mice die at birth due to diaphragm muscle defects, specifically lacking slow-twitch muscle fibers needed for respiratory function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This biological process annotation is strongly supported by the lethal phenotype in FGFRL1 knockout mice due to diaphragm malformation.|
|additional_reference_ids[?](https://w3id.org/ai4curation/gene_review/additional_reference_ids)|[FGFRL1-deep-research.md](FGFRL1-deep-research.md)|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042802|
|label[?](https://w3id.org/ai4curation/gene_review/label)|identical protein binding|
|IPI|[PMID:18061161](PMID:18061161)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 forms constitutive homodimers at the cell surface as demonstrated by FRET and co-precipitation experiments.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This molecular function annotation is well-supported by experimental evidence showing FGFRL1 forms homodimers.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|IDA|[PMID:18061161](PMID:18061161)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is another duplicate plasma membrane annotation. FGFRL1 localizes to the plasma membrane where it forms homodimers.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate cellular component annotation with strong experimental support.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008201|
|label[?](https://w3id.org/ai4curation/gene_review/label)|heparin binding|
|IMP|[PMID:18061161](PMID:18061161)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 has a basic region in its extracellular domain that binds heparin and heparan sulfate, which is required for both FGF binding and cell adhesion functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This molecular function annotation is well-supported by experimental evidence and represents an important functional property of FGFRL1.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0044291|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell-cell contact zone|
|IDA|[PMID:18061161](PMID:18061161)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 is enriched at cell-cell contact sites when overexpressed, consistent with its role in cell adhesion.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This cellular component annotation is supported by experimental localization studies showing FGFRL1 accumulates at sites of cell-cell contact.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0098742|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell-cell adhesion via plasma-membrane adhesion molecules|
|IMP|[PMID:18061161](PMID:18061161)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 promotes cell adhesion through its extracellular Ig-like domains, acting similarly to nectin-like cell adhesion molecules.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This biological process annotation accurately describes a key function of FGFRL1 beyond its role as a decoy receptor.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005794|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Golgi apparatus|
|IDA|[PMID:18061161](PMID:18061161)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 can be found in the Golgi apparatus as part of its trafficking through the secretory pathway before reaching the plasma membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This cellular component annotation reflects the trafficking route of FGFRL1 through the secretory pathway.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030133|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transport vesicle|
|IDA|[PMID:18061161](PMID:18061161)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|FGFRL1 is found in transport vesicles during its trafficking between cellular compartments, including endocytosis from the plasma membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This cellular component annotation reflects the dynamic trafficking of FGFRL1, which is rapidly internalized via its cytoplasmic trafficking motifs.|
|additional_reference_ids[?](https://w3id.org/ai4curation/gene_review/additional_reference_ids)|[FGFRL1-deep-research.md](FGFRL1-deep-research.md)|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|TAS|[Reactome:R-HSA-5654510](Reactome:R-HSA-5654510)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is another duplicate plasma membrane annotation. FGFRL1 is correctly localized to the plasma membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate cellular component annotation based on curated pathway information.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|TAS|[Reactome:R-HSA-5654511](Reactome:R-HSA-5654511)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is yet another duplicate plasma membrane annotation. FGFRL1 localizes to the plasma membrane where it binds FGF ligands.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate cellular component annotation based on curated pathway information.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005007|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor receptor activity|
|IDA|[PMID:12813049](PMID:12813049)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is another duplicate of the receptor activity annotation. FGFRL1 binds FGFs but lacks signal transduction capability.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Same reasoning as previous receptor activity annotations - FGFRL1 binds FGFs but does not have true receptor activity due to lack of kinase domain.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|IDA|[PMID:12813049](PMID:12813049)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is another duplicate plasma membrane annotation. FGFRL1 is correctly localized to the plasma membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate cellular component annotation with experimental support.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008201|
|label[?](https://w3id.org/ai4curation/gene_review/label)|heparin binding|
|IDA|[PMID:12813049](PMID:12813049)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a duplicate of the IMP heparin binding annotation above. FGFRL1 specifically binds heparin through its basic extracellular domain.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate molecular function annotation representing an important binding activity of FGFRL1.|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|anatomical_locations|
|---|---|---|---|---|---|
|FGF ligand sequestration and signaling modulation||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0017134|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor binding|
||||
|Cell-cell adhesion through constitutive homodimerization||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042802|
|label[?](https://w3id.org/ai4curation/gene_review/label)|identical protein binding|
||||
|Heparan sulfate-mediated cell surface interactions||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008201|
|label[?](https://w3id.org/ai4curation/gene_review/label)|heparin binding|
||||
|Regulation of cell proliferation through FGF signaling modulation||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0017134|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor binding|
||||
|Essential developmental processes regulation||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0017134|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast growth factor binding|
||||

## suggested_questions


|question|experts|
|---|---|
|What are the specific binding partners of FGFRL1's Ig3 domain that mediate cell-cell fusion?| Beat Trueb Thomas Rieckmann|
|How does FGFRL1's fusogenic activity contribute to muscle development and fiber type specification?| Beat Trueb|
|What is the physiological significance of FGFRL1 ectodomain shedding and when does it occur in vivo?| Beat Trueb Thomas Rieckmann|
|How does FGFRL1 interaction with Spred1/2 modulate FGF signaling in kidney development?||

## suggested_experiments


|hypothesis|description|experiment_type|
|---|---|---|
|FGFRL1 is required for proper muscle fiber formation across multiple muscle types, not just diaphragm|Create muscle-specific knockout of Fgfrl1 using Cre-lox system to determine if muscle fiber formation defects extend beyond the diaphragm|Genetic knockout|
|FGFRL1's fusogenic activity is essential for its developmental functions in muscle and bone|Generate knock-in mouse with fusion-defective FGFRL1 mutant (single amino acid change in Ig3 domain) to test if fusogenic activity is required for development|Gene targeting|
|FGFRL1 interacts with specific cell surface proteins to mediate cell-cell adhesion and fusion|Perform co-immunoprecipitation and mass spectrometry to identify FGFRL1-interacting proteins on adjacent cell surfaces|Protein interaction mapping|
|FGFRL1 shapes FGF morphogen gradients by sequestering ligands and affects their spatial distribution|Use fluorescent FGF ligands in organ culture to examine how FGFRL1 affects FGF gradient formation and persistence|Live imaging|
