
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|P04637|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|TP53|
|description[?](https://w3id.org/ai4curation/gene_review/description)|Tumor suppressor p53 is a multifunctional transcription factor that acts as a guardian of the genome, coordinating cellular responses to diverse stress signals including DNA damage, oxidative stress, hypoxia, and metabolic stress. It regulates cell fate decisions through transcriptional control of genes involved in cell cycle arrest, apoptosis, senescence, DNA repair, and metabolism. TP53 is the most frequently mutated gene in human cancers (~50%), underlining its critical role in preventing malignant transformation. Beyond its canonical tumor suppressive functions, p53 also regulates ferroptosis, autophagy, immune responses, and stemness.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|
|---|---|
|GO_REF:0000002|Gene Ontology annotation through association of InterPro records with GO terms.|
|GO_REF:0000107|Automatic transfer of experimentally verified manual GO annotation data to orthologs using Ensembl Compara.|
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.|
|PMID:10360174|Mutations in human ARF exon 2 disrupt its nucleolar localization and impair its ability to block nuclear export of MDM2 and p53.|
|PMID:10962037|Overexpression of MYC causes p53-dependent G2 arrest of normal fibroblasts.|
|PMID:14654789|A member of the Pyrin family, IFI16, is a novel BRCA1-associated protein involved in the p53-mediated apoptosis pathway.|
|PMID:14744935|Endoplasmic reticulum stress induces p53 cytoplasmic localization and prevents p53-dependent apoptosis by a pathway involving glycogen synthase kinase-3beta.|
|PMID:15710329|Human MUC1 oncoprotein regulates p53-responsive gene transcription in the genotoxic stress response.|
|PMID:16061649|Ckap2 regulates aneuploidy, cell cycling, and cell death in a p53-dependent manner.|
|PMID:17599062|CDIP, a novel pro-apoptotic gene, regulates TNFalpha-mediated apoptosis in a p53-dependent manner.|
|PMID:17996705|An acetylation switch in p53 mediates holo-TFIID recruitment.|
|PMID:30089260|p53 Regulates the Expression of LRP1 and Apoptosis through a Stress Intensity-Dependent MicroRNA Feedback Loop.|
|PMID:30514107|CELF1/p53 axis: a sustained antiproliferative signal leading to villus atrophy under total parenteral nutrition.|
|UniProt:P04637|UniProt entry for Human TP53|

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000976|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transcription cis-regulatory region binding|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 is a well-established DNA-binding transcription factor that directly binds to p53 response elements in target gene promoters. This is one of its core molecular functions, supported by extensive experimental evidence including the IDA annotations from PMID:15710329 and PMID:17996705.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This annotation accurately represents a core molecular function of p53. Multiple lines of experimental evidence confirm p53 directly binds cis-regulatory regions to regulate transcription.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000122|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of transcription by RNA polymerase II|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 functions as both a transcriptional activator and repressor. It negatively regulates transcription of genes like BCL2, survivin, and various cell cycle-promoting genes. This represents part of its core transcriptional regulatory function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Well-documented function of p53 as a transcriptional repressor, complementing its role as an activator. This is essential for its tumor suppressor function.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000423|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mitophagy|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 has complex roles in autophagy and mitophagy regulation, with both promoting and inhibiting functions depending on context. There is also a negative regulation of mitophagy annotation (GO:1901525), suggesting context-dependent regulation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 regulates mitophagy through multiple mechanisms including transcriptional control of mitophagy-related genes. This is part of its broader role in cellular stress response and metabolism.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001701|
|label[?](https://w3id.org/ai4curation/gene_review/label)|in utero embryonic development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 plays roles in embryonic development, though p53 knockout mice are viable, indicating it is not absolutely essential. Development-related functions are context-specific rather than core functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While p53 has developmental roles, these are not its primary functions. p53-null mice can complete embryonic development, though with increased cancer susceptibility.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001756|
|label[?](https://w3id.org/ai4curation/gene_review/label)|somitogenesis|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Somitogenesis is a specific developmental process. While p53 may influence developmental processes, this is a highly specific annotation that is not a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Very specific developmental process that is not central to p53 function. p53 knockout mice complete somitogenesis successfully.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001836|
|label[?](https://w3id.org/ai4curation/gene_review/label)|release of cytochrome c from mitochondria|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 promotes apoptosis through the mitochondrial pathway, including regulation of cytochrome c release via transcriptional targets like BAX and PUMA, as well as direct mitochondrial functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core component of p53-mediated intrinsic apoptotic pathway. p53 induces pro-apoptotic BCL2 family members that trigger cytochrome c release.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0002309|
|label[?](https://w3id.org/ai4curation/gene_review/label)|T cell proliferation involved in immune response|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 has emerging roles in immune regulation, but specific T cell proliferation is not a core function. This may relate to non-cell-autonomous tumor suppression through immune modulation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While p53 influences immune responses, direct regulation of T cell proliferation is not a primary function. This is a context-specific, non-core role.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0002326|
|label[?](https://w3id.org/ai4curation/gene_review/label)|B cell lineage commitment|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|B cell lineage commitment is a specific hematopoietic developmental process. While p53 may influence hematopoiesis, this is not a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Highly specific developmental process in the immune system. Not a primary function of p53 as a tumor suppressor.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0002360|
|label[?](https://w3id.org/ai4curation/gene_review/label)|T cell lineage commitment|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|T cell lineage commitment is a specific developmental process in hematopoiesis. Not a core p53 function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Specific immune system developmental process that is not central to p53 tumor suppressor function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0002931|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to ischemia|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 responds to hypoxic stress and ischemia as part of its cellular stress response functions. This involves HIF pathway interactions and metabolic regulation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Part of p53s broader cellular stress response function. Ischemia/hypoxia is a well-documented p53-activating stress signal.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006302|
|label[?](https://w3id.org/ai4curation/gene_review/label)|double-strand break repair|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 coordinates DNA damage response including regulation of DNA repair genes. It promotes both repair and apoptosis depending on damage severity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function of p53 as guardian of the genome. p53 regulates multiple DNA repair pathways through transcriptional control of repair genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006357|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of transcription by RNA polymerase II|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core molecular function of TP53 as a transcription factor. This is a general parent term for both positive and negative regulation of transcription.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Fundamental molecular function of p53. Encompasses both activation and repression of RNA pol II transcription, which are core p53 activities.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein import into nucleus|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This likely refers to p53s own nuclear import rather than it regulating import of other proteins. p53 nuclear localization is critical for its function but this annotation is unclear.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Ambiguous annotation - p53 undergoes nuclear import but does not regulate protein import into nucleus as a general function. This appears to be a misannotation.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006974|
|label[?](https://w3id.org/ai4curation/gene_review/label)|DNA damage response|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|DNA damage response is one of the most fundamental functions of TP53. It acts as a sensor and coordinator of cellular responses to DNA damage. Also supported by IDA evidence from PMID:14744935.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Absolutely core function - p53 is the guardian of the genome, activated by DNA damage to coordinate repair, cell cycle arrest, or apoptosis.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006979|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to oxidative stress|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 is activated by oxidative stress and regulates antioxidant responses, ROS metabolism, and cell fate decisions under oxidative conditions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core stress response function. p53 responds to oxidative stress and regulates redox balance through multiple target genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007179|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transforming growth factor beta receptor signaling pathway|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 has complex interactions with TGF-beta signaling, but this is not a core function. There is also annotation for negative regulation of this pathway (GO:0030512).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While p53 interacts with TGF-beta signaling, this is a context-specific crosstalk rather than a core p53 function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007346|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of mitotic cell cycle|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell cycle regulation is a fundamental function of TP53. It induces cell cycle arrest through p21/CDKN1A and other cell cycle regulators.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core tumor suppressor function. p53 induces cell cycle checkpoints in response to stress, preventing damaged cells from proliferating.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007369|
|label[?](https://w3id.org/ai4curation/gene_review/label)|gastrulation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Gastrulation is a specific early developmental process. p53 is not essential for gastrulation as p53-null mice complete this process.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Highly specific developmental process that is not a core p53 function. p53 knockout mice undergo normal gastrulation.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007405|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neuroblast proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Neuroblast proliferation is a specific neural developmental process. Not a core p53 function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Specific developmental process in nervous system. p53 may regulate proliferation broadly but neuroblast-specific regulation is not core.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007406|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of neuroblast proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 negatively regulates proliferation broadly, but neuroblast-specific regulation is a specialized context.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While p53 inhibits proliferation generally, neuroblast-specific regulation is not a core function but rather a tissue-specific manifestation.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007417|
|label[?](https://w3id.org/ai4curation/gene_review/label)|central nervous system development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CNS development is a broad developmental process. p53 has roles but is not essential as p53-null mice have functional CNS.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Developmental process that is not core to p53 tumor suppressor function. p53-null mice develop functional nervous systems.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007507|
|label[?](https://w3id.org/ai4curation/gene_review/label)|heart development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Heart development is a specific organ developmental process. Not essential as p53-null mice have normal hearts.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Organ-specific developmental process that is not a core p53 function. p53 knockout mice develop normal hearts.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008156|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of DNA replication|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 inhibits DNA replication as part of cell cycle arrest, preventing replication of damaged DNA.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function related to cell cycle control and preventing propagation of damaged DNA. Part of p53s tumor suppressor activity.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008283|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell population proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a very general term. p53 primarily NEGATIVELY regulates proliferation. The annotation GO:0008285 (negative regulation) is more accurate.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Too general and potentially misleading. p53 primarily inhibits rather than promotes proliferation. The negative regulation annotation is more appropriate.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009303|
|label[?](https://w3id.org/ai4curation/gene_review/label)|rRNA transcription|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 can regulate ribosomal biogenesis and rRNA transcription as part of metabolic control, typically inhibiting it under stress.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 regulates ribosome biogenesis and protein synthesis capacity, typically suppressing rRNA transcription under stress conditions.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009410|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to xenobiotic stimulus|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 responds to various xenobiotic stresses including drugs and toxins that cause DNA damage or cellular stress.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Part of p53s broad stress response function. Many xenobiotics activate p53 through DNA damage or other stress pathways.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009411|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to UV|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|UV radiation causes DNA damage that strongly activates p53. This is a well-characterized p53-activating stress.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core stress response function. UV-induced DNA damage is a classic p53 activator leading to repair, arrest, or apoptosis.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009651|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to salt stress|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Salt stress response is not a well-characterized p53 function in mammalian cells. This may be an over-annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|No strong evidence for salt stress as a p53-activating stimulus in mammalian cells. Likely computational over-annotation.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009792|
|label[?](https://w3id.org/ai4curation/gene_review/label)|embryo development ending in birth or egg hatching|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Very broad developmental term. p53-null mice are viable, so p53 is not essential for embryo development to birth.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While p53 has developmental roles, it is not essential for embryo development to birth. This is a non-core function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0010165|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to X-ray|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|X-ray radiation causes DNA damage that activates p53. Part of p53s DNA damage response function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - ionizing radiation like X-rays cause DNA damage that strongly activates p53.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0010332|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to gamma radiation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Gamma radiation causes DNA damage that strongly activates p53. Classic p53-activating genotoxic stress.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - gamma radiation is a well-established p53 activator through DNA damage signaling.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0010629|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of gene expression|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 functions as both activator and repressor of gene expression. It represses anti-apoptotic and proliferation genes.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core transcriptional function. p53 represses multiple genes including BCL2, survivin, and cell cycle genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0010659|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cardiac muscle cell apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell type-specific apoptosis. While p53 induces apoptosis broadly, cardiac-specific regulation is context-dependent.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific manifestation of p53s apoptotic function. Not a core function but a specialized context.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0010666|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of cardiac muscle cell apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Tissue-specific apoptotic regulation. p53 can induce apoptosis in cardiac cells but this is not a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Cell type-specific apoptotic function. The general apoptotic function is core, but cardiac-specific regulation is context-dependent.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0014009|
|label[?](https://w3id.org/ai4curation/gene_review/label)|glial cell proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell type-specific proliferation. Not a core p53 function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific proliferation regulation in nervous system. Not central to p53 tumor suppressor function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0019661|
|label[?](https://w3id.org/ai4curation/gene_review/label)|glucose catabolic process to lactate via pyruvate|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TP53 regulates glucose metabolism and suppresses glycolysis (Warburg effect). There is also annotation for negative regulation (GO:1904024).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 primarily SUPPRESSES glycolysis/Warburg effect, not promotes it. The negative regulation annotation is more accurate.|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0021549|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cerebellum development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Specific brain region development. Not essential as p53-null mice have functional cerebellum.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Region-specific developmental process that is not core to p53 function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030330|
|label[?](https://w3id.org/ai4curation/gene_review/label)|DNA damage response, signal transduction by p53 class mediator|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is the most specific and accurate term for p53s core DNA damage response function. p53 IS the p53 class mediator.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Absolutely core function - this term specifically describes p53s role as the central mediator of DNA damage response signaling.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030512|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of transforming growth factor beta receptor signaling pathway|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 has complex crosstalk with TGF-beta signaling but this is not a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Context-specific signaling crosstalk rather than core p53 function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0031571|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mitotic G1 DNA damage checkpoint signaling|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|G1 checkpoint is a critical p53 function. p53 induces p21/CDKN1A to arrest cells at G1/S transition after DNA damage.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core cell cycle checkpoint function. p53-mediated G1 arrest via p21 is one of the best-characterized p53 responses.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0033077|
|label[?](https://w3id.org/ai4curation/gene_review/label)|T cell differentiation in thymus|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Specific immune system developmental process. Not a core p53 function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific developmental process in immune system. Not central to p53 tumor suppressor function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0033554|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cellular response to stress|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Broad term encompassing p53s fundamental role as a stress sensor and responder. Covers DNA damage, oxidative stress, hypoxia, etc.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is the master cellular stress response coordinator, responding to diverse stress signals.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0034103|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of tissue remodeling|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Tissue remodeling is not a well-characterized p53 function. This appears to be an over-annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Not a core p53 function. May relate to indirect effects through apoptosis or senescence but not a primary role.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035264|
|label[?](https://w3id.org/ai4curation/gene_review/label)|multicellular organism growth|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Very broad developmental term. p53 influences growth through proliferation control but is not essential for organism growth.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Broad developmental process. p53-null mice show normal growth patterns despite cancer susceptibility.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035794|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of mitochondrial membrane permeability|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 promotes mitochondrial outer membrane permeabilization during apoptosis through BAX/BAK activation. Core apoptotic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core component of p53-mediated intrinsic apoptosis. p53 target genes like BAX and PUMA promote MOMP.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042127|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of cell population proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|General term for proliferation regulation. p53 primarily negatively regulates proliferation through cell cycle arrest.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is a key regulator of proliferation, primarily through negative regulation via cell cycle checkpoints.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043065|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Promotion of apoptosis is one of p53s most fundamental tumor suppressor functions. It induces apoptosis through multiple pathways.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Absolutely core function. p53 is a major pro-apoptotic transcription factor, inducing PUMA, BAX, NOXA, FAS, and other death genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043066|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 has context-dependent anti-apoptotic functions, though its primary role is pro-apoptotic. Can induce survival genes in some contexts.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While primarily pro-apoptotic, p53 can have anti-apoptotic functions in certain contexts, such as mild stress or through p21-mediated cell cycle arrest.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043504|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mitochondrial DNA repair|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 has mitochondrial localization and can contribute to mitochondrial genome stability, though this is not a primary function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While p53 can localize to mitochondria and influence mtDNA stability, this is not a core tumor suppressor function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043516|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of DNA damage response, signal transduction by p53 class mediator|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 both executes and regulates its own DNA damage response pathway. This is a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is the central regulator of its own pathway, with multiple feedback loops and regulatory mechanisms.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043523|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of neuron apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell type-specific apoptotic regulation. p53 can regulate neuronal apoptosis but this is context-specific.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific manifestation of p53s apoptotic function. Not a core function but relevant in neurodegeneration contexts.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043525|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of neuron apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Neuron-specific pro-apoptotic function. While p53 induces apoptosis broadly, neuron-specific regulation is context-dependent.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Cell type-specific apoptotic function relevant in neurodegeneration but not a core p53 function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045861|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of proteolysis|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 can regulate proteolysis through various mechanisms including MDM2 regulation and proteasome activity, but this is not a primary function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While p53 influences protein stability and degradation pathways, negative regulation of proteolysis is not a core function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045893|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of DNA-templated transcription|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core function of p53 as a transcriptional activator. This is a broad term covering p53s activation of target genes.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Fundamental molecular function. p53 is primarily a transcriptional activator of numerous target genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045930|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of mitotic cell cycle|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core function - p53 induces cell cycle arrest at multiple checkpoints to prevent proliferation of damaged cells.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Fundamental tumor suppressor function. p53 arrests cell cycle through p21 and other CDK inhibitors.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045944|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of transcription by RNA polymerase II|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core transcriptional activation function. Also supported by IGI and IDA evidence from other entries.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Fundamental molecular function - p53 is a sequence-specific transcriptional activator of numerous genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0048144|
|label[?](https://w3id.org/ai4curation/gene_review/label)|fibroblast proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell type-specific proliferation. p53 generally inhibits rather than promotes proliferation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 primarily inhibits proliferation. This positive proliferation annotation appears incorrect. The negative regulation annotation (GO:0048147) is accurate.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0048147|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of fibroblast proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 inhibits fibroblast proliferation. This is a cell type-specific manifestation of its anti-proliferative function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Cell type-specific anti-proliferative function. Consistent with p53s general growth suppressive role but not core.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0048568|
|label[?](https://w3id.org/ai4curation/gene_review/label)|embryonic organ development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Broad developmental term. p53 has developmental roles but is not essential for organ development.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Developmental process that is not core to p53 function. p53-null mice develop functional organs.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050821|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein stabilization|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 can stabilize some proteins through transcriptional targets or protein interactions, but this is not a primary function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Not a core p53 function. p53 itself is heavily regulated by stability but protein stabilization of other proteins is not primary.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0051276|
|label[?](https://w3id.org/ai4curation/gene_review/label)|chromosome organization|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Very broad term. p53 influences chromosome stability through DNA repair and checkpoints but does not directly organize chromosomes.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Too broad and indirect. p53 maintains genome stability but chromosome organization per se is not its function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0051402|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neuron apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Neuron-specific apoptosis. While p53 can induce neuronal apoptosis, this is context-specific.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Cell type-specific apoptotic process. Relevant in neurodegeneration but not a core p53 function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0051726|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of cell cycle|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell cycle regulation is a fundamental p53 function. Also supported by ISS evidence. This is a core tumor suppressor function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is a master regulator of cell cycle checkpoints at G1/S and G2/M transitions.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060253|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of glial cell proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell type-specific anti-proliferative function in nervous system.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific manifestation of p53s anti-proliferative function. Not a core function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060411|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cardiac septum morphogenesis|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Very specific cardiac developmental process. Not a core p53 function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Highly specific developmental process. p53-null mice have normal cardiac development.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0070242|
|label[?](https://w3id.org/ai4curation/gene_review/label)|thymocyte apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell type-specific apoptosis in immune system development.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific apoptotic process in immune system. Not a core p53 function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0070243|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of thymocyte apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Immune cell-specific apoptotic regulation during T cell development.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific process in immune system development. Not a core p53 function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0070266|
|label[?](https://w3id.org/ai4curation/gene_review/label)|necroptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 has emerging roles in regulating necroptosis, a form of programmed necrosis, though apoptosis is its primary death pathway.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Emerging function of p53 in non-apoptotic cell death. p53 can regulate necroptosis through various mechanisms.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0071480|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cellular response to gamma radiation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Gamma radiation induces DNA damage that strongly activates p53. Core stress response function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is activated by ionizing radiation-induced DNA damage to coordinate repair or cell death.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0071494|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cellular response to UV-C|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|UV-C radiation causes severe DNA damage that activates p53. Part of DNA damage response.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - UV-C induces DNA lesions that strongly activate p53-mediated responses.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0072089|
|label[?](https://w3id.org/ai4curation/gene_review/label)|stem cell proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 generally suppresses rather than promotes proliferation, including in stem cells.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 typically inhibits stem cell proliferation to maintain genomic stability. The negative regulation annotation (GO:2000647) is more accurate.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0072331|
|label[?](https://w3id.org/ai4curation/gene_review/label)|signal transduction by p53 class mediator|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This term specifically describes p53s signal transduction function. p53 IS the p53 class mediator. Absolutely core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - this term precisely describes p53s role as the central signal transducer in stress responses.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0072332|
|label[?](https://w3id.org/ai4curation/gene_review/label)|intrinsic apoptotic signaling pathway by p53 class mediator|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This specifically describes p53-mediated intrinsic apoptosis through mitochondrial pathway. Core apoptotic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is the master regulator of intrinsic apoptosis through BAX, PUMA, NOXA and other target genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0072593|
|label[?](https://w3id.org/ai4curation/gene_review/label)|reactive oxygen species metabolic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 regulates ROS through multiple mechanisms including antioxidant genes and metabolic regulation. Important for redox balance.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Important function - p53 regulates cellular redox state through multiple target genes affecting ROS production and scavenging.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1901525|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of mitophagy|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 has complex bidirectional effects on autophagy/mitophagy - can both promote and inhibit depending on context.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 can negatively regulate mitophagy in certain contexts while promoting it in others. This dual role is well-documented.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1902108|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of mitochondrial membrane permeability involved in apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core apoptotic function - p53 regulates mitochondrial outer membrane permeabilization through BAX/BAK activation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function in intrinsic apoptosis. p53 target genes like BAX and PUMA promote mitochondrial membrane permeabilization.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1902253|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of intrinsic apoptotic signaling pathway by p53 class mediator|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 both executes and regulates its own intrinsic apoptotic pathway. Core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is the master regulator of its own apoptotic pathway with multiple feedback mechanisms.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1903799|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of miRNA processing|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 has complex effects on miRNA processing - can both promote and inhibit. There is also positive regulation annotation (GO:1902895).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 regulates miRNA processing through interaction with Drosha complex and other mechanisms. Can have both positive and negative effects.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1904024|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of glucose catabolic process to lactate via pyruvate|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 suppresses glycolysis/Warburg effect as part of its metabolic regulatory and tumor suppressor functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Important metabolic function - p53 suppresses aerobic glycolysis (Warburg effect) that is characteristic of cancer cells.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1990144|
|label[?](https://w3id.org/ai4curation/gene_review/label)|intrinsic apoptotic signaling pathway in response to hypoxia|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 mediates apoptosis in response to severe hypoxia. Part of its stress response repertoire.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core stress response function - p53 induces apoptosis under severe hypoxic stress through HIF interactions and target genes.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:2000269|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of fibroblast apoptotic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell type-specific apoptotic regulation in fibroblasts.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Tissue-specific apoptotic function. The general apoptotic function is core, but fibroblast-specific regulation is context-dependent.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:2000378|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of reactive oxygen species metabolic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 regulates cellular redox balance and can suppress ROS through antioxidant gene expression.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Important function - p53 maintains redox homeostasis by inducing antioxidant genes like GPX1, SOD2, and others.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:2000647|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of stem cell proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 maintains stem cell genomic stability by limiting proliferation. Important for preventing cancer stem cells.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Important function - p53 restricts stem cell self-renewal to maintain genomic integrity and prevent transformation.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:2000772|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of cellular senescence|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Senescence induction is a core p53 tumor suppressor mechanism, alternative to apoptosis for preventing damaged cell proliferation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 induces senescence through p21 and other targets as an irreversible cell cycle arrest mechanism.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1902895|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of miRNA transcription|
|IMP|[PMID:30089260](PMID:30089260)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 transcriptionally regulates multiple miRNAs including miR-103 and miR-107 as shown in this paper. Important for stress response modulation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Well-documented function with strong experimental evidence. p53 regulates numerous miRNAs to modulate cellular responses.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008285|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of cell population proliferation|
|ISS|[PMID:30514107](PMID:30514107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core tumor suppressor function - p53 negatively regulates proliferation through multiple mechanisms including cell cycle arrest.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Fundamental tumor suppressor function supported by extensive evidence. p53 inhibits proliferation through p21 and other targets.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0051726|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of cell cycle|
|ISS|[PMID:30514107](PMID:30514107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cell cycle regulation is a core p53 function. Duplicate of the IEA annotation but with experimental support.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function with strong experimental support. p53 is a master cell cycle regulator.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1902749|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of cell cycle G2/M phase transition|
|IMP|[PMID:10962037](PMID:10962037)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 regulates the G2/M checkpoint, preventing mitosis of damaged cells. Core cell cycle control function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core cell cycle checkpoint function with experimental evidence. p53 induces G2 arrest through multiple mechanisms.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045944|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of transcription by RNA polymerase II|
|IGI|[PMID:16061649](PMID:16061649)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core transcriptional activation function with genetic interaction evidence. Duplicate annotation with different evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Fundamental molecular function with strong experimental support through genetic interactions.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000976|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transcription cis-regulatory region binding|
|IDA|[PMID:15710329](PMID:15710329)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Direct experimental evidence for p53 DNA binding activity. Core molecular function with IDA support.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core molecular function with direct experimental evidence of p53 binding to cis-regulatory regions.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042771|
|label[?](https://w3id.org/ai4curation/gene_review/label)|intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator|
|IDA|[PMID:14654789](PMID:14654789)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core apoptotic function - p53 mediates intrinsic apoptosis specifically in response to DNA damage. Strong IDA evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Absolutely core function with direct experimental evidence. p53 is THE p53 class mediator for DNA damage-induced apoptosis.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006974|
|label[?](https://w3id.org/ai4curation/gene_review/label)|DNA damage response|
|IDA|[PMID:14744935](PMID:14744935)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core function with direct experimental evidence. Duplicate of IEA annotation but with IDA support.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Fundamental function with strong experimental support. p53 is the guardian of the genome.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006983|
|label[?](https://w3id.org/ai4curation/gene_review/label)|ER overload response|
|IDA|[PMID:14744935](PMID:14744935)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 responds to ER stress. The cited paper shows ER stress affects p53 localization and function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 responds to ER stress as part of its cellular stress response repertoire. Supported by experimental evidence.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042771|
|label[?](https://w3id.org/ai4curation/gene_review/label)|intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator|
|IDA|[PMID:14744935](PMID:14744935)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation with different reference. Core apoptotic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53-mediated intrinsic apoptosis in response to DNA damage.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042981|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of apoptotic process|
|IDA|[PMID:14744935](PMID:14744935)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|General term for apoptosis regulation. p53 both promotes and inhibits apoptosis depending on context.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - p53 is a master regulator of apoptosis with strong experimental evidence.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016604|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nuclear body|
|IDA|[PMID:10360174](PMID:10360174)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 localizes to nuclear bodies including PML bodies and nucleoli under certain conditions. This is related to its regulatory functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 localizes to nuclear bodies like PML bodies, especially during stress responses and senescence. This subcellular localization is functionally relevant.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000976|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transcription cis-regulatory region binding|
|IDA|[PMID:17996705](PMID:17996705)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Additional IDA evidence for p53 DNA binding. Duplicate annotation with different experimental support.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core molecular function with direct experimental evidence. p53 binds to consensus sequences in target gene promoters.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045944|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of transcription by RNA polymerase II|
|IDA|[PMID:17599062](PMID:17599062)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Additional IDA evidence for transcriptional activation. p53 activates CDIP expression.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core transcriptional function with direct experimental evidence showing p53-dependent gene activation.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042149|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cellular response to glucose starvation|
|IDA|[PMID:14744935](PMID:14744935)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|p53 responds to metabolic stress including glucose deprivation. Part of its metabolic regulatory functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|p53 acts as a metabolic stress sensor, responding to glucose starvation through AMPK and other pathways.|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|
|---|---|---|---|---|
|p53 functions as a sequence-specific DNA-binding transcription factor that binds p53 response elements to regulate transcription of target genes involved in cell cycle arrest, apoptosis, DNA repair, and metabolism||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000976|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transcription cis-regulatory region binding|
|||
|p53 localizes to nuclear bodies and regulates both activation and repression of transcription, including senescence programs||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000976|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transcription cis-regulatory region binding|
|||

## suggested_questions


|question|
|---|
|How does p53 integrate diverse stress signals to determine cell fate decisions between survival, senescence, and apoptosis?|
|What determines the selectivity of p53 for different target gene promoters and how is this modulated by post-translational modifications?|
|How do p53 isoforms and mutant forms interact to modulate tumor suppressor function in heterozygous cancer cells?|

## suggested_experiments


|description|
|---|
|Single-cell time-lapse imaging to track p53 dynamics and correlate oscillation patterns with specific cell fate outcomes|
|ChIP-seq coupled with PRO-seq to map p53 binding and transcriptional outcomes under different stress conditions|
|Proximity labeling proteomics to identify context-specific p53 interactors that determine target gene selectivity|
