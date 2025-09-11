
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q96NL8|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|CFAP418|
|description[?](https://w3id.org/ai4curation/gene_review/description)|CFAP418 (C8orf37) is a ciliary scaffolding protein localized at the photoreceptor connecting cilium base that is essential for photoreceptor outer segment disc morphogenesis and organization. It forms a protein complex with FAM161A at the ciliary base, contributing to photoreceptor structural integrity and survival. Mutations cause retinal dystrophies including cone-rod dystrophy 16 (CORD16), retinitis pigmentosa 64 (RP64), and Bardet-Biedl syndrome 21 (BBS21).|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|findings|
|---|---|---|
|GO_REF:0000024|Manual transfer of experimentally-verified manual GO annotation data to orthologs by curator judgment of sequence similarity.||
|GO_REF:0000044|Gene Ontology annotation based on UniProtKB/Swiss-Prot Subcellular Location vocabulary mapping, accompanied by conservative changes to GO terms applied by UniProt.||
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.||
|PMID:22177090|Mutations in C8orf37, encoding a ciliary protein, are associated with autosomal-recessive retinal dystrophies with early macular involvement.||
|PMID:27173435|An organelle-specific protein landscape identifies novel diseases and molecular mechanisms.||
|PMID:36233334|Interactions between C8orf37 and FAM161A, Two Ciliary Proteins Essential for Photoreceptor Survival.||
|file:human/CFAP418/CFAP418-bioinformatics/RESULTS.md|CFAP418 Bioinformatics Analysis Results||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001917|
|label[?](https://w3id.org/ai4curation/gene_review/label)|photoreceptor inner segment|
|IEA|[GO_REF:0000044](GO_REF:0000044)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This IEA annotation is supported by experimental evidence. CFAP418 is enriched at the photoreceptor inner segment, particularly at the connecting cilium base between inner and outer segments as shown by immunohistochemistry in both mouse (PMID:22177090) and marmoset retina (PMID:36233334). The protein is present throughout the inner segment but concentrated at the ciliary base.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:22177090 and PMID:36233334 provide direct immunohistochemical evidence for CFAP418 localization at the photoreceptor inner segment, with particular enrichment at the connecting cilium base.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This automated annotation is overly general. While CFAP418 is cytoplasmic (non-membrane bound), the more specific localization to the ciliary base (GO:0097546) better represents its functional localization. The cytoplasm term doesn't capture the protein's specialized ciliary localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:27173435](PMID:27173435)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|While this IPI evidence is valid (PMID:27173435 is a large-scale protein interaction study), the generic 'protein binding' term is uninformative about CFAP418's actual molecular function. The protein functions as a ciliary scaffold that specifically binds FAM161A and potentially other ciliary proteins. A more specific term like 'scaffold protein binding' (GO:0097110) would better represent its molecular function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:27173435 provides broad protein interaction data, but CFAP418 specifically functions as a ciliary scaffolding protein. The more specific scaffold protein binding term better captures its molecular role.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:36233334](PMID:36233334)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This IPI evidence from PMID:36233334 specifically demonstrates CFAP418-FAM161A interaction via Y2H, co-IP, and proximity ligation assays. The N-terminus of CFAP418 (aa 1-75) interacts with FAM161A's UPF0564 domain (aa 341-517). However, 'protein binding' is too generic - this should be annotated with a more specific molecular function reflecting its scaffolding role at the ciliary base.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:36233334 provides strong experimental evidence for direct CFAP418-FAM161A interaction through multiple biochemical assays, demonstrating CFAP418's role as a ciliary scaffold protein. Bioinformatics analysis reveals extensive coiled-coil regions (positions 0-259, 273-315, 350-441) supporting its scaffolding function.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001917|
|label[?](https://w3id.org/ai4curation/gene_review/label)|photoreceptor inner segment|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This ISS annotation (inferred from sequence similarity) is correct and supported by direct experimental evidence. Multiple studies confirm CFAP418 localization to the photoreceptor inner segment, particularly at the connecting cilium base. This represents a core localization for the protein's function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008594|
|label[?](https://w3id.org/ai4curation/gene_review/label)|photoreceptor cell morphogenesis|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This ISS annotation is supported by strong experimental evidence. CFAP418 knockout mice show severely disorganized photoreceptor outer segment discs from early postnatal development, demonstrating its requirement for proper photoreceptor morphogenesis. The more specific process 'photoreceptor cell outer segment organization' (GO:0035845) would better capture its primary role in outer segment disc formation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:22177090 demonstrated that CFAP418 knockout mice exhibit severe photoreceptor outer segment defects with disorganized disc morphology, indicating a specific role in outer segment organization rather than general morphogenesis.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IDA|[PMID:22177090](PMID:22177090)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|While PMID:22177090 does show cytoplasmic localization via IDA, this is overly general. The same paper specifically demonstrates localization at the ciliary base using immunohistochemistry. The ciliary base annotation (GO:0097546) from the same paper is more informative and specific.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0097546|
|label[?](https://w3id.org/ai4curation/gene_review/label)|ciliary base|
|IDA|[PMID:22177090](PMID:22177090)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Excellent annotation with strong IDA evidence. PMID:22177090 shows via immunohistochemistry that CFAP418 localizes at the base of primary cilia in RPE cells and at the base of connecting cilia in mouse photoreceptors. This is a core localization essential for CFAP418's function in photoreceptor maintenance and represents its primary functional compartment.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:22177090 provides definitive immunohistochemical evidence for CFAP418 localization at ciliary bases in both RPE cells and photoreceptors, establishing this as the protein's primary functional compartment.|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|
|---|---|---|---|---|
|Scaffolds protein complexes at photoreceptor connecting cilium base to organize ciliary trafficking machinery||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0097110|
|label[?](https://w3id.org/ai4curation/gene_review/label)|scaffold protein binding|
|||
|Organizes FAM161A-containing protein complex essential for photoreceptor outer segment disc morphogenesis||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0097110|
|label[?](https://w3id.org/ai4curation/gene_review/label)|scaffold protein binding|
|||
|Mediates ciliary trafficking required for photoreceptor outer segment protein transport and maintenance||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0097110|
|label[?](https://w3id.org/ai4curation/gene_review/label)|scaffold protein binding|
|||

## suggested_questions


|question|
|---|
|How does CFAP418 function in ciliary assembly and what specific role does it play in axoneme structure and stability?|
|What are the molecular interactions between CFAP418 and other ciliary proteins that are essential for proper cilia function?|
|How do mutations in CFAP418 contribute to ciliopathy phenotypes and what are the downstream cellular consequences?|
|What determines the tissue-specific expression pattern of CFAP418 and why is it particularly important in certain cell types?|

## suggested_experiments


|description|
|---|
|Cryo-electron microscopy of cilia from CFAP418-deficient cells to identify structural abnormalities in the axoneme|
|Proximity labeling proteomics to identify the complete CFAP418 interactome in ciliated cells|
|Live-cell imaging of ciliary assembly and disassembly to study CFAP418 dynamics during the cell cycle|
|Functional complementation studies using CFAP418 orthologs from different species to identify conserved functional domains|
