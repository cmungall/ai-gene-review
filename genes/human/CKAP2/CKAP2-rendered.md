
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q8WWK9|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|CKAP2|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|TMAP|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|LB1|
|description[?](https://w3id.org/ai4curation/gene_review/description)|CKAP2 is a potent microtubule-associated protein that functions as the most powerful known microtubule growth factor and stabilizer. Essential for faithful chromosome segregation during mitosis, CKAP2 promotes microtubule nucleation (~100-fold enhancement), dramatically increases growth rates (~20%), and suppresses catastrophic depolymerization. Required for proper mitotic spindle organization, bipolar spindle assembly, and genomic stability. Cell cycle-regulated protein that peaks at G2/M phase and is degraded by APC/C during mitotic exit.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|findings|
|---|---|---|
|GO_REF:0000033|Annotation inferences using phylogenetic trees||
|GO_REF:0000043|Gene Ontology annotation based on UniProtKB/Swiss-Prot keyword mapping||
|GO_REF:0000044|Gene Ontology annotation based on UniProtKB/Swiss-Prot Subcellular Location vocabulary mapping, accompanied by conservative changes to GO terms applied by UniProt.||
|GO_REF:0000052|Gene Ontology annotation based on curation of immunofluorescence data||
|GO_REF:0000107|Automatic transfer of experimentally verified manual GO annotation data to orthologs using Ensembl Compara.||
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.||
|PMID:16061649|Ckap2 regulates aneuploidy, cell cycling, and cell death in a p53-dependent manner.||
|PMID:21399614|Novel asymmetrically localizing components of human centrosomes identified by complementary proteomics methods.||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0015630|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule cytoskeleton|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CKAP2 is directly involved with the microtubule cytoskeleton as a microtubule-binding protein that stabilizes and promotes microtubule growth. This is a core cellular component annotation supported by extensive experimental evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007026|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of microtubule depolymerization|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately reflects CKAP2s core function. Research shows CKAP2 strongly suppresses microtubule catastrophes and stabilizes microtubules against depolymerization. This represents a core molecular function of CKAP2.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000922|
|label[?](https://w3id.org/ai4curation/gene_review/label)|spindle pole|
|IEA|[GO_REF:0000044](GO_REF:0000044)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Experimental evidence strongly supports CKAP2 localization at spindle poles during mitosis. CKAP2 is required for maintaining focused spindle poles and concentrating microtubule minus-ends at spindle poles. This is a well-supported cellular component annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CKAP2 is indeed present in the cytoplasm, particularly during interphase when it is largely diffuse in the cytosol. However, this is a very general localization term that does not capture CKAP2s specific and functionally important localizations. This annotation is correct but not informative.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005819|
|label[?](https://w3id.org/ai4curation/gene_review/label)|spindle|
|IEA|[GO_REF:0000044](GO_REF:0000044)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CKAP2 is strongly enriched at the mitotic spindle during mitosis and is required for proper spindle assembly and maintenance. Multiple studies show CKAP2 localizes to spindle microtubules and is essential for spindle organization. This is a core cellular component annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005856|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoskeleton|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CKAP2 is indeed a cytoskeleton-associated protein, specifically a microtubule-associated protein (MAP). However, this term is too general - the more specific microtubule cytoskeleton term better captures CKAP2s function. This annotation is correct but lacks specificity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005874|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CKAP2 directly associates with microtubules as a microtubule-binding protein. This is supported by extensive experimental evidence showing CKAP2 binds to and stabilizes microtubules, particularly during mitosis. This is a core cellular component annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006915|
|label[?](https://w3id.org/ai4curation/gene_review/label)|apoptotic process|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:16061649 explicitly states that "In p53-competent cells, Ckap2 does not induce tetraploidy but activates p53-mediated cell cycle arrest and apoptosis." This demonstrates that CKAP2 overexpression can directly activate apoptosis in a p53-dependent manner. However, this is based on overexpression studies and may not represent normal physiological function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007026|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of microtubule depolymerization|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a duplicate of the IBA annotation above with the same GO term. The function is well-supported by experimental evidence showing CKAP2 suppresses microtubule catastrophes and stabilizes microtubules. This represents a core biological process for CKAP2.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0015630|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule cytoskeleton|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a duplicate of the IBA annotation above with the same GO term. CKAP2 is directly involved with the microtubule cytoskeleton as a microtubule-binding and stabilizing protein. This is a core cellular component annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005730|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nucleolus|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is based on immunofluorescence data (GO_REF:0000052), but neither PMID:16061649 nor PMID:21399614 mention nucleolar localization. PMID:16061649 shows CKAP2 is absent in G1 phase when nucleolar localization might be expected. This localization is not functionally significant compared to CKAP2s core mitotic functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005813|
|label[?](https://w3id.org/ai4curation/gene_review/label)|centrosome|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CKAP2 localizes to centrosomes during interphase and mitosis and is required for controlling centrosome duplication and ensuring proper centrosome function. This is a core cellular component annotation representing a functionally important subcellular localization for CKAP2.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005929|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cilium|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is based on immunofluorescence data (GO_REF:0000052), but neither PMID:16061649 nor PMID:21399614 mention ciliary localization or function. CKAP2s well-characterized functions are in mitotic spindle assembly and chromosome segregation. This appears to be an over-annotation possibly due to cross-reactivity or mislocalization artifacts.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0015630|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule cytoskeleton|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a duplicate annotation with the same GO term but based on immunofluorescence data. CKAP2 is indeed directly involved with the microtubule cytoskeleton as demonstrated by multiple experimental approaches. This is a core cellular component annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0036064|
|label[?](https://w3id.org/ai4curation/gene_review/label)|ciliary basal body|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is based on immunofluorescence data (GO_REF:0000052), but neither PMID:16061649 nor PMID:21399614 provide evidence for ciliary basal body localization or function. CKAP2s well-characterized functions are in mitotic spindle assembly and chromosome segregation, not ciliary processes. This appears to be an over-annotation possibly due to immunofluorescence artifacts.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0072686|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mitotic spindle|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CKAP2 is strongly enriched at the mitotic spindle and is essential for mitotic spindle assembly and maintenance. This is more specific than the general spindle term and accurately reflects CKAP2s core function during mitosis. This is a core cellular component annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0015630|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule cytoskeleton|
|IDA|[PMID:16061649](PMID:16061649)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:16061649 directly states that "Overexpressed Ckap2 colocalizes with and stabilizes microtubules" and shows CKAP2 localizes "to the spindle during mitosis" and "colocalizes with k-fibers during mitosis." This provides strong experimental evidence for microtubule cytoskeleton involvement.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000281|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mitotic cytokinesis|
|IGI|[PMID:16061649](PMID:16061649)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:16061649 states that overexpression of Ckap2 "suggests disturbed mitosis and cytokinesis" and shows that in p53-null cells, CKAP2 overexpression "induces tetraploidy with aberrant centrosome numbers, suggesting disturbed mitosis and cytokinesis." This suggests CKAP2 dysfunction affects cytokinesis, supporting this annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

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
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|After thoroughly reading PMID:16061649, there is absolutely no mention of CKAP2 involvement in transcriptional regulation or RNA polymerase II. The publication focuses entirely on CKAP2s role in microtubule dynamics, chromosome segregation, and cell cycle regulation. This annotation is not supported by the referenced publication.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005813|
|label[?](https://w3id.org/ai4curation/gene_review/label)|centrosome|
|IDA|[PMID:21399614](PMID:21399614)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:21399614 is a general centrosome proteomics study that identified CKAP2 among "126 known and 40 candidate centrosomal proteins" but provides no specific functional data about CKAP2. This supports centrosome localization but without detailed functional context. This is a core cellular component annotation based on proteomics identification.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|
|---|---|---|---|---|
|Promotes microtubule nucleation and polymerization while suppressing catastrophes||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008017|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule binding|
|||
|Organizes bipolar mitotic spindle assembly and maintains spindle pole integrity||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008017|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule binding|
|||
|Ensures accurate chromosome segregation and mitotic progression||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008017|
|label[?](https://w3id.org/ai4curation/gene_review/label)|microtubule binding|
|||

## suggested_questions


|question|
|---|
|How does CKAP2 coordinate with other cytoskeletal proteins to regulate microtubule dynamics during cell division?|
|What determines the cell cycle-dependent expression and phosphorylation of CKAP2?|
|How does CKAP2 contribute to proper chromosome segregation and what are the consequences of its dysfunction?|
|What role does CKAP2 play in non-mitotic cells and how does it affect microtubule organization in interphase?|

## suggested_experiments


|description|
|---|
|Live-cell imaging of fluorescently tagged CKAP2 to study its dynamics during mitosis and cell division|
|Cryo-electron tomography of mitotic spindles in CKAP2-depleted cells to visualize microtubule organization defects|
|Proteomics analysis to identify CKAP2 interacting partners and phosphorylation sites throughout the cell cycle|
|Single-molecule biophysics to study CKAP2 interactions with microtubules and effects on microtubule stability|
