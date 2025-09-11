
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q587I9|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|SFT2D3|
|description[?](https://w3id.org/ai4curation/gene_review/description)|SFT2D3 encodes vesicle transport protein SFT2C, an integral membrane component of the evolutionarily conserved SFT2 family. This tetra-span transmembrane protein (containing the pfam04178 Got1/Sft2 domain) functions as a vesicle fusion facilitator in the retrograde transport pathway from endosomes to Golgi apparatus. The protein localizes to Golgi membranes and endosomal compartments where it likely organizes or stabilizes SNARE complexes during vesicle fusion events. Based on strong experimental evidence from paralog SFT2D2 and yeast homologs Sft2p/Got1p, SFT2D3 operates as a SNARE-associated factor that promotes fusion of retrograde transport vesicles derived from endocytic compartments with the trans-Golgi network. The gene shows ubiquitous but low-level expression consistent with a fundamental housekeeping role in cellular trafficking machinery.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|findings|
|---|---|---|
|GO_REF:0000002|Gene Ontology annotation through association of InterPro records with GO terms.||
|GO_REF:0000043|Gene Ontology annotation based on UniProtKB/Swiss-Prot keyword mapping||
|GO_REF:0000117|Electronic Gene Ontology annotations created by ARBA machine learning models||
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.||
|UniProt:Q587I9|UniProtKB entry for human SFT2D3 (vesicle transport protein SFT2C)||
|file:human/SFT2D3/SFT2D3-deep-research.md|Deep research analysis of SFT2D3 function and literature||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016020|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This membrane annotation is accurate and appropriate for SFT2D3. The protein is clearly an integral membrane protein with 4 transmembrane domains as shown in UniProt topology data and confirmed by structural predictions. The generic \"membrane\" term is suitable as it correctly captures the protein localization without being overly specific about particular membrane types. Given SFT2D3 localizes to multiple membrane compartments (Golgi, endosomes), this broad term is acceptable.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016192|
|label[?](https://w3id.org/ai4curation/gene_review/label)|vesicle-mediated transport|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is accurate and well-supported by strong evidence. SFT2D3 is clearly involved in vesicle-mediated transport based on: (1) UniProt annotation stating involvement in fusion of retrograde transport vesicles; (2) experimental evidence from SFT2D2 paralog showing role in endosome-to-Golgi retrieval; (3) yeast homolog studies demonstrating vesicle fusion function; (4) evolutionary conservation of the SFT2 family in vesicle trafficking. This represents a core function and should be accepted as-is.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IEA|[GO_REF:0000117](GO_REF:0000117)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is technically correct but overly general and not informative for SFT2D3 function. While SFT2D3 is localized to cytoplasmic organelles (not nucleus or extracellular), this broad term fails to capture its specific subcellular localization. Based on functional evidence and paralog localization studies, SFT2D3 specifically localizes to Golgi apparatus and endosomal membranes. This should be modified to the more specific and functionally relevant "Golgi apparatus" (GO:0005794).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0012505|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endomembrane system|
|IEA|[GO_REF:0000117](GO_REF:0000117)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is accurate and well-supported. SFT2D3 is clearly a component of the endomembrane system as an integral membrane protein involved in trafficking between endosomes and Golgi. This term appropriately captures its localization to intracellular membrane-bound organelles within the secretory/endocytic pathway network. The term is at the right level of specificity - more general than specific organelles but more specific than just "membrane" or "cytoplasm". This should be accepted as representing accurate cellular localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0015031|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein transport|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is accurate but could be more specific. SFT2D3 is indeed involved in protein transport, particularly in the context of retrograde transport where membrane proteins and sorting receptors are retrieved from endosomes back to the Golgi. However, the evidence strongly points to a more specific role in "retrograde transport, endosome to Golgi" (GO:0042147). While the current term is not wrong, the more specific term better captures the precise biological process SFT2D3 participates in based on functional evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms


#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|
|---|---|---|---|---|
|Mediating vesicle fusion during retrograde transport from endosomes to trans-Golgi network||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000149|
|label[?](https://w3id.org/ai4curation/gene_review/label)|SNARE binding|
|||
|Organizing SNARE complex assembly and recycling at Golgi membranes||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000149|
|label[?](https://w3id.org/ai4curation/gene_review/label)|SNARE binding|
|||
|Anchoring as tetra-span membrane protein in Golgi and endosomal compartments||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|||

## suggested_questions


|question|
|---|
|What is the specific function of SFT2D3 in intracellular membrane trafficking and organelle biogenesis?|
|How does SFT2D3 interact with the SNARE machinery and other trafficking regulators?|
|What determines the subcellular localization of SFT2D3 and how does this relate to its function?|
|How is SFT2D3 expression and activity regulated in different cell types and conditions?|

## suggested_experiments


|description|
|---|
|Live-cell imaging using fluorescently tagged SFT2D3 to study its dynamics and localization in membrane trafficking|
|Proteomics approaches to identify SFT2D3 interacting partners in different subcellular compartments|
|Functional analysis using CRISPR knockout to determine the cellular processes dependent on SFT2D3|
|Electron microscopy analysis of cells lacking SFT2D3 to identify ultrastructural defects in organelle organization|
