
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q9ULD6|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|INTU|
|description[?](https://w3id.org/ai4curation/gene_review/description)|INTU (Inturned) is a scaffold protein that functions as a core component of the CPLANE (ciliogenesis and planar polarity effector) complex at basal bodies, where it recruits intraflagellar transport machinery, specifically IFT-A proteins. INTU also serves as an adaptor linking ciliary proteins (NPHP4) to actin-modifying proteins (DAAM1) to control the subapical actin network required for basal body docking and ciliary orientation. Essential for primary cilia assembly and Hedgehog signaling, INTU is mutated in ciliopathies including Oral-Facial-Digital syndrome XVII and Short-Rib Thoracic Dysplasia 20.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|
|---|---|
|GO_REF:0000002|Gene Ontology annotation through association of InterPro records with GO terms.|
|GO_REF:0000024|Manual transfer of experimentally-verified manual GO annotation data to orthologs by curator judgment of sequence similarity.|
|GO_REF:0000033|Annotation inferences using phylogenetic trees|
|GO_REF:0000043|Gene Ontology annotation based on UniProtKB/Swiss-Prot keyword mapping|
|GO_REF:0000044|Gene Ontology annotation based on UniProtKB/Swiss-Prot Subcellular Location vocabulary mapping, accompanied by conservative changes to GO terms applied by UniProt.|
|GO_REF:0000052|Gene Ontology annotation based on curation of immunofluorescence data|
|GO_REF:0000107|Automatic transfer of experimentally verified manual GO annotation data to orthologs using Ensembl Compara.|
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.|
|PMID:26644512|The polarity protein Inturned links NPHP4 to Daam1 to control the subapical actin network in multiciliated cells.|
|PMID:27158779|The ciliopathy-associated CPLANE proteins direct basal body recruitment of intraflagellar transport machinery.|
|PMID:27173435|An organelle-specific protein landscape identifies novel diseases and molecular mechanisms.|
|PMID:33961781|Dual proteome-scale networks reveal cell-specific remodeling of the human interactome.|

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IBA annotation for cytoplasmic localization is supported by experimental evidence. INTU has cytosolic fractions when not assembled at cilia and likely shuttles between cytosol and ciliary base. The deep research confirms cytosolic localization (GO_REF:0000052).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|GO_REF:0000052 provides direct immunofluorescence evidence for cytoplasmic localization, and PMID:27158779 shows INTU can exist in cytosolic pools when not assembled at ciliary structures.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060271|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cilium assembly|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core function strongly supported by multiple experimental studies. INTU is essential for ciliogenesis through its role in the CPLANE complex recruiting IFT-A machinery to basal bodies (PMID:27158779). Mouse knockouts lack primary cilia, and human mutations cause ciliopathies.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:27158779 demonstrated INTU is a core component of the CPLANE complex that recruits IFT-A proteins to basal bodies, with knockout mice showing complete absence of primary cilia and human mutations causing OFD syndrome XVII.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005929|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cilium|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|While INTU localizes primarily at the basal body/ciliary base rather than within the cilium proper, this broader cellular component term is acceptable as INTU is functionally associated with ciliary structures. More specific localization would be ciliary basal body (GO:0036064).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:26644512 and PMID:27158779 demonstrate that INTU specifically localizes to ciliary basal bodies rather than within the cilium itself. GO:0036064 more accurately captures this specific localization.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007399|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nervous system development|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|While INTU mutations can affect neural development (neural tube defects, developmental delay), this is a consequence of defective ciliogenesis/Hedgehog signaling rather than a direct role in nervous system development. The term is too broad for the specific neural tube patterning defects observed.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001736|
|label[?](https://w3id.org/ai4curation/gene_review/label)|establishment of planar polarity|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU is a planar cell polarity effector protein that controls ciliary orientation through its effects on basal body positioning and the subapical actin network. This is well-supported by experimental evidence showing INTU controls rotational polarity of cilia in multiciliated cells (PMID:26644512).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:26644512 demonstrated that INTU is essential for establishing planar cell polarity by linking NPHP4 to DAAM1 to control the subapical actin network required for proper ciliary orientation in multiciliated cells.|

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
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|No evidence supports a direct role in vesicle-mediated transport. This appears to be an incorrect automated annotation, possibly based on superficial similarity to IFT proteins. INTU functions in intraciliary transport, not vesicle transport.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060271|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cilium assembly|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate of the IBA annotation for cilium assembly. The automated annotation correctly identifies this core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate of the IBA annotation for cytoplasm. Correctly identifies cytoplasmic localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005814|
|label[?](https://w3id.org/ai4curation/gene_review/label)|centriole|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU localizes to the mother centriole/basal body area. While technically correct, the more specific term ciliary basal body (GO:0036064) better captures INTU localization in the context of ciliogenesis.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005856|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoskeleton|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|While INTU interacts with actin cytoskeleton components through DAAM1 and controls the subapical actin network, it is not itself a cytoskeletal protein. This term is too broad; INTU specifically regulates actin organization at the apical cortex.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009986|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell surface|
|IEA|[GO_REF:0000044](GO_REF:0000044)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Some INTU truncations show enhanced membrane association, and INTU may associate with membrane through predicted phosphatidylinositol binding. However, cell surface is too general; apical plasma membrane would be more accurate given INTU functions at the apical cortex.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030030|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell projection organization|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Cilia are cell projections, and INTU is essential for their organization. However, cilium assembly (GO:0060271) is more specific and informative for INTU function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042995|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell projection|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Too vague for a cellular component term. INTU localizes specifically to ciliary basal bodies, not broadly to cell projections.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:26644512](PMID:26644512)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper shows INTU binds NPHP4 and DAAM1, forming a ternary complex. While protein binding is correct, it is uninformative. INTU functions as a scaffold/adaptor protein linking ciliary and cytoskeletal proteins. A more specific molecular function term would be protein-protein adaptor activity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:26644512 demonstrated that INTU functions as an adaptor protein that specifically links the ciliary protein NPHP4 to the actin-regulating protein DAAM1, mediating communication between ciliary and cytoskeletal systems.|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:27173435](PMID:27173435)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Protein interaction study showing INTU interactions. Again, protein binding is too vague. INTU acts as a scaffold in protein complexes.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:33961781](PMID:33961781)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Another protein interaction network study. The generic protein binding term should be replaced with the more informative scaffold/adaptor function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007224|
|label[?](https://w3id.org/ai4curation/gene_review/label)|smoothened signaling pathway|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU affects Hedgehog signaling indirectly through its essential role in ciliogenesis. Without cilia, Hh signaling is disrupted. This is a valid annotation but represents an indirect effect.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007399|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nervous system development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate of earlier annotation. Too broad; neural tube development is more specific.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008589|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of smoothened signaling pathway|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU regulates Hedgehog signaling indirectly through its requirement for ciliogenesis. Valid but represents secondary effect of ciliary dysfunction.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0010839|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of keratinocyte proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Tissue-specific consequence of INTU loss affecting hair follicle development through disrupted Hh signaling. Too specific for a general annotation; represents a downstream developmental effect.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0021513|
|label[?](https://w3id.org/ai4curation/gene_review/label)|spinal cord dorsal/ventral patterning|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Neural patterning defects occur in INTU mutants due to disrupted Hedgehog signaling from lack of cilia. This is a downstream developmental consequence, not a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0021915|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neural tube development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Neural tube defects are observed in INTU mutants, supported by experimental evidence in multiple species. Valid developmental consequence of ciliary dysfunction.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030216|
|label[?](https://w3id.org/ai4curation/gene_review/label)|keratinocyte differentiation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU affects hair follicle differentiation through its role in ciliogenesis and Hh signaling. Tissue-specific developmental effect.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030278|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of ossification|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU mutations cause skeletal defects including delayed ossification through disrupted Indian Hedgehog signaling in growth plates. Valid but represents developmental consequence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0031069|
|label[?](https://w3id.org/ai4curation/gene_review/label)|hair follicle morphogenesis|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Experimental evidence shows INTU is required for hair follicle development through cilia-dependent Hh signaling. Tissue-specific developmental effect.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0033365|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein localization to organelle|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU recruits IFT-A proteins to basal bodies as part of the CPLANE complex. This is a core function but could be more specific as intraciliary transport.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035091|
|label[?](https://w3id.org/ai4curation/gene_review/label)|phosphatidylinositol binding|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Bioinformatically predicted phosphatidylinositol binding capacity, which could mediate membrane association at the ciliary base. Reasonable prediction but lacks experimental validation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035869|
|label[?](https://w3id.org/ai4curation/gene_review/label)|ciliary transition zone|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU localizes near but not exactly at the transition zone. It is primarily at the basal body with some extension toward the transition zone. Basal body is more accurate.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042733|
|label[?](https://w3id.org/ai4curation/gene_review/label)|embryonic digit morphogenesis|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU mutations cause polydactyly in humans and mice. Well-supported developmental consequence of disrupted Hh signaling during limb development.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0044458|
|label[?](https://w3id.org/ai4curation/gene_review/label)|motile cilium assembly|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU is required for assembly of both motile and primary cilia. Evidence from multiciliated cells shows INTU localizes to basal bodies of motile cilia and controls their polarization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045880|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of smoothened signaling pathway|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU enables Hedgehog signaling by building the cilia required for signal transduction. This is an indirect positive effect through ciliogenesis.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0051301|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell division|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|No evidence supports a direct role in cell division. Centrioles are involved in both ciliogenesis and cell division, but INTU functions specifically in the ciliary context.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0051782|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of cell division|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|No evidence for INTU regulating cell division. This appears to be an incorrect inference.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060173|
|label[?](https://w3id.org/ai4curation/gene_review/label)|limb development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU mutations cause limb defects including polydactyly and shortened limbs. Valid developmental process affected by ciliary dysfunction.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1905515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|non-motile cilium assembly|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU is essential for primary (non-motile) cilium assembly. Mouse knockouts lack primary cilia.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005829|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytosol|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Direct experimental evidence for cytosolic localization by immunofluorescence. INTU has cytosolic fractions when not assembled at cilia.|
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
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Direct experimental evidence for basal body localization. This is the primary and most specific localization for INTU, strongly supported by multiple studies.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001736|
|label[?](https://w3id.org/ai4curation/gene_review/label)|establishment of planar polarity|
|NAS|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper identifies INTU as part of the CPLANE complex controlling planar polarity through ciliary orientation. Well-supported core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005929|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cilium|
|NAS|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU localizes at the ciliary base/basal body rather than within the cilium proper. More specific term would be ciliary basal body.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0021915|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neural tube development|
|NAS|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper shows neural tube defects in INTU mutants. Valid developmental consequence of ciliary dysfunction.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042073|
|label[?](https://w3id.org/ai4curation/gene_review/label)|intraciliary transport|
|NAS|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper shows INTU recruits IFT-A machinery to basal bodies as part of the CPLANE complex. Core function directly supported by experimental evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:1902017|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of cilium assembly|
|NAS|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|INTU is not just a regulator but is essential for cilium assembly itself. The more specific term cilium assembly (GO:0060271) better captures its essential role.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035091|
|label[?](https://w3id.org/ai4curation/gene_review/label)|phosphatidylinositol binding|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Sequence similarity-based prediction for phosphatidylinositol binding, which could mediate membrane association. Reasonable but lacks direct experimental validation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042733|
|label[?](https://w3id.org/ai4curation/gene_review/label)|embryonic digit morphogenesis|
|IMP|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper shows polydactyly in INTU mutant mice. Direct experimental evidence for developmental defect.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043587|
|label[?](https://w3id.org/ai4curation/gene_review/label)|tongue morphogenesis|
|IMP|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper shows tongue defects (lobulated tongues, hamartomas) in INTU mutants, consistent with OFD syndrome. Direct experimental evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060021|
|label[?](https://w3id.org/ai4curation/gene_review/label)|roof of mouth development|
|IMP|[PMID:27158779](PMID:27158779)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper shows high-arched palate in INTU mutant mice, characteristic of OFD syndrome. Direct experimental evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0031514|
|label[?](https://w3id.org/ai4curation/gene_review/label)|motile cilium|
|IDA|[PMID:26644512](PMID:26644512)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper shows INTU localization at basal bodies of motile cilia in multiciliated cells. However, INTU is at the basal body, not in the motile cilium itself.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0036064|
|label[?](https://w3id.org/ai4curation/gene_review/label)|ciliary basal body|
|IDA|[PMID:26644512](PMID:26644512)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This paper directly shows INTU localization at basal bodies in multiciliated cells. Strong experimental evidence for this core localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

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
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Sequence similarity-based annotation for cytoplasmic localization, consistent with experimental evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007399|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nervous system development|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Sequence similarity-based annotation. Too broad; neural tube development is more specific.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008589|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of smoothened signaling pathway|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Sequence similarity-based annotation. INTU affects Hh signaling indirectly through ciliogenesis.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060173|
|label[?](https://w3id.org/ai4curation/gene_review/label)|limb development|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Sequence similarity-based annotation. Supported by polydactyly phenotypes in mutants.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0060271|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cilium assembly|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Sequence similarity-based annotation for core function. Well-supported by experimental evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|

## core_functions


|description|molecular_function|directly_involved_in|locations|
|---|---|---|---|
|Scaffolds CPLANE complex assembly at basal bodies to recruit IFT-A machinery for cilium biogenesis|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030674|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein-macromolecule adaptor activity|
|||
|Bridges NPHP4-DAAM1 interaction to organize subapical actin network required for basal body docking|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030674|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein-macromolecule adaptor activity|
|||
|Enables phosphatidylinositol binding to anchor protein complexes at ciliary membrane domains|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035091|
|label[?](https://w3id.org/ai4curation/gene_review/label)|phosphatidylinositol binding|
|||

## suggested_questions


|question|
|---|
|How does INTU regulate intraflagellar transport and what specific cargo does it help transport within cilia?|
|What are the molecular mechanisms by which INTU coordinates ciliary assembly with cell cycle progression?|
|How do mutations in INTU lead to left-right asymmetry defects and what role does it play in nodal cilia function?|
|What determines the specificity of INTU interactions with different intraflagellar transport complexes?|

## suggested_experiments


|description|
|---|
|Super-resolution microscopy to track INTU and IFT particle movements along the ciliary axoneme with nanometer precision|
|Biochemical reconstitution of IFT complexes containing INTU to study cargo loading and transport mechanisms in vitro|
|Developmental analysis of left-right patterning in INTU mutant embryos using whole-mount in situ hybridization|
|Proteomics analysis of INTU-associated complexes during different stages of ciliogenesis|
