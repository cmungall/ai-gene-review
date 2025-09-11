
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q9UQM7|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|CAMK2A|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|CAMKA|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|KIAA0968|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|CaMKIINalpha|
|description[?](https://w3id.org/ai4curation/gene_review/description)|CAMK2A encodes the alpha subunit of Ca2+/calmodulin-dependent protein kinase II (CaMKII), a multifunctional Ser/Thr protein kinase that serves as a molecular switch in synaptic plasticity and memory formation. CaMKIIα is highly enriched in the brain, particularly at excitatory synapses, where it responds to calcium influx through NMDA receptors. Upon Ca2+/calmodulin binding, the kinase undergoes autophosphorylation at Thr286, generating Ca2+-independent activity that persists after calcium levels return to baseline - effectively storing a molecular memory of synaptic activity. The kinase forms dodecameric holoenzymes that phosphorylate numerous synaptic substrates including glutamate receptors, thereby strengthening synaptic transmission during long-term potentiation. Mutations in CAMK2A cause intellectual disability (MRD53) and are associated with autism spectrum disorder, highlighting its critical role in cognitive development and function.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|findings|
|---|---|---|
|GO_REF:0000002|Gene Ontology annotation through association of InterPro records with GO terms.||
|GO_REF:0000024|Manual transfer of experimentally-verified manual GO annotation data to orthologs by curator judgment of sequence similarity.||
|GO_REF:0000033|Annotation inferences using phylogenetic trees||
|GO_REF:0000043|Gene Ontology annotation based on UniProtKB/Swiss-Prot keyword mapping||
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.||
|PMID:11972023|Requirement of Ca2+ and CaMKII for Stat1 Ser-727 phosphorylation in response to IFN-gamma.||
|PMID:16257975|The conserved Leu-724 residue is required for both serine phosphorylation and co-activator recruitment for Stat1-mediated transcription activation in response to interferon-gamma.||
|PMID:17052756|Bcl10 is phosphorylated on Ser138 by Ca2+/calmodulin-dependent protein kinase II.||
|PMID:19453375|Phosphorylation status of the NR2B subunit of NMDA receptor regulates its interaction with calcium/calmodulin-dependent protein kinase II.||
|PMID:20668654|Structure of the CaMKIIdelta/calmodulin complex reveals the molecular mechanism of CaMKII kinase activation.||
|PMID:22939624|Quantitative analysis of HSP90-client interactions reveals principles of substrate recognition.||
|PMID:25852190|Integrative analysis of kinase networks in TRAIL-induced apoptosis provides a source of potential targets for combination therapy.||
|PMID:27173435|An organelle-specific protein landscape identifies novel diseases and molecular mechanisms.||
|PMID:28130356|A Novel Human CAMK2A Mutation Disrupts Dendritic Morphology and Synaptic Transmission, and Causes ASD-Related Behaviors.||
|PMID:28753426|Methyltransferase SETD2-Mediated Methylation of STAT1 Is Critical for Interferon Antiviral Activity.||
|PMID:29426014|Network Analysis of UBE3A/E6AP-Associated Proteins Provides Connections to Several Distinct Cellular Processes.||
|PMID:31980649|Extensive rewiring of the EGFR network in colorectal cancer cells expressing transforming levels of KRAS(G13D).||
|PMID:32296183|A reference map of the human binary protein interactome.||
|PMID:32707033|Kinase Interaction Network Expands Functional and Disease Roles of Human Kinases.||
|PMID:32814053|Interactome Mapping Provides a Network of Neurodegenerative Disease Proteins and Uncovers Widespread Protein Aggregation in Affected Brains.||
|PMID:35568036|A family of conserved bacterial virulence factors dampens interferon responses by blocking calcium signaling.||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004683|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium/calmodulin-dependent protein kinase activity|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Core molecular function strongly supported by extensive biochemical and functional evidence. Multiple studies demonstrate direct Ca2+/calmodulin-dependent kinase activity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004672|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein kinase activity|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|General protein kinase activity is correct but too broad. The more specific Ca2+/calmodulin-dependent kinase activity (GO:0004683) better captures the core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004683|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium/calmodulin-dependent protein kinase activity|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation with different evidence code. The core function is captured by the IBA annotation above.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000166|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nucleotide binding|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|General nucleotide binding is implied by kinase activity but too broad. More specific ATP binding (GO:0005524) would be more appropriate, though still not core to CAMK2A function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004674|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein serine/threonine kinase activity|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Correct Ser/Thr kinase activity but the Ca2+/calmodulin-dependent aspect is essential to CAMK2A function. Better captured by GO:0004683.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:19453375](PMID:19453375)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|General protein binding annotation. PMID:19453375 shows interaction with NMDAR NR2B subunit - more specific annotations for NMDAR binding would be more informative.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:20668654](PMID:20668654)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:20668654 describes CaMKII/calmodulin complex structure. Calmodulin binding is already captured by the Ca2+/calmodulin-dependent kinase activity annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:22939624](PMID:22939624)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:22939624 is about HSP90 client interactions. While CAMK2A may interact with HSP90, this is not a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:25852190](PMID:25852190)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Proteomic study of TRAIL-induced apoptosis showing various interactions. Not specific to CAMK2A core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:27173435](PMID:27173435)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Organelle proteomics study. General protein-protein interaction, not specific to CAMK2A function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:29426014](PMID:29426014)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|UBE3A network analysis. Not central to CAMK2A core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:31980649](PMID:31980649)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Large-scale interactome study. General protein binding annotations from high-throughput studies are not specific to core CAMK2A functions.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:32296183](PMID:32296183)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Human reference interactome map. High-throughput data, not specific to CAMK2A core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:32707033](PMID:32707033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Kinase interaction network study. While informative for kinase networks, generic protein binding is not core.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:32814053](PMID:32814053)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Neurodegenerative disease protein interactome. Not specific to CAMK2A core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004683|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium/calmodulin-dependent protein kinase activity|
|IDA|[PMID:35568036](PMID:35568036)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Direct experimental evidence showing Ca2+/calmodulin-dependent kinase activity. Strong support for core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004683|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium/calmodulin-dependent protein kinase activity|
|IDA|[PMID:11972023](PMID:11972023)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Direct evidence of CaMKII activity in STAT1 phosphorylation. Core molecular function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004683|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium/calmodulin-dependent protein kinase activity|
|TAS|[PMID:11972023](PMID:11972023)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation with different evidence type. Core function already captured.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007259|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell surface receptor signaling pathway via JAK-STAT|
|IDA|[PMID:11972023](PMID:11972023)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CAMK2A phosphorylates STAT1 at Ser727 in response to IFN-gamma, participating in JAK-STAT signaling. While not the primary function, this is a well-documented role.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007259|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell surface receptor signaling pathway via JAK-STAT|
|IDA|[PMID:28753426](PMID:28753426)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation for JAK-STAT signaling. Already captured above.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035458|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cellular response to interferon-beta|
|IDA|[PMID:28753426](PMID:28753426)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|CAMK2A participates in interferon-beta response. While documented, this is not the core neuronal function of CAMK2A.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0046427|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of receptor signaling pathway via JAK-STAT|
|IDA|[PMID:16257975](PMID:16257975)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Positive regulation of JAK-STAT signaling through STAT1 phosphorylation. Related to interferon response, not core synaptic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0071346|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cellular response to type II interferon|
|IDA|[PMID:11972023](PMID:11972023)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Response to IFN-gamma (type II interferon) via STAT1 phosphorylation. Documented but not core neuronal function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004674|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein serine/threonine kinase activity|
|IDA|[PMID:28130356](PMID:28130356)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Direct evidence of Ser/Thr kinase activity. More specific Ca2+/calmodulin-dependent annotation is preferred.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:28130356](PMID:28130356)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:28130356 shows CAMK2A interactions with Shank3, NMDAR, and L-type Ca channels - all critical synaptic partners. Should have more specific annotations.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004674|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein serine/threonine kinase activity|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Inferred from sequence similarity. General Ser/Thr kinase activity is correct but Ca2+/calmodulin-dependence is key.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:17052756](PMID:17052756)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|PMID:17052756 shows CAMK2 phosphorylates Bcl10. While valid, this is not a core synaptic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|in_complex|
|---|---|---|---|---|---|
|Ca2+/calmodulin-activated protein kinase activity phosphorylating synaptic substrates||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004683|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium/calmodulin-dependent protein kinase activity|
||||
|Autophosphorylation at Thr286 generating Ca2+-independent kinase activity||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004683|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium/calmodulin-dependent protein kinase activity|
||||
|Assembly into dodecameric holoenzyme complex enabling cooperative activation||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005954|
|label[?](https://w3id.org/ai4curation/gene_review/label)|calcium- and calmodulin-dependent protein kinase complex|
|
|STAT1 Ser727 phosphorylation in response to interferon-gamma signaling||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004674|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein serine/threonine kinase activity|
||||

## suggested_questions


|question|
|---|
|How does autophosphorylation of CAMK2A create molecular memory and contribute to synaptic plasticity and learning?|
|What determines the subcellular localization of CAMK2A and how does this regulate its access to different substrates?|
|How do different splice variants of CAMK2A contribute to brain region-specific functions and neuronal plasticity?|
|What are the mechanisms by which CAMK2A integrates calcium signals with other signaling pathways during synaptic transmission?|

## suggested_experiments


|description|
|---|
|Two-photon calcium imaging combined with optogenetics to study CAMK2A activation dynamics in dendritic spines during synaptic plasticity|
|Cryo-EM structural analysis of CAMK2A holoenzymes in different activation states to understand autophosphorylation mechanisms|
|Single-molecule tracking of CAMK2A in live neurons to characterize its mobility and clustering at synaptic sites|
|Proteomics identification of context-specific CAMK2A substrates using chemical crosslinking and mass spectrometry|
