
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|P35498|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|SCN1A|
|description[?](https://w3id.org/ai4curation/gene_review/description)|SCN1A encodes the pore-forming α-subunit of Nav1.1, the most clinically important voltage-gated sodium channel in epilepsy genetics. As a 24-transmembrane domain protein, SCN1A directly mediates the depolarizing phase of action potentials through voltage-dependent conformational switching that allows selective Na+ influx along electrochemical gradients. Nav1.1 is particularly critical in GABAergic interneurons where it regulates inhibitory neuron excitability - loss-of-function mutations impair interneuron function leading to network hyperexcitability and seizures. SCN1A mutations cause >80% of Dravet syndrome cases and represent the most frequent target of epilepsy-related mutations, causing a phenotypic spectrum from febrile seizures (GEFS+) to severe developmental epileptic encephalopathy. The protein localizes to specialized membrane domains including axon initial segments, nodes of Ranvier, and neuronal cell bodies where it enables action potential initiation, propagation, and presynaptic membrane potential regulation. Understanding SCN1A pathophysiology has revealed that sodium channel blockers are contraindicated in SCN1A-related epilepsies as they worsen seizures by further impairing already compromised interneuron excitability.|

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
|GO_REF:0000052|Gene Ontology annotation based on curation of immunofluorescence data||
|GO_REF:0000107|Automatic transfer of experimentally verified manual GO annotation data to orthologs using Ensembl Compara.||
|GO_REF:0000108|Automatic assignment of GO terms using logical inference, based on on inter-ontology links.||
|GO_REF:0000117|Electronic Gene Ontology annotations created by ARBA machine learning models||
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.||
|PMID:10742094|Mutations of SCN1A, encoding a neuronal sodium channel, in two families with GEFS+2.||
|PMID:14672992|Epilepsy-associated dysfunction in the voltage-gated neuronal sodium channel SCN1A.||
|PMID:22150645|Pure haploinsufficiency for Dravet syndrome Na(V)1.1 (SCN1A) sodium channel truncating mutations.||
|PMID:27207958|Variants of Transient Receptor Potential Melastatin Member 4 in Childhood Atrioventricular Block.||
|UniProt:P35498|SCN1A UniProt functional annotation||
|clinical_literature|SCN1A clinical significance in epilepsy genetics||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001518|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel complex|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation correctly identifies SCN1A as part of the voltage-gated sodium channel complex. As the pore-forming α-subunit, SCN1A associates with auxiliary β-subunits (SCN1B, SCN2B, SCN3B, SCN4B) to form the functional Nav1.1 complex. IBA evidence from phylogenetic analysis provides strong support for this cellular component annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005248|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel activity|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is the core molecular function of SCN1A. The protein directly mediates voltage-gated sodium channel activity through conformational changes that allow selective Na+ passage. This annotation captures the essential enzymatic function supported by UniProt catalytic activity data (Na+(in) = Na+(out)) and extensive experimental evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0035725|
|label[?](https://w3id.org/ai4curation/gene_review/label)|sodium ion transmembrane transport|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This accurately describes the biological process mediated by SCN1A. The protein enables sodium ion transport across neuronal membranes, which is fundamental for action potential generation and propagation. IBA evidence from phylogenetic analysis strongly supports this core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0086002|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cardiac muscle cell action potential involved in contraction|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|While Nav1.1 sodium channels may be expressed in cardiac tissue, SCN1A is predominantly and specifically important in neuronal tissues, particularly GABAergic interneurons. Cardiac action potentials are primarily mediated by SCN5A (Nav1.5). This annotation represents over-annotation for SCN1A, which should be marked as peripheral to its core neuronal function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0099505|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of presynaptic membrane potential|
|IEA|[GO_REF:0000108](GO_REF:0000108)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is consistent with SCN1A role in presynaptic terminals where voltage-gated sodium channels regulate membrane potential and neurotransmitter release. However, IEA evidence alone is relatively weak, and this represents a more specific function than the core sodium channel activity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001518|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel complex|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate of earlier IBA annotation for the same GO term. This IEA annotation provides additional computational support for SCN1A being part of voltage-gated sodium channel complex, consistent with its role as the pore-forming α-subunit.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005216|
|label[?](https://w3id.org/ai4curation/gene_review/label)|monoatomic ion channel activity|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a broader molecular function term that encompasses sodium channel activity. While accurate, it is less informative than the specific voltage-gated sodium channel activity annotation. This represents a valid but general annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005248|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel activity|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Another annotation for the core molecular function of SCN1A with IEA evidence. This duplicates the IBA annotation but provides additional computational support for the essential voltage-gated sodium channel activity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005261|
|label[?](https://w3id.org/ai4curation/gene_review/label)|monoatomic cation channel activity|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a broader molecular function term that encompasses sodium channel activity. Sodium channels are indeed monoatomic cation channels, but this annotation is less specific than voltage-gated sodium channel activity. Still accurate but general.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Correct cellular component annotation. SCN1A is a multi-pass membrane protein embedded in the plasma membrane where it functions as a voltage-gated sodium channel. UniProt confirms this subcellular localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006811|
|label[?](https://w3id.org/ai4curation/gene_review/label)|monoatomic ion transport|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a broad biological process term that encompasses sodium ion transport. While accurate (sodium is a monoatomic ion), it is less informative than the more specific sodium ion transmembrane transport annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006814|
|label[?](https://w3id.org/ai4curation/gene_review/label)|sodium ion transport|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This accurately describes a core biological process mediated by SCN1A. The protein enables sodium ion transport, which is fundamental for neuronal excitability and action potential propagation. This is a key function supported by extensive evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016020|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a very general cellular component term. While SCN1A is indeed a membrane protein, this annotation provides minimal informative value compared to the more specific plasma membrane annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0055085|
|label[?](https://w3id.org/ai4curation/gene_review/label)|transmembrane transport|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a broad biological process term that encompasses the specific sodium ion transmembrane transport function. While accurate, it provides less informative value than more specific transport annotations.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0001508|
|label[?](https://w3id.org/ai4curation/gene_review/label)|action potential|
|IEA|[GO_REF:0000117](GO_REF:0000117)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation correctly identifies SCN1A role in action potentials. As the pore-forming subunit of Nav1.1, SCN1A directly mediates the depolarizing phase of action potentials by allowing Na+ influx. This is a core biological process supported by UniProt functional data.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005272|
|label[?](https://w3id.org/ai4curation/gene_review/label)|sodium channel activity|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation captures the core sodium channel activity of SCN1A. While slightly broader than voltage-gated sodium channel activity, it accurately describes the fundamental molecular function of the protein in enabling sodium ion passage through the membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0034220|
|label[?](https://w3id.org/ai4curation/gene_review/label)|monoatomic ion transmembrane transport|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a broader biological process term that encompasses sodium ion transmembrane transport. While accurate (sodium is a monoatomic ion), the more specific sodium ion transmembrane transport annotation is more informative for SCN1A function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0034702|
|label[?](https://w3id.org/ai4curation/gene_review/label)|monoatomic ion channel complex|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a broader cellular component term that encompasses voltage-gated sodium channel complex. SCN1A is indeed part of a monoatomic ion channel complex, but the more specific voltage-gated sodium channel complex annotation provides better functional specificity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007628|
|label[?](https://w3id.org/ai4curation/gene_review/label)|adult walking behavior|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is based on Ensembl orthology inference from mouse models. While SCN1A mutations can affect motor behavior through seizures and neurological impairment, adult walking behavior is a high-level phenotype not directly representative of the core molecular function. This is likely over-annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008340|
|label[?](https://w3id.org/ai4curation/gene_review/label)|determination of adult lifespan|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is based on Ensembl orthology inference from mouse models. While SCN1A mutations can be associated with sudden unexpected death in epilepsy (SUDEP), determination of adult lifespan is a high-level organismal phenotype not directly representative of the core sodium channel function. This represents over-annotation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0014704|
|label[?](https://w3id.org/ai4curation/gene_review/label)|intercalated disc|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a cardiac-specific cellular component annotation based on Ensembl orthology inference. While sodium channels may be present in cardiac tissue, SCN1A (Nav1.1) is predominantly neuronal, whereas SCN5A (Nav1.5) is the primary cardiac sodium channel. This annotation is likely inappropriate for SCN1A.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0019227|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neuronal action potential propagation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately captures a core biological process mediated by SCN1A. Nav1.1 channels are essential for neuronal action potential propagation, particularly in axon initial segments and nodes of Ranvier. This is well-supported by functional and localization data.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0019228|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neuronal action potential|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation correctly identifies SCN1A role in neuronal action potentials. As the pore-forming subunit of Nav1.1, SCN1A directly mediates the depolarizing phase of neuronal action potentials. This is a core neuronal function well-supported by evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0021675|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nerve development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is based on Ensembl orthology inference. While SCN1A may play roles during nervous system development, this is a developmental process annotation that is peripheral to the core adult function as a voltage-gated sodium channel. Evidence for direct developmental roles is limited.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030018|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Z disc|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a muscle-specific cellular component annotation based on Ensembl orthology inference. Z discs are structures in striated muscle. SCN1A (Nav1.1) is predominantly neuronal and this localization is not supported by evidence for the human protein.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030315|
|label[?](https://w3id.org/ai4curation/gene_review/label)|T-tubule|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a muscle-specific cellular component annotation based on Ensembl orthology inference. T-tubules are structures in muscle cells for excitation-contraction coupling. SCN1A (Nav1.1) is predominantly neuronal and this localization is inappropriate.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030424|
|label[?](https://w3id.org/ai4curation/gene_review/label)|axon|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately reflects SCN1A localization. Nav1.1 channels are enriched in axons, particularly at axon initial segments and nodes of Ranvier where they are crucial for action potential initiation and propagation. This cellular component annotation is well-supported.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0033268|
|label[?](https://w3id.org/ai4curation/gene_review/label)|node of Ranvier|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately reflects a key subcellular localization of SCN1A. Nav1.1 channels are highly concentrated at nodes of Ranvier where they mediate saltatory conduction and action potential propagation along myelinated axons. This is a well-established localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0034706|
|label[?](https://w3id.org/ai4curation/gene_review/label)|sodium channel complex|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation correctly identifies SCN1A as part of a sodium channel complex. As the pore-forming α-subunit, SCN1A associates with auxiliary β-subunits to form functional Nav1.1 sodium channel complexes. This cellular component annotation is accurate.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042391|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of membrane potential|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately captures a key biological process mediated by SCN1A. Voltage-gated sodium channels are fundamental regulators of membrane potential in excitable cells, controlling the depolarization phase of action potentials and overall neuronal excitability.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043025|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neuronal cell body|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation reflects SCN1A localization in neuronal cell bodies. While Nav1.1 channels are present in cell bodies, they are more highly concentrated and functionally important at axon initial segments and nodes of Ranvier. This represents a valid but less specific localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043194|
|label[?](https://w3id.org/ai4curation/gene_review/label)|axon initial segment|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately reflects a critical subcellular localization of SCN1A. Nav1.1 channels are highly concentrated at axon initial segments where they are essential for action potential initiation. This is one of the most functionally important localizations for SCN1A.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050884|
|label[?](https://w3id.org/ai4curation/gene_review/label)|neuromuscular process controlling posture|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is based on Ensembl orthology inference from mouse models. While SCN1A mutations can indirectly affect posture through seizures and neurological impairment, this high-level behavioral annotation is not directly representative of the core sodium channel function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050966|
|label[?](https://w3id.org/ai4curation/gene_review/label)|detection of mechanical stimulus involved in sensory perception of pain|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation has some support from UniProt which states that Nav1.1 contributes to sensory perception of mechanically-induced pain through controlling excitability of somatosensory neurons. However, this is a specialized function compared to the core neuronal excitability role.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0051649|
|label[?](https://w3id.org/ai4curation/gene_review/label)|establishment of localization in cell|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This is a very broad biological process annotation based on Ensembl orthology inference. While SCN1A protein must be localized to specific membrane compartments, this general localization process annotation provides minimal informative value about the core sodium channel function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0086010|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane depolarization during action potential|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately describes the core biological process mediated by SCN1A. Nav1.1 channels directly mediate membrane depolarization during the rising phase of action potentials through selective sodium influx. This is a fundamental function supported by extensive evidence.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005654|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nucleoplasm|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is problematic. SCN1A encodes a multi-pass membrane protein that functions in the plasma membrane. Nuclear localization is inconsistent with its known function as a voltage-gated sodium channel. This IDA evidence may be from immunofluorescence artifacts or cross-reactivity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is correct and supported by direct experimental evidence (IDA). SCN1A is a multi-pass membrane protein embedded in the plasma membrane where it functions as a voltage-gated sodium channel. This is the primary functional localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016604|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nuclear body|
|IDA|[GO_REF:0000052](GO_REF:0000052)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation is inconsistent with SCN1A function. SCN1A encodes a voltage-gated sodium channel that functions in the plasma membrane, not in nuclear bodies. This IDA evidence may be from immunofluorescence artifacts or antibody cross-reactivity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0099508|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated monoatomic ion channel activity involved in regulation of presynaptic membrane potential|
|NAS|[PMID:22150645](PMID:22150645)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation describes a specific molecular function of SCN1A at presynaptic terminals. The referenced study (PMID:22150645) demonstrates SCN1A role in regulating presynaptic membrane potential. While this is a specialized function, it is well-supported by experimental evidence and represents an important aspect of Nav1.1 function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0099508|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated monoatomic ion channel activity involved in regulation of presynaptic membrane potential|
|IDA|[PMID:22150645](PMID:22150645)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation for the same GO term with stronger IDA evidence from PMID:22150645. This direct experimental evidence supports SCN1A role in regulating presynaptic membrane potential, which is important for neurotransmitter release and synaptic transmission.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0099508|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated monoatomic ion channel activity involved in regulation of presynaptic membrane potential|
|IMP|[PMID:22150645](PMID:22150645)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Third annotation for the same GO term with IMP evidence from PMID:22150645. This mutant phenotype evidence further supports the role of SCN1A in presynaptic membrane potential regulation. The study examined SCN1A truncation mutants and their effects on channel function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005248|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel activity|
|IMP|[PMID:14672992](PMID:14672992)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation for the core molecular function has strong IMP (mutant phenotype) evidence from PMID:14672992. This represents high-quality experimental support for the voltage-gated sodium channel activity, complementing the IBA evidence for the same function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|IDA|[PMID:14672992](PMID:14672992)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation has strong IDA evidence from PMID:14672992 for plasma membrane localization. This complements other annotations for the same cellular component and provides experimental support for the primary functional localization of SCN1A.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0086010|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane depolarization during action potential|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation accurately describes the core biological process with ISS (sequence similarity) evidence. SCN1A mediates membrane depolarization during action potentials through selective sodium influx. This is a fundamental function well-conserved across species.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0086010|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane depolarization during action potential|
|IMP|[PMID:14672992](PMID:14672992)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation for the same core biological process with strong IMP evidence from PMID:14672992. This mutant phenotype evidence provides experimental support for SCN1A role in membrane depolarization during action potentials.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0086002|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cardiac muscle cell action potential involved in contraction|
|IMP|[PMID:27207958](PMID:27207958)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation has IMP evidence but is questionable for SCN1A. The referenced study (PMID:27207958) is about TRPM4 variants in childhood atrioventricular block, not SCN1A. SCN1A (Nav1.1) is predominantly neuronal, while SCN5A (Nav1.5) is the main cardiac sodium channel. This may be an annotation error.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005248|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel activity|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Another annotation for the core molecular function with ISS evidence based on sequence similarity. This provides additional computational support for the voltage-gated sodium channel activity, which is highly conserved across orthologs.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050966|
|label[?](https://w3id.org/ai4curation/gene_review/label)|detection of mechanical stimulus involved in sensory perception of pain|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation for the same specialized function with ISS evidence. This provides additional computational support for SCN1A role in pain perception through controlling somatosensory neuron excitability, as mentioned in UniProt.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030018|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Z disc|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate annotation for muscle-specific Z disc localization with ISS evidence. Like the previous annotation, this is inappropriate for SCN1A (Nav1.1) which is predominantly neuronal. Z discs are structures in striated muscle not relevant to Nav1.1 function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005248|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel activity|
|NAS|[PMID:10742094](PMID:10742094)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Another annotation for the core molecular function with NAS evidence from PMID:10742094, a study about SCN1A mutations in GEFS+2. This provides additional literature support for the voltage-gated sodium channel activity of SCN1A.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006814|
|label[?](https://w3id.org/ai4curation/gene_review/label)|sodium ion transport|
|NAS|[PMID:10742094](PMID:10742094)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|This annotation for the core biological process has NAS evidence from PMID:10742094. This provides additional literature support for SCN1A role in sodium ion transport, which is fundamental to its voltage-gated sodium channel function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|

## core_functions


|description|molecular_function|directly_involved_in|locations|
|---|---|---|---|
|Mediating voltage-gated sodium channel activity as the pore-forming Nav1.1 α-subunit that enables neuronal action potential generation and propagation|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005248|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel activity|
|||
|Regulating presynaptic membrane potential through specialized voltage-gated ion channel activity that controls neurotransmitter release|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0099508|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated monoatomic ion channel activity involved in regulation of presynaptic membrane potential|
|||
|Controlling GABAergic interneuron excitability to maintain excitation-inhibition balance in neural circuits|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005248|
|label[?](https://w3id.org/ai4curation/gene_review/label)|voltage-gated sodium channel activity|
|||

## suggested_questions


|question|
|---|
|How do different SCN1A mutations affect sodium channel gating properties and contribute to distinct epilepsy phenotypes?|
|What determines the brain region-specific effects of SCN1A dysfunction and why is inhibitory neuron function particularly affected?|
|How does SCN1A haploinsufficiency lead to the temperature-sensitive seizures characteristic of Dravet syndrome?|
|What are the developmental changes in SCN1A expression and function that contribute to age-dependent seizure patterns?|

## suggested_experiments


|description|
|---|
|Patch-clamp electrophysiology of SCN1A variants in different neuronal subtypes to correlate biophysical properties with clinical phenotypes|
|Organoid models of human brain development using patient-derived iPSCs to study SCN1A function in cortical circuit formation|
|Two-photon calcium imaging in brain slices to study how SCN1A mutations affect inhibitory circuit function and excitability|
|Cryo-EM structural analysis of SCN1A in different conformational states to understand mutation effects on channel structure|
