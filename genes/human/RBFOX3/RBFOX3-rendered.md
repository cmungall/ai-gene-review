
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|A6NFN3|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|RBFOX3|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|NeuN|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|Fox-1 homolog C|
|description[?](https://w3id.org/ai4curation/gene_review/description)|RBFOX3 (RNA Binding Protein Fox-1 Homolog 3) is a neuronal-specific RNA-binding protein (~314 amino acids) that regulates alternative splicing in the central nervous system. Also known as NeuN (Neuronal nuclei antigen), it serves as a definitive marker for post-mitotic neurons and is critical for neuronal differentiation, adult neurogenesis, and synaptic function. Contains a highly conserved RNA recognition  motif (RRM) that binds to UGCAUG sequences in pre-mRNA. Functions both as an activator and repressor of alternative  splicing in a position-dependent manner. Alternative splicing produces nuclear (~46 kDa) and cytoplasmic (~48 kDa)  isoforms with distinct functions. Essential for maintaining excitatory/inhibitory balance in neural circuits.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|findings|
|---|---|---|
|PMID:21747913|NeuN/Rbfox3 nuclear and cytoplasmic isoforms differentially regulate alternative splicing and nonsense-mediated decay of Rbfox2||
|PMID:23420872|Rbfox3-regulated alternative splicing of Numb promotes neuronal differentiation during development||
|file:human/RBFOX3/RBFOX3-deep-research.md|Deep Research Report: RBFOX3 comprehensive analysis||
|file:human/RBFOX3/RBFOX3-bioinformatics/rbfox3_analysis/RESULTS.md|RBFOX3 Bioinformatics Analysis: Sequence, Domain, and Family Comparative Study||
|GO_REF:0000033|Gene Ontology inferred from electronic annotation (IBA)||
|GO_REF:0000002|Gene Ontology inferred from electronic annotation based on InterPro||
|GO_REF:0000120|Gene Ontology inferred from electronic annotation (IEA)||
|GO_REF:0000043|Gene Ontology inferred from electronic annotation based on UniProtKB keywords||

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000381|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of alternative mRNA splicing, via spliceosome|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on protein family membership|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function of RBFOX3 - well-documented regulation of alternative splicing through position-dependent binding to (U)GCAUG motifs|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003676|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nucleic acid binding|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on RNA-binding domain|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Too general - RBFOX3 specifically binds RNA, not DNA. More specific RNA binding terms are available|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003723|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA binding|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on RNA recognition motif|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate and well-supported - RBFOX3 contains RRM domain and binds specific RNA sequences (PMID:21747913, file:human/RBFOX3/RBFOX3-bioinformatics/rbfox3_analysis/RESULTS.md)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007399|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nervous system development|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated and inferred annotation based on neuronal expression|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Well-supported - RBFOX3 is essential for neuronal differentiation and CNS development (PMID:23420872, file:human/RBFOX3/RBFOX3-deep-research.md)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043484|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of RNA splicing|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on splicing regulator function|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate but less specific than 'regulation of alternative mRNA splicing' - alternative splicing is RBFOX3's core function (PMID:21747913)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005634|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nucleus|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation and phylogenetic inference|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Incomplete - RBFOX3 has both nuclear and cytoplasmic isoforms with distinct functions (PMID:21747913). Nuclear localization is isoform-dependent|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation and phylogenetic inference|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Too general - RBFOX3 cytoplasmic isoforms have specific functions distinct from general cytoplasmic localization (PMID:21747913, file:human/RBFOX3/RBFOX3-deep-research.md)|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006397|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mRNA processing|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on RNA processing function|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Too general - RBFOX3 specifically regulates alternative splicing, not general mRNA processing (file:human/RBFOX3/RBFOX3-bioinformatics/rbfox3_analysis/RESULTS.md)|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008380|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA splicing|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Automated annotation based on splicing function|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Less specific than RBFOX3's actual function - it regulates alternative splicing, not general splicing (PMID:21747913)|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003729|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mRNA binding|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Phylogenetic inference based on ortholog functions|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Accurate and specific - RBFOX3 binds to specific mRNA sequences through its RRM domain (file:human/RBFOX3/RBFOX3-bioinformatics/rbfox3_analysis/RESULTS.md)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005634|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nucleus|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Phylogenetic inference of nuclear localization|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|RBFOX3 nuclear isoforms are well-documented and functionally important (PMID:21747913)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007399|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nervous system development|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Phylogenetic inference based on neuronal development role|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|RBFOX3 is essential for neuronal differentiation and CNS development (file:human/RBFOX3/RBFOX3-deep-research.md)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000381|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of alternative mRNA splicing, via spliceosome|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Phylogenetic inference of alternative splicing regulation|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Core function - RBFOX3 regulates alternative splicing through position-dependent binding (PMID:21747913, file:human/RBFOX3/RBFOX3-deep-research.md)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Phylogenetic inference of cytoplasmic localization|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|RBFOX3 cytoplasmic isoforms are well-documented with distinct functions (PMID:21747913)|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0032228|
|label[?](https://w3id.org/ai4curation/gene_review/label)|regulation of synaptic transmission, GABAergic|
|ISS|[file:human/RBFOX3/RBFOX3-deep-research.md](file:human/RBFOX3/RBFOX3-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|RBFOX3 regulates GABAergic synaptic transmission through alternative splicing of GABAergic pathway genes, maintaining excitatory/inhibitory balance.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Co-expression network analyses link RBFOX3 to GABAergic pathway genes. Loss of RBFOX3 disrupts excitatory/inhibitory balance in neural circuits.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0045773|
|label[?](https://w3id.org/ai4curation/gene_review/label)|positive regulation of axon extension|
|ISS|[file:human/RBFOX3/RBFOX3-deep-research.md](file:human/RBFOX3/RBFOX3-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|RBFOX3 promotes neurite outgrowth and axonal extension during neuronal differentiation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Knockdown of RBFOX3 impairs neurite outgrowth, while normal RBFOX3 expression promotes neuronal maturation including axon development.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042552|
|label[?](https://w3id.org/ai4curation/gene_review/label)|myelination|
|ISS|[file:human/RBFOX3/RBFOX3-deep-research.md](file:human/RBFOX3/RBFOX3-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|RBFOX3 regulates splicing of myelin-related transcripts critical for proper myelination.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|RBFOX3 alternative splicing targets include myelin basic protein and other myelination-related genes in the nervous system.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0048870|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cell motility|
|ISS|[file:human/RBFOX3/RBFOX3-deep-research.md](file:human/RBFOX3/RBFOX3-deep-research.md)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|RBFOX3 regulates neuronal migration through alternative splicing of cell motility genes including NUMB.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|NEW|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|RBFOX3-mediated splicing of NUMB affects cell fate decisions and neuronal migration during development.|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|anatomical_locations|substrates|
|---|---|---|---|---|---|---|
|Neuron-specific regulation of alternative splicing and synaptic function through sequence-specific RNA binding||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003729|
|label[?](https://w3id.org/ai4curation/gene_review/label)|mRNA binding|
|||||

## suggested_questions


|question|
|---|
|How do the nuclear and cytoplasmic isoforms of RBFOX3 coordinate to regulate both alternative splicing and mRNA localization in post-mitotic neurons?|
|What determines the specificity of RBFOX3 binding to UGCAUG motifs and how does binding position relative to exons influence splicing outcomes?|
|How does RBFOX3 expression and activity change during neuronal differentiation and maturation, and what role does this play in establishing neuronal identity?|
|What are the mechanisms by which RBFOX3 maintains the balance between excitatory and inhibitory neurotransmission at the molecular level?|

## suggested_experiments


|description|
|---|
|Single-cell RNA sequencing combined with RBFOX3 ChIP-seq to map tissue-specific and cell-type-specific alternative splicing programs in different brain regions|
|Super-resolution microscopy to visualize the subcellular localization and dynamics of RBFOX3 isoforms during neuronal activity and synaptic plasticity|
|CRISPR-Cas13 mediated knockdown of specific RBFOX3 isoforms to dissect their individual contributions to neuronal function and splicing regulation|
|Proteomics analysis using proximity labeling to identify RBFOX3 interacting partners in nuclear versus cytoplasmic compartments during neuronal development|
