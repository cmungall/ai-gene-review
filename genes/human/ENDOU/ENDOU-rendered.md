
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|P21128|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|ENDOU|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|PP11|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|Placental protein 11|
|aliases[?](https://w3id.org/ai4curation/gene_review/aliases)|PRSS26|
|description[?](https://w3id.org/ai4curation/gene_review/description)|Uridylate-specific endoribonuclease that cleaves single-stranded RNA at UU and GU motifs in a Mn2+-dependent manner. Originally misidentified as a serine protease (placental protein 11/PP11), ENDOU is now established as a poly(U)-specific ribonuclease with emerging roles in lipid homeostasis and immune regulation. The enzyme produces 2',3'-cyclic phosphate termini through a His-His-Lys catalytic triad mechanism. Secreted protein expressed predominantly in placenta (syncytiotrophoblast) and various tumor tissues, with newly discovered functions in downregulating lipolytic gene expression to maintain lipid storage during aging.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|
|---|---|
|GO_REF:0000002|Gene Ontology annotation through association of InterPro records with GO terms.|
|GO_REF:0000033|Annotation inferences using phylogenetic trees|
|GO_REF:0000043|Gene Ontology annotation based on UniProtKB/Swiss-Prot keyword mapping|
|GO_REF:0000108|Automatic assignment of GO terms using logical inference, based on on inter-ontology links.|
|GO_REF:0000117|Electronic Gene Ontology annotations created by ARBA machine learning models|
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.|
|PMID:1710108|Homology of placental protein 11 and pea seed albumin 2 with vitronectin.|
|PMID:18936097|The tumor marker human placental protein 11 is an endoribonuclease.|
|PMID:2350438|Cloning and expression of a cDNA encoding human placental protein 11, a putative serine protease with diagnostic significance as a tumor marker.|
|PMID:37803019|The endoribonuclease Arlr is required to maintain lipid homeostasis by downregulating lipolytic genes during aging.|
|PMID:6755403|Immunohistochemical detection of pregnancy-specific protein (SP1) and placenta-specific tissue proteins (PP5, PP10, PP11 and PP12) in ovarian adenocarcinomas.|

## existing_annotations


|term|evidence_type|original_reference_id|review|
|---|---|---|---|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004521|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA endonuclease activity|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Core enzymatic function strongly supported by IBA and confirmed experimentally. This phylogenetically-inferred annotation aligns perfectly with the demonstrated molecular function of ENDOU.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:18936097 definitively established ENDOU as an RNA endonuclease with Mn2+-dependent activity that cleaves single-stranded RNA at UU and GU dinucleotide motifs, producing 2',3'-cyclic phosphate termini. This core function is strongly supported by phylogenetic analysis.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007165|
|label[?](https://w3id.org/ai4curation/gene_review/label)|signal transduction|
|IEA|[GO_REF:0000108](GO_REF:0000108)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - This vague IEA annotation likely derives from obsolete protein family classifications. No direct evidence supports ENDOU involvement in signal transduction pathways. The protein functions as an RNA endonuclease that affects gene expression post-transcriptionally, not through signal transduction.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|No evidence from PMID:18936097, PMID:37803019, or other references supports signal transduction activity. ENDOU functions as an RNA-degrading enzyme, not a signaling molecule.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016192|
|label[?](https://w3id.org/ai4curation/gene_review/label)|vesicle-mediated transport|
|IEA|[GO_REF:0000108](GO_REF:0000108)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|MARK_AS_OVER_ANNOTATED - While ENDOU is a secreted protein that passes through the secretory pathway, this annotation implies active participation in vesicle transport mechanisms. The protein is a cargo, not a regulator of vesicle transport. This represents an over-interpretation of the secretory nature of the protein.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MARK_AS_OVER_ANNOTATED|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004521|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA endonuclease activity|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - This IEA annotation based on InterPro domain mapping correctly identifies the core molecular function, subsequently confirmed by experimental evidence (PMID:18936097).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:18936097 provided definitive experimental proof of RNA endonuclease activity through biochemical assays showing Mn2+-dependent cleavage of RNA substrates at specific dinucleotide sequences.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005044|
|label[?](https://w3id.org/ai4curation/gene_review/label)|scavenger receptor activity|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - No evidence supports scavenger receptor activity. This erroneous annotation likely stems from early misidentification of ENDOU as a serine protease or confusion with other placental proteins. ENDOU is an RNA endonuclease, not a receptor.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:18936097 definitively disproved receptor activity and established ENDOU as an RNA endonuclease. No binding or signaling assays support scavenger receptor function.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006955|
|label[?](https://w3id.org/ai4curation/gene_review/label)|immune response|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|MODIFY - While too broad, there is evidence for ENDOU involvement in B cell tolerance and activation-induced cell death (mouse studies). Should be refined to more specific immune process terms like "B cell tolerance induction" (GO:0002514) or "activation-induced cell death of T cells" (GO:0006924).|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Based on mouse studies showing ENDOU/EndoU role in immune cell regulation, but the broad immune response term should be replaced with more specific processes.|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030247|
|label[?](https://w3id.org/ai4curation/gene_review/label)|polysaccharide binding|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - No evidence for polysaccharide binding. ENDOU binds RNA (polynucleotide), not polysaccharides. This appears to be a misannotation, possibly from confusion between nucleotides and saccharides.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003723|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA binding|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct annotation supported by experimental evidence (PMID:18936097). ENDOU binds single-stranded RNA substrates, particularly poly(U) sequences, as part of its endonuclease function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:18936097 demonstrated RNA binding through electrophoretic mobility shift assays, showing specific binding to RNA substrates with a Kd of approximately 140 nM.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004518|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nuclease activity|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct but general parent term. More specific child term "RNA endonuclease activity" is preferred, but this annotation is not incorrect.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004519|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endonuclease activity|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct intermediate-level annotation. While "RNA endonuclease activity" is more specific, this parent term accurately describes the molecular function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004540|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA nuclease activity|
|IEA|[GO_REF:0000117](GO_REF:0000117)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct annotation. ENDOU is indeed an RNA nuclease, specifically an endoribonuclease. This is a valid parent term of the more specific "RNA endonuclease activity".|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005576|
|label[?](https://w3id.org/ai4curation/gene_review/label)|extracellular region|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct localization. ENDOU is a secreted protein with a signal peptide, functioning in the extracellular space as confirmed by multiple studies.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016787|
|label[?](https://w3id.org/ai4curation/gene_review/label)|hydrolase activity|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct high-level annotation. Nucleases are hydrolases that cleave phosphodiester bonds. While very general, this annotation is accurate.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016829|
|label[?](https://w3id.org/ai4curation/gene_review/label)|lyase activity|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - While ENDOU primarily functions as a hydrolase, it produces 2',3'-cyclic phosphate ends through a mechanism involving elimination (lyase-like activity). This dual classification reflects the complex chemistry of the reaction.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0046872|
|label[?](https://w3id.org/ai4curation/gene_review/label)|metal ion binding|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|MODIFY - Too general. Should be refined to "manganese ion binding" (GO:0030145) as ENDOU specifically requires Mn2+ for catalytic activity, as demonstrated experimentally.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:18936097 specifically demonstrated manganese-dependent catalytic activity, with Mn2+ being essential for RNA endonuclease function.|

#### proposed_replacement_terms

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004521|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA endonuclease activity|
|IGI|[PMID:37803019](PMID:37803019)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Strong genetic evidence from Drosophila Arlr (ENDOU ortholog) studies confirming RNA endonuclease activity. This study demonstrated that the EndoU domain is essential for function in lipid metabolism regulation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:37803019 provided strong genetic evidence using Drosophila Arlr mutants, demonstrating that the EndoU domain is required for normal lipid homeostasis through RNA endonuclease activity.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016441|
|label[?](https://w3id.org/ai4curation/gene_review/label)|post-transcriptional gene silencing|
|IGI|[PMID:37803019](PMID:37803019)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - This study demonstrated that Arlr/ENDOU degrades mRNAs of lipolytic genes, representing a form of post-transcriptional gene silencing. The RIP-seq data confirmed binding to and degradation of specific mRNA targets.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:37803019 used RIP-seq and functional assays to demonstrate that ENDOU/Arlr specifically binds to and degrades mRNAs of lipolytic genes, representing post-transcriptional gene silencing.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050995|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of lipid catabolic process|
|IGI|[PMID:37803019](PMID:37803019)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Excellent specific annotation based on recent high-quality research. The study clearly demonstrated that ENDOU/Arlr downregulates lipolytic genes (Lsd-1, regucalcin, yip2, CG5162) to maintain lipid homeostasis during aging. This represents a newly discovered core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:37803019 provided comprehensive evidence showing ENDOU/Arlr negatively regulates lipid catabolism by degrading specific lipolytic gene mRNAs, with loss-of-function causing accelerated lipid depletion during aging.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003723|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA binding|
|IDA|[PMID:18936097](PMID:18936097)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Direct experimental evidence showing ENDOU binds RNA substrates with Kd of 140 nM. Electrophoretic mobility shift assays clearly demonstrated RNA-protein complex formation.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004521|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA endonuclease activity|
|IDA|[PMID:18936097](PMID:18936097)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Definitive experimental proof that ENDOU is an RNA endonuclease. This pivotal study demonstrated Mn2+-dependent cleavage at UU/GU sites producing 2',3'-cyclic phosphate ends, definitively establishing the enzymatic function and correcting the prior misannotation as a serine protease.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008236|
|label[?](https://w3id.org/ai4curation/gene_review/label)|serine-type peptidase activity|
|IDA|[PMID:18936097](PMID:18936097)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - This annotation is INCORRECT and contradicts the findings of the cited paper. PMID:18936097 explicitly DISPROVED serine protease activity, showing no detectable protease activity in chromogenic assays. The paper established that ENDOU is an RNA endonuclease, NOT a protease. This represents a critical misannotation that must be removed.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|PMID:18936097 definitively disproved serine protease activity through negative chromogenic assays and established ENDOU as an RNA endonuclease. This annotation directly contradicts experimental evidence.|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030145|
|label[?](https://w3id.org/ai4curation/gene_review/label)|manganese ion binding|
|TAS|[PMID:18936097](PMID:18936097)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct and specific annotation. The study demonstrated that ENDOU activity is Mn2+-dependent, with manganese ions essential for the catalytic mechanism.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005886|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasma membrane|
|TAS|[PMID:1710108](PMID:1710108)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - ENDOU is a secreted protein, not a membrane-bound protein. While it may transiently associate with membranes during secretion, there is no evidence for stable plasma membrane localization. The protein lacks transmembrane domains.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008083|
|label[?](https://w3id.org/ai4curation/gene_review/label)|growth factor activity|
|NAS|[PMID:1710108](PMID:1710108)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - No evidence supports growth factor activity. This early annotation was based on sequence similarity to vitronectin, but ENDOU functions as an RNA endonuclease, not a growth factor. This represents an outdated misannotation from before the true function was discovered.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005576|
|label[?](https://w3id.org/ai4curation/gene_review/label)|extracellular region|
|TAS|[PMID:1710108](PMID:1710108)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct localization. ENDOU is indeed a secreted protein found in the extracellular space, consistent with its signal peptide and secretory pathway trafficking.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|TAS|[PMID:2350438](PMID:2350438)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - ENDOU is present in the cytoplasm of secretory cells prior to secretion. Immunohistochemistry shows cytoplasmic localization in syncytiotrophoblasts and other producing cells.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005737|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cytoplasm|
|IDA|[PMID:6755403](PMID:6755403)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Immunohistochemical evidence for cytoplasmic localization in ovarian carcinoma cells. Consistent with ENDOU being synthesized in the cytoplasm before secretion.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006508|
|label[?](https://w3id.org/ai4curation/gene_review/label)|proteolysis|
|IDA|[PMID:2350438](PMID:2350438)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - This annotation is based on the original misidentification of ENDOU as a serine protease. Later studies (PMID:18936097) definitively proved ENDOU has no protease activity. This outdated annotation must be removed.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007565|
|label[?](https://w3id.org/ai4curation/gene_review/label)|female pregnancy|
|IEP|[PMID:2350438](PMID:2350438)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|KEEP_AS_NON_CORE - ENDOU is highly expressed in placental syncytiotrophoblast during pregnancy. While this reflects tissue-specific expression rather than a core molecular function, the association with pregnancy is valid for this placental protein.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008236|
|label[?](https://w3id.org/ai4curation/gene_review/label)|serine-type peptidase activity|
|IDA|[PMID:2350438](PMID:2350438)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|REMOVE - Based on incorrect initial characterization. The 1990 paper assumed protease activity based on sequence similarity, but this was definitively disproven by PMID:18936097, which showed no protease activity and established ENDOU as an RNA endonuclease.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005615|
|label[?](https://w3id.org/ai4curation/gene_review/label)|extracellular space|
|TAS|[PMID:2350438](PMID:2350438)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ACCEPT - Correct annotation. ENDOU is a secreted protein found in the extracellular space, consistent with its signal peptide and lack of transmembrane domains.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|substrates|
|---|---|---|---|---|---|
|Cleaves single-stranded RNA at UU and GU dinucleotides through manganese-dependent endonuclease activity||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004521|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA endonuclease activity|
||||
|Degrades lipolytic gene mRNAs to suppress lipid catabolism and maintain lipid storage||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0004521|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA endonuclease activity|
||||
|Produces 2',3'-cyclic phosphate termini via His-His-Lys catalytic triad mechanism||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016829|
|label[?](https://w3id.org/ai4curation/gene_review/label)|lyase activity|
||||
|Binds RNA substrates with preference for uridine-rich sequences||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0003723|
|label[?](https://w3id.org/ai4curation/gene_review/label)|RNA binding|
||||
|Requires manganese ions as essential cofactor for catalytic activity||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030145|
|label[?](https://w3id.org/ai4curation/gene_review/label)|manganese ion binding|
||||

## suggested_questions


|question|
|---|
|How does ENDOU achieve specificity for polyuridine sequences and what determines its substrate selectivity?|
|What role does ENDOU play in placental development and how does it regulate trophoblast function?|
|How is ENDOU expression and activity regulated during pregnancy and in response to placental stress?|
|What are the downstream consequences of ENDOU-mediated RNA cleavage on gene expression and protein synthesis?|

## suggested_experiments


|description|
|---|
|RNA-seq analysis of ENDOU-deficient placental tissues to identify direct and indirect targets of ENDOU regulation|
|Biochemical characterization of ENDOU substrate specificity using synthetic RNA oligonucleotides and kinetic analyses|
|Single-cell RNA sequencing of placental cell types to understand ENDOU function in trophoblast differentiation|
|Structural biology approaches to determine the molecular basis of ENDOU endonuclease activity and substrate recognition|
