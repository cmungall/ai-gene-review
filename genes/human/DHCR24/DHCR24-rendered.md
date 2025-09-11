
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|Q15392|
|gene_symbol[?](https://w3id.org/ai4curation/gene_review/gene_symbol)|DHCR24|
|description[?](https://w3id.org/ai4curation/gene_review/description)|DHCR24 encodes Delta(24)-sterol reductase (also known as seladin-1), the terminal enzyme in cholesterol biosynthesis that catalyzes the reduction of the delta-24 double bond of sterol intermediates, primarily converting desmosterol to cholesterol. The enzyme is FAD-dependent, requires NADPH as a cofactor, and is localized to the endoplasmic reticulum membrane. DHCR24 also confers neuroprotection against oxidative stress and amyloid-beta toxicity. Mutations in DHCR24 cause desmosterolosis, a rare autosomal recessive disorder of cholesterol biosynthesis characterized by multiple congenital anomalies and elevated desmosterol levels.|

## taxon


|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|NCBITaxon:9606|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Homo sapiens|

## references


|id|title|
|---|---|
|GO_REF:0000002|Gene Ontology annotation through association of InterPro records with GO terms.|
|GO_REF:0000003|Gene Ontology annotation based on Enzyme Commission mapping|
|GO_REF:0000024|Manual transfer of experimentally-verified manual GO annotation data to orthologs by curator judgment of sequence similarity.|
|GO_REF:0000033|Annotation inferences using phylogenetic trees|
|GO_REF:0000043|Gene Ontology annotation based on UniProtKB/Swiss-Prot keyword mapping|
|GO_REF:0000044|Gene Ontology annotation based on UniProtKB/Swiss-Prot Subcellular Location vocabulary mapping, accompanied by conservative changes to GO terms applied by UniProt.|
|GO_REF:0000107|Automatic transfer of experimentally verified manual GO annotation data to orthologs using Ensembl Compara.|
|GO_REF:0000120|Combined Automated Annotation using Multiple IEA Methods.|
|PMID:11007892|The human DIMINUTO/DWARF1 homolog seladin-1 confers resistance to Alzheimer's disease-associated neurodegeneration and oxidative stress.|
|PMID:11519011|Mutations in the 3beta-hydroxysterol Delta24-reductase gene cause desmosterolosis, an autosomal recessive disorder of cholesterol biosynthesis.|
|PMID:12457401|Desmosterolosis presenting with multiple congenital anomalies and profound developmental delay.|
|PMID:15577914|Regulation of cellular response to oncogenic and oxidative stress by Seladin-1.|
|PMID:19946888|Defining the membrane proteome of NK cells.|
|PMID:25637936|The terminal enzymes of cholesterol synthesis, DHCR24 and DHCR7, interact physically and functionally.|
|PMID:30021884|Histone Interaction Landscapes Visualized by Crosslinking Mass Spectrometry in Intact Cell Nuclei.|
|PMID:35271311|OpenCell: Endogenous tagging for the cartography of human cellular organization.|
|Reactome:R-HSA-191273|Cholesterol biosynthesis|
|Reactome:R-HSA-196417|Reduction of desmosterol to cholesterol|
|Reactome:R-HSA-6807047|Cholesterol biosynthesis via desmosterol|
|Reactome:R-HSA-6807062|Cholesterol biosynthesis via lathosterol|
|Reactome:R-HSA-6807064|DHCR24 reduces ZYMOL to ZYMSTNL|
|Reactome:R-HSA-9755937|DHCR24 reduces LAN to 24,25-dhLAN|

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
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IBA annotation for cytoplasm is questionable. While DHCR24 may have some cytoplasmic presence, all experimental evidence from multiple studies consistently shows primary localization to the endoplasmic reticulum membrane with its catalytic domain facing the cytoplasm [PMID:11007892, PMID:22010141]. The cytoplasmic annotation likely reflects the cytoplasmic orientation of the catalytic domain rather than true cytoplasmic localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This annotation is misleading as it suggests cytoplasmic localization when the protein is actually an integral ER membrane protein with a cytoplasmic-facing catalytic domain. The UniProt record and experimental evidence clearly establish ER membrane as the primary localization [PMID:11007892, PMID:22010141].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008202|
|label[?](https://w3id.org/ai4curation/gene_review/label)|steroid metabolic process|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IBA annotation for steroid metabolic process is accurate. DHCR24 is a key enzyme in cholesterol biosynthesis, which is part of steroid metabolism. The annotation is at an appropriate level of specificity given that cholesterol is the precursor for all steroid hormones.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 clearly participates in steroid metabolism through its essential role in cholesterol biosynthesis. Cholesterol is both a steroid itself and the precursor for all steroid hormones. The IBA annotation is well-supported by experimental evidence [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000246|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Delta24(24-1) sterol reductase activity|
|IBA|[GO_REF:0000033](GO_REF:0000033)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IBA annotation for Delta24(24-1) sterol reductase activity is correct and represents a core molecular function. This term specifically describes the enzymatic activity of DHCR24 in reducing the delta-24 double bond in sterols, which is directly validated by experimental evidence [PMID:11519011].|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This is the core molecular function of DHCR24, experimentally validated through multiple studies. The enzyme specifically catalyzes the reduction of the delta-24 double bond in sterol intermediates, including desmosterol, lanosterol, and zymosterol [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050614|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Delta24-sterol reductase activity|
|IEA|[GO_REF:0000003](GO_REF:0000003)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for Delta24-sterol reductase activity based on EC number mapping. This is essentially the same function as GO:0000246 but with slightly different terminology. Both terms correctly describe the core enzymatic function of DHCR24.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This annotation accurately describes the enzymatic activity of DHCR24. The term is based on EC:1.3.1.72 mapping which is experimentally validated [PMID:11519011]. This is a duplicate of GO:0000246 with different terminology but both are correct.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000139|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Golgi membrane|
|IEA|[GO_REF:0000044](GO_REF:0000044)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for Golgi membrane based on UniProt subcellular location vocabulary. While some DHCR24 may be present in Golgi, the primary and functionally relevant localization is the ER membrane where cholesterol synthesis occurs [PMID:11007892].|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Minor Golgi localization may occur but is not the primary or functionally relevant site. The ER membrane is the established primary location for DHCR24 function in cholesterol biosynthesis. UniProt notes both localizations but emphasizes ER [PMID:11007892].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005789|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endoplasmic reticulum membrane|
|IEA|[GO_REF:0000044](GO_REF:0000044)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for endoplasmic reticulum membrane is accurate and represents the primary cellular localization. This is the established site where DHCR24 performs its enzymatic function in cholesterol biosynthesis, confirmed by multiple experimental studies [PMID:11007892, PMID:22010141].|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|ER membrane is the correct and primary localization for DHCR24. This is where the enzyme performs its function in cholesterol biosynthesis. Multiple experimental studies confirm this localization [PMID:11007892, PMID:22010141].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006629|
|label[?](https://w3id.org/ai4curation/gene_review/label)|lipid metabolic process|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for lipid metabolic process based on UniProt keyword mapping. This is accurate but very general. More specific terms like cholesterol biosynthetic process (GO:0006695) better describe the function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While accurate, this is a very broad parent term. DHCR24 does participate in lipid metabolism through its role in cholesterol biosynthesis. However, more specific child terms provide better functional description.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006694|
|label[?](https://w3id.org/ai4curation/gene_review/label)|steroid biosynthetic process|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for steroid biosynthetic process is accurate. DHCR24 is essential for cholesterol biosynthesis, and cholesterol is both a steroid and the precursor for all steroid hormones.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 directly participates in steroid biosynthesis through its essential role in producing cholesterol, which is a steroid molecule and the precursor for all steroid hormones [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006695|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol biosynthetic process|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for cholesterol biosynthetic process is accurate and represents a core function. DHCR24 is the terminal enzyme in cholesterol biosynthesis, converting desmosterol to cholesterol [PMID:11519011].|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This is a core biological process for DHCR24. The enzyme catalyzes the final step in cholesterol biosynthesis, converting desmosterol to cholesterol. This is experimentally validated [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008202|
|label[?](https://w3id.org/ai4curation/gene_review/label)|steroid metabolic process|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate IEA annotation for steroid metabolic process (also annotated with IBA evidence). The annotation is accurate as DHCR24 participates in cholesterol/steroid metabolism.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Duplicate but accurate annotation. DHCR24 participates in steroid metabolism through cholesterol biosynthesis. The IBA version of this annotation was already accepted.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008203|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol metabolic process|
|IEA|[GO_REF:0000120](GO_REF:0000120)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for cholesterol metabolic process is accurate. DHCR24 is directly involved in cholesterol metabolism as the terminal enzyme in cholesterol biosynthesis.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 is a key enzyme in cholesterol metabolism, specifically in the biosynthetic pathway. The annotation accurately captures this core function [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016126|
|label[?](https://w3id.org/ai4curation/gene_review/label)|sterol biosynthetic process|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for sterol biosynthetic process is accurate. DHCR24 is essential for sterol biosynthesis, specifically in the final step of converting sterol intermediates to cholesterol.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 is a key enzyme in sterol biosynthesis, catalyzing the reduction of the delta-24 double bond in various sterol intermediates including desmosterol, lanosterol, and zymosterol [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016491|
|label[?](https://w3id.org/ai4curation/gene_review/label)|oxidoreductase activity|
|IEA|[GO_REF:0000043](GO_REF:0000043)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for oxidoreductase activity is accurate but very general. DHCR24 is indeed an oxidoreductase, specifically a FAD-dependent oxidoreductase that uses NADPH. More specific terms like GO:0016628 provide better functional description.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 is a FAD-dependent oxidoreductase that catalyzes redox reactions using NADPH as electron donor. While accurate, more specific child terms better describe the function [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050660|
|label[?](https://w3id.org/ai4curation/gene_review/label)|flavin adenine dinucleotide binding|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for FAD binding based on InterPro domain mapping. This is accurate as DHCR24 contains a FAD-binding domain and FAD enhances its enzymatic activity [PMID:11519011].|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 is a FAD-dependent oxidoreductase with a conserved FAD-binding domain. Addition of FAD increases enzymatic activity twofold, suggesting noncovalent FAD binding [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0071949|
|label[?](https://w3id.org/ai4curation/gene_review/label)|FAD binding|
|IEA|[GO_REF:0000002](GO_REF:0000002)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Duplicate IEA annotation for FAD binding (same as GO:0050660). The annotation is accurate as DHCR24 is a FAD-dependent enzyme.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Duplicate of GO:0050660 but accurate. DHCR24 requires FAD as a cofactor for its oxidoreductase activity [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:30021884](PMID:30021884)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Generic protein binding annotation from a large-scale histone interaction study. This vague annotation provides no functional insight about DHCR24. The specific interaction with DHCR7 (PMID:25637936) is more informative.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Protein binding without specificity is uninformative. The curation guidelines explicitly state to avoid this vague term. The functionally relevant interaction with DHCR7 is captured in another annotation with PMID:25637936.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:35271311](PMID:35271311)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Generic protein binding annotation from OpenCell high-throughput study. This vague annotation provides no functional insight. Specific interactions like DHCR7 are more informative.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Non-specific protein binding annotation from high-throughput study lacks functional context. The curation guidelines explicitly state to avoid this vague term. The functionally relevant DHCR7 interaction is captured elsewhere.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006695|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol biosynthetic process|
|TAS|[Reactome:R-HSA-191273](Reactome:R-HSA-191273)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TAS annotation from Reactome for cholesterol biosynthetic process. This is a core function of DHCR24 as the terminal enzyme in cholesterol biosynthesis.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 is the terminal enzyme in cholesterol biosynthesis. This Reactome annotation accurately captures this core biological process [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0033489|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol biosynthetic process via desmosterol|
|TAS|[Reactome:R-HSA-6807047](Reactome:R-HSA-6807047)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TAS annotation from Reactome for cholesterol biosynthesis via desmosterol. This is the primary pathway where DHCR24 functions, converting desmosterol to cholesterol in the Bloch pathway.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This specifically describes the Bloch pathway where DHCR24 catalyzes the conversion of desmosterol to cholesterol. This is the primary route and core function of DHCR24 [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0033490|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol biosynthetic process via lathosterol|
|TAS|[Reactome:R-HSA-6807062](Reactome:R-HSA-6807062)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TAS annotation from Reactome for cholesterol biosynthesis via lathosterol (Kandutsch-Russell pathway). DHCR24 can act on various sterol intermediates including those in this alternate pathway.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 functions in both the Bloch and Kandutsch-Russell pathways of cholesterol biosynthesis. It can reduce the delta-24 double bond in various sterol intermediates [PMID:25637936].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0007265|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Ras protein signal transduction|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for Ras signaling transferred from orthologs. While DHCR24/seladin-1 mediates response to Ras-induced senescence [PMID:15577914], this is a stress response function rather than core Ras signaling.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24/seladin-1 responds to oncogenic Ras stress but is not a core component of Ras signaling. It mediates Ras-induced senescence through p53 interactions [PMID:15577914]. This is a secondary, non-core function.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008104|
|label[?](https://w3id.org/ai4curation/gene_review/label)|intracellular protein localization|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for intracellular protein localization from ortholog transfer. This is vague and lacks supporting evidence for DHCR24 having a role in localizing other proteins.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|No evidence that DHCR24 functions in protein localization. This appears to be an over-annotation from ortholog transfer without validation. The primary function is enzymatic in cholesterol biosynthesis.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008285|
|label[?](https://w3id.org/ai4curation/gene_review/label)|negative regulation of cell population proliferation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for negative regulation of proliferation. DHCR24/seladin-1 does mediate Ras-induced senescence and suppresses transformation [PMID:15577914], supporting an anti-proliferative role.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24/seladin-1 mediates oncogenic stress-induced senescence and suppresses cellular transformation [PMID:15577914]. This is a validated but non-core function distinct from its primary enzymatic role.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009725|
|label[?](https://w3id.org/ai4curation/gene_review/label)|response to hormone|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for response to hormone from ortholog transfer. This is too vague without specifying which hormones. Cholesterol biosynthesis is regulated by hormones but this annotation lacks specificity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Overly broad annotation without supporting evidence for specific hormone responses by DHCR24. While cholesterol biosynthesis is hormonally regulated, this vague annotation provides no functional insight.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009888|
|label[?](https://w3id.org/ai4curation/gene_review/label)|tissue development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for tissue development from ortholog transfer. While DHCR24 mutations cause developmental anomalies in desmosterolosis, this is an indirect effect of cholesterol deficiency.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 deficiency causes developmental defects in desmosterolosis patients [PMID:11519011], but this is secondary to cholesterol deficiency rather than a direct developmental role. Non-core but valid consequence.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016125|
|label[?](https://w3id.org/ai4curation/gene_review/label)|sterol metabolic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for sterol metabolic process is accurate. DHCR24 is directly involved in sterol metabolism through its enzymatic activity on various sterol intermediates.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 directly participates in sterol metabolism by catalyzing the reduction of the delta-24 double bond in multiple sterol intermediates [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0030539|
|label[?](https://w3id.org/ai4curation/gene_review/label)|male genitalia development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for male genitalia development from ortholog transfer. While desmosterolosis can include urogenital anomalies, this is too specific and indirect.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Over-specific developmental annotation. While desmosterolosis includes congenital anomalies, singling out male genitalia development is not supported as a specific function of DHCR24. This is an indirect effect of cholesterol deficiency.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0031639|
|label[?](https://w3id.org/ai4curation/gene_review/label)|plasminogen activation|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for plasminogen activation from ortholog transfer. No evidence supports DHCR24 involvement in plasminogen activation. This appears to be an erroneous transfer.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|No evidence linking DHCR24 to plasminogen activation. This appears to be an incorrect ortholog transfer. DHCR24 functions in cholesterol biosynthesis, not coagulation/fibrinolysis pathways.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042987|
|label[?](https://w3id.org/ai4curation/gene_review/label)|amyloid precursor protein catabolic process|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for APP catabolism from ortholog transfer. While DHCR24/seladin-1 protects against amyloid-beta toxicity [PMID:11007892], it does not directly catabolize APP.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24/seladin-1 protects against amyloid-beta toxicity but does not catabolize APP. The protective effect is through reducing oxidative stress and caspase activation, not APP processing [PMID:11007892].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043588|
|label[?](https://w3id.org/ai4curation/gene_review/label)|skin development|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for skin development from ortholog transfer. While cholesterol is important for skin barrier function, this is an indirect effect of DHCR24 enzymatic activity.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Cholesterol biosynthesis is important for skin development and barrier function. DHCR24 deficiency can affect skin, but this is secondary to its enzymatic role. Non-core but valid indirect function.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0061024|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane organization|
|IEA|[GO_REF:0000107](GO_REF:0000107)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IEA annotation for membrane organization from ortholog transfer. While cholesterol is crucial for membrane structure, DHCR24 does not directly organize membranes.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 produces cholesterol which affects membrane properties, but the enzyme itself does not organize membranes. This is an over-interpretation of its indirect effects through cholesterol production.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050614|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Delta24-sterol reductase activity|
|EXP|[PMID:11519011](PMID:11519011)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Experimental evidence for Delta24-sterol reductase activity from the key paper identifying DHCR24. This is the core molecular function, directly demonstrated through heterologous expression and enzyme activity measurements.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Direct experimental validation of DHCR24 enzymatic activity. The study used heterologous expression in yeast followed by enzyme activity measurements to confirm Delta24-sterol reductase activity [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050614|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Delta24-sterol reductase activity|
|TAS|[Reactome:R-HSA-9755937](Reactome:R-HSA-9755937)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TAS annotation from Reactome for Delta24-sterol reductase activity, specifically for lanosterol reduction. This is a validated core molecular function of DHCR24.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 catalyzes the reduction of the delta-24 double bond in multiple sterols including lanosterol. This Reactome annotation correctly captures this enzymatic activity [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0000246|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Delta24(24-1) sterol reductase activity|
|IMP|[PMID:11519011](PMID:11519011)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IMP evidence for Delta24(24-1) sterol reductase activity from mutant phenotype analysis. Mutations in DHCR24 cause desmosterolosis with accumulation of desmosterol, proving this enzymatic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Mutations in DHCR24 cause desmosterolosis with elevated desmosterol levels, directly demonstrating the enzyme's Delta24-sterol reductase activity through mutant phenotype [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005789|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endoplasmic reticulum membrane|
|NAS|[PMID:11007892](PMID:11007892)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|NAS evidence for ER membrane localization from the seladin-1 paper. Subcellular fractionation confirmed predominant ER localization. This is the primary and functionally relevant cellular location.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Experimental subcellular fractionation demonstrated that seladin-1/DHCR24 is predominantly localized to the ER membrane, where it performs its enzymatic function [PMID:11007892].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006695|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol biosynthetic process|
|IMP|[PMID:11519011](PMID:11519011)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IMP evidence for cholesterol biosynthetic process from mutant phenotype. DHCR24 mutations cause desmosterolosis with defective cholesterol biosynthesis, proving this core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Mutations in DHCR24 cause desmosterolosis, a cholesterol biosynthesis disorder, directly demonstrating the enzyme's essential role in cholesterol biosynthesis [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0033489|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol biosynthetic process via desmosterol|
|IMP|[PMID:11519011](PMID:11519011)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IMP evidence for cholesterol biosynthesis via desmosterol pathway. Patients with DHCR24 mutations accumulate desmosterol, proving this is the enzyme that converts desmosterol to cholesterol.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 deficiency causes desmosterol accumulation, directly demonstrating its role in the desmosterol-to-cholesterol conversion in the Bloch pathway [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005515|
|label[?](https://w3id.org/ai4curation/gene_review/label)|protein binding|
|IPI|[PMID:25637936](PMID:25637936)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IPI evidence for protein binding with DHCR7. While generic "protein binding" is usually uninformative, this specific DHCR7 interaction is functionally important for regulating DHCR7 activity in cholesterol synthesis.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|MODIFY|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|The interaction with DHCR7 is functionally significant, but should be annotated with a more specific term. The interaction regulates DHCR7 activity, suggesting "enzyme regulator activity" would be more informative [PMID:25637936].|

#### proposed_replacement_terms


#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005789|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endoplasmic reticulum membrane|
|TAS|[Reactome:R-HSA-196417](Reactome:R-HSA-196417)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|TAS annotation from Reactome for ER membrane localization. This is the correct primary localization where DHCR24 performs its enzymatic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|ER membrane is the established localization for DHCR24 where it catalyzes the conversion of desmosterol to cholesterol [PMID:11007892, PMID:22010141].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005789|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endoplasmic reticulum membrane|
|TAS|[Reactome:R-HSA-6807064](Reactome:R-HSA-6807064)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Another TAS annotation from Reactome for ER membrane. This is a duplicate but correct annotation for the primary cellular localization.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Duplicate but accurate ER membrane localization annotation from Reactome. Multiple evidence sources confirm this localization [PMID:11007892].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005789|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endoplasmic reticulum membrane|
|TAS|[Reactome:R-HSA-9755937](Reactome:R-HSA-9755937)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|Another TAS annotation from Reactome for ER membrane. Multiple Reactome pathways correctly place DHCR24 at the ER membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Another duplicate but accurate ER membrane annotation. The ER is where DHCR24 performs its enzymatic function in cholesterol biosynthesis [PMID:11007892].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016020|
|label[?](https://w3id.org/ai4curation/gene_review/label)|membrane|
|HDA|[PMID:19946888](PMID:19946888)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|HDA evidence for membrane localization from NK cell proteomics study. This is very general compared to the specific ER membrane localization established by other studies.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|While less specific than ER membrane annotations, this high-throughput data confirms DHCR24 is membrane-associated, consistent with its ER membrane localization [PMID:19946888].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0016628|
|label[?](https://w3id.org/ai4curation/gene_review/label)|oxidoreductase activity, acting on the CH-CH group of donors, NAD or NADP as acceptor|
|IDA|[PMID:11519011](PMID:11519011)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IDA evidence for specific oxidoreductase activity. This accurately describes DHCR24 mechanism - it reduces C-C double bonds (delta-24) using NADPH as electron donor.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Direct experimental demonstration that DHCR24 is an oxidoreductase acting on CH-CH groups (the delta-24 double bond) with NADPH as acceptor. This is more specific than general oxidoreductase activity [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0006695|
|label[?](https://w3id.org/ai4curation/gene_review/label)|cholesterol biosynthetic process|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ISS annotation for cholesterol biosynthetic process based on sequence similarity. This is a duplicate annotation but correctly identifies a core function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Sequence similarity-based annotation that correctly identifies DHCR24 role in cholesterol biosynthesis. This is validated by experimental evidence [PMID:11519011].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0019899|
|label[?](https://w3id.org/ai4curation/gene_review/label)|enzyme binding|
|IPI|[PMID:15577914](PMID:15577914)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IPI evidence for enzyme binding from the seladin-1 stress response paper. DHCR24/seladin-1 binds to Mdm2 (an E3 ubiquitin ligase) to regulate p53. This is a non-core stress response function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24/seladin-1 binds the E3 ubiquitin ligase Mdm2 during stress response, affecting p53 regulation. This is a validated but non-core function distinct from its primary enzymatic role [PMID:15577914].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0043588|
|label[?](https://w3id.org/ai4curation/gene_review/label)|skin development|
|ISS|[GO_REF:0000024](GO_REF:0000024)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|ISS annotation for skin development based on sequence similarity. While cholesterol is important for skin, this is an indirect effect of DHCR24 enzymatic function.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Cholesterol biosynthesis impacts skin development and barrier function. DHCR24 deficiency can cause skin abnormalities, but this is secondary to its enzymatic role. Non-core indirect function.|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005634|
|label[?](https://w3id.org/ai4curation/gene_review/label)|nucleus|
|IDA|[PMID:15577914](PMID:15577914)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IDA evidence for nuclear localization from stress response study. While some DHCR24/seladin-1 may translocate to nucleus during stress, the primary functional localization is ER membrane.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24/seladin-1 can localize to nucleus during stress response to interact with p53/Mdm2, but this is not the primary localization. Main function occurs at ER membrane [PMID:15577914, PMID:11007892].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0005783|
|label[?](https://w3id.org/ai4curation/gene_review/label)|endoplasmic reticulum|
|IDA|[PMID:11007892](PMID:11007892)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IDA evidence for ER localization from the seladin-1 paper. This is the correct primary localization, though ER membrane (GO:0005789) is more specific.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|ACCEPT|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|Direct experimental evidence for ER localization through subcellular fractionation. The more specific ER membrane term is also annotated [PMID:11007892].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0009888|
|label[?](https://w3id.org/ai4curation/gene_review/label)|tissue development|
|IMP|[PMID:12457401](PMID:12457401)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IMP evidence for tissue development from a desmosterolosis case report showing developmental delay and anomalies. This is an indirect effect of cholesterol deficiency rather than a direct developmental role.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|KEEP_AS_NON_CORE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|DHCR24 deficiency causes developmental anomalies in desmosterolosis due to cholesterol deficiency. This is an indirect consequence rather than a direct developmental function [PMID:12457401].|

#### supported_by

|
|
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0042605|
|label[?](https://w3id.org/ai4curation/gene_review/label)|peptide antigen binding|
|IPI|[PMID:15577914](PMID:15577914)|
|Slot|Value|
|---|---|
|summary[?](https://w3id.org/ai4curation/gene_review/summary)|IPI evidence for peptide antigen binding from stress response study. This annotation appears incorrect - the paper shows p53 and Mdm2 binding, not antigen binding.|
|action[?](https://w3id.org/ai4curation/gene_review/action)|REMOVE|
|reason[?](https://w3id.org/ai4curation/gene_review/reason)|This appears to be an annotation error. PMID:15577914 describes DHCR24/seladin-1 binding to p53 and Mdm2 proteins during stress response, not peptide antigen binding. No MHC or antigen presentation function is described.|

#### supported_by

|

## core_functions


|description|supported_by|molecular_function|directly_involved_in|locations|substrates|
|---|---|---|---|---|---|
|Catalyzes the reduction of the delta-24 double bond in sterol intermediates using NADPH and FAD as cofactors, converting desmosterol to cholesterol as the terminal step of cholesterol biosynthesis||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0050614|
|label[?](https://w3id.org/ai4curation/gene_review/label)|Delta24-sterol reductase activity|
||||
|Positively regulates DHCR7 oxidoreductase activity through direct protein interaction, enhancing 7-dehydrocholesterol reductase function in cholesterol biosynthesis||
|Slot|Value|
|---|---|
|id[?](https://w3id.org/ai4curation/gene_review/id)|GO:0008047|
|label[?](https://w3id.org/ai4curation/gene_review/label)|enzyme activator activity|
||||

## suggested_questions


|question|
|---|
|How does DHCR24 coordinate cholesterol biosynthesis with membrane homeostasis and cellular sterol requirements?|
|What are the regulatory mechanisms that control DHCR24 expression and activity in response to sterol levels?|
|How do mutations in DHCR24 lead to Smith-Lemli-Opitz syndrome and what are the developmental consequences of altered sterol metabolism?|
|What determines the subcellular localization of DHCR24 and how does this affect its function in sterol metabolism?|

## suggested_experiments


|description|
|---|
|Lipidomics analysis to characterize the complete sterol profile in DHCR24-deficient cells and tissues|
|Live-cell imaging using fluorescent sterol analogs to track cholesterol metabolism and membrane distribution in real-time|
|Cryo-EM structure determination of DHCR24 to understand the molecular basis of sterol reduction and enzyme specificity|
|Developmental analysis of DHCR24 mutant model organisms to study the role of cholesterol in embryogenesis and organogenesis|
