# Gpd3 (SPBC354.12) in *Schizosaccharomyces pombe*

Gpd3 (UniProt O43026) is a glyceraldehyde-3-phosphate dehydrogenase in fission yeast[\[1\]](https://thebiogrid.org/277491/summary/schizosaccharomyces-pombe/gpd3.html#:~:text=). It catalyzes the NAD⁺-dependent oxidation of glyceraldehyde-3-phosphate to 1,3-bisphosphoglycerate in glycolysis. BioGRID/GO annotations list “canonical glycolysis” and “glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity” for Gpd3[\[1\]](https://thebiogrid.org/277491/summary/schizosaccharomyces-pombe/gpd3.html#:~:text=), and place it in the cytosol. Gpd3 is one of two GAPDH isoforms in *S. pombe* (the other is Tdh1); Gpd3 is less abundantly expressed under rich-growth conditions.

## Enzymatic Role and Pathway Context

* **Glycolytic function:** Gpd3’s primary role is as a glycolytic enzyme[\[1\]](https://thebiogrid.org/277491/summary/schizosaccharomyces-pombe/gpd3.html#:~:text=). It likely forms a homotetramer with a Rossmann-fold NAD-binding domain, similar to other GAPDHs. Its catalytic cysteine residue (analogous to Cys-152 in Tdh1) attacks glyceraldehyde-3-P to form a thiohemiacetal, generating NADH and 1,3-bisphosphoglycerate.

* **Metabolic integration:** As part of the central carbon metabolism, Gpd3 works closely with triosephosphate isomerase, phosphoglycerate kinase, etc., to extract energy from glucose. It also influences the cellular NADH/NAD⁺ balance. STRING network data confirm Gpd3’s strong associations with glycolytic enzymes (e.g. tpi1, pgk1, fba1)[\[2\]](https://string-db.org/network/284812.O43026#:~:text=tpi1)[\[3\]](https://string-db.org/network/284812.O43026#:~:text=tdh1), reflecting its role in the glycolytic “metabolon.”

* **Isozyme redundancy:** *S. pombe* has two GAPDH genes (TDH1 and GPD3). Tdh1 is expressed at high levels and is the major glycolytic isozyme, whereas Gpd3 is expressed at a much lower level. For example, Morigasaki *et al.* note that a *tdh1*Δ strain remains viable on glucose because *gpd3* is present as a second GAPDH gene[\[4\]](https://pubmed.ncbi.nlm.nih.gov/18406331/#:~:text=protein%20and%20finally%20to%20the,sensitive%20cysteine%20residue%20may%20provide). Thus Gpd3 provides backup glycolytic activity. It is presumed that a double *tdh1Δ gpd3Δ* would be inviable (complete loss of GAPDH), though this double mutant has not been reported.

## Non-catalytic (“Moonlighting”) and Regulatory Functions

Evidence from fission yeast and other systems suggests Gpd3/GAPDH may have roles beyond metabolism:  
\- **Stress signaling:** GAPDH can act as a redox sensor. In *S. pombe*, Tdh1 (GAPDH1) was shown to bind components of a two-component MAPK cascade under H₂O₂ stress[\[4\]](https://pubmed.ncbi.nlm.nih.gov/18406331/#:~:text=protein%20and%20finally%20to%20the,sensitive%20cysteine%20residue%20may%20provide). Oxidation of Tdh1’s active-site Cys-152 upon peroxide exposure enhances its binding to the Mcs4 response regulator, which is required for downstream stress signaling[\[4\]](https://pubmed.ncbi.nlm.nih.gov/18406331/#:~:text=protein%20and%20finally%20to%20the,sensitive%20cysteine%20residue%20may%20provide). By analogy, Gpd3 likely has a corresponding cysteine and might similarly participate in redox signaling if activated.  
\- **Transcriptional cofactor:** Proteomic work revealed that GAPDH associates with RNA polymerase II subunits in *S. pombe*. Mitsuzawa *et al.* found GAPDH physically interacts with the Rpb7 subunit of Pol II; in fact, GAPDH was affinity-purified via an Rpb4/7 column[\[5\]](https://pubmed.ncbi.nlm.nih.gov/15620689/#:~:text=implicated%20in%20transcriptional%20activation%20in,interaction%20between%20actin%20and%20Rpb7). This suggests GAPDH can bind transcription machinery, potentially influencing gene expression. Although that study did not distinguish Tdh1 vs Gpd3, it indicates that *fission yeast* GAPDH can enter the nucleus and bind Pol II complexes. By extension, Gpd3 (the minor isozyme) may also engage in such interactions under certain conditions.  
\- **RNA processing:** Intriguingly, Gpd3 co-purified with the splicing factor Rbm10[\[6\]](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-021-00382-y/tables/1#:~:text=SPBC354). In a mass-spectrometry survey of Rbm10-associated proteins, the Rbm10 pulldown included SPBC354.12 (Gpd3, “glyceraldehyde 3-phosphate dehydrogenase Gpd3”)[\[6\]](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-021-00382-y/tables/1#:~:text=SPBC354). This implies a physical association between Gpd3 and the splicing apparatus, suggesting yet another nuclear role.  
\- **Cytoskeleton/membrane:** In budding yeast, GAPDH is known to bind actin and membranes; in *S. pombe*, actin was also found associated with Pol II[\[5\]](https://pubmed.ncbi.nlm.nih.gov/15620689/#:~:text=implicated%20in%20transcriptional%20activation%20in,interaction%20between%20actin%20and%20Rpb7). It is plausible Gpd3 interacts with the cytoskeleton or membranes under certain conditions, but direct evidence in fission yeast is lacking.

## Regulation, Genetic Interactions, and System-level Effects

* **Expression control:** Gpd3 transcription appears constitutive; it is not sharply induced by stress. In fact, Gpd3 mRNA has been used as a “housekeeping” control in expression studies (for example, Gpd3 was the normalization standard in a thiamine/glucose study[\[7\]](https://www.cellmolbiol.org/index.php/CMB/article/download/2679/1402/6665#:~:text=to%20gpd3%20gene,P%3C0%2C001)). Its promoter does not contain obvious stress-response elements known for *tdh1*. However, Gpd3 sits in a subtelomeric gene cluster repressed by chromatin factors: deletion of the Vtc4 polyphosphate polymerase (which normally contributes to heterochromatic silencing at telomeres) leads to de-repression of several subtelomeric genes including *gpd3*, and in vtc4Δ cells *gpd3* expression is down \~1.6-fold[\[8\]](https://www.micropublication.org/static/pdf/micropub-biology-000744.pdf#:~:text=SPBC354.12%20gpd3%20glyceraldehyde%203,1.29). This suggests that Gpd3 is normally repressed in that chromosomal region, possibly via TORC2/HDAC pathways.

* **Phenotypes and screens:** A genome-wide MBF reporter screen (tracking G1/S transcription) found that *gpd3* deletion modestly elevates MBF-driven reporter expression (YFP/FSC ratio 1.57)[\[9\]](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845923/#:~:text=SPBC354,56). This implies Gpd3 normally dampens G1/S transcription or that loss of Gpd3 triggers a compensatory cell-cycle response. No obvious growth or cell-cycle arrest was reported for gpd3Δ alone, consistent with its redundant role. However, gpd3Δ might exhibit subtle stress or metabolic phenotypes (e.g. altered NADH/NAD⁺ ratio) that remain to be characterized. In other genetic studies (e.g. cohesin function screens[\[10\]](https://research.bioinformatics.udel.edu/iptmnet/entry/O43026/#:~:text=K5%20%20Ubiquitination%20%20score1,Ubiquitination%20%20score1%20PomBase), CLS studies), Gpd3 has not emerged as a strong hit, but it may still contribute to broader networks.

* **Interaction networks:** The STRING database clusters Gpd3 with canonical glycolytic and NAD⁺/redox enzymes[\[2\]](https://string-db.org/network/284812.O43026#:~:text=tpi1)[\[3\]](https://string-db.org/network/284812.O43026#:~:text=tdh1), reaffirming its metabolic role. No high-confidence novel protein complexes are known for Gpd3 aside from the multi-enzyme glycolytic complexes.

## Comparative and Evolutionary Perspective

Gpd3 is a member of the universally conserved GAPDH family (PF00121). Multiple organisms have several GAPDH isoforms: for example, *S. cerevisiae* has three (Tdh1-3) and mammals have somatic vs sperm-specific GAPDH. In *S. pombe*, the two isoforms (Tdh1 and Gpd3) presumably arose by gene duplication. Sequence identity between Tdh1 and Gpd3 is high (likely \~70–80%), with conservation of catalytic residues. Orthologs of Gpd3 exist in other Schizosaccharomyces species (e.g. *S. japonicus* SJAG\_03828) and in a wide range of fungi. Phylogenetically, yeast GAPDHs form a clade separate from bacterial GAPDHs. The retention of a second GAPDH isozyme in *S. pombe* suggests adaptive significance: possibly one isoform (Tdh1) is optimized for bulk glycolysis, while Gpd3 may be specialized for regulation or backup under specific conditions. The split into glycolytic vs glycerol-3-PDH enzymes also differs: *S. cerevisiae* uses Gpd1/2 for glycerol metabolism, whereas *S. pombe* uses separate genes (gpd1/SPBC215.05 and gpd2) for glycerol-3-PDH, so there is no confusion between “GPD” for glycerol-PDH vs “TDH/GPD3” for glycolysis.

## Open Questions and Future Directions

Despite its clear glycolytic role, many aspects of Gpd3 remain unknown. Key questions include:

* **Nuclear roles:** Under what conditions does Gpd3 enter the nucleus or bind nuclear proteins? The Rbm10/Rpb7 associations suggest nuclear targeting – but is this condition-dependent (e.g. stress, cell cycle stage)? Does Gpd3 carry a nuclear localization signal or depend on other factors to relocate?

* **Crosstalk with Tdh1:** Does Gpd3 compensate when Tdh1 is inactivated (e.g. by oxidation) or absent? For instance, is Gpd3 upregulated when Tdh1 is mutated? Conversely, does overexpressing Gpd3 rescue any *tdh1* phenotypes?

* **Regulation by signaling pathways:** Is *gpd3* transcription or Gpd3 activity regulated by nutrient or stress-sensing kinases (e.g. TORC2, AMPK, Sty1)? The subtelomeric silencing hints at chromatin-level control; do glucose or oxidative signals modulate Gpd3 via post-translational modifications?

* **Physiological impact of loss:** What phenotypes does *gpd3Δ* display in detailed assays? Possible experiments: measure growth on various carbon sources, sensitivity to H₂O₂ or other stresses, and changes in glycolytic flux or NAD⁺ levels. Given the MBF screen result, does *gpd3Δ* alter cell-cycle timing or S-phase entry?

* **Protein interactions:** Are there additional, yet-unidentified binding partners of Gpd3 (e.g. metabolic enzymes, signaling proteins, RNA-binding proteins)? Does Gpd3 form mixed tetramers with Tdh1?

* **Evolutionary divergence:** Do the extra sequences (or lack thereof) between Tdh1 and Gpd3 confer specific regulatory features (e.g. unique phosphorylation sites)? For example, iPTMnet annotations list predicted phosphorylation/sumoylation sites in Gpd3[\[10\]](https://research.bioinformatics.udel.edu/iptmnet/entry/O43026/#:~:text=K5%20%20Ubiquitination%20%20score1,Ubiquitination%20%20score1%20PomBase)[\[11\]](https://research.bioinformatics.udel.edu/iptmnet/entry/O43026/#:~:text=K193%20%20Ubiquitination%20%20score1,PomBase%20UniProt%20%2018257517%2C%2030726745); are any of these functionally relevant?

**Experimental priorities:** To address these questions, useful approaches would be:  
\- *Localization studies:* Tag Gpd3 with GFP (or immunostain) to track subcellular distribution under different conditions (log vs stationary phase, \+/– stress).  
\- *Protein complex mapping:* Perform Gpd3 co-immunoprecipitation followed by mass spectrometry in wild-type vs stress conditions to identify novel partners (e.g. splicing or transcription factors). Reciprocal pulldowns of Rbm10 or Rpb7 could confirm Gpd3 involvement.  
\- *Mutational analysis:* Create a catalytically inactive Gpd3 (Cys-to-Ser) to see if it retains non-catalytic functions (e.g. stress signaling). Conversely, overexpress Gpd3 or Cys mutants to test effects on stress response and metabolism.  
\- *Gene expression and metabolism:* Quantify *gpd3* mRNA/protein in mutants of chromatin regulators (vtc4Δ, sir2Δ, etc.) and nutrient shifts. Use ^13C-labeling or enzymatic assays to compare glycolytic flux and NADH/NAD⁺ ratio in WT vs *gpd3Δ*.  
\- *Phenotype screens:* Systematically test *gpd3Δ* (and double *tdh1Δ gpd3Δ*) for growth on various media, lifespan assays, and synthetic interactions with stress-signaling mutants. Monitor MBF-target gene expression and cell-cycle profiles in *gpd3Δ*.

These studies would clarify Gpd3’s roles beyond glycolysis and how it integrates with cellular physiology.

**References:** Gpd3’s annotated glycolytic function is documented in PomBase/UniProt[\[1\]](https://thebiogrid.org/277491/summary/schizosaccharomyces-pombe/gpd3.html#:~:text=). Key findings on GAPDH moonlighting in *S. pombe* include Morigasaki *et al.* (2008), who showed GAPDH (Tdh1) is required for H₂O₂-activated MAPK signaling via its redox-sensitive Cys[\[4\]](https://pubmed.ncbi.nlm.nih.gov/18406331/#:~:text=protein%20and%20finally%20to%20the,sensitive%20cysteine%20residue%20may%20provide), and Mitsuzawa *et al.* (2005), who showed GAPDH binds RNA Pol II subunit Rpb7[\[5\]](https://pubmed.ncbi.nlm.nih.gov/15620689/#:~:text=implicated%20in%20transcriptional%20activation%20in,interaction%20between%20actin%20and%20Rpb7). Proteomic screens by Weigt *et al.* (2021) found Gpd3 in Rbm10 complexes[\[6\]](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-021-00382-y/tables/1#:~:text=SPBC354). Recent functional screens uncovered a *gpd3*Δ effect on G1/S transcription[\[9\]](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845923/#:~:text=SPBC354,56) and the linkage of *gpd3* to a subtelomeric silencing cluster[\[8\]](https://www.micropublication.org/static/pdf/micropub-biology-000744.pdf#:~:text=SPBC354.12%20gpd3%20glyceraldehyde%203,1.29). These sources underpin the views above.

---

[\[1\]](https://thebiogrid.org/277491/summary/schizosaccharomyces-pombe/gpd3.html#:~:text=) gpd3 (SPBC354.12) Result Summary | BioGRID

[https://thebiogrid.org/277491/summary/schizosaccharomyces-pombe/gpd3.html](https://thebiogrid.org/277491/summary/schizosaccharomyces-pombe/gpd3.html)

[\[2\]](https://string-db.org/network/284812.O43026#:~:text=tpi1) [\[3\]](https://string-db.org/network/284812.O43026#:~:text=tdh1) gpd3 protein (Schizosaccharomyces pombe) \- STRING interaction network

[https://string-db.org/network/284812.O43026](https://string-db.org/network/284812.O43026)

[\[4\]](https://pubmed.ncbi.nlm.nih.gov/18406331/#:~:text=protein%20and%20finally%20to%20the,sensitive%20cysteine%20residue%20may%20provide) Glycolytic enzyme GAPDH promotes peroxide stress signaling through multistep phosphorelay to a MAPK cascade \- PubMed

[https://pubmed.ncbi.nlm.nih.gov/18406331/](https://pubmed.ncbi.nlm.nih.gov/18406331/)

[\[5\]](https://pubmed.ncbi.nlm.nih.gov/15620689/#:~:text=implicated%20in%20transcriptional%20activation%20in,interaction%20between%20actin%20and%20Rpb7) Glyceraldehyde-3-phosphate dehydrogenase and actin associate with RNA polymerase II and interact with its Rpb7 subunit \- PubMed

[https://pubmed.ncbi.nlm.nih.gov/15620689/](https://pubmed.ncbi.nlm.nih.gov/15620689/)

[\[6\]](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-021-00382-y/tables/1#:~:text=SPBC354) Rbm10 facilitates heterochromatin assembly via the Clr6 HDAC complex | Epigenetics & Chromatin | Full Text

[https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-021-00382-y/tables/1](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-021-00382-y/tables/1)

[\[7\]](https://www.cellmolbiol.org/index.php/CMB/article/download/2679/1402/6665#:~:text=to%20gpd3%20gene,P%3C0%2C001) cellmolbiol.org

[https://www.cellmolbiol.org/index.php/CMB/article/download/2679/1402/6665](https://www.cellmolbiol.org/index.php/CMB/article/download/2679/1402/6665)

[\[8\]](https://www.micropublication.org/static/pdf/micropub-biology-000744.pdf#:~:text=SPBC354.12%20gpd3%20glyceraldehyde%203,1.29) micropublication.org

[https://www.micropublication.org/static/pdf/micropub-biology-000744.pdf](https://www.micropublication.org/static/pdf/micropub-biology-000744.pdf)

[\[9\]](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845923/#:~:text=SPBC354,56)  A functional genome-wide genetic screening identifies new pathways controlling the G1/S transcriptional wave \- PMC 

[https://pmc.ncbi.nlm.nih.gov/articles/PMC4845923/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845923/)

[\[10\]](https://research.bioinformatics.udel.edu/iptmnet/entry/O43026/#:~:text=K5%20%20Ubiquitination%20%20score1,Ubiquitination%20%20score1%20PomBase) [\[11\]](https://research.bioinformatics.udel.edu/iptmnet/entry/O43026/#:~:text=K193%20%20Ubiquitination%20%20score1,PomBase%20UniProt%20%2018257517%2C%2030726745) iPTMnet Report O43026 gpd3 

[https://research.bioinformatics.udel.edu/iptmnet/entry/O43026/](https://research.bioinformatics.udel.edu/iptmnet/entry/O43026/)