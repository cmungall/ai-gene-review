
# Class: Finding

A finding is a statement about a gene, which is supported by a reference. Similar to "comments" in uniprot

URI: [gene_review:Finding](https://w3id.org/ai4curation/gene_review/Finding)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Reference]++-%20findings%200..*>[Finding&#124;statement:string%20%3F;supporting_text:string%20%3F;reference_section_type:string%20%3F],[Reference])](https://yuml.me/diagram/nofunky;dir:TB/class/[Reference]++-%20findings%200..*>[Finding&#124;statement:string%20%3F;supporting_text:string%20%3F;reference_section_type:string%20%3F],[Reference])

## Referenced by Class

 *  **None** *[findings](findings.md)*  <sub>0..\*</sub>  **[Finding](Finding.md)**

## Attributes


### Own

 * [statement](statement.md)  <sub>0..1</sub>
     * Description: Concise statement describing an aspect of the gene
     * Range: [String](types/String.md)
 * [supporting_text](supporting_text.md)  <sub>0..1</sub>
     * Description: Supporting text from the publication. This should be exact substrings. Different substrings can be broken up by '...'s. These substrings will be checked against the actual text of the paper. If editorialization is necessary, put this in square brackets (this is not checked). For example, you can say '...[CFAP300 shows] transport within cilia is IFT dependent...'
     * Range: [String](types/String.md)
 * [reference_section_type](reference_section_type.md)  <sub>0..1</sub>
     * Description: Type of section in the reference (e.g., 'ABSTRACT', 'METHODS', 'RESULTS', 'DISCUSSION')
     * Range: [String](types/String.md)
