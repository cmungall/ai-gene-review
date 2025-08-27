
# Class: Finding



URI: [gene_review:Finding](https://w3id.org/ai4curation/gene_review/Finding)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Reference]++-%20findings%200..*>[Finding&#124;statement:string%20%3F;supporting_text:string%20%3F],[Reference])](https://yuml.me/diagram/nofunky;dir:TB/class/[Reference]++-%20findings%200..*>[Finding&#124;statement:string%20%3F;supporting_text:string%20%3F],[Reference])

## Referenced by Class

 *  **None** *[findings](findings.md)*  <sub>0..\*</sub>  **[Finding](Finding.md)**

## Attributes


### Own

 * [statement](statement.md)  <sub>0..1</sub>
     * Description: Concise statement describing an aspect of the gene
     * Range: [String](types/String.md)
 * [supporting_text](supporting_text.md)  <sub>0..1</sub>
     * Description: Supporting text from the publication. Can be interleaved text if there are multiple lines of evidence. Can also be annotated
     * Range: [String](types/String.md)
