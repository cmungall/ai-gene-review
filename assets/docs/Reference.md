
# Class: Reference



URI: [gene_review:Reference](https://w3id.org/ai4curation/gene_review/Reference)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Finding]<findings%200..*-++[Reference&#124;id:string;title:string%20%3F],[ExistingAnnotation]-%20original_reference_id%200..1>[Reference],[GeneReview]++-%20references%200..*>[Reference],[GeneReview],[Finding],[ExistingAnnotation])](https://yuml.me/diagram/nofunky;dir:TB/class/[Finding]<findings%200..*-++[Reference&#124;id:string;title:string%20%3F],[ExistingAnnotation]-%20original_reference_id%200..1>[Reference],[GeneReview]++-%20references%200..*>[Reference],[GeneReview],[Finding],[ExistingAnnotation])

## Referenced by Class

 *  **None** *[original_reference_id](original_reference_id.md)*  <sub>0..1</sub>  **[Reference](Reference.md)**
 *  **None** *[references](references.md)*  <sub>0..\*</sub>  **[Reference](Reference.md)**

## Attributes


### Own

 * [id](id.md)  <sub>1..1</sub>
     * Range: [String](types/String.md)
 * [title](title.md)  <sub>0..1</sub>
     * Description: Title of the entity
     * Range: [String](types/String.md)
 * [findings](findings.md)  <sub>0..\*</sub>
     * Range: [Finding](Finding.md)
