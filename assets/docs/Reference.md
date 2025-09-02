
# Class: Reference

A reference is a published text  that describes a finding or a method. References might be formal publications (where the ID is a PMID), or for methods, a GO_REF. Additionally, a reference to a local ad-hoc analysis or review can be made by using the `file:` prefix.

URI: [gene_review:Reference](https://w3id.org/ai4curation/gene_review/Reference)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Finding]<findings%200..*-++[Reference&#124;id:string;title:string;is_invalid:boolean%20%3F],[Review]-%20additional_reference_ids%200..*>[Reference],[ExistingAnnotation]-%20original_reference_id%200..1>[Reference],[SupportingTextInReference]-%20reference_id%200..1>[Reference],[GeneReview]++-%20references%200..*>[Reference],[SupportingTextInReference],[Review],[GeneReview],[Finding],[ExistingAnnotation])](https://yuml.me/diagram/nofunky;dir:TB/class/[Finding]<findings%200..*-++[Reference&#124;id:string;title:string;is_invalid:boolean%20%3F],[Review]-%20additional_reference_ids%200..*>[Reference],[ExistingAnnotation]-%20original_reference_id%200..1>[Reference],[SupportingTextInReference]-%20reference_id%200..1>[Reference],[GeneReview]++-%20references%200..*>[Reference],[SupportingTextInReference],[Review],[GeneReview],[Finding],[ExistingAnnotation])

## Referenced by Class

 *  **None** *[additional_reference_ids](additional_reference_ids.md)*  <sub>0..\*</sub>  **[Reference](Reference.md)**
 *  **None** *[original_reference_id](original_reference_id.md)*  <sub>0..1</sub>  **[Reference](Reference.md)**
 *  **None** *[reference_id](reference_id.md)*  <sub>0..1</sub>  **[Reference](Reference.md)**
 *  **None** *[references](references.md)*  <sub>0..\*</sub>  **[Reference](Reference.md)**

## Attributes


### Own

 * [id](id.md)  <sub>1..1</sub>
     * Range: [String](types/String.md)
 * [title](title.md)  <sub>1..1</sub>
     * Description: Title of the entity
     * Range: [String](types/String.md)
 * [findings](findings.md)  <sub>0..\*</sub>
     * Range: [Finding](Finding.md)
 * [is_invalid](is_invalid.md)  <sub>0..1</sub>
     * Description: Whether the reference is invalid (e.g., retracted or replaced)
     * Range: [Boolean](types/Boolean.md)
