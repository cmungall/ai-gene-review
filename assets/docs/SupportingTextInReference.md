
# Class: SupportingTextInReference

A supporting text in a reference.

URI: [gene_review:SupportingTextInReference](https://w3id.org/ai4curation/gene_review/SupportingTextInReference)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Reference]<reference_id%200..1-%20[SupportingTextInReference&#124;supporting_text:string%20%3F;full_text_unavailable:boolean%20%3F;reference_section_type:ManuscriptSection%20%3F],[CoreFunction]++-%20supported_by%200..*>[SupportingTextInReference],[ProposedOntologyTerm]++-%20supported_by%200..*>[SupportingTextInReference],[Review]++-%20supported_by%200..*>[SupportingTextInReference],[Review],[Reference],[ProposedOntologyTerm],[CoreFunction])](https://yuml.me/diagram/nofunky;dir:TB/class/[Reference]<reference_id%200..1-%20[SupportingTextInReference&#124;supporting_text:string%20%3F;full_text_unavailable:boolean%20%3F;reference_section_type:ManuscriptSection%20%3F],[CoreFunction]++-%20supported_by%200..*>[SupportingTextInReference],[ProposedOntologyTerm]++-%20supported_by%200..*>[SupportingTextInReference],[Review]++-%20supported_by%200..*>[SupportingTextInReference],[Review],[Reference],[ProposedOntologyTerm],[CoreFunction])

## Referenced by Class

 *  **None** *[➞supported_by](coreFunction__supported_by.md)*  <sub>0..\*</sub>  **[SupportingTextInReference](SupportingTextInReference.md)**
 *  **None** *[➞supported_by](proposedOntologyTerm__supported_by.md)*  <sub>0..\*</sub>  **[SupportingTextInReference](SupportingTextInReference.md)**
 *  **None** *[supported_by](supported_by.md)*  <sub>0..\*</sub>  **[SupportingTextInReference](SupportingTextInReference.md)**

## Attributes


### Own

 * [reference_id](reference_id.md)  <sub>0..1</sub>
     * Range: [Reference](Reference.md)
 * [supporting_text](supporting_text.md)  <sub>0..1</sub>
     * Description: Supporting text from the publication. This should be exact substrings. Different substrings can be broken up by '...'s. These substrings will be checked against the actual text of the paper. If editorialization is necessary, put this in square brackets (this is not checked). For example, you can say '...[CFAP300 shows] transport within cilia is IFT dependent...'
     * Range: [String](types/String.md)
 * [full_text_unavailable](full_text_unavailable.md)  <sub>0..1</sub>
     * Description: Whether the full text is unavailable
     * Range: [Boolean](types/Boolean.md)
 * [reference_section_type](reference_section_type.md)  <sub>0..1</sub>
     * Description: Type of section in the reference (e.g., 'ABSTRACT', 'METHODS', 'RESULTS', 'DISCUSSION')
     * Range: [ManuscriptSection](ManuscriptSection.md)
