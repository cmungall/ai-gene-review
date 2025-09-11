
# Class: ProposedOntologyTerm

A proposed new ontology term that should exist but doesn't currently

URI: [gene_review:ProposedOntologyTerm](https://w3id.org/ai4curation/gene_review/ProposedOntologyTerm)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[SupportingTextInReference],[SupportingTextInReference]<supported_by%200..*-++[ProposedOntologyTerm&#124;proposed_name:string;proposed_definition:string;justification:string%20%3F],[GeneReview]++-%20proposed_new_terms%200..*>[ProposedOntologyTerm],[GeneReview])](https://yuml.me/diagram/nofunky;dir:TB/class/[SupportingTextInReference],[SupportingTextInReference]<supported_by%200..*-++[ProposedOntologyTerm&#124;proposed_name:string;proposed_definition:string;justification:string%20%3F],[GeneReview]++-%20proposed_new_terms%200..*>[ProposedOntologyTerm],[GeneReview])

## Referenced by Class

 *  **None** *[proposed_new_terms](proposed_new_terms.md)*  <sub>0..\*</sub>  **[ProposedOntologyTerm](ProposedOntologyTerm.md)**

## Attributes


### Own

 * [➞proposed_name](proposedOntologyTerm__proposed_name.md)  <sub>1..1</sub>
     * Description: Proposed name for the new term
     * Range: [String](types/String.md)
 * [➞proposed_definition](proposedOntologyTerm__proposed_definition.md)  <sub>1..1</sub>
     * Description: Proposed definition for the new term
     * Range: [String](types/String.md)
 * [➞justification](proposedOntologyTerm__justification.md)  <sub>0..1</sub>
     * Description: Justification for why this term is needed
     * Range: [String](types/String.md)
 * [➞supported_by](proposedOntologyTerm__supported_by.md)  <sub>0..\*</sub>
     * Range: [SupportingTextInReference](SupportingTextInReference.md)
