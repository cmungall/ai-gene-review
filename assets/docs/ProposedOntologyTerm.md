
# Class: ProposedOntologyTerm

A proposed new ontology term that should exist but doesn't currently

URI: [gene_review:ProposedOntologyTerm](https://w3id.org/ai4curation/gene_review/ProposedOntologyTerm)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[GeneReview]++-%20proposed_new_terms%200..*>[ProposedOntologyTerm&#124;proposed_name:string;proposed_definition:string;justification:string%20%3F],[GeneReview])](https://yuml.me/diagram/nofunky;dir:TB/class/[GeneReview]++-%20proposed_new_terms%200..*>[ProposedOntologyTerm&#124;proposed_name:string;proposed_definition:string;justification:string%20%3F],[GeneReview])

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
