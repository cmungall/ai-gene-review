
# Class: Review

A review of an existing annotation.

URI: [gene_review:Review](https://w3id.org/ai4curation/gene_review/Review)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[SupportingTextInReference],[SupportingTextInReference]<supported_by%200..*-++[Review&#124;summary:string%20%3F;action:ActionEnum%20%3F;reason:string%20%3F],[Reference]<additional_reference_ids%200..*-%20[Review],[Term]<proposed_replacement_terms%200..*-++[Review],[ExistingAnnotation]++-%20review%200..1>[Review],[Reference],[ExistingAnnotation])](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[SupportingTextInReference],[SupportingTextInReference]<supported_by%200..*-++[Review&#124;summary:string%20%3F;action:ActionEnum%20%3F;reason:string%20%3F],[Reference]<additional_reference_ids%200..*-%20[Review],[Term]<proposed_replacement_terms%200..*-++[Review],[ExistingAnnotation]++-%20review%200..1>[Review],[Reference],[ExistingAnnotation])

## Referenced by Class

 *  **None** *[review](review.md)*  <sub>0..1</sub>  **[Review](Review.md)**

## Attributes


### Own

 * [summary](summary.md)  <sub>0..1</sub>
     * Description: Summary of the review
     * Range: [String](types/String.md)
 * [action](action.md)  <sub>0..1</sub>
     * Description: Action to be taken
     * Range: [ActionEnum](ActionEnum.md)
 * [reason](reason.md)  <sub>0..1</sub>
     * Description: Reason for the action
     * Range: [String](types/String.md)
 * [proposed_replacement_terms](proposed_replacement_terms.md)  <sub>0..\*</sub>
     * Description: Proposed replacement terms
     * Range: [Term](Term.md)
 * [additional_reference_ids](additional_reference_ids.md)  <sub>0..\*</sub>
     * Description: IDs of the references
     * Range: [Reference](Reference.md)
 * [supported_by](supported_by.md)  <sub>0..\*</sub>
     * Range: [SupportingTextInReference](SupportingTextInReference.md)
