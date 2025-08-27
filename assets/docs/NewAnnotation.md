
# Class: NewAnnotation



URI: [ai4curation:NewAnnotation](https://w3id.org/ai4curation/NewAnnotation)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Reference],[Reference]<original_reference_id%200..1-%20[NewAnnotation&#124;negated:boolean%20%3F;evidence_type:string%20%3F;reason:string%20%3F],[AnnotationExtension]<extensions%200..*-++[NewAnnotation],[Term]<term%200..1-++[NewAnnotation],[GeneReview]++-%20proposed_new_annotations%200..*>[NewAnnotation],[GeneReview],[AnnotationExtension])](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Reference],[Reference]<original_reference_id%200..1-%20[NewAnnotation&#124;negated:boolean%20%3F;evidence_type:string%20%3F;reason:string%20%3F],[AnnotationExtension]<extensions%200..*-++[NewAnnotation],[Term]<term%200..1-++[NewAnnotation],[GeneReview]++-%20proposed_new_annotations%200..*>[NewAnnotation],[GeneReview],[AnnotationExtension])

## Referenced by Class

 *  **None** *[proposed_new_annotations](proposed_new_annotations.md)*  <sub>0..\*</sub>  **[NewAnnotation](NewAnnotation.md)**

## Attributes


### Own

 * [term](term.md)  <sub>0..1</sub>
     * Description: Term to be annotated
     * Range: [Term](Term.md)
 * [extensions](extensions.md)  <sub>0..\*</sub>
     * Range: [AnnotationExtension](AnnotationExtension.md)
 * [negated](negated.md)  <sub>0..1</sub>
     * Description: Whether the term is negated
     * Range: [Boolean](types/Boolean.md)
 * [evidence_type](evidence_type.md)  <sub>0..1</sub>
     * Description: Evidence code (e.g., IDA, IBA, ISS, TAS)
     * Range: [String](types/String.md)
 * [original_reference_id](original_reference_id.md)  <sub>0..1</sub>
     * Description: ID of the original reference
     * Range: [Reference](Reference.md)
 * [reason](reason.md)  <sub>0..1</sub>
     * Description: Reason for the action
     * Range: [String](types/String.md)
