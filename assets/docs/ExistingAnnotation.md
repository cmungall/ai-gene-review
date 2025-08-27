
# Class: ExistingAnnotation



URI: [gene_review:ExistingAnnotation](https://w3id.org/ai4curation/gene_review/ExistingAnnotation)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Review],[Reference],[Review]<review%200..1-++[ExistingAnnotation&#124;negated:boolean%20%3F;evidence_type:string%20%3F],[Reference]<original_reference_id%200..1-%20[ExistingAnnotation],[AnnotationExtension]<extensions%200..*-++[ExistingAnnotation],[Term]<term%200..1-++[ExistingAnnotation],[GeneReview]++-%20existing_annotations%200..*>[ExistingAnnotation],[GeneReview],[AnnotationExtension])](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Review],[Reference],[Review]<review%200..1-++[ExistingAnnotation&#124;negated:boolean%20%3F;evidence_type:string%20%3F],[Reference]<original_reference_id%200..1-%20[ExistingAnnotation],[AnnotationExtension]<extensions%200..*-++[ExistingAnnotation],[Term]<term%200..1-++[ExistingAnnotation],[GeneReview]++-%20existing_annotations%200..*>[ExistingAnnotation],[GeneReview],[AnnotationExtension])

## Referenced by Class

 *  **None** *[existing_annotations](existing_annotations.md)*  <sub>0..\*</sub>  **[ExistingAnnotation](ExistingAnnotation.md)**

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
 * [review](review.md)  <sub>0..1</sub>
     * Description: Review of the gene
     * Range: [Review](Review.md)
