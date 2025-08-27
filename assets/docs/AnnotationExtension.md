
# Class: AnnotationExtension



URI: [gene_review:AnnotationExtension](https://w3id.org/ai4curation/gene_review/AnnotationExtension)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Term]<term%200..1-++[AnnotationExtension&#124;predicate:string%20%3F],[ExistingAnnotation]++-%20extensions%200..*>[AnnotationExtension],[ExistingAnnotation])](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Term]<term%200..1-++[AnnotationExtension&#124;predicate:string%20%3F],[ExistingAnnotation]++-%20extensions%200..*>[AnnotationExtension],[ExistingAnnotation])

## Referenced by Class

 *  **None** *[extensions](extensions.md)*  <sub>0..\*</sub>  **[AnnotationExtension](AnnotationExtension.md)**

## Attributes


### Own

 * [predicate](predicate.md)  <sub>0..1</sub>
     * Description: Predicate of the extension
     * Range: [String](types/String.md)
 * [term](term.md)  <sub>0..1</sub>
     * Description: Term to be annotated
     * Range: [Term](Term.md)
