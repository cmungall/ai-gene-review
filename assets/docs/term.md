
# Class: Term



URI: [gene_review:Term](https://w3id.org/ai4curation/gene_review/Term)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[CoreFunction]++-%20anatomical_locations%200..*>[Term&#124;id:string;label:string%20%3F],[CoreFunction]++-%20directly_involved_in%200..*>[Term],[CoreFunction]++-%20in_complex%200..1>[Term],[CoreFunction]++-%20locations%200..*>[Term],[CoreFunction]++-%20molecular_function%200..1>[Term],[CoreFunction]++-%20substrates%200..*>[Term],[Review]++-%20proposed_replacement_terms%200..*>[Term],[GeneReview]++-%20taxon%200..1>[Term],[ExistingAnnotation]++-%20term%200..1>[Term],[AnnotationExtension]++-%20term%200..1>[Term],[Review],[GeneReview],[ExistingAnnotation],[CoreFunction],[AnnotationExtension])](https://yuml.me/diagram/nofunky;dir:TB/class/[CoreFunction]++-%20anatomical_locations%200..*>[Term&#124;id:string;label:string%20%3F],[CoreFunction]++-%20directly_involved_in%200..*>[Term],[CoreFunction]++-%20in_complex%200..1>[Term],[CoreFunction]++-%20locations%200..*>[Term],[CoreFunction]++-%20molecular_function%200..1>[Term],[CoreFunction]++-%20substrates%200..*>[Term],[Review]++-%20proposed_replacement_terms%200..*>[Term],[GeneReview]++-%20taxon%200..1>[Term],[ExistingAnnotation]++-%20term%200..1>[Term],[AnnotationExtension]++-%20term%200..1>[Term],[Review],[GeneReview],[ExistingAnnotation],[CoreFunction],[AnnotationExtension])

## Referenced by Class

 *  **None** *[➞anatomical_locations](coreFunction__anatomical_locations.md)*  <sub>0..\*</sub>  **[Term](Term.md)**
 *  **None** *[➞directly_involved_in](coreFunction__directly_involved_in.md)*  <sub>0..\*</sub>  **[Term](Term.md)**
 *  **None** *[➞in_complex](coreFunction__in_complex.md)*  <sub>0..1</sub>  **[Term](Term.md)**
 *  **None** *[➞locations](coreFunction__locations.md)*  <sub>0..\*</sub>  **[Term](Term.md)**
 *  **None** *[➞molecular_function](coreFunction__molecular_function.md)*  <sub>0..1</sub>  **[Term](Term.md)**
 *  **None** *[➞substrates](coreFunction__substrates.md)*  <sub>0..\*</sub>  **[Term](Term.md)**
 *  **None** *[proposed_replacement_terms](proposed_replacement_terms.md)*  <sub>0..\*</sub>  **[Term](Term.md)**
 *  **None** *[taxon](taxon.md)*  <sub>0..1</sub>  **[Term](Term.md)**
 *  **None** *[term](term.md)*  <sub>0..1</sub>  **[Term](Term.md)**

## Attributes


### Own

 * [id](id.md)  <sub>1..1</sub>
     * Range: [String](types/String.md)
 * [label](label.md)  <sub>0..1</sub>
     * Description: Human readable name of the entity
     * Range: [String](types/String.md)
