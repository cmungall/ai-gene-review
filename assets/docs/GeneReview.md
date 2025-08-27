
# Class: GeneReview

Complete review for a gene

URI: [gene_review:GeneReview](https://w3id.org/ai4curation/gene_review/GeneReview)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Reference],[CoreFunction]<core_functions%200..*-++[GeneReview&#124;id:string;gene_symbol:string%20%3F;description:string%20%3F],[ExistingAnnotation]<existing_annotations%200..*-++[GeneReview],[Reference]<references%200..*-++[GeneReview],[Term]<taxon%200..1-++[GeneReview],[ExistingAnnotation],[CoreFunction])](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[Reference],[CoreFunction]<core_functions%200..*-++[GeneReview&#124;id:string;gene_symbol:string%20%3F;description:string%20%3F],[ExistingAnnotation]<existing_annotations%200..*-++[GeneReview],[Reference]<references%200..*-++[GeneReview],[Term]<taxon%200..1-++[GeneReview],[ExistingAnnotation],[CoreFunction])

## Attributes


### Own

 * [id](id.md)  <sub>1..1</sub>
     * Range: [String](types/String.md)
 * [gene_symbol](gene_symbol.md)  <sub>0..1</sub>
     * Description: Symbol of the gene
     * Range: [String](types/String.md)
 * [description](description.md)  <sub>0..1</sub>
     * Description: Description of the entity
     * Range: [String](types/String.md)
 * [taxon](taxon.md)  <sub>0..1</sub>
     * Range: [Term](Term.md)
 * [references](references.md)  <sub>0..\*</sub>
     * Range: [Reference](Reference.md)
 * [existing_annotations](existing_annotations.md)  <sub>0..\*</sub>
     * Range: [ExistingAnnotation](ExistingAnnotation.md)
 * [core_functions](core_functions.md)  <sub>0..\*</sub>
     * Range: [CoreFunction](CoreFunction.md)
