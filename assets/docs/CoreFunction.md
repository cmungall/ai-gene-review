
# Class: CoreFunction

A core function is a GO-CAM-like annotation of the core evolved functions of a gene. This is a synthesis of the reviewed core annotations, brought together into a unified GO-CAM-like representation.

URI: [gene_review:CoreFunction](https://w3id.org/ai4curation/gene_review/CoreFunction)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[SupportingTextInReference],[Term]<in_complex%200..1-++[CoreFunction&#124;description:string%20%3F],[Term]<substrates%200..*-++[CoreFunction],[Term]<anatomical_locations%200..*-++[CoreFunction],[Term]<locations%200..*-++[CoreFunction],[Term]<directly_involved_in%200..*-++[CoreFunction],[Term]<molecular_function%200..1-++[CoreFunction],[SupportingTextInReference]<supported_by%200..*-++[CoreFunction],[GeneReview]++-%20core_functions%200..*>[CoreFunction],[GeneReview])](https://yuml.me/diagram/nofunky;dir:TB/class/[Term],[SupportingTextInReference],[Term]<in_complex%200..1-++[CoreFunction&#124;description:string%20%3F],[Term]<substrates%200..*-++[CoreFunction],[Term]<anatomical_locations%200..*-++[CoreFunction],[Term]<locations%200..*-++[CoreFunction],[Term]<directly_involved_in%200..*-++[CoreFunction],[Term]<molecular_function%200..1-++[CoreFunction],[SupportingTextInReference]<supported_by%200..*-++[CoreFunction],[GeneReview]++-%20core_functions%200..*>[CoreFunction],[GeneReview])

## Referenced by Class

 *  **None** *[core_functions](core_functions.md)*  <sub>0..\*</sub>  **[CoreFunction](CoreFunction.md)**

## Attributes


### Own

 * [➞description](coreFunction__description.md)  <sub>0..1</sub>
     * Description: Description of the core function
     * Range: [String](types/String.md)
 * [➞supported_by](coreFunction__supported_by.md)  <sub>0..\*</sub>
     * Range: [SupportingTextInReference](SupportingTextInReference.md)
 * [➞molecular_function](coreFunction__molecular_function.md)  <sub>0..1</sub>
     * Range: [Term](Term.md)
 * [➞directly_involved_in](coreFunction__directly_involved_in.md)  <sub>0..\*</sub>
     * Range: [Term](Term.md)
 * [➞locations](coreFunction__locations.md)  <sub>0..\*</sub>
     * Range: [Term](Term.md)
 * [➞anatomical_locations](coreFunction__anatomical_locations.md)  <sub>0..\*</sub>
     * Range: [Term](Term.md)
 * [➞substrates](coreFunction__substrates.md)  <sub>0..\*</sub>
     * Range: [Term](Term.md)
 * [➞in_complex](coreFunction__in_complex.md)  <sub>0..1</sub>
     * Range: [Term](Term.md)
