
# gene_curation


**metamodel version:** 1.7.0

**version:** None


Schema for gene curation Top level entity is a GeneReview, which is about a single gene (and its equivalent swiss-prot entry). It contains a high level summary of the gene, plus a review of all existing annotations. It also contains a list of core functions, which are GO-CAM-like annotons describing the core evolved functions of the gene.


### Classes

 * [AnnotationExtension](AnnotationExtension.md)
 * [CoreFunction](CoreFunction.md) - A core function is a GO-CAM-like annotation of the core evolved functions of a gene. This is a synthesis of the reviewed core annotations, brought together into a unified GO-CAM-like representation.
 * [ExistingAnnotation](ExistingAnnotation.md) - An existing annotation from the GO database, plus a review of the annotation.
 * [Finding](Finding.md) - A finding is a statement about a gene, which is supported by a reference. Similar to "comments" in uniprot
 * [GeneReview](GeneReview.md) - Complete review for a gene
 * [ProposedOntologyTerm](ProposedOntologyTerm.md) - A proposed new ontology term that should exist but doesn't currently
 * [Reference](Reference.md) - A reference is a published text  that describes a finding or a method. References might be formal publications (where the ID is a PMID), or for methods, a GO_REF. Additionally, a reference to a local ad-hoc analysis or review can be made by using the `file:` prefix.
 * [Review](Review.md) - A review of an existing annotation.
 * [SupportingTextInReference](SupportingTextInReference.md) - A supporting text in a reference.
 * [Term](Term.md) - A term in a specific ontology

### Mixins


### Slots

 * [action](action.md) - Action to be taken
 * [additional_reference_ids](additional_reference_ids.md) - IDs of the references
 * [aliases](aliases.md)
 * [➞anatomical_locations](coreFunction__anatomical_locations.md)
 * [➞description](coreFunction__description.md) - Description of the core function
 * [➞directly_involved_in](coreFunction__directly_involved_in.md)
 * [➞in_complex](coreFunction__in_complex.md)
 * [➞locations](coreFunction__locations.md)
 * [➞molecular_function](coreFunction__molecular_function.md)
 * [➞substrates](coreFunction__substrates.md)
 * [➞supported_by](coreFunction__supported_by.md)
 * [core_functions](core_functions.md)
 * [description](description.md) - Description of the entity
 * [evidence_type](evidence_type.md) - Evidence code (e.g., IDA, IBA, ISS, TAS)
 * [existing_annotations](existing_annotations.md)
 * [extensions](extensions.md)
 * [findings](findings.md)
 * [gene_symbol](gene_symbol.md) - Symbol of the gene
 * [id](id.md)
 * [is_invalid](is_invalid.md) - Whether the reference is invalid (e.g., retracted or replaced)
 * [label](label.md) - Human readable name of the entity
 * [negated](negated.md) - Whether the term is negated
 * [ontology](ontology.md) - Ontology of the term. E.g `go`, `cl`, `hp`
 * [original_reference_id](original_reference_id.md) - ID of the original reference
 * [predicate](predicate.md) - Predicate of the extension
 * [➞justification](proposedOntologyTerm__justification.md) - Justification for why this term is needed
 * [➞proposed_definition](proposedOntologyTerm__proposed_definition.md) - Proposed definition for the new term
 * [➞proposed_name](proposedOntologyTerm__proposed_name.md) - Proposed name for the new term
 * [proposed_new_terms](proposed_new_terms.md) - Proposed new ontology terms that should exist but don't
 * [proposed_replacement_terms](proposed_replacement_terms.md) - Proposed replacement terms
 * [reason](reason.md) - Reason for the action
 * [reference_id](reference_id.md)
 * [reference_section_type](reference_section_type.md) - Type of section in the reference (e.g., 'ABSTRACT', 'METHODS', 'RESULTS', 'DISCUSSION')
 * [references](references.md)
 * [review](review.md) - Review of the gene
 * [statement](statement.md) - Concise statement describing an aspect of the gene
 * [summary](summary.md) - Summary of the review
 * [supported_by](supported_by.md)
 * [supporting_entities](supporting_entities.md) - IDs of the supporting entities
 * [supporting_text](supporting_text.md) - Supporting text from the publication. This should be exact substrings. Different substrings can be broken up by '...'s. These substrings will be checked against the actual text of the paper. If editorialization is necessary, put this in square brackets (this is not checked). For example, you can say '...[CFAP300 shows] transport within cilia is IFT dependent...'
 * [taxon](taxon.md)
 * [term](term.md) - Term to be annotated
     * [ExistingAnnotation➞term](ExistingAnnotation_term.md)
 * [title](title.md) - Title of the entity

### Enums

 * [ActionEnum](ActionEnum.md)
 * [EvidenceType](EvidenceType.md)
 * [GOBiologicalProcessEnum](GOBiologicalProcessEnum.md) - A biological process term in the GO ontology
 * [GOCellularLocationEnum](GOCellularLocationEnum.md) - A cellular location term in the GO ontology (excludes protein-containing complexes)
 * [GOMolecularActivityEnum](GOMolecularActivityEnum.md) - A molecular activity term in the GO ontology
 * [GOProteinContainingComplexEnum](GOProteinContainingComplexEnum.md) - A protein-containing complex term in the GO ontology
 * [GOTermEnum](GOTermEnum.md) - A term in the GO ontology

### Subsets


### Types


#### Built in

 * **Bool**
 * **Curie**
 * **Decimal**
 * **ElementIdentifier**
 * **NCName**
 * **NodeIdentifier**
 * **URI**
 * **URIorCURIE**
 * **XSDDate**
 * **XSDDateTime**
 * **XSDTime**
 * **float**
 * **int**
 * **str**

#### Defined

 * [Boolean](types/Boolean.md)  (**Bool**)  - A binary (true or false) value
 * [Curie](types/Curie.md)  (**Curie**)  - a compact URI
 * [Date](types/Date.md)  (**XSDDate**)  - a date (year, month and day) in an idealized calendar
 * [DateOrDatetime](types/DateOrDatetime.md)  (**str**)  - Either a date or a datetime
 * [Datetime](types/Datetime.md)  (**XSDDateTime**)  - The combination of a date and time
 * [Decimal](types/Decimal.md)  (**Decimal**)  - A real number with arbitrary precision that conforms to the xsd:decimal specification
 * [Double](types/Double.md)  (**float**)  - A real number that conforms to the xsd:double specification
 * [Float](types/Float.md)  (**float**)  - A real number that conforms to the xsd:float specification
 * [Integer](types/Integer.md)  (**int**)  - An integer
 * [Jsonpath](types/Jsonpath.md)  (**str**)  - A string encoding a JSON Path. The value of the string MUST conform to JSON Point syntax and SHOULD dereference to zero or more valid objects within the current instance document when encoded in tree form.
 * [Jsonpointer](types/Jsonpointer.md)  (**str**)  - A string encoding a JSON Pointer. The value of the string MUST conform to JSON Point syntax and SHOULD dereference to a valid object within the current instance document when encoded in tree form.
 * [Ncname](types/Ncname.md)  (**NCName**)  - Prefix part of CURIE
 * [Nodeidentifier](types/Nodeidentifier.md)  (**NodeIdentifier**)  - A URI, CURIE or BNODE that represents a node in a model.
 * [Objectidentifier](types/Objectidentifier.md)  (**ElementIdentifier**)  - A URI or CURIE that represents an object in the model.
 * [Sparqlpath](types/Sparqlpath.md)  (**str**)  - A string encoding a SPARQL Property Path. The value of the string MUST conform to SPARQL syntax and SHOULD dereference to zero or more valid objects within the current instance document when encoded as RDF.
 * [String](types/String.md)  (**str**)  - A character string
 * [Time](types/Time.md)  (**XSDTime**)  - A time object represents a (local) time of day, independent of any particular day
 * [Uri](types/Uri.md)  (**URI**)  - a complete URI
 * [Uriorcurie](types/Uriorcurie.md)  (**URIorCURIE**)  - a URI or a CURIE
