-- # Class: GeneReview Description: Complete review for a gene
--     * Slot: id
--     * Slot: gene_symbol Description: Symbol of the gene
--     * Slot: description Description: Description of the entity
--     * Slot: taxon_id
-- # Class: Term Description: A term in a specific ontology
--     * Slot: id
--     * Slot: label Description: Human readable name of the entity
--     * Slot: description Description: Description of the entity
--     * Slot: ontology Description: Ontology of the term. E.g `go`, `cl`, `hp`
--     * Slot: Review_id Description: Autocreated FK slot
--     * Slot: CoreFunction_id Description: Autocreated FK slot
-- # Class: Reference Description: A reference is a published text  that describes a finding or a method. References might be formal publications (where the ID is a PMID), or for methods, a GO_REF. Additionally, a reference to a local ad-hoc analysis or review can be made by using the `file:` prefix.
--     * Slot: id
--     * Slot: title Description: Title of the entity
--     * Slot: is_invalid Description: Whether the reference is invalid (e.g., retracted or replaced)
--     * Slot: GeneReview_id Description: Autocreated FK slot
-- # Class: Finding Description: A finding is a statement about a gene, which is supported by a reference. Similar to "comments" in uniprot
--     * Slot: id
--     * Slot: statement Description: Concise statement describing an aspect of the gene
--     * Slot: supporting_text Description: Supporting text from the publication. This should be exact substrings. Different substrings can be broken up by '...'s. These substrings will be checked against the actual text of the paper. If editorialization is necessary, put this in square brackets (this is not checked). For example, you can say '...[CFAP300 shows] transport within cilia is IFT dependent...'
--     * Slot: reference_section_type Description: Type of section in the reference (e.g., 'ABSTRACT', 'METHODS', 'RESULTS', 'DISCUSSION')
-- # Class: SupportingTextInReference Description: A supporting text in a reference.
--     * Slot: id
--     * Slot: reference_id
--     * Slot: supporting_text Description: Supporting text from the publication. This should be exact substrings. Different substrings can be broken up by '...'s. These substrings will be checked against the actual text of the paper. If editorialization is necessary, put this in square brackets (this is not checked). For example, you can say '...[CFAP300 shows] transport within cilia is IFT dependent...'
--     * Slot: reference_section_type Description: Type of section in the reference (e.g., 'ABSTRACT', 'METHODS', 'RESULTS', 'DISCUSSION')
-- # Class: ExistingAnnotation Description: An existing annotation from the GO database, plus a review of the annotation.
--     * Slot: id
--     * Slot: negated Description: Whether the term is negated
--     * Slot: evidence_type Description: Evidence code (e.g., IDA, IBA, ISS, TAS)
--     * Slot: original_reference_id Description: ID of the original reference
--     * Slot: term_id Description: Term to be annotated
--     * Slot: review_id Description: Review of the gene
-- # Class: Review Description: A review of an existing annotation.
--     * Slot: id
--     * Slot: summary Description: Summary of the review
--     * Slot: action Description: Action to be taken
--     * Slot: reason Description: Reason for the action
-- # Class: CoreFunction Description: A core function is a GO-CAM-like annotation of the core evolved functions of a gene. This is a synthesis of the reviewed core annotations, brought together into a unified GO-CAM-like representation.
--     * Slot: id
--     * Slot: description Description: Description of the core function
--     * Slot: molecular_function_id
--     * Slot: in_complex_id
-- # Class: AnnotationExtension
--     * Slot: id
--     * Slot: predicate Description: Predicate of the extension
--     * Slot: term_id Description: Term to be annotated
-- # Class: ProposedOntologyTerm Description: A proposed new ontology term that should exist but doesn't currently
--     * Slot: id
--     * Slot: proposed_name Description: Proposed name for the new term
--     * Slot: proposed_definition Description: Proposed definition for the new term
--     * Slot: justification Description: Justification for why this term is needed
-- # Class: GeneReview_aliases
--     * Slot: GeneReview_id Description: Autocreated FK slot
--     * Slot: aliases
-- # Class: GeneReview_existing_annotations
--     * Slot: GeneReview_id Description: Autocreated FK slot
--     * Slot: existing_annotations_id
-- # Class: GeneReview_core_functions
--     * Slot: GeneReview_id Description: Autocreated FK slot
--     * Slot: core_functions_id
-- # Class: GeneReview_proposed_new_terms
--     * Slot: GeneReview_id Description: Autocreated FK slot
--     * Slot: proposed_new_terms_id Description: Proposed new ontology terms that should exist but don't
-- # Class: Reference_findings
--     * Slot: Reference_id Description: Autocreated FK slot
--     * Slot: findings_id
-- # Class: ExistingAnnotation_extensions
--     * Slot: ExistingAnnotation_id Description: Autocreated FK slot
--     * Slot: extensions_id
-- # Class: ExistingAnnotation_supporting_entities
--     * Slot: ExistingAnnotation_id Description: Autocreated FK slot
--     * Slot: supporting_entities Description: IDs of the supporting entities
-- # Class: Review_additional_reference_ids
--     * Slot: Review_id Description: Autocreated FK slot
--     * Slot: additional_reference_ids_id Description: IDs of the references
-- # Class: Review_supported_by
--     * Slot: Review_id Description: Autocreated FK slot
--     * Slot: supported_by_id
-- # Class: CoreFunction_supported_by
--     * Slot: CoreFunction_id Description: Autocreated FK slot
--     * Slot: supported_by_id

CREATE TABLE "Term" (
	id TEXT NOT NULL,
	label TEXT NOT NULL,
	description TEXT,
	ontology TEXT,
	"Review_id" INTEGER,
	"CoreFunction_id" INTEGER,
	PRIMARY KEY (id),
	FOREIGN KEY("Review_id") REFERENCES "Review" (id),
	FOREIGN KEY("CoreFunction_id") REFERENCES "CoreFunction" (id)
);CREATE INDEX "ix_Term_id" ON "Term" (id);
CREATE TABLE "Finding" (
	id INTEGER NOT NULL,
	statement TEXT,
	supporting_text TEXT,
	reference_section_type TEXT,
	PRIMARY KEY (id)
);CREATE INDEX "ix_Finding_id" ON "Finding" (id);
CREATE TABLE "Review" (
	id INTEGER NOT NULL,
	summary TEXT,
	action VARCHAR(22),
	reason TEXT,
	PRIMARY KEY (id)
);CREATE INDEX "ix_Review_id" ON "Review" (id);
CREATE TABLE "CoreFunction" (
	id INTEGER NOT NULL,
	description TEXT,
	molecular_function_id TEXT,
	in_complex_id TEXT,
	PRIMARY KEY (id),
	FOREIGN KEY(molecular_function_id) REFERENCES "Term" (id),
	FOREIGN KEY(in_complex_id) REFERENCES "Term" (id)
);CREATE INDEX "ix_CoreFunction_id" ON "CoreFunction" (id);
CREATE TABLE "ProposedOntologyTerm" (
	id INTEGER NOT NULL,
	proposed_name TEXT NOT NULL,
	proposed_definition TEXT NOT NULL,
	justification TEXT,
	PRIMARY KEY (id)
);CREATE INDEX "ix_ProposedOntologyTerm_id" ON "ProposedOntologyTerm" (id);
CREATE TABLE "GeneReview" (
	id TEXT NOT NULL,
	gene_symbol TEXT NOT NULL,
	description TEXT,
	taxon_id TEXT,
	PRIMARY KEY (id),
	FOREIGN KEY(taxon_id) REFERENCES "Term" (id)
);CREATE INDEX "ix_GeneReview_id" ON "GeneReview" (id);
CREATE TABLE "AnnotationExtension" (
	id INTEGER NOT NULL,
	predicate TEXT,
	term_id TEXT,
	PRIMARY KEY (id),
	FOREIGN KEY(term_id) REFERENCES "Term" (id)
);CREATE INDEX "ix_AnnotationExtension_id" ON "AnnotationExtension" (id);
CREATE TABLE "Reference" (
	id TEXT NOT NULL,
	title TEXT NOT NULL,
	is_invalid BOOLEAN,
	"GeneReview_id" TEXT,
	PRIMARY KEY (id),
	FOREIGN KEY("GeneReview_id") REFERENCES "GeneReview" (id)
);CREATE INDEX "ix_Reference_id" ON "Reference" (id);
CREATE TABLE "GeneReview_aliases" (
	"GeneReview_id" TEXT,
	aliases TEXT,
	PRIMARY KEY ("GeneReview_id", aliases),
	FOREIGN KEY("GeneReview_id") REFERENCES "GeneReview" (id)
);CREATE INDEX "ix_GeneReview_aliases_aliases" ON "GeneReview_aliases" (aliases);CREATE INDEX "ix_GeneReview_aliases_GeneReview_id" ON "GeneReview_aliases" ("GeneReview_id");
CREATE TABLE "GeneReview_core_functions" (
	"GeneReview_id" TEXT,
	core_functions_id INTEGER,
	PRIMARY KEY ("GeneReview_id", core_functions_id),
	FOREIGN KEY("GeneReview_id") REFERENCES "GeneReview" (id),
	FOREIGN KEY(core_functions_id) REFERENCES "CoreFunction" (id)
);CREATE INDEX "ix_GeneReview_core_functions_core_functions_id" ON "GeneReview_core_functions" (core_functions_id);CREATE INDEX "ix_GeneReview_core_functions_GeneReview_id" ON "GeneReview_core_functions" ("GeneReview_id");
CREATE TABLE "GeneReview_proposed_new_terms" (
	"GeneReview_id" TEXT,
	proposed_new_terms_id INTEGER,
	PRIMARY KEY ("GeneReview_id", proposed_new_terms_id),
	FOREIGN KEY("GeneReview_id") REFERENCES "GeneReview" (id),
	FOREIGN KEY(proposed_new_terms_id) REFERENCES "ProposedOntologyTerm" (id)
);CREATE INDEX "ix_GeneReview_proposed_new_terms_proposed_new_terms_id" ON "GeneReview_proposed_new_terms" (proposed_new_terms_id);CREATE INDEX "ix_GeneReview_proposed_new_terms_GeneReview_id" ON "GeneReview_proposed_new_terms" ("GeneReview_id");
CREATE TABLE "SupportingTextInReference" (
	id INTEGER NOT NULL,
	reference_id TEXT,
	supporting_text TEXT,
	reference_section_type TEXT,
	PRIMARY KEY (id),
	FOREIGN KEY(reference_id) REFERENCES "Reference" (id)
);CREATE INDEX "ix_SupportingTextInReference_id" ON "SupportingTextInReference" (id);
CREATE TABLE "ExistingAnnotation" (
	id INTEGER NOT NULL,
	negated BOOLEAN,
	evidence_type TEXT,
	original_reference_id TEXT,
	term_id TEXT,
	review_id INTEGER,
	PRIMARY KEY (id),
	FOREIGN KEY(original_reference_id) REFERENCES "Reference" (id),
	FOREIGN KEY(term_id) REFERENCES "Term" (id),
	FOREIGN KEY(review_id) REFERENCES "Review" (id)
);CREATE INDEX "ix_ExistingAnnotation_id" ON "ExistingAnnotation" (id);
CREATE TABLE "Reference_findings" (
	"Reference_id" TEXT,
	findings_id INTEGER,
	PRIMARY KEY ("Reference_id", findings_id),
	FOREIGN KEY("Reference_id") REFERENCES "Reference" (id),
	FOREIGN KEY(findings_id) REFERENCES "Finding" (id)
);CREATE INDEX "ix_Reference_findings_findings_id" ON "Reference_findings" (findings_id);CREATE INDEX "ix_Reference_findings_Reference_id" ON "Reference_findings" ("Reference_id");
CREATE TABLE "Review_additional_reference_ids" (
	"Review_id" INTEGER,
	additional_reference_ids_id TEXT,
	PRIMARY KEY ("Review_id", additional_reference_ids_id),
	FOREIGN KEY("Review_id") REFERENCES "Review" (id),
	FOREIGN KEY(additional_reference_ids_id) REFERENCES "Reference" (id)
);CREATE INDEX "ix_Review_additional_reference_ids_Review_id" ON "Review_additional_reference_ids" ("Review_id");CREATE INDEX "ix_Review_additional_reference_ids_additional_reference_ids_id" ON "Review_additional_reference_ids" (additional_reference_ids_id);
CREATE TABLE "GeneReview_existing_annotations" (
	"GeneReview_id" TEXT,
	existing_annotations_id INTEGER,
	PRIMARY KEY ("GeneReview_id", existing_annotations_id),
	FOREIGN KEY("GeneReview_id") REFERENCES "GeneReview" (id),
	FOREIGN KEY(existing_annotations_id) REFERENCES "ExistingAnnotation" (id)
);CREATE INDEX "ix_GeneReview_existing_annotations_existing_annotations_id" ON "GeneReview_existing_annotations" (existing_annotations_id);CREATE INDEX "ix_GeneReview_existing_annotations_GeneReview_id" ON "GeneReview_existing_annotations" ("GeneReview_id");
CREATE TABLE "ExistingAnnotation_extensions" (
	"ExistingAnnotation_id" INTEGER,
	extensions_id INTEGER,
	PRIMARY KEY ("ExistingAnnotation_id", extensions_id),
	FOREIGN KEY("ExistingAnnotation_id") REFERENCES "ExistingAnnotation" (id),
	FOREIGN KEY(extensions_id) REFERENCES "AnnotationExtension" (id)
);CREATE INDEX "ix_ExistingAnnotation_extensions_extensions_id" ON "ExistingAnnotation_extensions" (extensions_id);CREATE INDEX "ix_ExistingAnnotation_extensions_ExistingAnnotation_id" ON "ExistingAnnotation_extensions" ("ExistingAnnotation_id");
CREATE TABLE "ExistingAnnotation_supporting_entities" (
	"ExistingAnnotation_id" INTEGER,
	supporting_entities TEXT,
	PRIMARY KEY ("ExistingAnnotation_id", supporting_entities),
	FOREIGN KEY("ExistingAnnotation_id") REFERENCES "ExistingAnnotation" (id)
);CREATE INDEX "ix_ExistingAnnotation_supporting_entities_ExistingAnnotation_id" ON "ExistingAnnotation_supporting_entities" ("ExistingAnnotation_id");CREATE INDEX "ix_ExistingAnnotation_supporting_entities_supporting_entities" ON "ExistingAnnotation_supporting_entities" (supporting_entities);
CREATE TABLE "Review_supported_by" (
	"Review_id" INTEGER,
	supported_by_id INTEGER,
	PRIMARY KEY ("Review_id", supported_by_id),
	FOREIGN KEY("Review_id") REFERENCES "Review" (id),
	FOREIGN KEY(supported_by_id) REFERENCES "SupportingTextInReference" (id)
);CREATE INDEX "ix_Review_supported_by_Review_id" ON "Review_supported_by" ("Review_id");CREATE INDEX "ix_Review_supported_by_supported_by_id" ON "Review_supported_by" (supported_by_id);
CREATE TABLE "CoreFunction_supported_by" (
	"CoreFunction_id" INTEGER,
	supported_by_id INTEGER,
	PRIMARY KEY ("CoreFunction_id", supported_by_id),
	FOREIGN KEY("CoreFunction_id") REFERENCES "CoreFunction" (id),
	FOREIGN KEY(supported_by_id) REFERENCES "SupportingTextInReference" (id)
);CREATE INDEX "ix_CoreFunction_supported_by_supported_by_id" ON "CoreFunction_supported_by" (supported_by_id);CREATE INDEX "ix_CoreFunction_supported_by_CoreFunction_id" ON "CoreFunction_supported_by" ("CoreFunction_id");
