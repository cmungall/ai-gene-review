-- # Class: GeneReview Description: Complete review for a gene
--     * Slot: id
--     * Slot: gene_symbol Description: Symbol of the gene
--     * Slot: description Description: Description of the entity
--     * Slot: taxon_id
-- # Class: Term
--     * Slot: id
--     * Slot: label Description: Human readable name of the entity
--     * Slot: Review_id Description: Autocreated FK slot
--     * Slot: CoreFunction_id Description: Autocreated FK slot
-- # Class: Reference
--     * Slot: id
--     * Slot: title Description: Title of the entity
--     * Slot: GeneReview_id Description: Autocreated FK slot
-- # Class: Finding
--     * Slot: id
--     * Slot: statement Description: Concise statement describing an aspect of the gene
--     * Slot: supporting_text Description: Supporting text from the publication. Can be interleaved text if there are multiple lines of evidence. Can also be annotated
-- # Class: ExistingAnnotation
--     * Slot: id
--     * Slot: negated Description: Whether the term is negated
--     * Slot: evidence_type Description: Evidence code (e.g., IDA, IBA, ISS, TAS)
--     * Slot: original_reference_id Description: ID of the original reference
--     * Slot: term_id Description: Term to be annotated
--     * Slot: review_id Description: Review of the gene
-- # Class: Review
--     * Slot: id
--     * Slot: summary Description: Summary of the review
--     * Slot: action Description: Action to be taken
--     * Slot: reason Description: Reason for the action
-- # Class: CoreFunction
--     * Slot: id
--     * Slot: description Description: Description of the core function
--     * Slot: molecular_function_id
--     * Slot: in_complex_id
-- # Class: AnnotationExtension
--     * Slot: id
--     * Slot: predicate Description: Predicate of the extension
--     * Slot: term_id Description: Term to be annotated
-- # Class: GeneReview_existing_annotations
--     * Slot: GeneReview_id Description: Autocreated FK slot
--     * Slot: existing_annotations_id
-- # Class: GeneReview_core_functions
--     * Slot: GeneReview_id Description: Autocreated FK slot
--     * Slot: core_functions_id
-- # Class: Reference_findings
--     * Slot: Reference_id Description: Autocreated FK slot
--     * Slot: findings_id
-- # Class: ExistingAnnotation_extensions
--     * Slot: ExistingAnnotation_id Description: Autocreated FK slot
--     * Slot: extensions_id

CREATE TABLE "Term" (
	id TEXT NOT NULL,
	label TEXT,
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
CREATE TABLE "GeneReview" (
	id TEXT NOT NULL,
	gene_symbol TEXT,
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
	title TEXT,
	"GeneReview_id" TEXT,
	PRIMARY KEY (id),
	FOREIGN KEY("GeneReview_id") REFERENCES "GeneReview" (id)
);CREATE INDEX "ix_Reference_id" ON "Reference" (id);
CREATE TABLE "GeneReview_core_functions" (
	"GeneReview_id" TEXT,
	core_functions_id INTEGER,
	PRIMARY KEY ("GeneReview_id", core_functions_id),
	FOREIGN KEY("GeneReview_id") REFERENCES "GeneReview" (id),
	FOREIGN KEY(core_functions_id) REFERENCES "CoreFunction" (id)
);CREATE INDEX "ix_GeneReview_core_functions_core_functions_id" ON "GeneReview_core_functions" (core_functions_id);CREATE INDEX "ix_GeneReview_core_functions_GeneReview_id" ON "GeneReview_core_functions" ("GeneReview_id");
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
);CREATE INDEX "ix_ExistingAnnotation_extensions_ExistingAnnotation_id" ON "ExistingAnnotation_extensions" ("ExistingAnnotation_id");CREATE INDEX "ix_ExistingAnnotation_extensions_extensions_id" ON "ExistingAnnotation_extensions" (extensions_id);
