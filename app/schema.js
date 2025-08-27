window.searchSchema = {
  "name": "GeneAnnotationReviewSchema",
  "description": "Schema for displaying gene annotation review data in linkml-browser",
  "title": "Gene Annotation Review Browser",
  "id": "https://example.org/gene-annotation-review",
  "searchableFields": [
    "gene_symbol",
    "protein_id",
    "term_label",
    "term_id"
  ],
  "searchPlaceholder": "Search annotations...",
  "facets": [
    {
      "field": "review.action",
      "label": "Review Action",
      "type": "string",
      "sortBy": "alphabetical"
    },
    {
      "field": "gene_symbol",
      "label": "Gene Symbol",
      "type": "string",
      "sortBy": "alphabetical"
    },
    {
      "field": "taxon_label",
      "label": "Species",
      "type": "string",
      "sortBy": "count"
    },
    {
      "field": "term_ontology",
      "label": "Ontology Aspect",
      "type": "string",
      "sortBy": "alphabetical"
    },
    {
      "field": "evidence_type",
      "label": "Evidence Code",
      "type": "string",
      "sortBy": "alphabetical"
    },
    
    {
      "field": "negated",
      "label": "Negated",
      "type": "boolean",
      "sortBy": "alphabetical"
    }
  ],
  "displayFields": [
    {
      "field": "gene_symbol",
      "label": "Gene",
      "type": "string"
    },
    {
      "field": "protein_id",
      "label": "Protein",
      "type": "string"
    },
    {
      "field": "term_label",
      "label": "GO Term",
      "type": "string"
    },
    {
      "field": "term_id",
      "label": "GO ID",
      "type": "curie"
    },
    {
      "field": "evidence_type",
      "label": "Evidence",
      "type": "string"
    },
    {
      "field": "original_reference_id",
      "label": "Ref",
      "type": "curie"
    },
    {
      "field": "original_reference_title",
      "label": "Title",
      "type": "string"
    },
    {
      "field": "review.summary",
      "label": "Summary",
      "type": "string"
    },
    {
      "field": "review.action",
      "label": "Action",
      "type": "string"
    },
    {
      "field": "review.reason",
      "label": "Reason",
      "type": "string"
    },
    {
      "field": "review.proposed_replacement_terms",
      "label": "Replacement Terms",
      "type": "string"
    }
  ]
};
window.dispatchEvent(new Event('searchDataReady'));
