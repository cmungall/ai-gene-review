"""JSON exporter for gene review data."""

import json
import yaml
from pathlib import Path
from typing import Dict, List, Any, Union, Optional

from ai_gene_review.datamodel.gene_review_model import GeneReview, ExistingAnnotation


class JSONExporter:
    """Export gene review data to JSON format, flattening for linkml-browser."""

    def __init__(self):
        """Initialize the JSON exporter."""
        pass

    def export_from_file(self, file_path: Union[str, Path]) -> List[Dict[str, Any]]:
        """
        Export existing_annotations from a single gene review file to list of dicts.

        Args:
            file_path: Path to the gene review YAML file

        Returns:
            List of dictionaries with one dict per existing annotation
        """
        file_path = Path(file_path)

        with open(file_path, "r") as f:
            data = yaml.safe_load(f)

        # Parse into GeneReview model
        gene_review = GeneReview.model_validate(data)

        return self._flatten_existing_annotations(gene_review)

    def export_from_files(
        self, file_paths: List[Union[str, Path]]
    ) -> List[Dict[str, Any]]:
        """
        Export existing_annotations from multiple gene review files to list of dicts.

        Args:
            file_paths: List of paths to gene review YAML files

        Returns:
            Combined list with one dict per existing annotation from all files
        """
        all_annotations = []

        for file_path in file_paths:
            try:
                annotations = self.export_from_file(file_path)
                all_annotations.extend(annotations)
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue

        return all_annotations

    def _flatten_existing_annotations(
        self, gene_review: GeneReview
    ) -> List[Dict[str, Any]]:
        """
        Flatten existing_annotations from a GeneReview into list of dicts.

        Args:
            gene_review: The parsed GeneReview object

        Returns:
            List of flattened annotation dictionaries
        """
        if not gene_review.existing_annotations:
            return []

        annotations = []

        for annotation in gene_review.existing_annotations:
            flat_annotation = self._flatten_annotation(annotation, gene_review)
            annotations.append(flat_annotation)

        return annotations

    def _flatten_annotation(
        self, annotation: ExistingAnnotation, gene_review: GeneReview
    ) -> Dict[str, Any]:
        """
        Flatten a single existing annotation into a dictionary.

        Args:
            annotation: The ExistingAnnotation object
            gene_review: The parent GeneReview object

        Returns:
            Flattened dictionary with all relevant fields
        """
        row: Dict[str, Any] = {}

        # Create a unique ID for this annotation
        if gene_review.id and annotation.term and annotation.term.id:
            row["id"] = (
                f"{gene_review.id}_{annotation.term.id}_{annotation.evidence_type or 'NA'}"
            )
        else:
            row["id"] = f"{gene_review.id}_annotation_{id(annotation)}"

        # Parent object information
        row["protein_id"] = gene_review.id
        row["gene_symbol"] = gene_review.gene_symbol
        row["taxon_id"] = gene_review.taxon.id if gene_review.taxon else None
        row["taxon_label"] = gene_review.taxon.label if gene_review.taxon else None

        # Term information
        if annotation.term:
            row["term_id"] = annotation.term.id
            row["term_label"] = annotation.term.label
            row["term_description"] = annotation.term.description
            row["term_ontology"] = annotation.term.ontology
        else:
            row["term_id"] = None
            row["term_label"] = None
            row["term_description"] = None
            row["term_ontology"] = None

        # Evidence information
        row["evidence_type"] = annotation.evidence_type
        row["negated"] = annotation.negated or False

        # Original reference information
        if annotation.original_reference_id:
            ref_id = annotation.original_reference_id
            row["original_reference_id"] = ref_id

            # Find the matching reference in the parent's references list
            original_reference_title = self._find_reference_title(gene_review, ref_id)
            row["original_reference_title"] = original_reference_title
        else:
            row["original_reference_id"] = None
            row["original_reference_title"] = None

        # Supporting entities
        if annotation.supporting_entities:
            row["supporting_entities"] = "; ".join(annotation.supporting_entities)
        else:
            row["supporting_entities"] = None

        # Extensions (as string for simplicity)
        if annotation.extensions:
            extension_strs = []
            for ext in annotation.extensions:
                if ext.predicate and ext.term:
                    ext_str = f"{ext.predicate}({ext.term.id}:{ext.term.label})"
                    extension_strs.append(ext_str)
            row["extensions"] = "; ".join(extension_strs) if extension_strs else None
        else:
            row["extensions"] = None

        # Review information (flattened with dot notation for compatibility)
        if annotation.review:
            review = annotation.review
            row["review.summary"] = review.summary
            row["review.action"] = review.action if review.action else None
            row["review.reason"] = review.reason
            
            # Flatten supported_by field
            if review.supported_by:
                supporting_texts = []
                reference_ids = []
                for sup in review.supported_by:
                    if sup.supporting_text:
                        supporting_texts.append(sup.supporting_text)
                    if sup.reference_id:
                        reference_ids.append(sup.reference_id)
                row["review.supporting_text"] = " | ".join(supporting_texts) if supporting_texts else None
                row["review.supporting_reference_ids"] = "; ".join(reference_ids) if reference_ids else None
            else:
                row["review.supporting_text"] = None
                row["review.supporting_reference_ids"] = None

            # Proposed replacement terms
            if review.proposed_replacement_terms:
                replacements = []
                for term in review.proposed_replacement_terms:
                    replacements.append(f"{term.id}:{term.label}")
                row["review.proposed_replacement_terms"] = "; ".join(replacements)
            else:
                row["review.proposed_replacement_terms"] = None

            # Additional reference IDs
            if review.additional_reference_ids:
                row["review.additional_reference_ids"] = "; ".join(
                    review.additional_reference_ids
                )
            else:
                row["review.additional_reference_ids"] = None
        else:
            row["review.summary"] = None
            row["review.action"] = None
            row["review.reason"] = None
            row["review.supporting_text"] = None
            row["review.proposed_replacement_terms"] = None
            row["review.additional_reference_ids"] = None

        return row

    def _find_reference_title(
        self, gene_review: GeneReview, ref_id: str
    ) -> Optional[str]:
        """
        Find the title of a reference by its ID in the gene review's reference list.

        Args:
            gene_review: The GeneReview object
            ref_id: The reference ID to find

        Returns:
            The reference title if found, None otherwise
        """
        if not gene_review.references:
            return None

        for ref in gene_review.references:
            if ref.id == ref_id:
                return ref.title

        return None

    def export_to_json(
        self, file_paths: List[Union[str, Path]], output_path: Union[str, Path]
    ) -> None:
        """
        Export existing_annotations from multiple files to JSON.

        Args:
            file_paths: List of paths to gene review YAML files
            output_path: Path for the output JSON file
        """
        annotations = self.export_from_files(file_paths)

        output_path = Path(output_path)
        with open(output_path, "w") as f:
            json.dump(annotations, f, indent=2, default=str)
