"""Test validation of GO_REF references in annotations."""

from pathlib import Path
import tempfile
import yaml

from ai_gene_review.validation import validate_gene_review, ValidationSeverity


def test_valid_goref_reference():
    """Test that annotations with valid GO_REF reference IDs pass validation."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST_GOREF",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene for GO_REF validation",
        "references": [
            {"id": "PMID:12345", "title": "Test paper 1"},
            {
                "id": "GO_REF:0000033",
                "title": "Annotation inferred from sequence orthology",
            },
            {
                "id": "GO_REF:0000002",
                "title": "Gene Ontology annotation through association of InterPro records with GO terms.",
            },
        ],
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IEA",
                "original_reference_id": "GO_REF:0000033",
                "review": {"summary": "Valid GO_REF reference", "action": "ACCEPT"},
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "IEA",
                "original_reference_id": "GO_REF:0000002",
                "review": {
                    "summary": "Another valid GO_REF reference",
                    "action": "ACCEPT",
                },
            },
            {
                "term": {"id": "GO:0005634", "label": "nucleus"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:12345",
                "review": {"summary": "Regular PMID reference", "action": "ACCEPT"},
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        # Disable best practices to avoid term validation errors
        report = validate_gene_review(temp_path, check_best_practices=False)

        # Should not have any errors about reference IDs
        ref_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "non-existent reference ID" in issue.message
        ]
        assert len(ref_errors) == 0, f"Unexpected reference errors: {ref_errors}"
        assert report.is_valid
    finally:
        temp_path.unlink()


def test_invalid_goref_reference():
    """Test that annotations with invalid GO_REF reference IDs are caught."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST_GOREF_BAD",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene with invalid GO_REF reference",
        "references": [
            {
                "id": "GO_REF:0000033",
                "title": "Annotation inferred from sequence orthology",
            }
        ],
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IEA",
                "original_reference_id": "GO_REF:0000033",
                "review": {"summary": "Valid reference", "action": "ACCEPT"},
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "IEA",
                "original_reference_id": "GO_REF:9999999",  # This GO_REF doesn't exist in references
                "review": {"summary": "Invalid GO_REF reference", "action": "ACCEPT"},
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        # Enable best practices to check references
        report = validate_gene_review(temp_path, check_best_practices=True)

        # Should have an error about the invalid GO_REF reference ID
        ref_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "GO_REF:9999999" in issue.message
        ]
        assert len(ref_errors) == 1, (
            f"Expected 1 reference error, got {len(ref_errors)}: {ref_errors}"
        )
        assert "non-existent reference ID" in ref_errors[0].message
    finally:
        temp_path.unlink()


def test_mixed_reference_types():
    """Test that both PMID and GO_REF references can coexist."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST_MIXED",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene with mixed reference types",
        "references": [
            {"id": "PMID:12345", "title": "A research paper"},
            {
                "id": "GO_REF:0000033",
                "title": "Annotation inferred from sequence orthology",
            },
            {"id": "PMID:67890", "title": "Another paper"},
            {
                "id": "GO_REF:0000120",
                "title": "Combined Automated Annotation using Multiple IEA Methods",
            },
        ],
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:12345",
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "IEA",
                "original_reference_id": "GO_REF:0000033",
            },
            {
                "term": {"id": "GO:0005634", "label": "nucleus"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:67890",
            },
            {
                "term": {"id": "GO:0008150", "label": "biological_process"},
                "evidence_type": "IEA",
                "original_reference_id": "GO_REF:0000120",
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        # Disable best practices to avoid term validation errors
        report = validate_gene_review(temp_path, check_best_practices=False)

        # Should not have any errors about reference IDs
        ref_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "non-existent reference ID" in issue.message
        ]
        assert len(ref_errors) == 0, f"Unexpected reference errors: {ref_errors}"
        assert report.is_valid
    finally:
        temp_path.unlink()
