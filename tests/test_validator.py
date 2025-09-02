"""Tests for the gene review validator with updated API."""

import tempfile
from pathlib import Path

import pytest
import yaml

from ai_gene_review.validation import (
    validate_gene_review,
    validate_multiple_files,
    ValidationSeverity,
)


def test_validate_valid_file():
    """Test validation of a valid gene review file."""
    # This uses the example file we copied to tests/input
    yaml_file = Path("tests/input/CFAP300-ai-review.yaml")

    if yaml_file.exists():
        report = validate_gene_review(yaml_file)
        assert report.is_valid, (
            f"Valid file should pass validation, but got errors: {report.issues}"
        )
        assert report.error_count == 0


def test_validate_invalid_file():
    """Test validation of an invalid gene review file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        # Write invalid YAML (missing required fields)
        invalid_data = {
            "id": "TEST123",
            # Missing required fields: gene_symbol
        }
        yaml.dump(invalid_data, f)
        temp_file = Path(f.name)

    try:
        report = validate_gene_review(temp_file)
        assert not report.is_valid, "Invalid file should fail validation"
        assert report.error_count > 0, "Should have validation errors"
    finally:
        temp_file.unlink()


def test_validate_nonexistent_file():
    """Test validation of a non-existent file."""
    report = validate_gene_review("nonexistent.yaml")
    assert not report.is_valid
    assert any("not found" in issue.message.lower() for issue in report.issues)


def test_validate_invalid_taxon_id():
    """Test validation with invalid taxon ID (should be string)."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        invalid_data = {
            "id": "Q456",
            "gene_symbol": "GENE2",
            "description": "Test gene",
            "taxon": {
                "id": 9606,  # Should be string like "NCBITaxon:9606"
                "label": "Homo sapiens",
            },
        }
        yaml.dump(invalid_data, f)
        temp_file = Path(f.name)

    try:
        report = validate_gene_review(temp_file)
        assert not report.is_valid, "Should fail with numeric taxon ID"
        assert any("not of type 'string'" in issue.message for issue in report.issues)
    finally:
        temp_file.unlink()


def test_validate_multiple_files():
    """Test batch validation of multiple files."""
    # Create some test files
    test_files = []

    # Create a valid file
    with tempfile.NamedTemporaryFile(mode="w", suffix="_valid.yaml", delete=False) as f:
        valid_data = {
            "id": "Q123",
            "gene_symbol": "GENE1",
            "description": "Test gene 1",
        }
        yaml.dump(valid_data, f)
        test_files.append(Path(f.name))

    # Create an invalid file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix="_invalid.yaml", delete=False
    ) as f:
        invalid_data = {
            "id": "Q456",
            # Missing gene_symbol
        }
        yaml.dump(invalid_data, f)
        test_files.append(Path(f.name))

    try:
        batch_report = validate_multiple_files(test_files)

        # Check that we got results for both files
        assert len(batch_report.reports) == 2

        # One should be valid, one invalid
        assert batch_report.valid_files == 1, "Should have one valid file"
        assert batch_report.invalid_files == 1, "Should have one invalid file"

    finally:
        for f in test_files:
            f.unlink()


def test_validation_summary():
    """Test the validation summary function."""
    # Create a batch report
    from ai_gene_review.validation import BatchValidationReport, ValidationReport

    batch = BatchValidationReport()

    # Add valid reports
    batch.reports.append(ValidationReport(file_path=Path("file1.yaml"), is_valid=True))
    batch.reports.append(ValidationReport(file_path=Path("file3.yaml"), is_valid=True))

    # Add invalid report
    invalid_report = ValidationReport(file_path=Path("file2.yaml"), is_valid=False)
    invalid_report.add_issue(
        ValidationSeverity.ERROR,
        "Missing required field: gene_symbol",
        path="gene_symbol",
    )
    batch.reports.append(invalid_report)

    summary = batch.summary()

    # Check that summary contains expected information
    assert "3" in summary  # 3 total files
    assert "Valid: 2" in summary  # 2 valid files
    assert "Invalid: 1" in summary  # 1 invalid file


@pytest.mark.parametrize(
    "gene_data,should_be_valid",
    [
        # Valid minimal structure
        ({"id": "Q123", "gene_symbol": "GENE1", "description": "Test gene"}, True),
        # Valid with taxon
        (
            {
                "id": "Q456",
                "gene_symbol": "GENE2",
                "description": "Test gene",
                "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
            },
            True,
        ),
        # Valid with references
        (
            {
                "id": "Q789",
                "gene_symbol": "GENE3",
                "description": "Test gene",
                "references": [{"id": "PMID:12345", "title": "Test paper"}],
            },
            True,
        ),
        # Invalid - missing gene_symbol
        ({"id": "Q111", "description": "Test gene"}, False),
        # Invalid - missing id
        ({"gene_symbol": "GENE1", "description": "Test gene"}, False),
    ],
)
def test_various_gene_review_structures(gene_data, should_be_valid):
    """Test validation of various gene review structures."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(gene_data, f)
        temp_file = Path(f.name)

    try:
        # Disable best practices to avoid issues with term/publication validation
        report = validate_gene_review(temp_file, check_best_practices=False)
        if should_be_valid:
            assert report.is_valid, f"Should be valid but got errors: {report.issues}"
        else:
            assert not report.is_valid, "Should be invalid but passed validation"
    finally:
        temp_file.unlink()
