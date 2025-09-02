"""Test validation of file: references."""

import tempfile
from pathlib import Path

import pytest
import yaml

from ai_gene_review.validation import validate_gene_review, ValidationSeverity


def test_valid_file_reference(tmp_path):
    """Test that file: references to existing files pass validation."""
    # Create a test file structure
    genes_dir = tmp_path / "genes"
    worm_dir = genes_dir / "worm"
    lrx1_dir = worm_dir / "lrx-1"
    bio_dir = lrx1_dir / "bioinformatics"
    bio_dir.mkdir(parents=True)

    # Create the referenced file
    results_file = bio_dir / "RESULTS.md"
    results_file.write_text("# Results\nSome analysis results here.")

    # Create test YAML data with file reference
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST1",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene with file reference",
        "references": [
            {
                "id": "file:worm/lrx-1/bioinformatics/RESULTS.md",
                "title": "Bioinformatics analysis results",
            }
        ],
    }

    yaml_file = tmp_path / "test.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(data, f)

    # Change to the tmp directory to simulate running from project root
    import os

    original_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)

        report = validate_gene_review(yaml_file)

        # Should not have errors about the file reference
        file_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "File reference" in issue.message
        ]
        assert len(file_errors) == 0, f"Unexpected file reference errors: {file_errors}"

    finally:
        os.chdir(original_cwd)


def test_invalid_file_reference_nonexistent(tmp_path):
    """Test that file: references to non-existent files are caught."""
    # Create a minimal genes directory structure
    genes_dir = tmp_path / "genes"
    genes_dir.mkdir()

    # Create test YAML data with reference to non-existent file
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST2",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene with invalid file reference",
        "references": [
            {
                "id": "file:worm/lrx-1/bioinformatics/NONEXISTENT.md",
                "title": "Non-existent file",
            }
        ],
    }

    yaml_file = tmp_path / "test.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(data, f)

    # Change to the tmp directory to simulate running from project root
    import os

    original_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)

        report = validate_gene_review(yaml_file)

        # Should have an error about the non-existent file
        file_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "non-existent file" in issue.message
        ]
        assert len(file_errors) == 1
        assert "worm/lrx-1/bioinformatics/NONEXISTENT.md" in file_errors[0].message

    finally:
        os.chdir(original_cwd)


def test_invalid_file_reference_directory(tmp_path):
    """Test that file: references to directories (not files) are caught."""
    # Create a test directory structure
    genes_dir = tmp_path / "genes"
    worm_dir = genes_dir / "worm"
    lrx1_dir = worm_dir / "lrx-1"
    bio_dir = lrx1_dir / "bioinformatics"
    bio_dir.mkdir(parents=True)

    # Create test YAML data with reference to a directory
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST3",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene with directory reference",
        "references": [
            {
                "id": "file:worm/lrx-1/bioinformatics",  # This is a directory!
                "title": "Bioinformatics directory",
            }
        ],
    }

    yaml_file = tmp_path / "test.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(data, f)

    # Change to the tmp directory to simulate running from project root
    import os

    original_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)

        report = validate_gene_review(yaml_file)

        # Should have an error about pointing to a directory
        file_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "directory, not a file" in issue.message
        ]
        assert len(file_errors) == 1
        assert "worm/lrx-1/bioinformatics" in file_errors[0].message

    finally:
        os.chdir(original_cwd)


def test_file_prefix_not_validated_as_ontology():
    """Test that file: prefix is excluded from ontology validation."""
    from ai_gene_review.validation import TermValidator

    validator = TermValidator()

    data = {
        "references": [{"id": "file:some/path/to/file.md", "title": "A file reference"}]
    }

    results = validator.validate_terms_in_data(data)

    # Should not try to validate file: as an ontology term
    assert len(results) == 0


def test_mixed_reference_types():
    """Test validation with mixed PMID and file references."""
    # Use actual CFAP300 structure for realistic test
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST4",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene with mixed references",
        "references": [
            {"id": "PMID:12345", "title": "A published paper"},
            {"id": "file:human/TEST4/analysis/results.md", "title": "Analysis results"},
            {"id": "DOI:10.1234/test", "doi": "10.1234/test"},
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        # This will fail the file reference check since we're not in a proper directory
        # structure, but it should still validate the YAML structure
        report = validate_gene_review(temp_path)

        # Check that different reference types are handled
        # PMID and DOI should not cause ontology validation errors
        ontology_errors = [
            issue
            for issue in report.issues
            if "Unrecognized identifier prefix" in issue.message
        ]
        # Should not have errors for PMID, DOI, or file prefixes
        for error in ontology_errors:
            assert "PMID" not in error.message
            assert "DOI" not in error.message
            assert "file" not in error.message

    finally:
        temp_path.unlink()


@pytest.mark.integration
def test_real_file_reference():
    """Test file reference validation with real gene directory structure."""
    # Check if we have the worm/lrx-1 directory
    lrx1_path = Path("genes/worm/lrx-1")
    if not lrx1_path.exists():
        pytest.skip("worm/lrx-1 gene directory not found")

    # Create a test results file
    bio_dir = lrx1_path / "bioinformatics"
    bio_dir.mkdir(exist_ok=True)
    results_file = bio_dir / "TEST_RESULTS.md"
    results_file.write_text("# Test Results\nThis is a test file for validation.")

    try:
        # Create a test YAML with valid file reference
        data = {
            "id": "Q12345",
            "gene_symbol": "LRX1_TEST",
            "taxon": {"id": "NCBITaxon:6239", "label": "Caenorhabditis elegans"},
            "description": "Test gene with real file reference",
            "references": [
                {
                    "id": "file:worm/lrx-1/bioinformatics/TEST_RESULTS.md",
                    "title": "Test bioinformatics results",
                }
            ],
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(data, f)
            temp_path = Path(f.name)

        try:
            report = validate_gene_review(temp_path)

            # Should not have errors about the file reference
            file_errors = [
                issue
                for issue in report.issues
                if issue.severity == ValidationSeverity.ERROR
                and "File reference" in issue.message
            ]
            assert len(file_errors) == 0

        finally:
            temp_path.unlink()

    finally:
        # Clean up test file
        if results_file.exists():
            results_file.unlink()
