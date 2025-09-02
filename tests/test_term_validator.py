"""Tests for ontology term validation."""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import yaml

from ai_gene_review.validation import (
    TermValidator,
    TermValidationResult,
    validate_gene_review,
    ValidationSeverity,
)


def test_term_format_validation():
    """Test validation of term ID formats."""
    validator = TermValidator()

    # Valid formats
    assert validator._is_valid_term_format("GO:0001234")
    assert validator._is_valid_term_format("HP:0000001")
    assert validator._is_valid_term_format("NCBITaxon:9606")
    assert validator._is_valid_term_format("UBERON:0001234")

    # Invalid formats
    assert not validator._is_valid_term_format("not_a_term")
    assert not validator._is_valid_term_format("12345")
    assert not validator._is_valid_term_format("GO")
    assert not validator._is_valid_term_format("")


def test_is_ontology_term():
    """Test detection of ontology term IDs."""
    validator = TermValidator()

    # Valid ontology terms
    assert validator._is_ontology_term("GO:0001234")
    assert validator._is_ontology_term("HP:0000001")
    assert validator._is_ontology_term("NCBITaxon:9606")

    # Not ontology terms
    assert not validator._is_ontology_term("not_a_term")
    assert not validator._is_ontology_term(123)
    assert not validator._is_ontology_term(None)
    assert not validator._is_ontology_term("")


def test_validate_term_with_mock():
    """Test term validation with mocked API responses."""
    validator = TermValidator()

    # Mock the API call
    with patch.object(validator, "_get_term_label") as mock_get_label:
        mock_get_label.return_value = "protein binding"

        # Valid term with correct label
        result = validator.validate_term("GO:0005515", "protein binding")
        assert result.is_valid
        assert result.correct_label == "protein binding"
        assert result.error_message is None

        # Valid term with wrong label
        result = validator.validate_term("GO:0005515", "wrong label")
        assert not result.is_valid
        assert result.correct_label == "protein binding"
        assert "Label mismatch" in result.error_message
        assert "expected 'protein binding'" in result.error_message

        # Valid term with no label provided
        result = validator.validate_term("GO:0005515")
        assert result.is_valid
        assert result.correct_label == "protein binding"
        assert result.provided_label is None


def test_validate_term_not_found():
    """Test validation when term is not found."""
    validator = TermValidator()

    with patch.object(validator, "_get_term_label") as mock_get_label:
        mock_get_label.return_value = None

        result = validator.validate_term("GO:9999999", "some label")
        assert not result.is_valid
        assert result.correct_label is None
        assert "not found in ontology" in result.error_message


def test_validate_terms_in_data():
    """Test recursive validation of terms in data structures."""
    validator = TermValidator()

    # Mock the API calls
    with patch.object(validator, "_get_term_label") as mock_get_label:

        def mock_label(term_id):
            labels = {
                "GO:0005515": "protein binding",
                "GO:0005737": "cytoplasm",
                "NCBITaxon:9606": "Homo sapiens",
            }
            return labels.get(term_id)

        mock_get_label.side_effect = mock_label

        data = {
            "taxon": {
                "id": "NCBITaxon:9606",
                "label": "Homo sapiens",  # Correct
            },
            "existing_annotations": [
                {
                    "term": {
                        "id": "GO:0005515",
                        "label": "protein binding",  # Correct
                    }
                },
                {
                    "term": {
                        "id": "GO:0005737",
                        "label": "wrong label",  # Wrong!
                    }
                },
            ],
        }

        results = validator.validate_terms_in_data(data)

        # Should find 3 terms
        assert len(results) == 3

        # Check validation results
        valid_results = [r for r in results if r.is_valid]
        invalid_results = [r for r in results if not r.is_valid]

        assert len(valid_results) == 2  # NCBITaxon and first GO term
        assert len(invalid_results) == 1  # Second GO term with wrong label

        # Check the invalid result
        invalid = invalid_results[0]
        assert invalid.term_id == "GO:0005737"
        assert invalid.provided_label == "wrong label"
        assert invalid.correct_label == "cytoplasm"


@pytest.mark.integration
def test_validate_term_real_api():
    """Test term validation with real OLS API.

    This test requires network access and is marked as integration.
    """
    validator = TermValidator()

    # Test a well-known GO term
    result = validator.validate_term("GO:0005515", "protein binding")
    assert result.is_valid
    assert result.correct_label == "protein binding"

    # Test with wrong label
    result = validator.validate_term("GO:0005515", "wrong label")
    assert not result.is_valid
    assert result.correct_label == "protein binding"

    # Test NCBITaxon term
    result = validator.validate_term("NCBITaxon:9606", "Homo sapiens")
    assert result.is_valid


def test_integration_with_main_validator():
    """Test that term validation is integrated into the main validator."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {
            "id": "NCBITaxon:9606",
            "label": "Wrong Species Name",  # Should be "Homo sapiens"
        },
        "description": "Test gene for validation",
        "existing_annotations": [
            {
                "term": {
                    "id": "GO:0005515",
                    "label": "wrong binding",  # Should be "protein binding"
                },
                "evidence_type": "IPI",
                "review": {"summary": "Test", "action": "ACCEPT"},
            }
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        # Mock the term validator
        with patch(
            "ai_gene_review.validation.validator.TermValidator"
        ) as MockTermValidator:
            mock_validator = Mock()
            MockTermValidator.return_value = mock_validator

            # Set up mock results
            mock_results = [
                TermValidationResult(
                    term_id="NCBITaxon:9606",
                    provided_label="Wrong Species Name",
                    correct_label="Homo sapiens",
                    is_valid=False,
                    error_message="Label mismatch for NCBITaxon:9606: got 'Wrong Species Name', expected 'Homo sapiens'",
                    path="taxon",
                ),
                TermValidationResult(
                    term_id="GO:0005515",
                    provided_label="wrong binding",
                    correct_label="protein binding",
                    is_valid=False,
                    error_message="Label mismatch for GO:0005515: got 'wrong binding', expected 'protein binding'",
                    path="existing_annotations[0].term",
                ),
            ]
            mock_validator.validate_terms_in_data.return_value = mock_results

            report = validate_gene_review(temp_path)

            # Should have errors for the mismatched labels
            term_errors = [
                issue
                for issue in report.issues
                if issue.severity == ValidationSeverity.ERROR
                and "Label mismatch" in issue.message
            ]

            assert len(term_errors) == 2
            assert not report.is_valid

    finally:
        temp_path.unlink()


def test_cache_behavior():
    """Test that caching works correctly."""
    validator = TermValidator(use_cache=True)

    # Pre-populate cache to test caching behavior
    validator._label_cache["GO:0005515"] = "protein binding"

    # First call should use cache
    label1 = validator._get_term_label("GO:0005515")
    assert label1 == "protein binding"

    # Second call should also use cache
    label2 = validator._get_term_label("GO:0005515")
    assert label2 == "protein binding"

    # Check cache contains the term
    assert "GO:0005515" in validator._label_cache
    assert validator._label_cache["GO:0005515"] == "protein binding"

    # Test with a new term not in cache - mock the adapter call
    with patch.object(validator, "_get_ontology_adapter") as mock_get_adapter:
        mock_adapter = Mock()
        mock_adapter.label.return_value = "cytoplasm"
        mock_get_adapter.return_value = mock_adapter

        # This should call the adapter since it's not cached
        label3 = validator._get_term_label("GO:0005737")
        assert label3 == "cytoplasm"
        assert mock_adapter.label.call_count == 1

        # Now it should be cached
        assert "GO:0005737" in validator._label_cache
        assert validator._label_cache["GO:0005737"] == "cytoplasm"
