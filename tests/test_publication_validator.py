"""Tests for publication validation."""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import yaml

from ai_gene_review.validation import (
    PublicationValidator,
    PublicationValidationResult,
    validate_yaml_file_publications,
    validate_gene_review,
)


def test_pmid_extraction():
    """Test that PMIDs are extracted correctly."""
    validator = PublicationValidator()

    # Test with mock to avoid actual validation
    with patch.object(validator, "_get_publication_title") as mock_get:
        mock_get.return_value = ("Test Title", False)

        result = validator.validate_publication("PMID:12345")
        assert result.pmid == "12345"

        result = validator.validate_publication("12345")
        assert result.pmid == "12345"

        result = validator.validate_publication("pmid:12345")
        assert result.pmid == "12345"


def test_title_matching():
    """Test title matching logic."""
    validator = PublicationValidator()

    # Exact match
    assert validator._titles_match("Test Article Title", "Test Article Title")

    # Case insensitive
    assert validator._titles_match("test article title", "TEST ARTICLE TITLE")

    # Punctuation differences
    assert validator._titles_match("Test Article: A Study", "Test Article A Study")

    # Extra whitespace
    assert validator._titles_match("Test  Article   Title", "Test Article Title")

    # Different titles
    assert not validator._titles_match("Test Article", "Different Article")


def test_substring_matching():
    """Test detection of truncated titles."""
    validator = PublicationValidator()

    # Truncated title (prefix)
    assert validator._is_substring_match(
        "C11orf70 Mutations Disrupting",
        "C11orf70 Mutations Disrupting the Intraflagellar Transport-Dependent Assembly",
    )

    # Too short to be considered valid
    assert not validator._is_substring_match(
        "C11orf70", "C11orf70 Mutations Disrupting the Intraflagellar Transport"
    )

    # Not a match
    assert not validator._is_substring_match(
        "Different Title", "C11orf70 Mutations Disrupting the Intraflagellar Transport"
    )


def test_validate_publication_with_cache(tmp_path):
    """Test validation using cached publication."""
    # Create a cached publication file
    cache_dir = tmp_path / "publications"
    cache_dir.mkdir()

    cached_pub = cache_dir / "PMID_12345.md"
    cached_pub.write_text("""---
pmid: '12345'
title: Cached Test Article
authors:
- Smith J
journal: Test Journal
year: '2024'
---

# Cached Test Article

**Authors:** Smith J
**Journal:** Test Journal (2024)

## Abstract

Test abstract.
""")

    validator = PublicationValidator(cache_dir=cache_dir)

    # Valid title
    result = validator.validate_publication("PMID:12345", "Cached Test Article")
    assert result.is_valid
    assert result.correct_title == "Cached Test Article"
    assert result.from_cache

    # Invalid title
    result = validator.validate_publication("PMID:12345", "Wrong Title")
    assert not result.is_valid
    assert result.correct_title == "Cached Test Article"
    assert "Title mismatch" in result.error_message

    # No title provided - just check existence
    result = validator.validate_publication("PMID:12345")
    assert result.is_valid
    assert result.correct_title == "Cached Test Article"


def test_validate_publication_with_mock_fetch():
    """Test validation with mocked PubMed fetch."""
    validator = PublicationValidator(auto_fetch=True)

    # Mock the fetch function
    with patch(
        "ai_gene_review.validation.publication_validator.fetch_pubmed_data"
    ) as mock_fetch:
        mock_pub = Mock()
        mock_pub.title = "Fetched Article Title"
        mock_fetch.return_value = mock_pub

        # Use a PMID that won't be cached (very high number)
        test_pmid = "PMID:999999999"

        # Valid title
        result = validator.validate_publication(test_pmid, "Fetched Article Title")
        assert result.is_valid
        assert result.correct_title == "Fetched Article Title"
        assert not result.from_cache  # Fetched, not cached

        # Invalid title
        result = validator.validate_publication(test_pmid, "Wrong Title")
        assert not result.is_valid
        assert result.correct_title == "Fetched Article Title"


def test_validate_publication_not_found():
    """Test validation when publication is not found."""
    validator = PublicationValidator(auto_fetch=False)

    result = validator.validate_publication("PMID:99999999", "Some Title")
    assert not result.is_valid
    assert result.correct_title is None
    assert "Could not verify" in result.error_message


def test_validate_publications_in_data():
    """Test recursive validation in data structures."""
    validator = PublicationValidator()

    # Mock the title fetching
    with patch.object(validator, "_get_publication_title") as mock_get:

        def mock_titles(pmid):
            titles = {
                "12345": ("First Article", True),
                "67890": ("Second Article", True),
                "11111": ("Third Article", True),
            }
            return titles.get(pmid, (None, False))

        mock_get.side_effect = mock_titles

        data = {
            "references": [
                {
                    "id": "PMID:12345",
                    "title": "First Article",  # Correct
                },
                {
                    "id": "PMID:67890",
                    "title": "Wrong Title",  # Wrong
                },
                {
                    "id": "DOI:10.1234/test",  # Not a PMID
                    "title": "Article with DOI",
                },
            ],
            "existing_annotations": [
                {
                    "original_reference_id": "PMID:11111"  # No title to check
                }
            ],
        }

        results = validator.validate_publications_in_data(data)

        # Should find 3 PMIDs
        assert len(results) == 3

        # Check individual results
        valid_results = [r for r in results if r.is_valid]
        invalid_results = [r for r in results if not r.is_valid]

        assert len(valid_results) == 2  # First ref and original_reference_id
        assert len(invalid_results) == 1  # Second ref with wrong title

        # Check the invalid one
        invalid = invalid_results[0]
        assert invalid.pmid == "67890"
        assert invalid.provided_title == "Wrong Title"
        assert invalid.correct_title == "Second Article"


def test_integration_with_main_validator(tmp_path):
    """Test that publication validation is integrated into main validator."""
    # Create a cached publication
    cache_dir = tmp_path / "publications"
    cache_dir.mkdir()

    cached_pub = cache_dir / "PMID_12345.md"
    cached_pub.write_text("""---
pmid: '12345'
title: Correct Publication Title
authors:
- Author A
journal: Test Journal
year: '2024'
---

# Correct Publication Title
""")

    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene",
        "references": [
            {
                "id": "PMID:12345",
                "title": "Wrong Title Here",  # This should trigger an error
            }
        ],
    }

    yaml_file = tmp_path / "test.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(data, f)

    # Mock the PublicationValidator to use our temp cache
    with patch(
        "ai_gene_review.validation.validator.PublicationValidator"
    ) as MockPubValidator:
        mock_validator = PublicationValidator(cache_dir=cache_dir)
        MockPubValidator.return_value = mock_validator

        report = validate_gene_review(yaml_file)

        # Should have an error for the title mismatch
        pub_errors = [
            issue for issue in report.issues if "Title mismatch" in issue.message
        ]

        assert len(pub_errors) > 0
        assert not report.is_valid


@pytest.mark.integration
def test_validate_real_publication():
    """Test validation with real cached CFAP300 publications.

    This test assumes the CFAP300 publications have been cached.
    """
    # Check if publications are cached
    cache_dir = Path("publications")
    if not cache_dir.exists() or not (cache_dir / "PMID_29727692.md").exists():
        pytest.skip("CFAP300 publications not cached")

    validator = PublicationValidator(cache_dir=cache_dir)

    # Test with correct title (abbreviated but should match as substring)
    result = validator.validate_publication(
        "PMID:29727692",
        "C11orf70 Mutations Disrupting the Intraflagellar Transport-Dependent Assembly of Multiple Axonemal Dyneins Cause Primary Ciliary Dyskinesia",
    )
    assert result.is_valid
    assert result.from_cache

    # Test with wrong title
    result = validator.validate_publication("PMID:29727692", "Completely Wrong Title")
    assert not result.is_valid
    assert "Title mismatch" in result.error_message


def test_validate_yaml_file_publications(tmp_path):
    """Test the convenience function for validating a YAML file."""
    # Create test YAML file
    data = {"references": [{"id": "PMID:12345", "title": "Test Article"}]}

    yaml_file = tmp_path / "test.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(data, f)

    # Mock the validator
    with patch(
        "ai_gene_review.validation.publication_validator.PublicationValidator"
    ) as MockValidator:
        mock_validator = Mock()
        mock_validator.validate_publications_in_data.return_value = [
            PublicationValidationResult(
                pmid="12345",
                provided_title="Test Article",
                correct_title="Test Article",
                is_valid=True,
                from_cache=False,
            )
        ]
        MockValidator.return_value = mock_validator

        all_valid, results = validate_yaml_file_publications(yaml_file)

        assert all_valid
        assert len(results) == 1
        assert results[0].is_valid
