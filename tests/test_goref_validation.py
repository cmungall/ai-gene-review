"""Tests for GO_REF reference validation and title fetching."""

from pathlib import Path
from unittest.mock import Mock, patch

import yaml

from ai_gene_review.validation.goa_validator import GOAValidator
from ai_gene_review.validation.term_validator import TermValidator


def test_goref_allowed_in_term_validator():
    """Test that GO_REF references are properly excluded from term validation."""
    validator = TermValidator()

    # GO_REF should be in excluded prefixes
    assert "GO_REF" in validator.EXCLUDED_PREFIXES

    # Test data with GO_REF reference
    data = {
        "references": [
            {"id": "GO_REF:0000033", "title": "Some title"},
            {"id": "PMID:12345", "title": "Another title"},
        ]
    }

    # Validate - should not report errors for GO_REF
    results = validator.validate_terms_in_data(data)

    # Should have no validation results since both are excluded prefixes
    assert len(results) == 0


def test_goref_validation_in_yaml_file():
    """Test that YAML files with GO_REF references pass validation."""
    test_file = Path("tests/input/test_goref.yaml")

    # Load and validate the file
    with open(test_file) as f:
        data = yaml.safe_load(f)

    validator = TermValidator()
    results = validator.validate_terms_in_data(data)

    # Check that no errors are reported for GO_REF references
    invalid_results = [r for r in results if not r.is_valid]
    goref_errors = [r for r in invalid_results if "GO_REF" in r.term_id]

    assert len(goref_errors) == 0, (
        f"Unexpected GO_REF validation errors: {goref_errors}"
    )


def test_get_goref_title_with_mock():
    """Test fetching GO_REF titles from GO site."""
    validator = GOAValidator()

    # Mock the GO refs data
    mock_goref_data = [
        {
            "id": "GO_REF:0000033",
            "title": "Annotation inferred from sequence orthology",
            "authors": ["GO Consortium"],
        },
        {
            "id": "GO_REF:0000002",
            "title": "Gene Ontology annotation through association of InterPro records with GO terms.",
            "authors": ["InterPro2GO"],
        },
    ]

    with patch("requests.get") as mock_get:
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = yaml.dump(mock_goref_data)
        mock_get.return_value = mock_response

        # Test fetching titles
        cache = {}

        title1 = validator._get_go_ref_title("GO_REF:0000033", cache)
        assert title1 == "Annotation inferred from sequence orthology"
        assert "GO_REF:0000033" in cache

        title2 = validator._get_go_ref_title("GO_REF:0000002", cache)
        assert (
            title2
            == "Gene Ontology annotation through association of InterPro records with GO terms."
        )
        assert "GO_REF:0000002" in cache

        # Test with unknown GO_REF
        title3 = validator._get_go_ref_title("GO_REF:9999999", cache)
        assert title3 == "TODO: Fetch title"
        assert "GO_REF:9999999" in cache


def test_get_goref_title_network_failure():
    """Test GO_REF title fetching handles network failures gracefully."""
    validator = GOAValidator()

    with patch("requests.get") as mock_get:
        # Simulate network failure
        mock_get.side_effect = Exception("Network error")

        cache = {}
        title = validator._get_go_ref_title("GO_REF:0000033", cache)

        # Should return placeholder on failure
        assert title == "TODO: Fetch title"
        assert cache["GO_REF:0000033"] == "TODO: Fetch title"


def test_seed_with_goref_references():
    """Test seeding handles both PMID and GO_REF references properly."""
    validator = GOAValidator()

    # Create test YAML data
    test_data = {"gene_symbol": "TEST1", "existing_annotations": []}

    # Create temporary test file with proper naming
    import tempfile

    temp_dir = Path(tempfile.mkdtemp())
    test_file = temp_dir / "test-ai-review.yaml"
    with open(test_file, "w") as f:
        yaml.dump(test_data, f)

    # Create a dummy GOA file
    goa_file = temp_dir / "test-goa.tsv"
    goa_file.touch()  # Just create empty file since we mock parse_goa_file

    try:
        # Mock GOA annotations with mixed references
        mock_annotations = [
            Mock(
                go_id="GO:0005515",
                go_term="protein binding",
                evidence_code="IEA",
                reference="GO_REF:0000033",
            ),
            Mock(
                go_id="GO:0005634",
                go_term="nucleus",
                evidence_code="IEA",
                reference="PMID:12345678",
            ),
        ]

        with patch.object(validator, "parse_goa_file", return_value=mock_annotations):
            with patch.object(validator, "_get_go_ref_title") as mock_goref:
                mock_goref.return_value = "Test GO_REF title"

                with patch(
                    "ai_gene_review.etl.publication.get_cached_title"
                ) as mock_pmid:
                    mock_pmid.return_value = "Test PMID title"

                    # Seed with title fetching - provide GOA file explicitly
                    added, output_file, refs_added = validator.seed_missing_annotations(
                        test_file, goa_file=goa_file, fetch_titles=True
                    )

                    # Check results
                    assert added == 2  # Both annotations should be added

                    # Verify the seeded data
                    with open(output_file) as f:
                        result = yaml.safe_load(f)

                    # Check annotations were added
                    assert len(result["existing_annotations"]) == 2

                    # Check references were added - at least some should be present
                    assert "references" in result
                    assert len(result["references"]) > 0

                    # Verify the actual references that were added
                    ref_ids = [r["id"] for r in result["references"]]
                    print(f"Added references: {ref_ids}")

                    # The seed function collects unique references from GOA annotations
                    # Both PMID and GO_REF should be collected
                    unique_refs = set(["GO_REF:0000033", "PMID:12345678"])
                    added_refs = set(ref_ids)
                    assert added_refs == unique_refs, (
                        f"Expected {unique_refs}, got {added_refs}"
                    )

    finally:
        # Clean up
        import shutil

        shutil.rmtree(temp_dir)


def test_goa_validator_accepts_goref():
    """Test that GOA validator properly handles GO_REF references in existing annotations."""
    validator = GOAValidator()

    # Test YAML with GO_REF references
    test_file = Path("tests/input/test_goref.yaml")

    # Mock GOA data matching the test file
    mock_annotations = [
        Mock(
            go_id="GO:0005515",
            go_term="protein binding",
            evidence_code="IEA",
            reference="GO_REF:0000033",
        ),
        Mock(
            go_id="GO:0005634",
            go_term="nucleus",
            evidence_code="IEA",
            reference="GO_REF:0000002",
        ),
    ]

    # Create a temporary dummy GOA file
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix="-goa.tsv", delete=False) as f:
        goa_file = Path(f.name)
        # Write a minimal header
        f.write(
            "DB\tDB_Object_ID\tDB_Object_Symbol\tQualifier\tGO_ID\tGO_Term\tGO_Aspect\tEvidenceType\tEvidence\tReference\n"
        )

    try:
        with patch.object(validator, "parse_goa_file", return_value=mock_annotations):
            result = validator.validate_against_goa(test_file, goa_file)

            # Should be valid since GO_REF references match
            assert result.is_valid
            assert len(result.missing_in_yaml) == 0
            assert len(result.missing_in_goa) == 0
    finally:
        # Clean up
        goa_file.unlink()
