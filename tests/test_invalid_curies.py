"""Test validation of invalid/made-up CURIEs."""

import tempfile
from pathlib import Path

import yaml

from ai_gene_review.validation import (
    TermValidator,
    validate_gene_review,
    ValidationSeverity,
)


def test_invalid_curie_detection():
    """Test that made-up CURIEs are detected as invalid."""
    validator = TermValidator()

    data = {
        "core_functions": [
            {
                "substrates": [
                    {"id": "CFAP300:DNAAF2_complex", "label": "CFAP300-DNAAF2 complex"},
                    {
                        "id": "MYPROTEIN:binding_site",
                        "label": "My protein binding site",
                    },
                ]
            }
        ]
    }

    results = validator.validate_terms_in_data(data)

    # Should find 2 invalid CURIEs
    invalid_results = [r for r in results if not r.is_valid]
    assert len(invalid_results) == 2

    # Check error messages
    for result in invalid_results:
        assert "Unrecognized identifier prefix" in result.error_message
        assert "Valid ontology prefixes are:" in result.error_message

    # Check specific invalid IDs
    invalid_ids = {r.term_id for r in invalid_results}
    assert "CFAP300:DNAAF2_complex" in invalid_ids
    assert "MYPROTEIN:binding_site" in invalid_ids


def test_excluded_prefixes_not_validated():
    """Test that known non-ontology IDs are not validated (except TEMP needs description)."""
    validator = TermValidator()

    data = {
        "references": [
            {"id": "PMID:12345", "title": "A paper"},
            {"id": "DOI:10.1234/test", "doi": "10.1234/test"},
        ],
        "database_refs": [
            {"id": "UniProt:P12345", "label": "Some protein"},
            {"id": "PDB:1ABC", "label": "Crystal structure"},
        ],
    }

    results = validator.validate_terms_in_data(data)

    # Should not validate any of these (all are excluded prefixes)
    assert len(results) == 0


def test_valid_ontology_terms_still_validated():
    """Test that valid ontology terms are still validated."""
    validator = TermValidator()

    # Mock the label lookup
    validator._label_cache = {
        "GO:0005515": "protein binding",
        "CHEBI:12345": "some chemical",
    }

    data = {
        "annotations": [
            {
                "id": "GO:0005515",
                "label": "protein binding",  # Correct
            },
            {
                "id": "CHEBI:12345",
                "label": "wrong label",  # Wrong
            },
        ]
    }

    results = validator.validate_terms_in_data(data)

    # Should validate both ontology terms
    assert len(results) == 2

    valid = [r for r in results if r.is_valid]
    invalid = [r for r in results if not r.is_valid]

    assert len(valid) == 1
    assert len(invalid) == 1

    assert valid[0].term_id == "GO:0005515"
    assert invalid[0].term_id == "CHEBI:12345"


def test_integration_with_main_validator():
    """Test that invalid CURIE detection works in main validator."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene",
        "core_functions": [
            {
                "description": "Test function",
                "substrates": [
                    {"id": "MADEUP:complex_123", "label": "Made-up complex"}
                ],
            }
        ],
        "references": [
            {
                "id": "PMID:12345",  # Should be ignored
                "title": "Test paper",
            }
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        report = validate_gene_review(temp_path)

        # Should have error for the made-up CURIE
        curie_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "Unrecognized identifier prefix 'MADEUP'" in issue.message
        ]

        assert len(curie_errors) == 1
        assert not report.is_valid

        # Should NOT have errors for PMID
        pmid_errors = [
            issue
            for issue in report.issues
            if "PMID" in issue.message and "Unrecognized" in issue.message
        ]
        assert len(pmid_errors) == 0

    finally:
        temp_path.unlink()


def test_temp_identifier_requires_description():
    """Test that TEMP identifiers require a description field."""
    validator = TermValidator()

    # TEMP without description - should fail
    data_no_desc = {
        "substrates": [
            {
                "id": "TEMP:placeholder_1",
                "label": "Temporary complex",
                # No description - should fail
            }
        ]
    }

    results = validator.validate_terms_in_data(data_no_desc)
    assert len(results) == 1
    assert not results[0].is_valid
    assert "requires a 'description' field" in results[0].error_message

    # TEMP with description - should pass (no validation errors)
    data_with_desc = {
        "substrates": [
            {
                "id": "TEMP:placeholder_1",
                "label": "Temporary complex",
                "description": "This is a placeholder for a protein complex to be characterized",
            }
        ]
    }

    results = validator.validate_terms_in_data(data_with_desc)
    assert len(results) == 0  # No errors, TEMP is excluded when it has description


def test_temp_multiple_items():
    """Test TEMP validation with multiple items, some with and some without descriptions."""
    validator = TermValidator()

    data = {
        "items": [
            {
                "id": "TEMP:item1",
                "label": "Item 1",
                "description": "Properly documented temporary item",
            },
            {
                "id": "TEMP:item2",
                "label": "Item 2",
                # Missing description
            },
            {
                "id": "GO:0005515",
                "label": "protein binding",
                # Regular ontology term, no description needed
            },
            {
                "id": "TEMP:item3",
                "label": "Item 3",
                # Missing description
            },
        ]
    }

    # Mock GO term validation
    validator._label_cache["GO:0005515"] = "protein binding"

    results = validator.validate_terms_in_data(data)

    # Should have 2 errors for TEMP items without descriptions
    # and 1 successful validation for GO term
    errors = [r for r in results if not r.is_valid]
    valid = [r for r in results if r.is_valid]

    assert len(errors) == 2  # Two TEMP items without descriptions
    assert len(valid) == 1  # GO term validates successfully

    # Check that the errors are for the right items
    error_ids = {r.term_id for r in errors}
    assert "TEMP:item2" in error_ids
    assert "TEMP:item3" in error_ids


def test_mixed_valid_and_invalid_curies():
    """Test a mix of valid ontology terms, excluded IDs, and invalid CURIEs."""
    validator = TermValidator()

    # Mock some valid terms
    validator._label_cache = {
        "GO:0005515": "protein binding",
        "NCBITaxon:9606": "Homo sapiens",
    }

    data = {
        "taxon": {
            "id": "NCBITaxon:9606",  # Valid ontology term
            "label": "Homo sapiens",
        },
        "annotations": [
            {
                "id": "GO:0005515",  # Valid ontology term
                "label": "protein binding",
            }
        ],
        "references": [
            {
                "id": "PMID:12345",  # Excluded prefix - should be ignored
                "title": "Paper",
            }
        ],
        "complexes": [
            {
                "id": "COMPLEX:made_up_1",  # Invalid CURIE
                "label": "Complex 1",
            },
            {
                "id": "GENE:fake_2",  # Invalid CURIE
                "label": "Gene product",
            },
        ],
    }

    results = validator.validate_terms_in_data(data)

    # Should have 2 valid terms and 2 invalid CURIEs
    # PMIDs should be ignored
    assert len(results) == 4

    valid = [r for r in results if r.is_valid]
    invalid = [r for r in results if not r.is_valid]

    assert len(valid) == 2  # GO and NCBITaxon terms
    assert len(invalid) == 2  # COMPLEX and GENE CURIEs

    # Check the invalid ones are the made-up CURIEs
    invalid_ids = {r.term_id for r in invalid}
    assert "COMPLEX:made_up_1" in invalid_ids
    assert "GENE:fake_2" in invalid_ids
