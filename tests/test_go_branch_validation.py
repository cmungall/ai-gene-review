"""Test validation of GO term branch constraints in core_functions."""

import tempfile
from pathlib import Path

import pytest
import yaml

from ai_gene_review.validation import validate_gene_review, ValidationSeverity


def test_molecular_function_must_be_mf_branch():
    """Test that molecular_function field requires MF branch GO terms."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST1",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene for GO branch validation",
        "core_functions": [
            {
                "description": "Test function",
                "molecular_function": {
                    "id": "GO:0008150",  # biological_process root - wrong branch!
                    "label": "biological_process",
                },
            }
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        report = validate_gene_review(temp_path)

        # Should have error about wrong branch
        branch_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "biological_process branch but should be in the molecular_function branch"
            in issue.message
        ]
        assert len(branch_errors) == 1
        assert not report.is_valid
    finally:
        temp_path.unlink()


def test_directly_involved_in_must_be_bp_branch():
    """Test that directly_involved_in field requires BP branch GO terms."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST2",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene for GO branch validation",
        "core_functions": [
            {
                "description": "Test function",
                "molecular_function": {
                    "id": "GO:0003674",  # molecular_function root - correct
                    "label": "molecular_function",
                },
                "directly_involved_in": [
                    {
                        "id": "GO:0003674",  # molecular_function root - wrong branch!
                        "label": "molecular_function",
                    }
                ],
            }
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        report = validate_gene_review(temp_path)

        # Should have error about wrong branch for directly_involved_in
        branch_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "molecular_function branch but should be in the biological_process branch"
            in issue.message
        ]
        assert len(branch_errors) == 1
        assert not report.is_valid
    finally:
        temp_path.unlink()


def test_locations_must_be_cc_branch():
    """Test that locations field requires CC branch GO terms."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST3",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene for GO branch validation",
        "core_functions": [
            {
                "description": "Test function",
                "molecular_function": {
                    "id": "GO:0003674",  # molecular_function root - correct
                    "label": "molecular_function",
                },
                "locations": [
                    {
                        "id": "GO:0008150",  # biological_process root - wrong branch!
                        "label": "biological_process",
                    }
                ],
            }
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        report = validate_gene_review(temp_path)

        # Should have error about wrong branch for locations
        branch_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "biological_process branch but should be in the cellular_component branch"
            in issue.message
        ]
        assert len(branch_errors) == 1
        assert not report.is_valid
    finally:
        temp_path.unlink()


def test_correct_go_branches_pass():
    """Test that correctly placed GO terms pass validation."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST4",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene with correct GO branches",
        "core_functions": [
            {
                "description": "Test function with correct branches",
                "molecular_function": {
                    "id": "GO:0003723",  # RNA binding - MF branch
                    "label": "RNA binding",
                },
                "directly_involved_in": [
                    {
                        "id": "GO:0006397",  # mRNA processing - BP branch
                        "label": "mRNA processing",
                    }
                ],
                "locations": [
                    {
                        "id": "GO:0005634",  # nucleus - CC branch
                        "label": "nucleus",
                    }
                ],
            }
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_path = Path(f.name)

    try:
        report = validate_gene_review(temp_path)

        # Should not have any GO branch errors
        branch_errors = [
            issue
            for issue in report.issues
            if issue.severity == ValidationSeverity.ERROR
            and "branch but should be in" in issue.message
        ]
        assert len(branch_errors) == 0

        # May have other validation issues (like term label checks) but no branch issues
        # Check that if it's invalid, it's not due to branch errors
        if not report.is_valid:
            # Make sure failures are not branch-related
            for issue in report.issues:
                if issue.severity == ValidationSeverity.ERROR:
                    assert "branch but should be in" not in issue.message
    finally:
        temp_path.unlink()


@pytest.mark.integration
def test_real_gene_branch_validation():
    """Test GO branch validation on a real gene file that should pass."""
    cfap300_file = Path("genes/human/CFAP300/CFAP300-ai-review.yaml")
    if not cfap300_file.exists():
        pytest.skip("CFAP300 file not found")

    report = validate_gene_review(cfap300_file)

    # CFAP300 should have correct branches now
    # molecular_function: GO:0030674 (protein-macromolecule adaptor activity) is in MF
    # directly_involved_in items are in BP
    # locations items are in CC

    branch_errors = [
        issue
        for issue in report.issues
        if issue.severity == ValidationSeverity.ERROR
        and "branch but should be in" in issue.message
    ]

    # Should not have any branch errors
    assert len(branch_errors) == 0, f"Unexpected branch errors: {branch_errors}"
