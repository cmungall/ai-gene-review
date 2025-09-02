"""Test the core functions validation rules."""

import tempfile
from pathlib import Path
import yaml

from ai_gene_review.validation import validate_gene_review
from ai_gene_review.validation.validation_report import ValidationSeverity


def test_core_function_without_support_or_accepted_term():
    """Test that core functions must either come from ACCEPTED annotations or have supported_by."""
    
    # Create a test YAML with a core function that is neither from ACCEPTED nor has supported_by
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "A test gene for validation",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005524", "label": "ATP binding"},
                "evidence_type": "IEA",
                "review": {
                    "action": "REMOVE",  # Not ACCEPTED
                    "reason": "Too general"
                }
            }
        ],
        "core_functions": [
            {
                "molecular_function": {
                    "id": "GO:0004672",  # protein kinase activity - not from ACCEPTED
                    "label": "protein kinase activity"
                },
                "description": "Test function"
                # No supported_by field
            }
        ]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(data, f)
        test_file = Path(f.name)
    
    try:
        report = validate_gene_review(test_file, check_goa=False)
        
        # Should have an error about the core function
        errors = [issue for issue in report.issues if issue.severity == ValidationSeverity.ERROR]
        assert any(
            "GO:0004672" in issue.message and "not from an ACCEPTED annotation" in issue.message
            for issue in errors
        ), "Should report error for core function not from ACCEPTED annotation and lacking supported_by"
    finally:
        test_file.unlink()


def test_core_function_from_accepted_annotation():
    """Test that core functions from ACCEPTED annotations are valid."""
    
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "A test gene for validation",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005524", "label": "ATP binding"},
                "evidence_type": "IEA",
                "review": {
                    "action": "ACCEPT",  # ACCEPTED
                    "reason": "Well supported"
                }
            }
        ],
        "core_functions": [
            {
                "molecular_function": {
                    "id": "GO:0005524",  # Same as ACCEPTED annotation
                    "label": "ATP binding"
                },
                "description": "Binds ATP"
                # No supported_by needed since it's from ACCEPTED
            }
        ]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(data, f)
        test_file = Path(f.name)
    
    try:
        report = validate_gene_review(test_file, check_goa=False)
        
        # Should not have errors about this core function
        errors = [issue for issue in report.issues if issue.severity == ValidationSeverity.ERROR]
        assert not any(
            "GO:0005524" in issue.message and "not from an ACCEPTED annotation" in issue.message
            for issue in errors
        ), "Should not report error for core function from ACCEPTED annotation"
    finally:
        test_file.unlink()


def test_core_function_with_supported_by():
    """Test that core functions with supported_by are valid even if not from ACCEPTED."""
    
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "A test gene for validation",
        "references": [
            {"id": "PMID:12345", "title": "Test paper"}
        ],
        "existing_annotations": [],  # No existing annotations
        "core_functions": [
            {
                "molecular_function": {
                    "id": "GO:0004672",
                    "label": "protein kinase activity"
                },
                "description": "Novel kinase function",
                "supported_by": [
                    {
                        "reference_id": "PMID:12345",
                        "supporting_text": "The protein shows kinase activity"
                    }
                ]
            }
        ]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(data, f)
        test_file = Path(f.name)
    
    try:
        report = validate_gene_review(test_file, check_goa=False, check_supporting_text=False)
        
        # Should not have errors about this core function
        errors = [issue for issue in report.issues if issue.severity == ValidationSeverity.ERROR]
        assert not any(
            "GO:0004672" in issue.message and "not from an ACCEPTED annotation" in issue.message
            for issue in errors
        ), "Should not report error for core function with supported_by"
    finally:
        test_file.unlink()


def test_core_function_from_proposed_replacement():
    """Test that core functions from proposed_replacement_terms in MODIFY actions are valid."""
    
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "A test gene for validation",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005524", "label": "ATP binding"},
                "evidence_type": "IEA",
                "review": {
                    "action": "MODIFY",
                    "reason": "More specific term available",
                    "proposed_replacement_terms": [
                        {"id": "GO:0004672", "label": "protein kinase activity"}
                    ]
                }
            }
        ],
        "core_functions": [
            {
                "molecular_function": {
                    "id": "GO:0004672",  # From proposed_replacement_terms
                    "label": "protein kinase activity"
                },
                "description": "Kinase function"
            }
        ]
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(data, f)
        test_file = Path(f.name)
    
    try:
        report = validate_gene_review(test_file, check_goa=False)
        
        # Should not have errors about this core function
        errors = [issue for issue in report.issues if issue.severity == ValidationSeverity.ERROR]
        assert not any(
            "GO:0004672" in issue.message and "not from an ACCEPTED annotation" in issue.message
            for issue in errors
        ), "Should not report error for core function from proposed_replacement_terms"
    finally:
        test_file.unlink()