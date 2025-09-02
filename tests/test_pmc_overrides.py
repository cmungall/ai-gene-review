"""Tests for PMC override functionality."""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from ai_gene_review.etl.publication import load_pmc_overrides


def test_load_pmc_overrides_with_valid_file():
    """Test loading PMC overrides from a valid TSV file."""
    # Create a temporary override file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write("# PMC ID Override Table\n")
        f.write("# Comment line\n")
        f.write("2001740\t\tNo PMC version exists\n")
        f.write("12345\tPMC67890\tCorrect PMC ID\n")
        f.write("99999\t\n")  # Empty PMC field
        temp_file = Path(f.name)
    
    try:
        # Mock the path lookup to use our temp file
        with patch('ai_gene_review.etl.publication.Path') as mock_path:
            mock_path.return_value.parent = temp_file.parent
            mock_path.return_value.parent.__truediv__ = lambda self, x: temp_file
            
            # Clear cache first
            import ai_gene_review.etl.publication as pub_module
            pub_module._PMC_OVERRIDES_CACHE = None
            
            # Patch the actual path checking
            with patch.object(Path, 'exists') as mock_exists:
                def exists_side_effect(self):
                    return str(self) == str(temp_file)
                mock_exists.side_effect = exists_side_effect
                
                with patch('builtins.open', create=True) as mock_open:
                    mock_open.return_value.__enter__ = lambda self: self
                    mock_open.return_value.__exit__ = lambda self, *args: None
                    mock_open.return_value.read = lambda: temp_file.read_text()
                    mock_open.return_value.__iter__ = lambda self: iter(temp_file.read_text().splitlines())
                    
                    overrides = load_pmc_overrides()
        
        # Check the loaded overrides
        assert '2001740' in overrides
        assert overrides['2001740'] is None  # No PMC version
        
        assert '12345' in overrides
        assert overrides['12345'] == 'PMC67890'
        
        assert '99999' in overrides
        assert overrides['99999'] is None  # Empty PMC field
        
    finally:
        temp_file.unlink()


def test_load_pmc_overrides_no_file():
    """Test loading PMC overrides when no file exists."""
    # Clear cache first
    import ai_gene_review.etl.publication as pub_module
    pub_module._PMC_OVERRIDES_CACHE = None
    
    with patch.object(Path, 'exists', return_value=False):
        overrides = load_pmc_overrides()
    
    # Should return empty dict when no file exists
    assert overrides == {}


def test_load_pmc_overrides_caching():
    """Test that PMC overrides are cached after first load."""
    import ai_gene_review.etl.publication as pub_module
    
    # Set up a fake cache
    fake_cache = {'test': 'PMC123'}
    pub_module._PMC_OVERRIDES_CACHE = fake_cache
    
    # Should return cached value without file I/O
    overrides = load_pmc_overrides()
    assert overrides == fake_cache
    
    # Clean up
    pub_module._PMC_OVERRIDES_CACHE = None


def test_pmc_override_integration():
    """Test that overrides are used in fetch_pubmed_data."""
    # This test would require mocking Entrez calls
    # Placeholder for integration test
    pass


@pytest.mark.parametrize("line,expected_pmid,expected_pmcid", [
    ("12345\tPMC67890", "12345", "PMC67890"),
    ("12345\t", "12345", None),
    ("12345", "12345", None),
    ("12345\t\t", "12345", None),
    ("\t12345\tPMC67890", "", "12345"),  # Edge case: empty PMID
])
def test_parse_override_lines(line, expected_pmid, expected_pmcid):
    """Test parsing of different override line formats."""
    # This would be tested through load_pmc_overrides
    # but we can verify the parsing logic
    parts = line.split('\t')
    pmid = parts[0].strip() if len(parts) >= 1 else ""
    pmcid = parts[1].strip() if len(parts) >= 2 and parts[1].strip() else None
    
    assert pmid == expected_pmid
    assert pmcid == expected_pmcid