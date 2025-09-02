"""Tests for publication ETL functionality."""

import tempfile
from pathlib import Path

import pytest
import yaml

from ai_gene_review.etl.publication import (
    Publication,
    cache_publication,
    extract_pmid,
    extract_pmids_from_yaml,
    fetch_pubmed_data,
)


def test_extract_pmid():
    """Test PMID extraction from various formats."""
    assert extract_pmid("12345678") == "12345678"
    assert extract_pmid("PMID:12345678") == "12345678"
    assert extract_pmid("pmid: 12345678") == "12345678"
    assert extract_pmid("PMID12345678") == "12345678"
    assert extract_pmid(" PMID: 12345678 ") == "12345678"
    assert extract_pmid("pmid:12345678") == "12345678"


def test_publication_to_markdown():
    """Test Publication markdown generation."""
    pub = Publication(
        pmid="12345",
        title="Test Article",
        authors=["Smith J", "Doe J"],
        journal="Test Journal",
        year="2024",
        abstract="This is a test abstract.",
        doi="10.1234/test",
    )

    markdown = pub.to_markdown()

    # Check frontmatter
    assert "pmid: '12345'" in markdown
    assert "title: Test Article" in markdown
    assert "authors:" in markdown
    assert "- Smith J" in markdown
    assert "- Doe J" in markdown
    assert "journal: Test Journal" in markdown
    assert "year: '2024'" in markdown
    assert "doi: 10.1234/test" in markdown

    # Check body
    assert "# Test Article" in markdown
    assert "**Authors:** Smith J, Doe J" in markdown
    assert "**Journal:** Test Journal (2024)" in markdown
    assert "**DOI:** [10.1234/test](https://doi.org/10.1234/test)" in markdown
    assert "## Abstract" in markdown
    assert "This is a test abstract." in markdown


def test_publication_to_frontmatter_dict():
    """Test Publication frontmatter dictionary generation."""
    pub = Publication(
        pmid="12345",
        title="Test Article",
        authors=["Smith J"],
        journal="Test Journal",
        year="2024",
        abstract="Abstract text",
        pmcid="PMC123456",
        doi="10.1234/test",
        keywords=["test", "publication"],
    )

    frontmatter = pub.to_frontmatter_dict()

    assert frontmatter["pmid"] == "12345"
    assert frontmatter["title"] == "Test Article"
    assert frontmatter["authors"] == ["Smith J"]
    assert frontmatter["journal"] == "Test Journal"
    assert frontmatter["year"] == "2024"
    assert frontmatter["pmcid"] == "PMC123456"
    assert frontmatter["doi"] == "10.1234/test"
    assert frontmatter["keywords"] == ["test", "publication"]

    # Abstract should not be in frontmatter
    assert "abstract" not in frontmatter


def test_extract_pmids_from_yaml():
    """Test extracting PMIDs from a gene review YAML file."""
    data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "references": [
            {"id": "PMID:12345", "title": "Paper 1"},
            {"id": "PMID:67890", "title": "Paper 2"},
            {
                "id": "DOI:10.1234/test",
                "title": "Paper with DOI only",
            },  # Should be ignored
        ],
        "existing_annotations": [
            {"term": {"id": "GO:0005515"}, "original_reference_id": "PMID:11111"},
            {
                "term": {"id": "GO:0005737"},
                "original_reference_id": "PMID:12345",  # Duplicate, should be deduplicated
            },
            {
                "term": {"id": "GO:0005856"},
                # No original_reference_id
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_file = Path(f.name)

    try:
        pmids = extract_pmids_from_yaml(temp_file)

        # Should extract unique PMIDs
        assert len(pmids) == 3
        assert set(pmids) == {"12345", "67890", "11111"}
    finally:
        temp_file.unlink()


def test_extract_pmids_from_empty_yaml():
    """Test extracting PMIDs from an empty or minimal YAML file."""
    data = {"id": "Q12345", "gene_symbol": "TEST", "description": "Test gene"}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        temp_file = Path(f.name)

    try:
        pmids = extract_pmids_from_yaml(temp_file)
        assert pmids == []
    finally:
        temp_file.unlink()


@pytest.mark.integration
def test_fetch_pubmed_data():
    """Test fetching real data from PubMed.

    This test requires network access and is marked as integration.
    """
    # Use a known PMID from CFAP300 literature
    pmid = "29727692"

    pub = fetch_pubmed_data(pmid)

    assert pub is not None
    assert pub.pmid == pmid
    assert "C11orf70" in pub.title or "CFAP300" in pub.title
    assert len(pub.authors) > 0
    assert pub.year == "2018"
    assert len(pub.abstract) > 100  # Should have substantial abstract
    assert pub.pmcid == "PMC5986720"  # This paper has PMC full text
    assert pub.doi is not None


@pytest.mark.integration
def test_cache_publication():
    """Test caching a publication to filesystem.

    This test requires network access and is marked as integration.
    """
    pmid = "29727692"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        # First cache should succeed
        success = cache_publication(pmid, output_dir, force=False)
        assert success

        # Check file was created
        expected_file = output_dir / f"PMID_{pmid}.md"
        assert expected_file.exists()

        # Check content
        content = expected_file.read_text()
        assert f"pmid: '{pmid}'" in content
        assert "# C11orf70" in content or "# CFAP300" in content
        assert "## Abstract" in content

        # Second cache without force should still succeed (already cached)
        success = cache_publication(pmid, output_dir, force=False)
        assert success

        # Cache with force should re-download
        original_mtime = expected_file.stat().st_mtime
        import time

        time.sleep(0.1)  # Ensure time difference
        success = cache_publication(pmid, output_dir, force=True)
        assert success
        new_mtime = expected_file.stat().st_mtime
        assert new_mtime > original_mtime  # File was updated


def test_cache_publication_invalid_pmid():
    """Test caching with an invalid PMID."""
    pmid = "99999999999"  # Unlikely to exist

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        success = cache_publication(pmid, output_dir)
        assert not success

        # No file should be created
        expected_file = output_dir / f"PMID_{pmid}.md"
        assert not expected_file.exists()
