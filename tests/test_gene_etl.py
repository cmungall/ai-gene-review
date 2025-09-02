"""Tests for the gene ETL module."""

import tempfile
from pathlib import Path
import pytest
from ai_gene_review.etl.gene import fetch_gene_data


@pytest.mark.integration
@pytest.mark.parametrize(
    "organism,gene,uniprot_id",
    [
        ("human", "CFAP300", "Q9BRQ4"),
        ("human", "TP53", "P04637"),
        ("human", "BRCA1", "P38398"),
    ],
)
def test_fetch_gene_data_creates_correct_directory_structure(
    organism, gene, uniprot_id
):
    """Test that fetch_gene_data creates the expected directory structure and fetches real data.

    This integration test verifies:
    - Correct directory structure is created
    - UniProt data is fetched and saved
    - GOA data is fetched and saved
    - Files contain expected content
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        base_path = Path(tmpdir)

        # Call the function with real data
        fetch_gene_data((organism, gene), uniprot_id=uniprot_id, base_path=base_path)

        # Check directory structure
        gene_dir = base_path / "genes" / organism / gene
        assert gene_dir.exists(), f"Directory {gene_dir} was not created"

        # Check UniProt file
        uniprot_file = gene_dir / f"{gene}-uniprot.txt"
        assert uniprot_file.exists(), f"UniProt file {uniprot_file} was not created"

        uniprot_content = uniprot_file.read_text()
        assert len(uniprot_content) > 100, "UniProt file seems too small"
        assert uniprot_id in uniprot_content, (
            f"UniProt ID {uniprot_id} not found in content"
        )

        # Check GOA file
        goa_file = gene_dir / f"{gene}-goa.tsv"
        assert goa_file.exists(), f"GOA file {goa_file} was not created"

        goa_content = goa_file.read_text()
        assert len(goa_content) > 50, "GOA file seems too small"
        # GOA files should have tab-separated headers
        assert "\t" in goa_content, "GOA file should be tab-separated"


@pytest.mark.integration
def test_fetch_gene_data_without_uniprot_id():
    """Test that fetch_gene_data can resolve gene names to UniProt IDs.

    This test uses real UniProt API to resolve gene names.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        base_path = Path(tmpdir)

        # Call without providing UniProt ID - should resolve automatically
        fetch_gene_data(("human", "CFAP300"), base_path=base_path)

        # Check that files were created
        gene_dir = base_path / "genes" / "human" / "CFAP300"
        assert gene_dir.exists()

        uniprot_file = gene_dir / "CFAP300-uniprot.txt"
        assert uniprot_file.exists()

        # Verify it found the correct UniProt entry
        uniprot_content = uniprot_file.read_text()
        assert "Q9BRQ4" in uniprot_content or "CFAP300" in uniprot_content


@pytest.mark.integration
def test_fetch_gene_data_handles_invalid_gene():
    """Test graceful handling when gene doesn't exist."""
    with tempfile.TemporaryDirectory() as tmpdir:
        base_path = Path(tmpdir)

        # Should raise an exception for non-existent gene
        with pytest.raises(ValueError, match="Could not find.*UniProt ID"):
            fetch_gene_data(("human", "NONEXISTENTGENE12345"), base_path=base_path)


@pytest.mark.integration
def test_fetch_gene_data_handles_invalid_uniprot_id():
    """Test graceful handling when UniProt ID is invalid."""
    with tempfile.TemporaryDirectory() as tmpdir:
        base_path = Path(tmpdir)

        # Should raise an exception for invalid UniProt ID
        with pytest.raises(ValueError, match="Failed to fetch UniProt data"):
            fetch_gene_data(
                ("human", "FAKEGENE"), uniprot_id="INVALID123", base_path=base_path
            )


@pytest.mark.integration
@pytest.mark.parametrize(
    "organism,gene",
    [
        ("mouse", "Trp53"),  # Mouse p53
        ("yeast", "CDC28"),  # Yeast cell division control protein
    ],
)
def test_fetch_gene_data_different_organisms(organism, gene):
    """Test fetching genes from different organisms."""
    with tempfile.TemporaryDirectory() as tmpdir:
        base_path = Path(tmpdir)

        # This will attempt to resolve the gene name for the organism
        fetch_gene_data((organism, gene), base_path=base_path)

        # Check that files were created
        gene_dir = base_path / "genes" / organism / gene
        assert gene_dir.exists()

        uniprot_file = gene_dir / f"{gene}-uniprot.txt"
        assert uniprot_file.exists()
        assert uniprot_file.read_text()  # Should have content
