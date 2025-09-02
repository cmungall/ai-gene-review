"""Test organism code support in gene fetching.

This module tests the ability to use UniProt organism codes (like PSEPK, ECOLI)
in addition to common names and taxonomy IDs.
"""

from unittest.mock import patch, MagicMock
from ai_gene_review.etl.gene import (
    resolve_organism_code_to_taxon,
    uniprot_by_gene_taxon,
    get_organism_name_from_uniprot,
    expand_organism_name,
)


class TestOrganismCodeResolution:
    """Test organism code to taxon ID resolution."""

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_resolve_organism_code_to_taxon(self, mock_get):
        """Test resolving UniProt organism codes to taxonomy IDs."""
        # Mock response for PSEPK
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [
                {
                    "taxonId": 160488,
                    "scientificName": "Pseudomonas putida KT2440",
                    "mnemonic": "PSEPK",
                }
            ]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        # Test PSEPK
        taxon_id = resolve_organism_code_to_taxon("PSEPK")
        assert taxon_id == "160488"

        # Verify API call
        mock_get.assert_called_with(
            "https://rest.uniprot.org/taxonomy/stream",
            params={"query": "mnemonic:PSEPK", "format": "json", "size": "1"},
            timeout=10,
        )

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_resolve_organism_code_not_found(self, mock_get):
        """Test handling of unknown organism codes."""
        # Mock empty response
        mock_response = MagicMock()
        mock_response.json.return_value = {"results": []}
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        # Test unknown code
        taxon_id = resolve_organism_code_to_taxon("XXXXX")
        assert taxon_id is None

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_resolve_organism_code_api_error(self, mock_get):
        """Test handling of API errors."""
        # Mock API error
        mock_get.side_effect = Exception("API Error")

        # Should return None on error
        taxon_id = resolve_organism_code_to_taxon("PSEPK")
        assert taxon_id is None


class TestOrganismNameFetching:
    """Test fetching organism names from UniProt."""

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_get_organism_name_from_uniprot(self, mock_get):
        """Test fetching organism name from UniProt ID."""
        # Mock response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "organism": {"scientificName": "Pseudomonas putida (strain KT2440)"}
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        # Test fetching name
        name = get_organism_name_from_uniprot("Q88JH0")
        assert name == "Pseudomonas putida (strain KT2440)"

        # Verify API call
        mock_get.assert_called_with(
            "https://rest.uniprot.org/uniprotkb/Q88JH0?format=json&fields=organism_name",
            timeout=10,
        )

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_get_organism_name_api_error(self, mock_get):
        """Test handling of API errors when fetching organism name."""
        mock_get.side_effect = Exception("API Error")

        # Should return None on error
        name = get_organism_name_from_uniprot("Q88JH0")
        assert name is None


class TestUniprotLookupWithOrganismCodes:
    """Test UniProt lookups with different organism identifiers."""

    @patch("ai_gene_review.etl.gene.requests.get")
    @patch("ai_gene_review.etl.gene.resolve_organism_code_to_taxon")
    def test_uniprot_lookup_with_organism_code(self, mock_resolve, mock_get):
        """Test UniProt lookup using organism code."""
        # Mock organism code resolution
        mock_resolve.return_value = "160488"

        # Mock UniProt API response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [
                {
                    "primaryAccession": "Q88JH0",
                    "uniProtkbId": "PEDH_PSEPK",
                    "genes": [{"geneName": {"value": "pedH"}}],
                    "organism": {"scientificName": "Pseudomonas putida KT2440"},
                }
            ]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        # Test lookup
        results = uniprot_by_gene_taxon("pedH", "PSEPK", limit=1)

        assert len(results) == 1
        assert results[0]["accession"] == "Q88JH0"
        assert results[0]["gene"] == "pedH"
        assert results[0]["organism"] == "Pseudomonas putida KT2440"

        # Verify organism code was resolved
        mock_resolve.assert_called_once_with("PSEPK")

        # Verify correct query was built
        expected_query = (
            "(gene_exact:pedH) AND (organism_id:160488) AND (reviewed:true)"
        )
        call_args = mock_get.call_args[1]["params"]
        assert call_args["query"] == expected_query

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_uniprot_lookup_with_common_name(self, mock_get):
        """Test UniProt lookup using common organism name."""
        # Mock UniProt API response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [
                {
                    "primaryAccession": "P04637",
                    "uniProtkbId": "P53_HUMAN",
                    "genes": [{"geneName": {"value": "TP53"}}],
                    "organism": {"scientificName": "Homo sapiens"},
                }
            ]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        # Test lookup with common name
        results = uniprot_by_gene_taxon("TP53", "human", limit=1)

        assert len(results) == 1
        assert results[0]["accession"] == "P04637"
        assert results[0]["gene"] == "TP53"
        assert results[0]["organism"] == "Homo sapiens"

        # Verify correct query was built (should use taxon ID for human)
        expected_query = "(gene_exact:TP53) AND (organism_id:9606) AND (reviewed:true)"
        call_args = mock_get.call_args[1]["params"]
        assert call_args["query"] == expected_query

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_uniprot_lookup_with_taxon_id(self, mock_get):
        """Test UniProt lookup using NCBI taxon ID."""
        # Mock UniProt API response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [
                {
                    "primaryAccession": "P00722",
                    "uniProtkbId": "BGAL_ECOLI",
                    "genes": [{"geneName": {"value": "lacZ"}}],
                    "organism": {"scientificName": "Escherichia coli K12"},
                }
            ]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        # Test lookup with taxon ID
        results = uniprot_by_gene_taxon("lacZ", "83333", limit=1)

        assert len(results) == 1
        assert results[0]["accession"] == "P00722"
        assert results[0]["gene"] == "lacZ"

        # Verify correct query was built
        expected_query = "(gene_exact:lacZ) AND (organism_id:83333) AND (reviewed:true)"
        call_args = mock_get.call_args[1]["params"]
        assert call_args["query"] == expected_query


class TestOrganismNameExpansion:
    """Test organism name expansion functionality."""

    def test_expand_common_names(self):
        """Test expansion of common organism names."""
        assert expand_organism_name("human") == "Homo sapiens"
        assert expand_organism_name("mouse") == "Mus musculus"
        assert expand_organism_name("yeast") == "Saccharomyces cerevisiae"
        assert expand_organism_name("fly") == "Drosophila melanogaster"
        assert expand_organism_name("worm") == "Caenorhabditis elegans"

    def test_preserve_scientific_names(self):
        """Test that scientific names are preserved."""
        assert expand_organism_name("Homo sapiens") == "Homo sapiens"
        assert expand_organism_name("Escherichia coli") == "Escherichia coli"

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_expand_organism_code(self, mock_get):
        """Test expansion of UniProt organism codes."""
        # Mock response for PSEPK
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [
                {
                    "scientificName": "Pseudomonas putida (strain KT2440)",
                    "taxonId": 160488,
                    "mnemonic": "PSEPK",
                }
            ]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        result = expand_organism_name("PSEPK")
        assert result == "Pseudomonas putida (strain KT2440)"

        # Verify API call
        mock_get.assert_called_with(
            "https://rest.uniprot.org/taxonomy/stream",
            params={"query": "mnemonic:PSEPK", "format": "json", "size": "1"},
            timeout=10,
        )

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_expand_taxon_id(self, mock_get):
        """Test expansion of taxon IDs."""
        # Mock response for taxon ID
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [
                {
                    "scientificName": "Escherichia coli str. K-12 substr. MG1655",
                    "taxonId": 511145,
                }
            ]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        result = expand_organism_name("511145")
        assert result == "Escherichia coli str. K-12 substr. MG1655"

        # Verify API call
        mock_get.assert_called_with(
            "https://rest.uniprot.org/taxonomy/stream",
            params={"query": "id:511145", "format": "json", "size": "1"},
            timeout=10,
        )

    @patch("ai_gene_review.etl.gene.requests.get")
    def test_expand_unknown_code(self, mock_get):
        """Test handling of unknown organism codes."""
        # Mock empty response
        mock_response = MagicMock()
        mock_response.json.return_value = {"results": []}
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        # Should return original if not found
        result = expand_organism_name("XXXXX")
        assert result == "XXXXX"

    def test_case_insensitive_common_names(self):
        """Test that common names are case-insensitive."""
        assert expand_organism_name("Human") == "Homo sapiens"
        assert expand_organism_name("HUMAN") == "Homo sapiens"
        assert expand_organism_name("ecoli") == "Escherichia coli"
        assert expand_organism_name("E.coli") == "Escherichia coli"
