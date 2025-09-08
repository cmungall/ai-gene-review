"""Tests for supporting_text validator."""

import pytest
from pathlib import Path
import tempfile
import yaml
from ai_gene_review.validation.supporting_text_validator import (
    SupportingTextValidator,
    validate_supporting_text_in_file,
)


class TestSupportingTextValidator:
    """Test suite for SupportingTextValidator."""

    @pytest.fixture
    def validator(self):
        """Create a validator instance."""
        return SupportingTextValidator()

    @pytest.fixture
    def temp_publications_dir(self):
        """Create a temporary directory with mock publications."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create mock publication
            pub_file = tmpdir_path / "PMID_12345678.md"
            pub_content = """
---
pmid: '12345678'
title: Test Publication
---

# Test Publication

JAK1 is a non-receptor tyrosine kinase that plays a critical role in cytokine signaling.
It contains four major domains including a FERM domain and a kinase domain.
JAK1 is essential for Type I and Type II interferon signaling.
"""
            pub_file.write_text(pub_content)

            yield tmpdir_path

    def test_basic_functionality(self, validator):
        """Test basic functionality of supporting text validator."""
        data = {
            "existing_annotations": [
                {
                    "term": {"id": "GO:0000001", "label": "test"},
                    "evidence_type": "IEA",
                    "original_reference_id": "PMID:12345678",
                    "review": {
                        "summary": "Test annotation",
                        "action": "ACCEPT",
                        "reason": "Good annotation (PMID:12345678)",
                        "supported_by": [
                            {
                                "reference_id": "PMID:12345678",
                                "supporting_text": "This text should be in the publication",
                            }
                        ],
                    },
                }
            ]
        }

        report = validator.validate_data(data)

        # Should have found one annotation with supporting_text (in supported_by)
        assert report.annotations_with_supporting_text == 1

    @pytest.mark.parametrize(
        "input_text,expected",
        [
            ("  This is   a TEST   string!  ", "this is a test string!"),
            ("JAK1-mediated (STAT) activation", "jak1-mediated stat activation"),
            ("", ""),
            ("UPPER CASE TEXT", "upper case text"),
            ("  multiple     spaces    ", "multiple spaces"),
        ],
    )
    def test_normalize_text(self, validator, input_text, expected):
        """Test text normalization with various inputs."""
        normalized = validator.normalize_text(input_text)
        assert normalized == expected

    @pytest.mark.parametrize(
        "reference,expected_pmid",
        [
            ("PMID:12345678", "12345678"),
            ("GO_REF:0000033", None),
            ("PMID:99999999", "99999999"),
            ("", None),
            ("Invalid", None),
            ("pmid:12345", "12345"),  # lowercase pmid
        ],
    )
    def test_extract_pmid(self, validator, reference, expected_pmid):
        """Test PMID extraction from various reference formats."""
        pmid = validator.extract_pmid_from_reference(reference)
        assert pmid == expected_pmid

    @pytest.mark.parametrize(
        "text,publication,should_find,min_score",
        [
            # Exact match
            (
                "JAK1 is essential for cytokine signaling",
                "JAK1 is essential for cytokine signaling. It phosphorylates STAT proteins.",
                True,
                1.0,
            ),
            # Partial match with high similarity (will find because threshold is 0.85 by default)
            (
                "JAK1 is critical for cytokine signaling",
                "JAK1 is essential for cytokine signaling. It phosphorylates STAT proteins.",
                True,
                0.85,
            ),  # Should find with high similarity
            # No match
            (
                "This text does not exist",
                "JAK1 is essential for cytokine signaling.",
                False,
                0.0,
            ),
            # Case insensitive match
            (
                "jak1 is ESSENTIAL for cytokine signaling",
                "JAK1 is essential for cytokine signaling.",
                True,
                1.0,
            ),
        ],
    )
    def test_find_text_in_publication(
        self, validator, text, publication, should_find, min_score
    ):
        """Test finding text in publication content."""
        found, score, match = validator.find_text_in_publication(text, publication)
        assert found == should_find
        if should_find:
            assert score >= min_score

    @pytest.mark.parametrize(
        "text,publication,should_find,min_score",
        [
            # Single bracketed editorial note
            (
                "JAK1 is a [non-receptor] tyrosine kinase",
                "JAK1 is a tyrosine kinase that binds receptors.",
                True,
                1.0,
            ),
            # Multiple bracketed sections
            (
                "[The protein] CFAP300 shows [direct] interaction with DNAAF2 [cytoplasmic] assembly factor",
                "CFAP300 shows interaction with DNAAF2 assembly factor in cells.",
                True,
                1.0,
            ),
            # Brackets at beginning and end
            (
                "[According to studies,] CFAP300 is essential for dynein assembly [in cilia]",
                "CFAP300 is essential for dynein assembly throughout the cell.",
                True,
                1.0,
            ),
            # Empty brackets should be handled
            (
                "The protein [] functions in [] the cytoplasm",
                "The protein functions in the cytoplasm normally.",
                True,
                1.0,
            ),
            # Nested brackets (though not recommended)
            (
                "CFAP300 [a protein [also known as C11orf70]] binds DNAAF2",
                "CFAP300 binds DNAAF2 directly.",
                True,
                1.0,
            ),
            # Text that doesn't match even after removing brackets
            (
                "[Editorial note:] This protein does not exist",
                "CFAP300 is a real protein with important functions.",
                False,
                0.0,
            ),
        ],
    )
    def test_square_brackets_ignored(
        self, validator, text, publication, should_find, min_score
    ):
        """Test that text in square brackets is properly ignored during matching."""
        found, score, match = validator.find_text_in_publication(text, publication)
        assert found == should_find
        if should_find:
            assert score >= min_score

    @pytest.mark.parametrize(
        "text,publication,should_find,min_score",
        [
            # Simple two-part quote
            (
                "protein functions in cytoplasm ... helps with assembly",
                "The protein functions in cytoplasm and does other things. It helps with assembly of motors.",
                True,
                1.0,
            ),
            # Three-part quote
            (
                "CFAP300 localizes ... interacts with DNAAF2 ... essential for assembly",
                "CFAP300 localizes to cytoplasm. In cells, it interacts with DNAAF2 directly. Studies show it is essential for assembly.",
                True,
                1.0,
            ),
            # Multi-part with one part not found
            (
                "protein is important ... this text doesn't exist ... helps assembly",
                "The protein is important for function. Other things happen. It helps assembly.",
                False,
                0.0,
            ),
            # Empty parts should be ignored
            (
                "protein functions ... ... ... helps assembly",
                "The protein functions normally. It helps assembly.",
                True,
                1.0,
            ),
        ],
    )
    def test_ellipsis_multi_part_quotes(
        self, validator, text, publication, should_find, min_score
    ):
        """Test multi-part quotes separated by ellipsis."""
        found, score, match = validator.find_text_in_publication(text, publication)
        assert found == should_find
        if should_find:
            assert score >= min_score

    @pytest.mark.parametrize(
        "text,publication,should_find,min_score",
        [
            # Brackets and ellipsis combined
            (
                "[The] protein functions [mainly] in cytoplasm ... [It] helps with assembly [of motors]",
                "The protein functions in cytoplasm primarily. Later, it helps with assembly process.",
                True,
                1.0,
            ),
            # Complex real-world example
            (
                "C11orf70 [also known as CFAP300] is essential for ... assembly of dynein arms ... [and] mutations cause [severe] PCD",
                "C11orf70 is essential for proper function. The assembly of dynein arms is critical. Additionally, mutations cause PCD.",
                True,
                1.0,
            ),
            # Multiple brackets within ellipsis parts - note: commas remain after bracket removal
            (
                "[According to recent] studies, CFAP300 binds [directly to] DNAAF2 ... [Furthermore,] expression increases [significantly] during ciliogenesis",
                "Recent studies, CFAP300 binds DNAAF2 in cells. Expression increases during ciliogenesis.",
                True,
                0.85,
            ),  # Lower score expected due to remaining punctuation
            # Interleaved brackets and ellipsis
            (
                "The [important] protein ... [which is] CFAP300 ... shows [strong] interaction",
                "The protein is critical. CFAP300 has many roles. It shows interaction with partners.",
                True,
                1.0,
            ),
            # Edge case: brackets containing ellipsis (should remove entire bracket content)
            (
                "CFAP300 [... unclear mechanism ...] binds to DNAAF2",
                "CFAP300 binds to DNAAF2 directly.",
                True,
                1.0,
            ),
        ],
    )
    def test_complex_brackets_and_ellipsis(
        self, validator, text, publication, should_find, min_score
    ):
        """Test complex combinations of brackets and ellipsis."""
        found, score, match = validator.find_text_in_publication(text, publication)
        assert found == should_find
        if should_find:
            assert score >= min_score

    @pytest.mark.parametrize(
        "text,expected_count",
        [
            (
                "This is a long first sentence. This is a long second sentence! This is a long third sentence? Fourth short",
                3,
            ),  # Only sentences > 20 chars
            (
                "Line one\nLine two\nLine three that is longer than twenty characters.",
                1,
            ),  # Only long sentences
            ("Short.", 0),  # Too short
            ("", 0),
            ("One very long sentence that definitely exceeds twenty characters.", 1),
        ],
    )
    def test_split_sentences(self, validator, text, expected_count):
        """Test sentence splitting with various inputs."""
        sentences = validator.split_into_sentences(text)
        assert len(sentences) == expected_count

    def test_validate_with_mock_publication(self, temp_publications_dir):
        """Test validation with a mock publication file."""
        # Create test YAML with matching and non-matching supporting_text
        yaml_data = {
            "id": "P23458",
            "gene_symbol": "JAK1",
            "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
            "description": "Test gene",
            "existing_annotations": [
                {
                    "term": {"id": "GO:0000001", "label": "test1"},
                    "evidence_type": "IEA",
                    "original_reference_id": "PMID:12345678",
                    "review": {
                        "summary": "Test",
                        "action": "ACCEPT",
                        "reason": "Good annotation (PMID:12345678)",
                        "supported_by": [
                            {
                                "reference_id": "PMID:12345678",
                                "supporting_text": "JAK1 is essential for Type I and Type II interferon signaling",
                            }
                        ],
                    },
                },
                {
                    "term": {"id": "GO:0000002", "label": "test2"},
                    "evidence_type": "IEA",
                    "original_reference_id": "PMID:12345678",
                    "review": {
                        "summary": "Test",
                        "action": "ACCEPT",
                        "reason": "Another annotation (PMID:12345678)",
                        "supported_by": [
                            {
                                "reference_id": "PMID:12345678",
                                "supporting_text": "This text does not appear in the publication at all",
                            }
                        ],
                    },
                },
            ],
        }

        yaml_file = temp_publications_dir / "test.yaml"
        with open(yaml_file, "w") as f:
            yaml.dump(yaml_data, f)

        # Validate
        validator = SupportingTextValidator(publications_dir=temp_publications_dir)
        report = validator.validate_file(yaml_file)

        # Check results
        assert report.total_annotations == 2
        assert report.annotations_with_supporting_text == 2
        assert report.valid_supporting_texts == 1  # First one should be valid
        assert report.invalid_supporting_texts == 1  # Second one should be invalid

        # Check specific results
        valid_results = [r for r in report.results if r.is_valid]
        invalid_results = [r for r in report.results if not r.is_valid]

        assert len(valid_results) == 1
        assert len(invalid_results) == 1

        valid_result = valid_results[0]
        assert valid_result.found_in_publication
        assert valid_result.similarity_score >= 0.85

        invalid_result = invalid_results[0]
        assert not invalid_result.found_in_publication
        assert invalid_result.error_message is not None

    @pytest.mark.parametrize(
        "yaml_annotations,expected_valid,expected_invalid",
        [
            # All valid
            (
                [
                    {
                        "term": {"id": "GO:0000001", "label": "test"},
                        "evidence_type": "IEA",
                        "original_reference_id": "PMID:12345678",
                        "review": {
                            "summary": "Test",
                            "action": "ACCEPT",
                            "reason": "Test",
                            "supported_by": [
                                {
                                    "reference_id": "PMID:12345678",
                                    "supporting_text": "JAK1 is essential for Type I and Type II interferon signaling",
                                }
                            ],
                        },
                    }
                ],
                1,
                0,
            ),
            # All invalid
            (
                [
                    {
                        "term": {"id": "GO:0000001", "label": "test"},
                        "evidence_type": "IEA",
                        "original_reference_id": "PMID:12345678",
                        "review": {
                            "summary": "Test",
                            "action": "ACCEPT",
                            "reason": "Test",
                            "supported_by": [
                                {
                                    "reference_id": "PMID:12345678",
                                    "supporting_text": "This text is not in the publication",
                                }
                            ],
                        },
                    }
                ],
                0,
                1,
            ),
            # Mixed
            (
                [
                    {
                        "term": {"id": "GO:0000001", "label": "test1"},
                        "evidence_type": "IEA",
                        "original_reference_id": "PMID:12345678",
                        "review": {
                            "summary": "Test",
                            "action": "ACCEPT",
                            "reason": "Test",
                            "supported_by": [
                                {
                                    "reference_id": "PMID:12345678",
                                    "supporting_text": "JAK1 is a non-receptor tyrosine kinase",
                                }
                            ],
                        },
                    },
                    {
                        "term": {"id": "GO:0000002", "label": "test2"},
                        "evidence_type": "IEA",
                        "original_reference_id": "PMID:12345678",
                        "review": {
                            "summary": "Test",
                            "action": "ACCEPT",
                            "reason": "Test",
                            "supported_by": [
                                {
                                    "reference_id": "PMID:12345678",
                                    "supporting_text": "Not in publication",
                                }
                            ],
                        },
                    },
                ],
                1,
                1,
            ),
        ],
    )
    def test_validation_scenarios(
        self, temp_publications_dir, yaml_annotations, expected_valid, expected_invalid
    ):
        """Test various validation scenarios."""
        yaml_data = {
            "id": "P23458",
            "gene_symbol": "JAK1",
            "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
            "description": "Test gene",
            "existing_annotations": yaml_annotations,
        }

        yaml_file = temp_publications_dir / "test.yaml"
        with open(yaml_file, "w") as f:
            yaml.dump(yaml_data, f)

        validator = SupportingTextValidator(publications_dir=temp_publications_dir)
        report = validator.validate_file(yaml_file)

        assert report.valid_supporting_texts == expected_valid
        assert report.invalid_supporting_texts == expected_invalid

    def test_findings_validation(self, temp_publications_dir):
        """Test validation of supporting_text in references.findings."""
        yaml_data = {
            "id": "P23458",
            "gene_symbol": "JAK1",
            "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
            "description": "Test gene",
            "references": [
                {
                    "id": "PMID:12345678",
                    "title": "Test Publication",
                    "findings": [
                        {
                            "statement": "JAK1 has FERM domain",
                            "supporting_text": "It contains four major domains including a FERM domain and a kinase domain",
                        },
                        {
                            "statement": "JAK1 is involved in signaling",
                            "supporting_text": "This text does not exist in the publication",
                        },
                    ],
                }
            ],
            "existing_annotations": [],
        }

        yaml_file = temp_publications_dir / "test.yaml"
        with open(yaml_file, "w") as f:
            yaml.dump(yaml_data, f)

        validator = SupportingTextValidator(publications_dir=temp_publications_dir)
        report = validator.validate_file(yaml_file)

        # Should find 2 findings with supporting_text
        assert report.annotations_with_supporting_text == 2
        assert report.valid_supporting_texts == 1  # First finding is valid
        assert report.invalid_supporting_texts == 1  # Second finding is invalid

        # Check specific results
        valid_results = [r for r in report.results if r.is_valid]
        invalid_results = [r for r in report.results if not r.is_valid]

        assert len(valid_results) == 1
        assert len(invalid_results) == 1

        # Check paths are correct
        assert "references[0].findings[0]" in valid_results[0].annotation_path
        assert "references[0].findings[1]" in invalid_results[0].annotation_path

    def test_mixed_findings_and_annotations(self, temp_publications_dir):
        """Test validation with both findings and annotations having supporting_text."""
        yaml_data = {
            "id": "P23458",
            "gene_symbol": "JAK1",
            "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
            "description": "Test gene",
            "references": [
                {
                    "id": "PMID:12345678",
                    "title": "Test Publication",
                    "findings": [
                        {
                            "statement": "JAK1 has FERM domain",
                            "supporting_text": "It contains four major domains including a FERM domain",
                        }
                    ],
                }
            ],
            "existing_annotations": [
                {
                    "term": {"id": "GO:0000001", "label": "test"},
                    "evidence_type": "IEA",
                    "original_reference_id": "PMID:12345678",
                    "review": {
                        "action": "ACCEPT",
                        "reason": "Test",
                        "supported_by": [
                            {
                                "reference_id": "PMID:12345678",
                                "supporting_text": "JAK1 is essential for Type I and Type II interferon signaling",
                            }
                        ],
                    },
                }
            ],
        }

        yaml_file = temp_publications_dir / "test.yaml"
        with open(yaml_file, "w") as f:
            yaml.dump(yaml_data, f)

        validator = SupportingTextValidator(publications_dir=temp_publications_dir)
        report = validator.validate_file(yaml_file)

        # Should have 2 total (1 finding + 1 annotation)
        assert report.total_annotations == 2
        assert report.annotations_with_supporting_text == 2
        assert report.valid_supporting_texts == 2  # Both should be valid
        assert report.invalid_supporting_texts == 0


    def test_uniprot_supporting_text_validation(self, temp_publications_dir):
        """Test validation of supporting_text from UniProt files."""
        # Create a mock UniProt file for gene A0B297
        gene_dir = temp_publications_dir.parent / "genes" / "BURCH" / "A0B297"
        gene_dir.mkdir(parents=True, exist_ok=True)
        
        uniprot_file = gene_dir / "A0B297-uniprot.txt"
        uniprot_content = """ID   RGMG2_BURCH             Reviewed;         512 AA.
AC   A0B297;
CC   -!- FUNCTION: Part of an ABC transporter complex involved in carbohydrate
CC       import. Could be involved in ribose, galactose and/or methyl
CC       galactoside import. Responsible for energy coupling to the transport
CC       system. {ECO:0000255|HAMAP-Rule:MF_01717}.
CC   -!- CATALYTIC ACTIVITY:
CC       Reaction=D-ribose(out) + ATP + H2O = D-ribose(in) + ADP + phosphate +
CC         H(+); Xref=Rhea:RHEA:29903, ChEBI:CHEBI:15377, ChEBI:CHEBI:15378,
CC         ChEBI:CHEBI:30616, ChEBI:CHEBI:43474, ChEBI:CHEBI:47013,
CC         ChEBI:CHEBI:456216; EC=7.5.2.7; Evidence={ECO:0000255|HAMAP-
CC         Rule:MF_01717};
CC   -!- SUBCELLULAR LOCATION: Cell inner membrane {ECO:0000255|HAMAP-
CC       Rule:MF_01717}; Peripheral membrane protein {ECO:0000255|HAMAP-
CC       Rule:MF_01717}.
CC   -!- SIMILARITY: Belongs to the ABC transporter superfamily. Carbohydrate
CC       importer 2 (CUT2) (TC 3.A.1.2) family. {ECO:0000255|HAMAP-
CC       Rule:MF_01717}.
"""
        uniprot_file.write_text(uniprot_content)
        
        # Create test YAML with UniProt reference and findings
        yaml_data = {
            "id": "A0B297",
            "gene_symbol": "A0B297",
            "taxon": {"id": "NCBITaxon:331272", "label": "Burkholderia cenocepacia"},
            "description": "Test gene",
            "references": [
                {
                    "id": "uniprot:A0B297",
                    "title": "UniProt Entry",
                    "findings": [
                        {
                            "statement": "Catalyzes ATP-coupled sugar transport",
                            "supporting_text": "EC=7.5.2.11 (D-galactose transport) and EC=7.5.2.7 (D-ribose transport)"
                        },
                        {
                            "statement": "Peripheral membrane protein at inner membrane",
                            "supporting_text": "SUBCELLULAR LOCATION: Cell inner membrane; Peripheral membrane protein"
                        },
                        {
                            "statement": "Member of CUT2 family",
                            "supporting_text": "Belongs to the ABC transporter superfamily. Carbohydrate importer 2 (CUT2) (TC 3.A.1.2) family"
                        },
                        {
                            "statement": "TEMPORARY TEST STATEMENT",
                            "supporting_text": "TEMPORARY TEST SUPPORTING TEXT"  # This should fail validation
                        }
                    ]
                }
            ]
        }
        
        yaml_file = temp_publications_dir / "test_uniprot.yaml"
        with open(yaml_file, "w") as f:
            yaml.dump(yaml_data, f)
        
        # Validate
        validator = SupportingTextValidator(publications_dir=temp_publications_dir, gene_dir=gene_dir.parent.parent)
        report = validator.validate_file(yaml_file)
        
        # Check results
        assert report.total_annotations == 4
        assert report.annotations_with_supporting_text == 4
        # With stricter validation, short texts require exact matches
        # "EC=7.5.2.11 (D-galactose transport) and EC=7.5.2.7 (D-ribose transport)" - not found exactly
        # "SUBCELLULAR LOCATION: Cell inner membrane; Peripheral membrane protein" - not found exactly  
        # "Belongs to the ABC transporter superfamily. Carbohydrate importer 2 (CUT2) (TC 3.A.1.2) family" - found
        # "TEMPORARY TEST SUPPORTING TEXT" - not found
        assert report.valid_supporting_texts == 1  # Only the third one is valid
        assert report.invalid_supporting_texts == 3  # Three are invalid due to strict matching
        assert not report.is_valid  # Overall should be invalid
        
        # Check specific error messages
        invalid_results = [r for r in report.results if not r.is_valid]
        assert len(invalid_results) == 3  # Three invalid results
        # Check that the TEMPORARY TEST one is among them
        temp_test_results = [r for r in invalid_results if "TEMPORARY TEST" in r.supporting_text]
        assert len(temp_test_results) == 1


@pytest.mark.integration
class TestIntegration:
    """Integration tests that use real files."""

    def test_validate_real_file(self):
        """Test validation on a real gene review file if it exists."""
        jak1_file = Path("genes/human/JAK1/JAK1-ai-review.yaml")
        if not jak1_file.exists():
            pytest.skip("JAK1 file not found")

        report = validate_supporting_text_in_file(jak1_file)

        # Should have some results
        assert report.total_annotations > 0

        # Validate the report structure
        assert hasattr(report, "validation_rate")
        assert hasattr(report, "accuracy_rate")
        assert 0 <= report.validation_rate <= 100
        assert 0 <= report.accuracy_rate <= 100
