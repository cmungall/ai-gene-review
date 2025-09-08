"""Tests for GOA annotation validator."""

import pytest
from pathlib import Path
import tempfile
import yaml

from ai_gene_review.validation.goa_validator import GOAValidator, GOAAnnotation


@pytest.fixture
def sample_goa_tsv():
    """Create a sample GOA TSV file for testing."""
    # Add header row matching actual GOA format
    content = """GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE
UniProtKB	Q9BRQ4	CFAP300		GO:0005515	protein binding	MF	ECO:0000353	IPI	PMID:29727692	UniProtKB:Q9NVR5	NCBITaxon:9606	Homo sapiens	UniProt	Cilia- and flagella-associated protein 300	20180515
UniProtKB	Q9BRQ4	CFAP300		GO:0005737	cytoplasm	CC	ECO:0000250	ISS	PMID:29727692		NCBITaxon:9606	Homo sapiens	UniProt	Cilia- and flagella-associated protein 300	20180515
UniProtKB	Q9BRQ4	CFAP300		GO:0005856	cytoskeleton	CC	ECO:0007322	IEA	PMID:29727692		NCBITaxon:9606	Homo sapiens	UniProt	Cilia- and flagella-associated protein 300	20180515
UniProtKB	Q9BRQ4	CFAP300		GO:0031514	motile cilium	CC	ECO:0000250	ISS	PMID:29727692		NCBITaxon:9606	Homo sapiens	UniProt	Cilia- and flagella-associated protein 300	20180515
UniProtKB	Q9BRQ4	CFAP300		GO:0042995	cell projection	CC	ECO:0007322	IEA	PMID:29727692		NCBITaxon:9606	Homo sapiens	UniProt	Cilia- and flagella-associated protein 300	20180515"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(content)
        return Path(f.name)


@pytest.fixture
def sample_yaml_matching():
    """Create a YAML file with annotations matching the GOA file."""
    data = {
        "id": "Q9BRQ4",
        "gene_symbol": "CFAP300",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0005856", "label": "cytoskeleton"},
                "evidence_type": "IEA",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0031514", "label": "motile cilium"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0042995", "label": "cell projection"},
                "evidence_type": "IEA",
                "original_reference_id": "PMID:29727692",
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        return Path(f.name)


@pytest.fixture
def sample_yaml_missing():
    """Create a YAML file missing some annotations from GOA."""
    data = {
        "id": "Q9BRQ4",
        "gene_symbol": "CFAP300",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:29727692",
            },
            # Missing GO:0005856, GO:0031514, GO:0042995
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        return Path(f.name)


@pytest.fixture
def sample_yaml_extra():
    """Create a YAML file with extra annotations not in GOA."""
    data = {
        "id": "Q9BRQ4",
        "gene_symbol": "CFAP300",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0005856", "label": "cytoskeleton"},
                "evidence_type": "IEA",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0031514", "label": "motile cilium"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0042995", "label": "cell projection"},
                "evidence_type": "IEA",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0099999", "label": "fake annotation"},
                "evidence_type": "IMP",
                "original_reference_id": "PMID:12345678",
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        return Path(f.name)


@pytest.fixture
def sample_yaml_label_mismatch():
    """Create a YAML file with mismatched labels."""
    data = {
        "id": "Q9BRQ4",
        "gene_symbol": "CFAP300",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene",
        "existing_annotations": [
            {
                "term": {
                    "id": "GO:0005515",
                    "label": "wrong label",
                },  # Should be 'protein binding'
                "evidence_type": "IPI",
                "original_reference_id": "PMID:29727692",
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:29727692",
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(data, f)
        return Path(f.name)


def test_parse_goa_annotation():
    """Test parsing a single GOA annotation from TSV row."""
    # Row format matching actual GOA file structure:
    # 0: DB, 1: ID, 2: Symbol, 3: Qualifier, 4: GO Term, 5: GO Name
    # 6: GO Aspect, 7: ECO ID, 8: GO Evidence Code, 9: Reference, etc.
    row = [
        "UniProtKB",
        "Q9BRQ4",
        "CFAP300",
        "",
        "GO:0005515",
        "protein binding",
        "MF",
        "ECO:0000353",
        "IPI",
        "PMID:29727692",
        "UniProtKB:Q9NVR5",
        "NCBITaxon:9606",
        "Homo sapiens",
        "UniProt",
        "GO_REF:0000024",
        "20180515",
    ]

    ann = GOAAnnotation.from_tsv_row(row)

    assert ann.database == "UniProtKB"
    assert ann.db_object_id == "Q9BRQ4"
    assert ann.db_object_symbol == "CFAP300"
    assert ann.go_id == "GO:0005515"
    assert ann.go_term == "protein binding"
    assert ann.go_aspect == "MF"
    assert ann.evidence_code == "IPI"  # From index 8
    assert ann.evidence_type == "ECO:0000353"  # From index 7
    assert ann.reference == "PMID:29727692"
    assert ann.taxon_id == "NCBITaxon:9606"


def test_parse_goa_file(sample_goa_tsv):
    """Test parsing a GOA TSV file."""
    validator = GOAValidator()
    annotations = validator.parse_goa_file(sample_goa_tsv)

    assert len(annotations) == 5
    assert annotations[0].go_id == "GO:0005515"
    assert annotations[1].go_id == "GO:0005737"
    assert annotations[2].go_id == "GO:0005856"
    assert annotations[3].go_id == "GO:0031514"
    assert annotations[4].go_id == "GO:0042995"

    # Clean up
    sample_goa_tsv.unlink()


def test_validate_matching_annotations(sample_goa_tsv, sample_yaml_matching):
    """Test validation when YAML and GOA match perfectly."""
    validator = GOAValidator()
    result = validator.validate_against_goa(sample_yaml_matching, sample_goa_tsv)

    assert result.is_valid
    assert len(result.missing_in_yaml) == 0
    assert len(result.missing_in_goa) == 0
    assert len(result.mismatched_labels) == 0
    assert len(result.mismatched_evidence) == 0

    # Clean up
    sample_goa_tsv.unlink()
    sample_yaml_matching.unlink()


def test_validate_missing_in_yaml(sample_goa_tsv, sample_yaml_missing):
    """Test validation when YAML is missing annotations from GOA."""
    validator = GOAValidator(strict_tuple_matching=False)  # Use legacy mode
    result = validator.validate_against_goa(sample_yaml_missing, sample_goa_tsv)

    assert not result.is_valid
    assert len(result.missing_in_yaml) == 3

    missing_go_ids = {ann.go_id for ann in result.missing_in_yaml}
    assert "GO:0005856" in missing_go_ids
    assert "GO:0031514" in missing_go_ids
    assert "GO:0042995" in missing_go_ids

    assert len(result.missing_in_goa) == 0
    assert len(result.mismatched_labels) == 0

    # Clean up
    sample_goa_tsv.unlink()
    sample_yaml_missing.unlink()


def test_validate_extra_in_yaml(sample_goa_tsv, sample_yaml_extra):
    """Test validation when YAML has extra annotations not in GOA."""
    validator = GOAValidator()
    result = validator.validate_against_goa(sample_yaml_extra, sample_goa_tsv)

    assert not result.is_valid
    assert len(result.missing_in_yaml) == 0
    assert len(result.missing_in_goa) == 1

    extra_ann = result.missing_in_goa[0]
    assert extra_ann["term"]["id"] == "GO:0099999"
    assert extra_ann["term"]["label"] == "fake annotation"

    # Clean up
    sample_goa_tsv.unlink()
    sample_yaml_extra.unlink()


def test_validate_label_mismatch(sample_goa_tsv, sample_yaml_label_mismatch):
    """Test validation when YAML has mismatched labels.
    
    Note: Label validation is now handled by term_validator.py against the ontology,
    not by GOA validator. GOA files may use synonyms instead of primary labels.
    This test verifies that annotations missing from YAML are still properly detected.
    """
    validator = GOAValidator()
    result = validator.validate_against_goa(sample_yaml_label_mismatch, sample_goa_tsv)

    # The YAML only has 2 annotations while GOA has 5, so 3 are missing
    # This is the critical check - ensuring all GOA annotations are covered
    assert not result.is_valid  # Should be invalid due to missing annotations
    assert len(result.missing_in_yaml) == 3  # GO:0005856, GO:0031514, GO:0042995 are missing
    
    # Verify the specific missing annotations
    missing_go_ids = {ann.go_id for ann in result.missing_in_yaml}
    assert "GO:0005856" in missing_go_ids  # cytoskeleton
    assert "GO:0031514" in missing_go_ids  # motile cilium  
    assert "GO:0042995" in missing_go_ids  # cell projection
    
    # Label mismatches are not validated by GOA validator anymore
    assert len(result.mismatched_labels) == 0  # Labels are validated against ontology, not GOA

    # Clean up
    sample_goa_tsv.unlink()
    sample_yaml_label_mismatch.unlink()


def test_validation_summary(sample_goa_tsv, sample_yaml_missing):
    """Test getting a human-readable summary of validation results."""
    validator = GOAValidator(strict_tuple_matching=False)  # Use legacy mode
    result = validator.validate_against_goa(sample_yaml_missing, sample_goa_tsv)

    summary = validator.get_summary(result)

    assert "✗ Annotations do not match GOA file" in summary
    assert "3 annotations in GOA but not in YAML" in summary
    assert "GO:0005856" in summary
    assert "cytoskeleton" in summary

    # Clean up
    sample_goa_tsv.unlink()
    sample_yaml_missing.unlink()


def test_derive_goa_path():
    """Test deriving GOA file path from YAML file path."""
    validator = GOAValidator()

    # Create a test YAML file with expected naming pattern
    yaml_path = Path("genes/human/CFAP300/CFAP300-ai-review.yaml")

    # We can't actually validate without the files existing,
    # but we can test the path derivation logic
    with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as yaml_file:
        yaml_file.write(b"id: test\n")
        yaml_path = Path(yaml_file.name)

    # The validator should handle missing GOA file gracefully
    result = validator.validate_against_goa(yaml_path)
    assert not result.is_valid
    assert result.error_message is not None

    # Clean up
    yaml_path.unlink()


def test_empty_goa_file():
    """Test handling of empty GOA file."""
    validator = GOAValidator()

    # Create GOA file with only header (no data rows)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as goa_file:
        goa_file.write(
            "GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE\n"
        )
        goa_path = Path(goa_file.name)

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False
    ) as yaml_file:
        yaml.dump({"id": "test", "existing_annotations": []}, yaml_file)
        yaml_path = Path(yaml_file.name)

    result = validator.validate_against_goa(yaml_path, goa_path)

    assert result.is_valid
    assert len(result.missing_in_yaml) == 0
    assert len(result.missing_in_goa) == 0

    # Clean up
    goa_path.unlink()
    yaml_path.unlink()


def test_evidence_type_mismatch():
    """Test validation when evidence types don't match."""
    # Create GOA with one evidence type (with header)
    goa_content = """GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE
UniProtKB	Q9BRQ4	CFAP300		GO:0005515	protein binding	MF	ECO:0000353	IPI	PMID:29727692	UniProtKB:Q9NVR5	NCBITaxon:9606	Homo sapiens	UniProt	Cilia- and flagella-associated protein 300	20180515"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(goa_content)
        goa_path = Path(f.name)

    # Create YAML with different evidence type
    yaml_data = {
        "id": "Q9BRQ4",
        "gene_symbol": "CFAP300",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "ISS",  # Different from IPI in GOA
                "original_reference_id": "PMID:29727692",
            }
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump(yaml_data, f)
        yaml_path = Path(f.name)

    # Test with legacy mode - evidence mismatch is not a hard failure
    validator = GOAValidator(strict_tuple_matching=False)
    result = validator.validate_against_goa(yaml_path, goa_path)
    assert result.is_valid  # Evidence mismatch doesn't fail in legacy mode

    # Test with strict tuple matching - should fail
    validator_strict = GOAValidator(strict_tuple_matching=True)
    result_strict = validator_strict.validate_against_goa(yaml_path, goa_path)
    assert not result_strict.is_valid  # Should fail with strict tuple matching
    assert len(result.mismatched_evidence) == 1

    go_id, yaml_ev, goa_ev = result.mismatched_evidence[0]
    assert go_id == "GO:0005515"
    assert yaml_ev == "ISS"
    assert "IPI" in goa_ev

    # In strict mode, evidence mismatch is a failure
    validator.strict_mode = True
    result = validator.validate_against_goa(yaml_path, goa_path)
    assert not result.is_valid

    # Clean up
    goa_path.unlink()
    yaml_path.unlink()


@pytest.mark.parametrize(
    "row_length,expected_fields",
    [
        (0, 0),  # Empty row
        (5, 5),  # Partial row
        (16, 16),  # Full row
        (20, 16),  # Extra columns (should handle gracefully)
    ],
)
def test_parse_goa_row_edge_cases(row_length, expected_fields):
    """Test parsing GOA rows with various lengths."""
    row = ["field" + str(i) for i in range(row_length)]
    ann = GOAAnnotation.from_tsv_row(row)

    # Should handle any row length gracefully
    assert ann is not None

    # Check that fields are populated or empty as expected
    if row_length > 4:
        assert ann.go_id == "field4"
    else:
        assert ann.go_id == ""


def test_seed_missing_annotations():
    """Test seeding missing annotations from GOA into YAML."""
    validator = GOAValidator()

    # Create a GOA file with annotations
    goa_content = """GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE
UniProtKB	Q12345	TEST		GO:0001234	test function	MF	ECO:0000353	IPI	PMID:12345		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515
UniProtKB	Q12345	TEST		GO:0005678	test process	BP	ECO:0000250	ISS	PMID:67890		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515
UniProtKB	Q12345	TEST		GO:0009999	test location	CC	ECO:0007322	IEA	PMID:99999		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as goa_file:
        goa_file.write(goa_content)
        goa_path = Path(goa_file.name)

    # Create a YAML file with only one annotation (missing two)
    yaml_data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "description": "Test gene",
        "existing_annotations": [
            {
                "term": {"id": "GO:0001234", "label": "test function"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:12345",
                "review": {"summary": "Already reviewed", "action": "ACCEPT"},
            }
        ],
    }

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False
    ) as yaml_file:
        yaml.dump(yaml_data, yaml_file)
        yaml_path = Path(yaml_file.name)

    try:
        # Seed missing annotations
        added_count, output_path, refs_added = validator.seed_missing_annotations(
            yaml_path, goa_path
        )

        # Should have added 2 annotations
        assert added_count == 2
        assert output_path == yaml_path

        # Load and verify the updated YAML
        with open(yaml_path, "r") as f:
            updated_data = yaml.safe_load(f)

        # Should now have 3 annotations
        assert len(updated_data["existing_annotations"]) == 3

        # Check that original annotation is unchanged
        original_ann = updated_data["existing_annotations"][0]
        assert original_ann["term"]["id"] == "GO:0001234"
        assert original_ann["review"]["summary"] == "Already reviewed"

        # Check that new annotations were added with TODO review sections
        new_anns = updated_data["existing_annotations"][1:]
        go_ids = {ann["term"]["id"] for ann in new_anns}
        assert "GO:0005678" in go_ids
        assert "GO:0009999" in go_ids

        # New annotations should have TODO review sections
        for ann in new_anns:
            assert "review" in ann
            assert ann["review"]["summary"] == "TODO: Review this GOA annotation"
            assert ann["review"]["action"] == "PENDING"

    finally:
        # Clean up
        goa_path.unlink()
        yaml_path.unlink()


def test_new_action_validation():
    """Test validation of annotations with action=NEW."""
    validator = GOAValidator()
    
    # Create GOA file with existing annotations
    goa_content = """GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE
UniProtKB	Q12345	TEST		GO:0005515	protein binding	MF	ECO:0000353	IPI	PMID:12345		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515
UniProtKB	Q12345	TEST		GO:0005737	cytoplasm	CC	ECO:0000250	ISS	PMID:67890		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515"""
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as goa_file:
        goa_file.write(goa_content)
        goa_path = Path(goa_file.name)
    
    # Test 1: NEW annotation that does NOT exist in GOA (valid)
    yaml_data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:12345",
                "review": {"summary": "Existing annotation", "action": "ACCEPT"},
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:67890",
                "review": {"summary": "Existing annotation", "action": "ACCEPT"},
            },
            {
                "term": {"id": "GO:0099999", "label": "new function"},
                "evidence_type": "IMP",
                "original_reference_id": "PMID:99999",
                "review": {"summary": "Proposed new annotation", "action": "NEW"},
            },
        ],
    }
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as yaml_file:
        yaml.dump(yaml_data, yaml_file)
        yaml_path = Path(yaml_file.name)
    
    # Should be valid - NEW annotation doesn't exist in GOA
    result = validator.validate_against_goa(yaml_path, goa_path)
    assert result.is_valid, f"NEW annotation for non-existing GO term should be valid: {result.error_message}"
    
    # Clean up
    yaml_path.unlink()
    
    # Test 2: NEW annotation that DOES exist in GOA (invalid)
    yaml_data_invalid = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:67890",
                "review": {"summary": "This already exists", "action": "NEW"},
            },
        ],
    }
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as yaml_file:
        yaml.dump(yaml_data_invalid, yaml_file)
        yaml_path = Path(yaml_file.name)
    
    # Should be invalid - NEW annotation exists in GOA
    result = validator.validate_against_goa(yaml_path, goa_path)
    assert not result.is_valid, "NEW annotation for existing GO term should be invalid"
    assert "action=NEW exists in GOA" in result.error_message
    
    # Clean up
    yaml_path.unlink()
    goa_path.unlink()


def test_new_action_exclusion_from_validation():
    """Test that NEW annotations are properly excluded from GOA completeness checks."""
    validator = GOAValidator()
    
    # Create GOA file with annotations
    goa_content = """GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE
UniProtKB	Q12345	TEST		GO:0005515	protein binding	MF	ECO:0000353	IPI	PMID:12345		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515
UniProtKB	Q12345	TEST		GO:0005737	cytoplasm	CC	ECO:0000250	ISS	PMID:67890		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515"""
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as goa_file:
        goa_file.write(goa_content)
        goa_path = Path(goa_file.name)
    
    # YAML has all GOA annotations plus a NEW one
    yaml_data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:12345",
                "review": {"summary": "Existing", "action": "ACCEPT"},
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:67890",
                "review": {"summary": "Existing", "action": "ACCEPT"},
            },
            {
                "term": {"id": "GO:0099999", "label": "new function"},
                "evidence_type": "IMP",
                "original_reference_id": "PMID:99999",
                "review": {"summary": "This is a new annotation", "action": "NEW"},
            },
        ],
    }
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as yaml_file:
        yaml.dump(yaml_data, yaml_file)
        yaml_path = Path(yaml_file.name)
    
    # Should be valid - all GOA annotations are present, NEW annotation is ignored
    result = validator.validate_against_goa(yaml_path, goa_path)
    assert result.is_valid, f"Should be valid when all GOA annotations are present and NEW is properly excluded: {result.error_message}"
    assert len(result.missing_in_yaml) == 0, "No GOA annotations should be missing"
    assert len(result.missing_in_goa) == 0, "NEW annotations should not be reported as missing in GOA"
    
    # Clean up
    yaml_path.unlink()
    
    # Now test with a NEW annotation that incorrectly matches GOA
    yaml_data_invalid = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IPI",
                "original_reference_id": "PMID:12345",
                "review": {"summary": "Existing", "action": "ACCEPT"},
            },
            {
                "term": {"id": "GO:0005737", "label": "cytoplasm"},
                "evidence_type": "ISS",
                "original_reference_id": "PMID:67890",
                "review": {"summary": "Should not be NEW", "action": "NEW"},  # This is wrong!
            },
        ],
    }
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as yaml_file:
        yaml.dump(yaml_data_invalid, yaml_file)
        yaml_path = Path(yaml_file.name)
    
    # Should be invalid - NEW annotation exists in GOA
    result = validator.validate_against_goa(yaml_path, goa_path)
    assert not result.is_valid, "Should be invalid when NEW annotation exists in GOA"
    assert "action=NEW exists in GOA" in str(result.error_message)
    
    # Clean up
    yaml_path.unlink()
    goa_path.unlink()


def test_new_action_legacy_mode():
    """Test validation of NEW annotations in legacy mode."""
    validator = GOAValidator(strict_tuple_matching=False)
    
    # Create GOA file
    goa_content = """GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE
UniProtKB	Q12345	TEST		GO:0005515	protein binding	MF	ECO:0000353	IPI	PMID:12345		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515"""
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as goa_file:
        goa_file.write(goa_content)
        goa_path = Path(goa_file.name)
    
    # Test with NEW annotation for existing GO ID
    yaml_data = {
        "id": "Q12345",
        "gene_symbol": "TEST",
        "existing_annotations": [
            {
                "term": {"id": "GO:0005515", "label": "protein binding"},
                "evidence_type": "IMP",  # Different evidence type
                "original_reference_id": "PMID:99999",  # Different reference
                "review": {"summary": "Different evidence", "action": "NEW"},
            },
        ],
    }
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as yaml_file:
        yaml.dump(yaml_data, yaml_file)
        yaml_path = Path(yaml_file.name)
    
    # Should be invalid even in legacy mode - NEW annotation with same GO ID exists in GOA
    result = validator.validate_against_goa(yaml_path, goa_path)
    assert not result.is_valid, "NEW annotation with existing GO ID should be invalid even in legacy mode"
    assert "action=NEW exists in GOA" in result.error_message
    
    # Clean up
    yaml_path.unlink()
    goa_path.unlink()


def test_seed_creates_new_file():
    """Test that seeding can create a new YAML file if it doesn't exist."""
    validator = GOAValidator()

    # Create a GOA file
    goa_content = """GENE PRODUCT DB	GENE PRODUCT ID	SYMBOL	QUALIFIER	GO TERM	GO NAME	GO ASPECT	ECO ID	GO EVIDENCE CODE	REFERENCE	WITH/FROM	TAXON ID	TAXON NAME	ASSIGNED BY	GENE NAME	DATE
UniProtKB	Q12345	TEST		GO:0001234	test function	MF	ECO:0000353	IPI	PMID:12345		NCBITaxon:9606	Homo sapiens	UniProt	Test protein	20180515"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as goa_file:
        goa_file.write(goa_content)
        goa_path = Path(goa_file.name)

    # Use a non-existent YAML path
    yaml_path = Path(tempfile.gettempdir()) / "test-new-file.yaml"

    try:
        # Seed should create the file
        added_count, output_path, refs_added = validator.seed_missing_annotations(
            yaml_path, goa_path
        )

        assert added_count == 1
        assert output_path.exists()

        # Verify the created file
        with open(output_path, "r") as f:
            data = yaml.safe_load(f)

        assert len(data["existing_annotations"]) == 1
        assert data["existing_annotations"][0]["term"]["id"] == "GO:0001234"

    finally:
        # Clean up
        goa_path.unlink()
        if yaml_path.exists():
            yaml_path.unlink()
