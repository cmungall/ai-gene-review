"""Test strict GOA validation where all tuples must match exactly."""

from pathlib import Path
import tempfile
import yaml


from ai_gene_review.validation.goa_validator import GOAValidator


def create_test_goa_file(annotations: list) -> Path:
    """Create a test GOA TSV file with given annotations."""
    with tempfile.NamedTemporaryFile(mode="w", suffix="-goa.tsv", delete=False) as f:
        # Write header
        f.write(
            "DB\tDB_Object_ID\tDB_Object_Symbol\tQualifier\tGO_ID\tGO_Term\tGO_Aspect\tEvidence_Type\tEvidence_Code\tReference\tDate\tAssignedBy\n"
        )

        # Write annotations
        for ann in annotations:
            f.write(
                f"UniProtKB\t{ann['id']}\t{ann['symbol']}\t{ann.get('qualifier', '')}\t"
                f"{ann['go_id']}\t{ann['go_term']}\t{ann.get('aspect', 'P')}\t"
                f"{ann['evidence']}\t{ann['evidence']}\t{ann['reference']}\t"
                f"20240101\tGOA\n"
            )

        return Path(f.name)


def create_test_yaml_file(gene_symbol: str, annotations: list) -> Path:
    """Create a test YAML file with given annotations."""
    data = {
        "gene_symbol": gene_symbol,
        "organism": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        "existing_annotations": annotations,
    }

    with tempfile.NamedTemporaryFile(
        mode="w", suffix="-ai-review.yaml", delete=False
    ) as f:
        yaml.dump(data, f)
        return Path(f.name)


def test_exact_tuple_match_passes():
    """Test that exact tuple matches pass validation."""
    # Create GOA annotations
    goa_anns = [
        {
            "id": "P12345",
            "symbol": "TEST1",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IPI",
            "reference": "PMID:12345",
        },
        {
            "id": "P12345",
            "symbol": "TEST1",
            "go_id": "GO:0005634",
            "go_term": "nucleus",
            "evidence": "IEA",
            "reference": "GO_REF:0000033",
        },
    ]
    goa_file = create_test_goa_file(goa_anns)

    # Create matching YAML annotations
    yaml_anns = [
        {
            "term": {"id": "GO:0005515", "label": "protein binding"},
            "evidence_type": "IPI",
            "original_reference_id": "PMID:12345",
        },
        {
            "term": {"id": "GO:0005634", "label": "nucleus"},
            "evidence_type": "IEA",
            "original_reference_id": "GO_REF:0000033",
        },
    ]
    yaml_file = create_test_yaml_file("TEST1", yaml_anns)

    try:
        validator = GOAValidator()
        result = validator.validate_against_goa(yaml_file, goa_file)

        assert result.is_valid, (
            f"Should be valid but got: {validator.get_summary(result)}"
        )
        assert len(result.missing_in_goa) == 0
        assert len(result.missing_in_yaml) == 0
    finally:
        goa_file.unlink()
        yaml_file.unlink()


def test_different_evidence_type_fails():
    """Test that same GO term with different evidence type fails."""
    # Create GOA annotations
    goa_anns = [
        {
            "id": "P12345",
            "symbol": "TEST2",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IPI",
            "reference": "PMID:12345",
        }
    ]
    goa_file = create_test_goa_file(goa_anns)

    # Create YAML with different evidence type
    yaml_anns = [
        {
            "term": {"id": "GO:0005515", "label": "protein binding"},
            "evidence_type": "IDA",  # Different evidence type
            "original_reference_id": "PMID:12345",
        }
    ]
    yaml_file = create_test_yaml_file("TEST2", yaml_anns)

    try:
        validator = GOAValidator()
        result = validator.validate_against_goa(yaml_file, goa_file)

        assert not result.is_valid
        assert len(result.missing_in_goa) == 1
        # The annotation with IDA evidence is not in GOA
        missing = result.missing_in_goa[0]
        assert missing["evidence_type"] == "IDA"
    finally:
        goa_file.unlink()
        yaml_file.unlink()


def test_different_reference_fails():
    """Test that same GO term with different reference fails."""
    # Create GOA annotations
    goa_anns = [
        {
            "id": "P12345",
            "symbol": "TEST3",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IPI",
            "reference": "PMID:12345",
        }
    ]
    goa_file = create_test_goa_file(goa_anns)

    # Create YAML with different reference
    yaml_anns = [
        {
            "term": {"id": "GO:0005515", "label": "protein binding"},
            "evidence_type": "IPI",
            "original_reference_id": "PMID:67890",  # Different reference
        }
    ]
    yaml_file = create_test_yaml_file("TEST3", yaml_anns)

    try:
        validator = GOAValidator()
        result = validator.validate_against_goa(yaml_file, goa_file)

        assert not result.is_valid
        assert len(result.missing_in_goa) == 1
        # The annotation with different reference is not in GOA
        missing = result.missing_in_goa[0]
        assert missing["original_reference_id"] == "PMID:67890"
    finally:
        goa_file.unlink()
        yaml_file.unlink()


def test_extra_annotation_in_yaml_fails():
    """Test that extra annotations in YAML that aren't in GOA fail."""
    # Create GOA annotations
    goa_anns = [
        {
            "id": "P12345",
            "symbol": "TEST4",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IPI",
            "reference": "PMID:12345",
        }
    ]
    goa_file = create_test_goa_file(goa_anns)

    # Create YAML with extra annotation
    yaml_anns = [
        {
            "term": {"id": "GO:0005515", "label": "protein binding"},
            "evidence_type": "IPI",
            "original_reference_id": "PMID:12345",
        },
        {
            "term": {"id": "GO:0005634", "label": "nucleus"},
            "evidence_type": "IDA",
            "original_reference_id": "PMID:99999",
        },
    ]
    yaml_file = create_test_yaml_file("TEST4", yaml_anns)

    try:
        validator = GOAValidator()
        result = validator.validate_against_goa(yaml_file, goa_file)

        assert not result.is_valid
        assert len(result.missing_in_goa) == 1
        # The extra annotation is not in GOA
        missing = result.missing_in_goa[0]
        assert missing["term"]["id"] == "GO:0005634"
    finally:
        goa_file.unlink()
        yaml_file.unlink()


def test_multiple_annotations_same_go_term():
    """Test that multiple annotations for same GO term with different evidence/refs work."""
    # Create GOA with multiple annotations for same GO term
    goa_anns = [
        {
            "id": "P12345",
            "symbol": "TEST5",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IPI",
            "reference": "PMID:12345",
        },
        {
            "id": "P12345",
            "symbol": "TEST5",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IDA",
            "reference": "PMID:67890",
        },
        {
            "id": "P12345",
            "symbol": "TEST5",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IEA",
            "reference": "GO_REF:0000033",
        },
    ]
    goa_file = create_test_goa_file(goa_anns)

    # Create YAML with matching annotations
    yaml_anns = [
        {
            "term": {"id": "GO:0005515", "label": "protein binding"},
            "evidence_type": "IPI",
            "original_reference_id": "PMID:12345",
        },
        {
            "term": {"id": "GO:0005515", "label": "protein binding"},
            "evidence_type": "IEA",
            "original_reference_id": "GO_REF:0000033",
        },
    ]
    yaml_file = create_test_yaml_file("TEST5", yaml_anns)

    try:
        validator = GOAValidator()
        result = validator.validate_against_goa(yaml_file, goa_file)

        # Should be valid - both YAML annotations have exact matches in GOA
        assert result.is_valid, (
            f"Should be valid but got: {validator.get_summary(result)}"
        )
        assert len(result.missing_in_goa) == 0
        # Note: missing_in_yaml will have one entry (the IDA annotation not in YAML)
        # but that's OK - we only require YAML annotations to be in GOA, not vice versa
    finally:
        goa_file.unlink()
        yaml_file.unlink()


def test_strict_mode_fails_on_missing_from_yaml():
    """Test that strict mode fails when GOA has annotations not in YAML."""
    # Create GOA annotations
    goa_anns = [
        {
            "id": "P12345",
            "symbol": "TEST6",
            "go_id": "GO:0005515",
            "go_term": "protein binding",
            "evidence": "IPI",
            "reference": "PMID:12345",
        },
        {
            "id": "P12345",
            "symbol": "TEST6",
            "go_id": "GO:0005634",
            "go_term": "nucleus",
            "evidence": "IDA",
            "reference": "PMID:67890",
        },
    ]
    goa_file = create_test_goa_file(goa_anns)

    # Create YAML with only one annotation
    yaml_anns = [
        {
            "term": {"id": "GO:0005515", "label": "protein binding"},
            "evidence_type": "IPI",
            "original_reference_id": "PMID:12345",
        }
    ]
    yaml_file = create_test_yaml_file("TEST6", yaml_anns)

    try:
        # Test with strict mode
        validator = GOAValidator(strict_mode=True)
        result = validator.validate_against_goa(yaml_file, goa_file)

        assert not result.is_valid
        assert len(result.missing_in_yaml) == 1
    finally:
        goa_file.unlink()
        yaml_file.unlink()
