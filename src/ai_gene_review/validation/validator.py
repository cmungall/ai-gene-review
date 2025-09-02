"""Validator for gene review YAML files using LinkML schema.

This module provides functionality to validate gene review YAML files
against the LinkML schema defined in schema/gene_review.yaml.

Example:
    >>> from ai_gene_review.validation import validate_gene_review
    >>> from pathlib import Path
    >>> # Example usage with a gene review file
    >>> # report = validate_gene_review(Path("genes/human/GENE/GENE-ai-review.yaml"))
    >>> # if report.is_valid:
    >>> #     print("Valid!")
    >>> # else:
    >>> #     print(f"Validation errors: {report.error_count}")
"""

from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import yaml
import re
import traceback
from linkml_runtime.utils.schemaview import SchemaView  # type: ignore[import-untyped]
from linkml.validators.jsonschemavalidator import JsonSchemaDataValidator  # type: ignore[import-untyped]

from ai_gene_review.validation.validation_report import (
    ValidationReport,
    ValidationSeverity,
    BatchValidationReport,
)
from ai_gene_review.validation.term_validator import TermValidator
from ai_gene_review.validation.publication_validator import PublicationValidator
from ai_gene_review.validation.goa_validator import GOAValidator
from ai_gene_review.validation.supporting_text_validator import SupportingTextValidator


def get_schema_path() -> Path:
    """Get the path to the LinkML schema file.

    Returns:
        Path to the gene_review.yaml schema file
    """
    return Path(__file__).parent.parent / "schema" / "gene_review.yaml"


def load_schema() -> SchemaView:
    """Load the LinkML schema.

    Returns:
        SchemaView object for the gene review schema

    Raises:
        FileNotFoundError: If schema file is not found
    """
    schema_path = get_schema_path()
    if not schema_path.exists():
        raise FileNotFoundError(f"Schema file not found: {schema_path}")

    return SchemaView(str(schema_path))


def validate_gene_review(
    yaml_file: Path | str,
    schema_path: Optional[Path | str] = None,
    check_best_practices: bool = True,
    check_goa: bool = True,
    check_supporting_text: bool = True,
) -> ValidationReport:
    """Validate a gene review YAML file against the LinkML schema.

    Args:
        yaml_file: Path to the YAML file to validate
        schema_path: Optional path to schema file (uses default if not provided)
        check_best_practices: Whether to check for best practices (soft failures)
        check_goa: Whether to validate against GOA file (enabled by default)
        check_supporting_text: Whether to validate supporting_text against cached publications

    Returns:
        ValidationReport with detailed validation results

    Example:
        >>> # Create a test file
        >>> import tempfile, yaml
        >>> from pathlib import Path
        >>> data = {
        ...     "id": "Q12345",
        ...     "gene_symbol": "TEST",
        ...     "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        ...     "description": "A test gene for demonstration purposes"
        ... }
        >>> with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        ...     yaml.dump(data, f)
        ...     test_file = Path(f.name)
        >>> report = validate_gene_review(test_file)
        >>> print("Valid!" if report.is_valid else f"Errors: {report.error_count}")
        Valid!
        >>> test_file.unlink()  # Clean up
    """
    yaml_file = Path(yaml_file)
    report = ValidationReport(file_path=yaml_file, is_valid=True)

    # Check if file exists
    if not yaml_file.exists():
        report.add_issue(
            ValidationSeverity.ERROR, f"File not found: {yaml_file}", path=None
        )
        return report

    try:
        # Load schema
        if schema_path:
            schema_path_obj = Path(schema_path)
            if not schema_path_obj.exists():
                report.add_issue(
                    ValidationSeverity.ERROR,
                    f"Schema file not found: {schema_path_obj}",
                    path=None,
                )
                return report
            SchemaView(str(schema_path_obj))  # Validate schema can be loaded
            schema_path = schema_path_obj  # Update schema_path to Path object
        else:
            load_schema()  # Validate schema can be loaded

        # Load YAML data
        with open(yaml_file, "r") as f:
            data = yaml.safe_load(f)

        # Create validator
        validator = JsonSchemaDataValidator(str(schema_path or get_schema_path()))

        # Validate against schema
        linkml_report = validator.validate_dict(data, target_class="GeneReview")

        # Process LinkML validation results
        if linkml_report and linkml_report.results:
            for result in linkml_report.results:
                # Extract error message and path
                message = result.message if hasattr(result, "message") else str(result)

                # Try to extract JSON path from message
                path_match = re.search(r"in \$\.(.+)$", message)
                path = path_match.group(1) if path_match else None

                # Schema validation errors are hard failures
                report.add_issue(ValidationSeverity.ERROR, message, path=path)

        # Check best practices if enabled and no hard errors
        if check_best_practices and not report.has_errors:
            check_best_practices_rules(
                data, report, yaml_file if check_goa else None, check_supporting_text
            )

        return report

    except yaml.YAMLError as e:
        report.add_issue(
            ValidationSeverity.ERROR, f"YAML parsing error: {str(e)}", path=None
        )
        return report
    except FileNotFoundError as e:
        report.add_issue(ValidationSeverity.ERROR, f"File error: {str(e)}", path=None)
        return report
    except ImportError as e:
        report.add_issue(
            ValidationSeverity.ERROR,
            f"Import error (check dependencies): {str(e)}",
            path=None,
        )
        return report
    except Exception as e:
        # For debugging: include the exception type
        error_type = type(e).__name__
        tb = traceback.format_exc()
        # Only show the last line of traceback to avoid clutter
        last_tb_line = tb.strip().split("\n")[-1] if tb else ""
        report.add_issue(
            ValidationSeverity.ERROR,
            f"Validation error ({error_type}): {str(e)} - {last_tb_line}",
            path=None,
        )
        return report


def check_best_practices_rules(
    data: Dict[str, Any],
    report: ValidationReport,
    yaml_file: Optional[Path | str] = None,
    check_supporting_text: bool = True,
) -> None:
    """Check for best practices and add soft failures (warnings).

    Args:
        data: The parsed YAML data
        report: ValidationReport to add warnings to
        yaml_file: Path to YAML file for GOA validation (if enabled)
        check_supporting_text: Whether to validate supporting_text against publications
    """
    # Validate ontology terms
    term_validator = TermValidator()
    term_results = term_validator.validate_terms_in_data(data)

    for result in term_results:
        if not result.is_valid:
            report.add_issue(
                ValidationSeverity.ERROR,
                result.error_message or f"Invalid term: {result.term_id}",
                path=result.path,
                suggestion=f"Use correct label: '{result.correct_label}'"
                if result.correct_label
                else None,
            )
        elif result.is_obsolete:
            report.add_issue(
                ValidationSeverity.WARNING,
                f"Term {result.term_id} is obsolete",
                path=result.path,
                suggestion="Consider using a non-obsolete term",
            )

    # Validate GO branch constraints in core_functions
    core_func_results = term_validator.validate_terms_in_core_functions(data)
    for result in core_func_results:
        if not result.is_valid:
            report.add_issue(
                ValidationSeverity.ERROR,
                result.error_message or f"Invalid term: {result.term_id}",
                path=result.path,
                suggestion=f"Use correct label: '{result.correct_label}'"
                if result.correct_label
                else None,
            )

    # Validate publication references
    pub_validator = PublicationValidator()
    pub_results = pub_validator.validate_publications_in_data(data)

    for pub_result in pub_results:
        if not pub_result.is_valid:
            # Determine severity based on the error
            if pub_result.correct_title is None:
                # Could not find publication - this is a warning since it might be a network issue
                severity = ValidationSeverity.WARNING
                suggestion = "Check if PMID is correct or fetch the publication"
            else:
                # Title mismatch - this is an error
                severity = ValidationSeverity.ERROR
                suggestion = f"Use correct title: '{pub_result.correct_title}'"

            report.add_issue(
                severity,
                pub_result.error_message
                or f"Invalid publication: PMID:{pub_result.pmid}",
                path=pub_result.path,
                suggestion=suggestion,
            )
    # Check for TODO in description
    if "description" in data and "TODO" in str(data["description"]):
        report.add_issue(
            ValidationSeverity.WARNING,
            "Description contains TODO placeholder",
            path="description",
            suggestion="Complete the gene description",
        )

    # Check for missing optional but recommended fields
    if "aliases" not in data:
        report.add_issue(
            ValidationSeverity.INFO,
            "No aliases provided for the gene",
            path="aliases",
            suggestion="Consider adding known gene aliases if available",
        )

    # Check for missing references
    if "references" not in data or not data["references"]:
        report.add_issue(
            ValidationSeverity.WARNING,
            "No references provided",
            path="references",
            suggestion="Add relevant literature references",
        )

    # Check that all original_reference_ids in annotations point to valid references
    if "existing_annotations" in data and data["existing_annotations"]:
        # Build set of valid reference IDs
        valid_ref_ids = set()
        if "references" in data and data["references"]:
            for ref in data["references"]:
                if isinstance(ref, dict) and "id" in ref:
                    valid_ref_ids.add(ref["id"])

        # Check each annotation's original_reference_id
        for i, annotation in enumerate(data["existing_annotations"]):
            if isinstance(annotation, dict) and "original_reference_id" in annotation:
                ref_id = annotation["original_reference_id"]
                if ref_id and ref_id not in valid_ref_ids:
                    report.add_issue(
                        ValidationSeverity.ERROR,
                        f"Annotation references non-existent reference ID: {ref_id}",
                        path=f"existing_annotations[{i}].original_reference_id",
                        suggestion=f"Add a reference with ID '{ref_id}' to the references section or correct the reference ID",
                    )

    # Check file: references
    if "references" in data and data["references"]:
        for i, ref in enumerate(data["references"]):
            if isinstance(ref, dict) and "id" in ref:
                ref_id = ref["id"]
                if ref_id.startswith("file:"):
                    # Extract the file path after 'file:'
                    file_path_str = ref_id[5:]  # Remove 'file:' prefix

                    # Check if the file exists relative to the genes directory
                    # Assuming validation is run from project root
                    base_path = Path("genes")
                    full_path = base_path / file_path_str

                    if not full_path.exists():
                        report.add_issue(
                            ValidationSeverity.ERROR,
                            f"File reference points to non-existent file: {file_path_str}",
                            path=f"references[{i}].id",
                            suggestion=f"Ensure the file exists at: genes/{file_path_str}",
                        )
                    elif not full_path.is_file():
                        report.add_issue(
                            ValidationSeverity.ERROR,
                            f"File reference points to a directory, not a file: {file_path_str}",
                            path=f"references[{i}].id",
                            suggestion="File references must point to actual files, not directories",
                        )

    # Check for PENDING annotations (placeholder state that should be resolved)
    if "existing_annotations" in data and data["existing_annotations"]:
        pending_count = 0
        pending_examples: List[str] = []
        for i, annotation in enumerate(data["existing_annotations"]):
            if isinstance(annotation, dict) and "review" in annotation:
                review = annotation.get("review", {})
                if isinstance(review, dict) and review.get("action") == "PENDING":
                    pending_count += 1
                    if len(pending_examples) < 3:  # Collect first 3 examples
                        term = annotation.get("term", {})
                        term_label = term.get("label", "unknown") if isinstance(term, dict) else "unknown"
                        pending_examples.append(f"{term_label}")
        
        if pending_count > 0:
            examples_text = ", ".join(pending_examples)
            if pending_count > 3:
                examples_text += f" ... and {pending_count - 3} more"
            report.add_issue(
                ValidationSeverity.WARNING,
                f"Review incomplete: {pending_count} annotations marked as PENDING ({examples_text})",
                path="existing_annotations",
                suggestion="Complete the review by resolving PENDING annotations",
            )
    
    # Check for missing core functions
    if "core_functions" not in data or not data["core_functions"]:
        report.add_issue(
            ValidationSeverity.WARNING,
            "No core functions defined",
            path="core_functions",
            suggestion="Define the core molecular functions of the gene",
        )
    else:
        # Validate that core function terms are properly supported
        # Build set of accepted GO terms from existing_annotations
        accepted_terms = set()
        proposed_replacement_terms = set()
        
        # First pass: collect all ACCEPTED terms and proposed replacements
        if "existing_annotations" in data and data["existing_annotations"]:
            for annotation in data["existing_annotations"]:
                if isinstance(annotation, dict):
                    review = annotation.get("review", {})
                    if isinstance(review, dict):
                        action = review.get("action")
                        
                        # Collect ACCEPTED terms
                        if action == "ACCEPT":
                            term = annotation.get("term", {})
                            if isinstance(term, dict) and "id" in term:
                                accepted_terms.add(term["id"])
                        
                        # Collect proposed replacement terms from MODIFY actions
                        if action == "MODIFY":
                            replacements = review.get("proposed_replacement_terms", [])
                            if replacements:
                                for replacement in replacements:
                                    if isinstance(replacement, dict) and "id" in replacement:
                                        proposed_replacement_terms.add(replacement["id"])
        
        # Second pass: identify terms that should NOT be used
        # A term should not be used ONLY if:
        # 1. It has at least one MODIFY action AND
        # 2. It has NO ACCEPT actions AND  
        # 3. It's not a proposed replacement for something else
        modified_terms = set()
        if "existing_annotations" in data and data["existing_annotations"]:
            for annotation in data["existing_annotations"]:
                if isinstance(annotation, dict):
                    review = annotation.get("review", {})
                    if isinstance(review, dict):
                        action = review.get("action")
                        
                        if action == "MODIFY":
                            term = annotation.get("term", {})
                            if isinstance(term, dict) and "id" in term:
                                term_id = term["id"]
                                # Only mark as "should not use" if it's not accepted elsewhere
                                # and not a proposed replacement
                                if (term_id not in accepted_terms and 
                                    term_id not in proposed_replacement_terms):
                                    modified_terms.add(term_id)
        
        # Now check each core function
        for i, core_func in enumerate(data["core_functions"]):
            if isinstance(core_func, dict):
                # Get the molecular_function term
                mol_func = core_func.get("molecular_function", {})
                if isinstance(mol_func, dict) and "id" in mol_func:
                    term_id = mol_func["id"]
                    term_label = mol_func.get("label", term_id)
                    
                    # Check if this term is from accepted annotations or proposed replacements
                    is_from_accepted = term_id in accepted_terms or term_id in proposed_replacement_terms
                    
                    # Check if it has supported_by references
                    supported_by = core_func.get("supported_by", [])
                    has_support = bool(supported_by)
                    
                    # If neither condition is met, it's an error
                    if not is_from_accepted and not has_support:
                        report.add_issue(
                            ValidationSeverity.ERROR,
                            f"Core function term {term_id} ({term_label}) is not from an ACCEPTED annotation and lacks supported_by references",
                            path=f"core_functions[{i}].molecular_function",
                            suggestion="Either use a term from ACCEPTED existing annotations/proposed replacements, or add supported_by references",
                        )
                
                # Check locations field (should be CC terms)
                locations = core_func.get("locations", [])
                for j, location in enumerate(locations):
                    if isinstance(location, dict) and "id" in location:
                        loc_id = location["id"]
                        loc_label = location.get("label", loc_id)
                        
                        # First check if this is a term marked for modification
                        if loc_id in modified_terms:
                            report.add_issue(
                                ValidationSeverity.ERROR,
                                f"Location term {loc_id} ({loc_label}) was marked for MODIFY and should not be used directly",
                                path=f"core_functions[{i}].locations[{j}]",
                                suggestion="Use the proposed replacement term instead of the original term marked for modification",
                            )
                            continue
                        
                        # Check if this term is from accepted annotations or proposed replacements
                        is_from_accepted = loc_id in accepted_terms or loc_id in proposed_replacement_terms
                        
                        # Locations don't have their own supported_by, they rely on the core function's supported_by
                        # So check if the core function has supported_by
                        supported_by = core_func.get("supported_by", [])
                        has_support = bool(supported_by)
                        
                        # If neither condition is met, it's an error
                        if not is_from_accepted and not has_support:
                            report.add_issue(
                                ValidationSeverity.ERROR,
                                f"Location term {loc_id} ({loc_label}) is not from an ACCEPTED annotation and the core function lacks supported_by references",
                                path=f"core_functions[{i}].locations[{j}]",
                                suggestion="Either use a location from ACCEPTED existing annotations/proposed replacements, or add supported_by references to the core function",
                            )
                
                # Check directly_involved_in field (should be BP terms)
                directly_involved = core_func.get("directly_involved_in", [])
                for j, process in enumerate(directly_involved):
                    if isinstance(process, dict) and "id" in process:
                        proc_id = process["id"]
                        proc_label = process.get("label", proc_id)
                        
                        # Check if this term is from accepted annotations or proposed replacements
                        is_from_accepted = proc_id in accepted_terms or proc_id in proposed_replacement_terms
                        
                        # Check if the core function has supported_by
                        supported_by = core_func.get("supported_by", [])
                        has_support = bool(supported_by)
                        
                        # If neither condition is met, it's an error
                        if not is_from_accepted and not has_support:
                            report.add_issue(
                                ValidationSeverity.ERROR,
                                f"Process term {proc_id} ({proc_label}) is not from an ACCEPTED annotation and the core function lacks supported_by references",
                                path=f"core_functions[{i}].directly_involved_in[{j}]",
                                suggestion="Either use a process from ACCEPTED existing annotations/proposed replacements, or add supported_by references to the core function",
                            )
                
                # Check in_complex field
                in_complex = core_func.get("in_complex")
                if isinstance(in_complex, dict) and "id" in in_complex:
                    complex_id = in_complex["id"]
                    complex_label = in_complex.get("label", complex_id)
                    
                    # Check if this term is from accepted annotations or proposed replacements
                    is_from_accepted = complex_id in accepted_terms or complex_id in proposed_replacement_terms
                    
                    # Check if the core function has supported_by
                    supported_by = core_func.get("supported_by", [])
                    has_support = bool(supported_by)
                    
                    # If neither condition is met, it's an error
                    if not is_from_accepted and not has_support:
                        report.add_issue(
                            ValidationSeverity.ERROR,
                            f"Complex term {complex_id} ({complex_label}) is not from an ACCEPTED annotation and the core function lacks supported_by references",
                            path=f"core_functions[{i}].in_complex",
                            suggestion="Either use a complex from ACCEPTED existing annotations/proposed replacements, or add supported_by references to the core function",
                        )

    # Check for ACCEPT annotations with PMIDs lacking supported_by
    # Only warn if the publication file exists (full text is available)
    if "existing_annotations" in data and data["existing_annotations"]:
        for i, annotation in enumerate(data["existing_annotations"]):
            if isinstance(annotation, dict):
                # Check if this annotation has a review with ACCEPT action
                review = annotation.get("review", {})
                if isinstance(review, dict) and review.get("action") == "ACCEPT":
                    # Check if it has a PMID reference
                    ref_id = annotation.get("original_reference_id", "")
                    if ref_id and ref_id.startswith("PMID:"):
                        # Check if supported_by is missing or empty
                        supported_by = review.get("supported_by", [])
                        if not supported_by:
                            # Check if publication file exists and has full text available
                            pmid_number = ref_id.replace("PMID:", "")
                            # Look for publications directory relative to the YAML file or in project root
                            if yaml_file is not None:
                                # Look relative to the YAML file's project root
                                project_root = yaml_file.parent
                                while project_root.parent != project_root and not (project_root / "publications").exists():
                                    project_root = project_root.parent
                                pub_file = project_root / "publications" / f"PMID_{pmid_number}.md"
                            else:
                                pub_file = Path("publications") / f"PMID_{pmid_number}.md"
                            
                            # Check if full text is available in the publication file
                            full_text_available = False
                            if pub_file.exists():
                                try:
                                    import yaml as yaml_lib
                                    with open(pub_file, 'r') as f:
                                        content = f.read()
                                        # Extract frontmatter between --- markers
                                        if content.startswith('---'):
                                            end_marker = content.find('---', 3)
                                            if end_marker != -1:
                                                frontmatter = content[3:end_marker]
                                                pub_data = yaml_lib.safe_load(frontmatter)
                                                full_text_available = pub_data.get('full_text_available', False)
                                except Exception:
                                    # If we can't parse the file, assume no full text
                                    pass
                            
                            # Only warn if full text is available
                            if full_text_available:
                                report.add_issue(
                                    ValidationSeverity.WARNING,
                                    f"ACCEPT annotation with {ref_id} lacks supported_by references",
                                    path=f"existing_annotations[{i}].review.supported_by",
                                    suggestion=f"Add supported_by with reference_id and supporting_text from {ref_id}",
                                )

    # Check for incomplete taxon information
    if "taxon" in data:
        taxon = data["taxon"]
        if isinstance(taxon, dict):
            if "id" in taxon and not taxon["id"].startswith("NCBITaxon:"):
                report.add_issue(
                    ValidationSeverity.INFO,
                    "Taxon ID should use NCBITaxon prefix",
                    path="taxon.id",
                    suggestion=f"Use 'NCBITaxon:{taxon['id']}' format",
                )

    # Check description length
    if "description" in data and len(str(data["description"])) < 20:
        report.add_issue(
            ValidationSeverity.WARNING,
            "Description is very short",
            path="description",
            suggestion="Provide a more detailed description of the gene",
        )

    # Validate against GOA file if enabled
    if yaml_file is not None and "existing_annotations" in data:
        goa_validator = GOAValidator()
        # yaml_file is guaranteed to be Path here (converted at line 92)
        assert isinstance(yaml_file, Path)
        goa_result = goa_validator.validate_against_goa(yaml_file)

        if not goa_result.is_valid:
            # Report GOA validation issues
            if goa_result.error_message:
                report.add_issue(
                    ValidationSeverity.WARNING,
                    f"GOA validation: {goa_result.error_message}",
                    path="existing_annotations",
                )

            if goa_result.missing_in_yaml:
                # Show first few missing annotations for debugging
                missing_count = len(goa_result.missing_in_yaml)
                examples = []
                for ann in goa_result.missing_in_yaml[:5]:  # Show first 5
                    examples.append(
                        f"{ann.go_id} ({ann.go_term}) - {ann.evidence_code} - {ann.reference}"
                    )

                if missing_count <= 5:
                    detail = ": " + "; ".join(examples)
                else:
                    detail = (
                        " (showing first 5): "
                        + "; ".join(examples)
                        + f" ... and {missing_count - 5} more"
                    )

                report.add_issue(
                    ValidationSeverity.ERROR,
                    f"Missing {missing_count} annotations from GOA{detail}",
                    path="existing_annotations",
                    suggestion="Review GOA file and add missing annotations or document why they were excluded",
                )

            if goa_result.missing_in_goa:
                # Show first few annotations not in GOA for debugging
                missing_count = len(goa_result.missing_in_goa)
                examples = []
                for ann_dict in goa_result.missing_in_goa[:5]:  # Show first 5
                    # ann_dict is a Dict from the YAML, not a GOAAnnotation
                    go_id = ann_dict.get("term", {}).get("id", "unknown")
                    go_label = ann_dict.get("term", {}).get("label", "unknown")
                    evidence = ann_dict.get("evidence_type", "unknown")
                    ref = ann_dict.get("original_reference_id", "unknown")
                    examples.append(f"{go_id} ({go_label}) - {evidence} - {ref}")

                if missing_count <= 5:
                    detail = ": " + "; ".join(examples)
                else:
                    detail = (
                        " (showing first 5): "
                        + "; ".join(examples)
                        + f" ... and {missing_count - 5} more"
                    )

                report.add_issue(
                    ValidationSeverity.ERROR,
                    f"Found {missing_count} annotations not in GOA{detail}",
                    path="existing_annotations",
                    suggestion="Verify these annotations are correct or remove if not supported by GOA",
                )

            if goa_result.mismatched_labels:
                report.add_issue(
                    ValidationSeverity.ERROR,
                    f"Found {len(goa_result.mismatched_labels)} GO term label mismatches with GOA file",
                    path="existing_annotations",
                    suggestion="Update GO term labels to match GOA file",
                )

            if goa_result.mismatched_evidence:
                # Evidence mismatches are warnings, not errors
                report.add_issue(
                    ValidationSeverity.WARNING,
                    f"Found {len(goa_result.mismatched_evidence)} evidence type mismatches with GOA file",
                    path="existing_annotations",
                    suggestion="Consider updating evidence types to match GOA file",
                )

    # Validate supporting_text against cached publications if enabled
    if check_supporting_text and "existing_annotations" in data:
        st_validator = SupportingTextValidator()
        st_report = st_validator.validate_data(data)

        if not st_report.is_valid:
            # Report invalid supporting texts
            for st_result in st_report.results:
                if not st_result.is_valid:
                    severity = ValidationSeverity.WARNING

                    # Check if this is a PMID finding without supporting_text
                    if "PMID reference" in (
                        st_result.error_message or ""
                    ) and "without supporting_text" in (st_result.error_message or ""):
                        # Always WARNING for missing supporting_text in PMID findings
                        severity = ValidationSeverity.WARNING
                    elif st_result.similarity_score < 0.5:
                        # Very low similarity suggests incorrect reference
                        severity = ValidationSeverity.ERROR

                    report.add_issue(
                        severity,
                        st_result.error_message
                        or "Supporting text not found in referenced publication",
                        path=st_result.annotation_path,
                        suggestion=st_result.suggested_fix,
                    )

        # Report coverage statistics as info
        if st_report.total_annotations > 0:
            coverage = st_report.validation_rate
            if coverage < 50:
                report.add_issue(
                    ValidationSeverity.INFO,
                    f"Only {coverage:.1f}% of annotations have supporting_text",
                    path="existing_annotations",
                    suggestion="Consider adding supporting_text to more annotations for better documentation",
                )


def validate_multiple_files(
    yaml_files: List[Path | str],
    schema_path: Optional[Path | str] = None,
    check_best_practices: bool = True,
    check_goa: bool = True,
    check_supporting_text: bool = True,
) -> BatchValidationReport:
    """Validate multiple gene review YAML files.

    Args:
        yaml_files: List of paths to YAML files
        schema_path: Optional path to schema file
        check_best_practices: Whether to check for best practices
        check_goa: Whether to validate against GOA files
        check_supporting_text: Whether to validate supporting_text against cached publications

    Returns:
        BatchValidationReport with results for all files

    Example:
        >>> # Create test files
        >>> import tempfile, yaml
        >>> from pathlib import Path
        >>> files = []
        >>> for i in range(2):
        ...     data = {
        ...         "id": f"Q{i}",
        ...         "gene_symbol": f"TEST{i}",
        ...         "taxon": {"id": "NCBITaxon:9606", "label": "Homo sapiens"},
        ...         "description": f"Test gene {i}"
        ...     }
        ...     with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        ...         yaml.dump(data, f)
        ...         files.append(Path(f.name))
        >>> batch_report = validate_multiple_files(files)
        >>> print(f"Valid: {batch_report.valid_files}, Invalid: {batch_report.invalid_files}")
        Valid: 2, Invalid: 0
        >>> for f in files:
        ...     f.unlink()  # Clean up
    """
    batch_report = BatchValidationReport()

    for yaml_file in yaml_files:
        yaml_file = Path(yaml_file)
        report = validate_gene_review(
            yaml_file,
            schema_path,
            check_best_practices,
            check_goa,
            check_supporting_text,
        )
        batch_report.reports.append(report)

    return batch_report


def get_validation_summary(results: Dict[str, Tuple[bool, List[str]]]) -> str:
    """Legacy function for backward compatibility.

    Args:
        results: Dictionary from old validate_multiple_files

    Returns:
        Human-readable summary string
    """
    # Convert old format to new format
    batch_report = BatchValidationReport()

    for file_path, (is_valid, errors) in results.items():
        report = ValidationReport(file_path=Path(file_path), is_valid=is_valid)
        for error in errors:
            report.add_issue(ValidationSeverity.ERROR, error)
        batch_report.reports.append(report)

    return batch_report.summary()
