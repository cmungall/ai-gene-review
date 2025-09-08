"""Validator for checking existing_annotations against GOA source files."""

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests
import yaml


@dataclass
class GOAAnnotation:
    """Represents a single GOA annotation from the TSV file."""

    database: str
    db_object_id: str
    db_object_symbol: str
    qualifier: str
    go_id: str
    go_term: str
    go_aspect: str
    evidence_code: str
    evidence_type: str
    reference: str
    with_from: str
    taxon_id: str
    taxon_label: str
    assigned_by: str
    annotation_extension: str
    date: str

    @classmethod
    def from_tsv_row(cls, row: list) -> "GOAAnnotation":
        """Parse a GOA annotation from a TSV row."""
        # GOA format columns:
        # 0: GENE PRODUCT DB, 1: GENE PRODUCT ID, 2: SYMBOL, 3: QUALIFIER
        # 4: GO TERM, 5: GO NAME, 6: GO ASPECT, 7: ECO ID, 8: GO EVIDENCE CODE
        # 9: REFERENCE, 10: WITH/FROM, 11: TAXON ID, 12: TAXON NAME
        # 13: ASSIGNED BY, 14: GENE NAME, 15: DATE
        return cls(
            database=row[0] if len(row) > 0 else "",
            db_object_id=row[1] if len(row) > 1 else "",
            db_object_symbol=row[2] if len(row) > 2 else "",
            qualifier=row[3] if len(row) > 3 else "",
            go_id=row[4] if len(row) > 4 else "",
            go_term=row[5] if len(row) > 5 else "",
            go_aspect=row[6] if len(row) > 6 else "",
            evidence_code=row[8] if len(row) > 8 else "",  # GO EVIDENCE CODE at index 8
            evidence_type=row[7] if len(row) > 7 else "",  # ECO ID at index 7
            reference=row[9] if len(row) > 9 else "",
            with_from=row[10] if len(row) > 10 else "",
            taxon_id=row[11] if len(row) > 11 else "",
            taxon_label=row[12] if len(row) > 12 else "",
            assigned_by=row[13] if len(row) > 13 else "",
            annotation_extension=row[14] if len(row) > 14 else "",  # Actually GENE NAME
            date=row[15] if len(row) > 15 else "",
        )


@dataclass
class GOAValidationResult:
    """Result of validating existing_annotations against GOA file."""

    is_valid: bool
    missing_in_yaml: List[GOAAnnotation] = field(default_factory=list)
    missing_in_goa: List[Dict] = field(default_factory=list)
    mismatched_labels: List[Tuple[str, str, str]] = field(
        default_factory=list
    )  # (go_id, yaml_label, goa_label)
    mismatched_evidence: List[Tuple[str, str, str]] = field(
        default_factory=list
    )  # (go_id, yaml_evidence, goa_evidence)
    error_message: Optional[str] = None


@dataclass
class GOAValidator:
    """Validates that existing_annotations match the GOA source file.

    This validator ensures consistency between the curated existing_annotations
    in the YAML file and the actual annotations in the GOA TSV file.

    Example:
        >>> validator = GOAValidator()
        >>> result = validator.validate_against_goa(
        ...     yaml_file=Path("genes/human/CFAP300/CFAP300-ai-review.yaml"),
        ...     goa_file=Path("genes/human/CFAP300/CFAP300-goa.tsv")
        ... )
        >>> if result.is_valid:
        ...     print("Annotations match!")
        ... else:
        ...     print(f"Found {len(result.missing_in_yaml)} annotations missing from YAML")
        Annotations match!
    """

    # Whether to check for exact matches (strict mode)
    strict_mode: bool = False
    # Whether to check complete tuples (term, evidence, reference) vs just GO terms
    strict_tuple_matching: bool = True

    def parse_goa_file(self, goa_file: Path) -> List[GOAAnnotation]:
        """Parse a GOA TSV file into annotation objects.

        Args:
            goa_file: Path to the GOA TSV file

        Returns:
            List of GOAAnnotation objects
        """
        annotations: List[GOAAnnotation] = []

        if not goa_file.exists():
            return annotations

        try:
            with open(goa_file, "r", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                header_skipped = False
                for row in reader:
                    # Skip empty rows or comments
                    if not row or (row[0] and row[0].startswith("#")):
                        continue

                    # Skip header row (first non-comment row)
                    if not header_skipped:
                        # Check for common header patterns
                        if (
                            row[0] in ["DB", "GENE PRODUCT DB", "Database"]
                            or "DB" in row[0]
                            or "Database" in row[0]
                            or row[0].startswith("DB")
                        ):
                            header_skipped = True
                            continue
                        header_skipped = True

                    annotation = GOAAnnotation.from_tsv_row(row)
                    annotations.append(annotation)
        except Exception as e:
            print(f"Error parsing GOA file {goa_file}: {e}")

        return annotations

    def parse_yaml_annotations(self, yaml_file: Path) -> List[Dict]:
        """Extract existing_annotations from YAML file.

        Args:
            yaml_file: Path to the gene review YAML file

        Returns:
            List of annotation dictionaries from the YAML
        """
        if not yaml_file.exists():
            return []

        try:
            with open(yaml_file, "r") as f:
                data = yaml.safe_load(f)
                return data.get("existing_annotations", [])
        except Exception as e:
            print(f"Error parsing YAML file {yaml_file}: {e}")
            return []

    def validate_against_goa(
        self, yaml_file: Path, goa_file: Optional[Path] = None
    ) -> GOAValidationResult:
        """Validate existing_annotations against GOA source file.

        Args:
            yaml_file: Path to the gene review YAML file
            goa_file: Path to GOA file (if None, derives from yaml_file path)

        Returns:
            GOAValidationResult with validation details
        """
        # Derive GOA file path if not provided
        if goa_file is None:
            # Assume standard structure: gene-name-ai-review.yaml -> gene-name-goa.tsv
            yaml_stem = yaml_file.stem  # e.g., "CFAP300-ai-review"
            if yaml_stem.endswith("-ai-review"):
                gene_name = yaml_stem[:-10]  # Remove "-ai-review"
                goa_file = yaml_file.parent / f"{gene_name}-goa.tsv"
            else:
                return GOAValidationResult(
                    is_valid=False,
                    error_message=f"Could not derive GOA file path from {yaml_file}",
                )

        # Check if GOA file exists
        if not goa_file.exists():
            return GOAValidationResult(
                is_valid=False, error_message=f"GOA file not found: {goa_file}"
            )

        # Parse both files
        goa_annotations = self.parse_goa_file(goa_file)
        yaml_annotations = self.parse_yaml_annotations(yaml_file)

        # Initialize result
        result = GOAValidationResult(is_valid=True)

        if self.strict_tuple_matching:
            # Strict mode: check complete tuples (go_id, evidence_code, reference)
            return self._validate_strict_tuples(
                goa_annotations, yaml_annotations, result
            )
        else:
            # Legacy mode: check by GO ID only
            return self._validate_by_go_id(goa_annotations, yaml_annotations, result)

    def _validate_strict_tuples(
        self,
        goa_annotations: List[GOAAnnotation],
        yaml_annotations: List[Dict],
        result: GOAValidationResult,
    ) -> GOAValidationResult:
        """Validate using strict tuple matching (go_id, evidence_code, reference)."""
        # Create lookup structure for GOA using complete tuples
        # Key: (go_id, evidence_code, reference)
        goa_tuples = set()
        goa_by_tuple: Dict[Tuple[str, str, str], GOAAnnotation] = {}
        for ann in goa_annotations:
            tuple_key = (ann.go_id, ann.evidence_code, ann.reference)
            goa_tuples.add(tuple_key)
            goa_by_tuple[tuple_key] = ann

        # Check each YAML annotation against GOA
        for yaml_ann in yaml_annotations:
            if not isinstance(yaml_ann, dict):
                continue

            # Check if this annotation has action=NEW
            review = yaml_ann.get("review", {})
            if isinstance(review, dict) and review.get("action") == "NEW":
                # NEW annotations should NOT exist in GOA
                term = yaml_ann.get("term", {})
                if isinstance(term, dict):
                    go_id = term.get("id", "")
                    evidence_type = yaml_ann.get("evidence_type", "")
                    original_ref = yaml_ann.get("original_reference_id", "")
                    
                    # Check if this annotation exists in GOA
                    yaml_tuple = (go_id, evidence_type, original_ref)
                    if yaml_tuple in goa_tuples:
                        # This is an error - NEW annotation should not exist in GOA
                        result.is_valid = False
                        result.error_message = f"Annotation with action=NEW exists in GOA: {go_id} ({evidence_type}, {original_ref})"
                # Skip further validation for NEW annotations
                continue

            # Extract the tuple from YAML annotation
            term = yaml_ann.get("term", {})
            if not isinstance(term, dict):
                continue

            go_id = term.get("id", "")
            evidence_type = yaml_ann.get("evidence_type", "")
            original_ref = yaml_ann.get("original_reference_id", "")

            if not go_id:
                continue

            # Create the tuple to check
            yaml_tuple = (go_id, evidence_type, original_ref)

            # Check if this exact tuple exists in GOA
            if yaml_tuple not in goa_tuples:
                # This exact annotation is not in GOA
                result.missing_in_goa.append(yaml_ann)
                result.is_valid = False

                # Check if it's a partial match (same GO term but different evidence/ref)
                partial_matches = [
                    goa_ann for key, goa_ann in goa_by_tuple.items() if key[0] == go_id
                ]

                if partial_matches:
                    # There are annotations for this GO term, but with different evidence/ref
                    # Check for evidence mismatch
                    goa_evidence_for_term = {
                        ann.evidence_code for ann in partial_matches
                    }
                    if evidence_type and evidence_type not in goa_evidence_for_term:
                        result.mismatched_evidence.append(
                            (
                                go_id,
                                evidence_type,
                                ", ".join(sorted(goa_evidence_for_term)),
                            )
                        )

                    # Check for reference mismatch
                    goa_refs_for_term = {ann.reference for ann in partial_matches}
                    if original_ref and original_ref not in goa_refs_for_term:
                        # This is a reference mismatch - add to result if we track it
                        pass  # We could add a mismatched_references list if needed

            # Skip label consistency check - we validate against ontology, not GOA
            # GOA files may use synonyms instead of primary labels
            # Label validation is handled by term_validator.py against the ontology

        # Check for annotations in GOA but not in YAML
        # Build set of YAML annotations, excluding NEW annotations (which shouldn't be in GOA)
        yaml_tuples = set()
        for yaml_ann in yaml_annotations:
            if isinstance(yaml_ann, dict):
                # Skip NEW annotations - they should NOT be in GOA
                review = yaml_ann.get("review", {})
                if isinstance(review, dict) and review.get("action") == "NEW":
                    continue
                    
                term = yaml_ann.get("term", {})
                if isinstance(term, dict):
                    go_id = term.get("id", "")
                    evidence_type = yaml_ann.get("evidence_type", "")
                    original_ref = yaml_ann.get("original_reference_id", "")
                    if go_id:
                        yaml_tuples.add((go_id, evidence_type, original_ref))

        for goa_tuple in goa_tuples:
            if goa_tuple not in yaml_tuples:
                goa_ann = goa_by_tuple[goa_tuple]
                result.missing_in_yaml.append(goa_ann)
                # Missing GOA annotations in YAML is always a validation failure
                result.is_valid = False

        return result

    def _validate_by_go_id(
        self,
        goa_annotations: List[GOAAnnotation],
        yaml_annotations: List[Dict],
        result: GOAValidationResult,
    ) -> GOAValidationResult:
        """Legacy validation: check by GO ID only (less strict)."""
        # Create lookup structures
        goa_by_go_id: Dict[str, List[GOAAnnotation]] = {}
        for ann in goa_annotations:
            go_id = ann.go_id
            if go_id not in goa_by_go_id:
                goa_by_go_id[go_id] = []
            goa_by_go_id[go_id].append(ann)

        yaml_by_go_id: Dict[str, List[Dict]] = {}
        for yaml_ann in yaml_annotations:
            if isinstance(yaml_ann, dict) and "term" in yaml_ann:
                # Check if this annotation has action=NEW
                review = yaml_ann.get("review", {})
                if isinstance(review, dict) and review.get("action") == "NEW":
                    # NEW annotations should NOT exist in GOA
                    term = yaml_ann["term"]
                    if isinstance(term, dict) and "id" in term:
                        go_id = term["id"]
                        if go_id in goa_by_go_id:
                            # This is an error - NEW annotation should not exist in GOA
                            result.is_valid = False
                            result.error_message = f"Annotation with action=NEW exists in GOA: {go_id}"
                    # Skip adding NEW annotations to yaml_by_go_id
                    continue
                    
                term = yaml_ann["term"]
                if isinstance(term, dict) and "id" in term:
                    go_id = term["id"]
                    if go_id not in yaml_by_go_id:
                        yaml_by_go_id[go_id] = []
                    yaml_by_go_id[go_id].append(yaml_ann)

        # Check for annotations in GOA but not in YAML
        for go_id, goa_anns in goa_by_go_id.items():
            if go_id not in yaml_by_go_id:
                # This GO term is in GOA but not in YAML
                for goa_ann in goa_anns:
                    result.missing_in_yaml.append(goa_ann)
                    result.is_valid = False

        # Check for annotations in YAML but not in GOA
        for go_id, yaml_anns in yaml_by_go_id.items():
            if go_id not in goa_by_go_id:
                # This GO term is in YAML but not in GOA
                for yaml_ann in yaml_anns:
                    result.missing_in_goa.append(yaml_ann)
                    result.is_valid = False
            else:
                # Check for mismatches in annotations that exist in both
                goa_anns = goa_by_go_id[go_id]

                for yaml_ann in yaml_anns:
                    # Skip label mismatch check - we validate against ontology, not GOA
                    # GOA files may use synonyms instead of primary labels
                    # Label validation is handled by term_validator.py against the ontology
                    pass

                    # Check evidence type mismatch
                    yaml_evidence = yaml_ann.get("evidence_type", "")
                    goa_evidence_codes = {ann.evidence_code for ann in goa_anns}

                    if (
                        yaml_evidence
                        and goa_evidence_codes
                        and yaml_evidence not in goa_evidence_codes
                    ):
                        # Evidence mismatch
                        result.mismatched_evidence.append(
                            (
                                go_id,
                                yaml_evidence,
                                ", ".join(sorted(goa_evidence_codes)),
                            )
                        )
                        # In legacy mode, evidence mismatch doesn't fail validation
                        # unless strict_mode is enabled
                        if self.strict_mode:
                            result.is_valid = False

        return result

    def get_summary(self, result: GOAValidationResult) -> str:
        """Get a human-readable summary of validation results.

        Args:
            result: The validation result

        Returns:
            Summary string
        """
        if result.is_valid:
            return "✓ All existing_annotations match GOA file"

        lines = ["✗ Annotations do not match GOA file:"]

        if result.error_message:
            lines.append(f"  Error: {result.error_message}")

        if result.missing_in_yaml:
            lines.append(
                f"  - {len(result.missing_in_yaml)} annotations in GOA but not in YAML:"
            )
            for ann in result.missing_in_yaml[:3]:  # Show first 3
                lines.append(f"    • {ann.go_id} ({ann.go_term}) - {ann.evidence_type}")
            if len(result.missing_in_yaml) > 3:
                lines.append(f"    ... and {len(result.missing_in_yaml) - 3} more")

        if result.missing_in_goa:
            lines.append(
                f"  - {len(result.missing_in_goa)} annotations in YAML but not in GOA (exact tuple mismatch):"
            )
            for yaml_ann in result.missing_in_goa[:3]:  # Show first 3
                go_id = yaml_ann.get("term", {}).get("id", "unknown")
                label = yaml_ann.get("term", {}).get("label", "unknown")
                evidence = yaml_ann.get("evidence_type", "unknown")
                ref = yaml_ann.get("original_reference_id", "unknown")
                lines.append(f"    • {go_id} ({label}) - {evidence} - {ref}")
            if len(result.missing_in_goa) > 3:
                lines.append(f"    ... and {len(result.missing_in_goa) - 3} more")

        if result.mismatched_labels:
            lines.append(f"  - {len(result.mismatched_labels)} label mismatches:")
            for go_id, yaml_label, goa_label in result.mismatched_labels[:3]:
                lines.append(f"    • {go_id}: YAML='{yaml_label}' vs GOA='{goa_label}'")
            if len(result.mismatched_labels) > 3:
                lines.append(f"    ... and {len(result.mismatched_labels) - 3} more")

        if result.mismatched_evidence:
            lines.append(
                f"  - {len(result.mismatched_evidence)} evidence type mismatches:"
            )
            for go_id, yaml_ev, goa_ev in result.mismatched_evidence[:3]:
                lines.append(f"    • {go_id}: YAML='{yaml_ev}' vs GOA='{goa_ev}'")
            if len(result.mismatched_evidence) > 3:
                lines.append(f"    ... and {len(result.mismatched_evidence) - 3} more")

        return "\n".join(lines)

    def seed_missing_annotations(
        self,
        yaml_file: Path,
        goa_file: Optional[Path] = None,
        output_file: Optional[Path] = None,
        fetch_titles: bool = False,
    ) -> Tuple[int, Path, int]:
        """Seed missing annotations from GOA file into YAML.

        This function adds any annotations present in the GOA file but missing
        from the YAML file. It does NOT overwrite existing annotations.
        The review section is left empty for AI to fill in.

        Args:
            yaml_file: Path to the gene review YAML file
            goa_file: Path to GOA file (if None, derives from yaml_file path)
            output_file: Output path (if None, overwrites input file)
            fetch_titles: If True, fetch actual titles from PubMed (may be slow)

        Returns:
            Tuple of (number of annotations added, output file path, number of references added)
        """
        # Derive GOA file path if not provided
        if goa_file is None:
            yaml_stem = yaml_file.stem
            if yaml_stem.endswith("-ai-review"):
                gene_name = yaml_stem[:-10]
                goa_file = yaml_file.parent / f"{gene_name}-goa.tsv"
            else:
                raise ValueError(f"Could not derive GOA file path from {yaml_file}")

        # Check if GOA file exists
        if not goa_file.exists():
            raise FileNotFoundError(f"GOA file not found: {goa_file}")

        # Parse GOA annotations
        goa_annotations = self.parse_goa_file(goa_file)

        # Load the full YAML data (not just annotations)
        if yaml_file.exists():
            with open(yaml_file, "r") as f:
                yaml_data = yaml.safe_load(f) or {}
        else:
            # Create minimal structure if file doesn't exist
            yaml_data = {
                "id": "",
                "gene_symbol": "",
                "taxon": {"id": "", "label": ""},
                "description": "TODO: Add gene description",
            }

        # Get existing annotations
        existing_annotations = yaml_data.get("existing_annotations", [])

        # Build set of existing tuples (GO ID, evidence_type, reference)
        existing_tuples = set()
        for ann in existing_annotations:
            if isinstance(ann, dict) and "term" in ann:
                term = ann["term"]
                if isinstance(term, dict) and "id" in term:
                    go_id = term["id"]
                    evidence = ann.get("evidence_type", "")
                    ref = ann.get("original_reference_id", "")
                    existing_tuples.add((go_id, evidence, ref))

        # Add missing annotations from GOA
        added_count = 0
        seen_tuples = set()  # Track which tuples we've already added
        pmids_to_add = set()  # Collect PMIDs and GO_REFs to add to references

        # First pass - collect ALL references from GOA (not just missing annotations)
        for goa_ann in goa_annotations:
            if goa_ann.reference and (
                goa_ann.reference.startswith("PMID:")
                or goa_ann.reference.startswith("GO_REF:")
                or goa_ann.reference.startswith("Reactome:")
            ):
                pmids_to_add.add(goa_ann.reference)

        # Second pass - add missing annotations based on complete tuple
        for goa_ann in goa_annotations:
            go_id = goa_ann.go_id
            evidence = goa_ann.evidence_code
            reference = goa_ann.reference
            tuple_key = (go_id, evidence, reference)

            # Skip if this exact tuple already exists in YAML or was already added
            if tuple_key in existing_tuples or tuple_key in seen_tuples:
                continue

            # Create new annotation entry (stub for review)
            new_annotation = {
                "term": {"id": go_id, "label": goa_ann.go_term},
                "evidence_type": evidence,
                "original_reference_id": reference,
            }

            # Add optional review section placeholder
            # This signals that the annotation needs review
            new_annotation["review"] = {
                "summary": "TODO: Review this GOA annotation",
                "action": "PENDING",
            }

            existing_annotations.append(new_annotation)
            seen_tuples.add(tuple_key)
            added_count += 1

        # Update the YAML data
        yaml_data["existing_annotations"] = existing_annotations

        # Seed references section if we have PMIDs and it doesn't exist
        if pmids_to_add and "references" not in yaml_data:
            yaml_data["references"] = []

        # Get existing reference IDs
        existing_refs = yaml_data.get("references", [])
        existing_ref_ids = set()
        for ref in existing_refs:
            if isinstance(ref, dict) and "id" in ref:
                existing_ref_ids.add(ref["id"])

        # Add missing PMIDs/GO_REFs to references
        refs_added = 0
        max_fetches = 500  # High limit since cache makes it fast

        # Import publication fetcher (only if we need it)
        if pmids_to_add and fetch_titles:
            from ai_gene_review.etl.publication import (
                cache_publication,
                get_cached_title,
            )

        # Cache for GO_REF and Reactome titles
        go_ref_titles: Dict[str, str] = {}
        reactome_titles: Dict[str, str] = {}

        pmids_list = sorted(pmids_to_add)  # Sort for consistency
        for i, pmid in enumerate(pmids_list):
            if pmid not in existing_ref_ids:
                # Try to fetch the actual title (but limit to avoid timeouts)
                title = "TODO: Fetch title"  # Default if fetch fails or skipped

                if fetch_titles and i < max_fetches:
                    if pmid.startswith("PMID:"):
                        pmid_num = pmid[5:]  # Remove 'PMID:' prefix
                        try:
                            # Use publications directory for cache
                            cache_dir = Path("publications")

                            # First try to get just the title from cache (fast)
                            cached_title = get_cached_title(pmid_num, cache_dir)
                            if cached_title:
                                title = cached_title
                            else:
                                # Cache the full publication to publications/PMID_*.md
                                # This is the same as what fetch-gene does
                                if cache_publication(pmid_num, cache_dir):
                                    # Try to get the title from the newly cached file
                                    cached_title = get_cached_title(pmid_num, cache_dir)
                                    if cached_title:
                                        title = cached_title
                        except Exception as e:
                            # If fetch fails, keep the TODO placeholder
                            print(f"Error fetching {pmid}: {e}")
                            pass
                    elif pmid.startswith("GO_REF:"):
                        # Fetch GO_REF title from GO site
                        title = self._get_go_ref_title(pmid, go_ref_titles)
                    elif pmid.startswith("Reactome:"):
                        # Fetch Reactome pathway title
                        title = self._get_reactome_title(pmid, reactome_titles)

                new_ref = {
                    "id": pmid,
                    "title": title,
                    "findings": [],  # Leave empty for AI to fill
                }
                existing_refs.append(new_ref)
                refs_added += 1

        if refs_added > 0:
            yaml_data["references"] = existing_refs
            print(f"    Also seeded {refs_added} references from GOA")

        # Determine output path
        if output_file is None:
            output_file = yaml_file

        # Write the updated YAML
        with open(output_file, "w") as f:
            yaml.dump(
                yaml_data,
                f,
                default_flow_style=False,
                sort_keys=False,
                allow_unicode=True,
            )

        return added_count, output_file, refs_added

    def _get_go_ref_title(self, go_ref_id: str, cache: Dict[str, str]) -> str:
        """Fetch GO_REF title from GO site YAML file.

        Args:
            go_ref_id: The GO_REF identifier (e.g., 'GO_REF:0000033')
            cache: Dictionary to cache fetched titles

        Returns:
            The title of the GO_REF, or a TODO placeholder if not found
        """
        # Check cache first
        if go_ref_id in cache:
            return cache[go_ref_id]

        # Extract the numeric part (e.g., '0000033' from 'GO_REF:0000033')
        ref_num = go_ref_id.replace("GO_REF:", "")

        try:
            # Fetch the GO refs YAML file if not cached
            if not hasattr(self, "_go_refs_data"):
                url = "https://raw.githubusercontent.com/geneontology/go-site/refs/heads/master/metadata/gorefs.yaml"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    self._go_refs_data = yaml.safe_load(response.text)
                else:
                    self._go_refs_data = []

            # Search for the matching GO_REF
            for ref in self._go_refs_data:
                if isinstance(ref, dict):
                    ref_id = ref.get("id", "")
                    # Match either format: 'GO_REF:0000033' or just '0000033'
                    if ref_id == go_ref_id or ref_id == ref_num:
                        title = ref.get("title", "TODO: Fetch title")
                        cache[go_ref_id] = title
                        return title
        except Exception:
            # If fetch fails, return TODO placeholder
            pass

        # Not found, cache and return placeholder
        placeholder = "TODO: Fetch title"
        cache[go_ref_id] = placeholder
        return placeholder

    def _get_reactome_title(self, reactome_id: str, cache: Dict[str, str]) -> str:
        """Fetch Reactome pathway title from Reactome API or cache.

        Args:
            reactome_id: The Reactome identifier (e.g., 'Reactome:R-HSA-6785807')
            cache: Dictionary to cache fetched titles

        Returns:
            The title of the Reactome pathway, or a TODO placeholder if not found
        """
        # Check cache first
        if reactome_id in cache:
            return cache[reactome_id]

        # Import Reactome fetcher
        from ai_gene_review.etl.reactome import fetch_reactome_data, get_cached_title

        # Extract the ID part (e.g., 'R-HSA-6785807' from 'Reactome:R-HSA-6785807')
        clean_id = reactome_id
        if ":" in clean_id:
            clean_id = clean_id.split(":", 1)[1].strip()

        try:
            # Use reactome directory for cache
            cache_dir = Path("reactome")

            # First try to get just the title from cache (fast)
            cached_title = get_cached_title(clean_id, cache_dir)
            if cached_title:
                title = cached_title
            else:
                # Fetch full pathway data and cache it
                pathway = fetch_reactome_data(
                    clean_id, use_cache=True, cache_dir=cache_dir
                )
                if pathway and pathway.display_name:
                    title = pathway.display_name
                else:
                    title = "TODO: Fetch title"
        except Exception:
            # If fetch fails, keep the TODO placeholder
            title = "TODO: Fetch title"

        # Cache and return
        cache[reactome_id] = title
        return title
