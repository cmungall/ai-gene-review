"""Validator for ontology terms to prevent hallucination."""

import pickle
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import requests
import yaml
from oaklib import get_adapter  # type: ignore[import-untyped]
from oaklib.datamodels.vocabulary import IS_A  # type: ignore[import-untyped]


@dataclass
class TermValidationResult:
    """Result of validating a single term."""

    term_id: str
    provided_label: Optional[str]
    correct_label: Optional[str]
    is_valid: bool
    is_obsolete: bool = False
    error_message: Optional[str] = None
    path: Optional[str] = None


@dataclass
class TermValidator:
    """Validates ontology term ID/label pairs against actual ontology data.

    This validator helps prevent hallucination by ensuring that when an
    ontology term ID is provided with a label, the label matches exactly
    what's in the ontology.

    Example:
        >>> validator = TermValidator()
        >>> result = validator.validate_term("GO:0005515", "protein binding")
        >>> result.is_valid
        True
        >>> result = validator.validate_term("GO:0005515", "wrong label")
        >>> result.is_valid
        False
    """

    # Cache for term labels to avoid repeated lookups
    _label_cache: Dict[str, Optional[str]] = field(default_factory=dict)
    _obsolete_cache: Dict[str, bool] = field(default_factory=dict)
    _go_adapter: Optional[Any] = field(default=None, init=False)
    _go_branch_cache: Dict[str, Set[str]] = field(default_factory=dict)
    _ontology_adapters: Dict[str, Any] = field(default_factory=dict)

    # OLS API base URL (can be overridden for testing)
    ols_base_url: str = "https://www.ebi.ac.uk/ols/api"

    # Timeout for API requests
    timeout: int = 10

    # Whether to use cache
    use_cache: bool = True

    # Known ontology prefixes and their OLS ontology IDs
    ONTOLOGY_MAP: Dict[str, str] = field(
        default_factory=lambda: {
            "GO": "go",  # Gene Ontology
            "HP": "hp",  # Human Phenotype Ontology
            "MONDO": "mondo",  # Mondo Disease Ontology
            "CL": "cl",  # Cell Ontology
            "UBERON": "uberon",  # Uberon Anatomy Ontology
            "CHEBI": "chebi",  # Chemical Entities of Biological Interest
            "PR": "pr",  # Protein Ontology
            "SO": "so",  # Sequence Ontology
            "PATO": "pato",  # Phenotype And Trait Ontology
            "NCBITaxon": "ncbitaxon",  # NCBI Taxonomy
        }
    )

    # Known non-ontology ID prefixes that should be ignored
    EXCLUDED_PREFIXES: Set[str] = field(
        default_factory=lambda: {
            "PMID",  # PubMed ID
            "PMC",  # PubMed Central ID
            "DOI",  # Digital Object Identifier
            "ISBN",  # Book identifier
            "ISSN",  # Journal identifier
            "UniProt",  # UniProt accession
            "UniProtKB",  # UniProt Knowledgebase accession
            "RefSeq",  # RefSeq accession
            "ENSEMBL",  # Ensembl ID
            "EC",  # Enzyme Commission number
            "PDB",  # Protein Data Bank
            "OMIM",  # Online Mendelian Inheritance in Man
            "TEMP",  # Temporary/placeholder identifiers
            "file",  # File path references
            "GO_REF",  # Gene Ontology reference collection
            "Reactome",  # Reactome pathway database
        }
    )

    # GO branch roots
    GO_MF_ROOT = "GO:0003674"  # molecular_function
    GO_BP_ROOT = "GO:0008150"  # biological_process
    GO_CC_ROOT = "GO:0005575"  # cellular_component

    def validate_term(
        self,
        term_id: str,
        provided_label: Optional[str] = None,
        path: Optional[str] = None,
    ) -> TermValidationResult:
        """Validate a single ontology term.

        Args:
            term_id: The ontology term ID (e.g., "GO:0005515")
            provided_label: The label provided by the user (optional)
            path: Path in the document where this term appears (for error reporting)

        Returns:
            TermValidationResult with validation details
        """
        # Check if term ID format is valid
        if not self._is_valid_term_format(term_id):
            return TermValidationResult(
                term_id=term_id,
                provided_label=provided_label,
                correct_label=None,
                is_valid=False,
                error_message=f"Invalid term ID format: {term_id}",
                path=path,
            )

        # Get the correct label from the ontology
        correct_label = self._get_term_label(term_id)

        if correct_label is None:
            # Could not verify term (might not exist or network error)
            return TermValidationResult(
                term_id=term_id,
                provided_label=provided_label,
                correct_label=None,
                is_valid=False,
                error_message=f"Term {term_id} not found in ontology",
                path=path,
            )

        # Check if term is obsolete
        is_obsolete = self._is_term_obsolete(term_id)

        # If no label was provided, we just validate the ID exists
        if provided_label is None:
            return TermValidationResult(
                term_id=term_id,
                provided_label=None,
                correct_label=correct_label,
                is_valid=True,
                is_obsolete=is_obsolete,
                path=path,
            )

        # Check if provided label matches correct label
        is_valid = provided_label == correct_label
        error_message = None
        if not is_valid:
            error_message = f"Label mismatch for {term_id}: got '{provided_label}', expected '{correct_label}'"

        return TermValidationResult(
            term_id=term_id,
            provided_label=provided_label,
            correct_label=correct_label,
            is_valid=is_valid,
            is_obsolete=is_obsolete,
            error_message=error_message,
            path=path,
        )

    def validate_terms_in_data(
        self, data: Any, path: str = ""
    ) -> List[TermValidationResult]:
        """Recursively validate all ontology terms in a data structure.

        Args:
            data: The data structure to validate (dict, list, or primitive)
            path: Current path in the data structure (for error reporting)

        Returns:
            List of validation results for all terms found
        """
        results = []

        if isinstance(data, dict):
            # Check if this dict represents a term (has id field)
            if "id" in data and isinstance(data["id"], str):
                term_id = data["id"]

                # Check if it looks like a CURIE (contains colon)
                if ":" in term_id:
                    prefix = term_id.split(":")[0]

                    # Special handling for TEMP prefix - requires description
                    if prefix == "TEMP":
                        # Check if description is provided
                        if not data.get("description"):
                            result = TermValidationResult(
                                term_id=term_id,
                                provided_label=data.get("label"),
                                correct_label=None,
                                is_valid=False,
                                error_message=f"TEMP identifier {term_id} requires a 'description' field to explain what this placeholder represents",
                                path=path,
                            )
                            results.append(result)
                    # Skip other known non-ontology identifiers
                    elif prefix in self.EXCLUDED_PREFIXES:
                        pass  # Skip validation for these
                    # If it looks like an ontology term, validate it
                    elif self._is_ontology_term(term_id):
                        result = self.validate_term(term_id, data.get("label"), path)
                        results.append(result)
                    else:
                        # It's a CURIE but not a recognized ontology term or excluded prefix
                        # This could be a made-up identifier
                        # Report as invalid
                        result = TermValidationResult(
                            term_id=term_id,
                            provided_label=data.get("label"),
                            correct_label=None,
                            is_valid=False,
                            error_message=f"Unrecognized identifier prefix '{prefix}' in {term_id}. Valid ontology prefixes are: {', '.join(sorted(self.ONTOLOGY_MAP.keys()))}",
                            path=path,
                        )
                        results.append(result)

            # Recurse into dict values
            for key, value in data.items():
                new_path = f"{path}.{key}" if path else key
                results.extend(self.validate_terms_in_data(value, new_path))

        elif isinstance(data, list):
            # Recurse into list items
            for i, item in enumerate(data):
                new_path = f"{path}[{i}]"
                results.extend(self.validate_terms_in_data(item, new_path))

        return results

    def _is_valid_term_format(self, term_id: str) -> bool:
        """Check if term ID has valid format.

        Examples:
            >>> validator = TermValidator()
            >>> validator._is_valid_term_format("GO:0001234")
            True
            >>> validator._is_valid_term_format("HP:0000001")
            True
            >>> validator._is_valid_term_format("NCBITaxon:9606")
            True
            >>> validator._is_valid_term_format("not_a_term")
            False
        """
        # Check for standard ontology ID patterns
        # Format: PREFIX:NUMBERS or PREFIX_NUMBERS
        pattern = r"^[A-Za-z]+[_:]?\d+$"
        return bool(re.match(pattern, term_id))

    def _is_ontology_term(self, value: Any) -> bool:
        """Check if value looks like an ontology term ID.

        Examples:
            >>> validator = TermValidator()
            >>> validator._is_ontology_term("GO:0001234")
            True
            >>> validator._is_ontology_term("not_a_term")
            False
            >>> validator._is_ontology_term(123)
            False
        """
        if not isinstance(value, str):
            return False

        # Check if it matches a known ontology prefix
        for prefix in self.ONTOLOGY_MAP:
            if value.startswith(f"{prefix}:") or value.startswith(f"{prefix}_"):
                return True
        return False

    def _get_term_label(self, term_id: str) -> Optional[str]:
        """Get the correct label for a term using OAK or cache.

        Args:
            term_id: The ontology term ID

        Returns:
            The correct label, or None if not found
        """
        # Check cache first
        if self.use_cache and term_id in self._label_cache:
            return self._label_cache[term_id]

        # Determine the ontology from the prefix
        prefix = term_id.split(":")[0] if ":" in term_id else term_id.split("_")[0]
        ontology = self.ONTOLOGY_MAP.get(prefix)

        if not ontology:
            return None

        try:
            # Use OAK adapter for the ontology
            adapter = self._get_ontology_adapter(ontology)
            if adapter is None:
                # Fallback to OLS API if adapter can't be loaded
                return self._get_term_label_from_ols(term_id, ontology)

            # Get the label using OAK
            label = adapter.label(term_id)

            # Cache the result
            if self.use_cache and label:
                self._label_cache[term_id] = label

            return label

        except Exception as e:
            # Log error but don't crash
            print(f"Error fetching term {term_id} from OAK: {e}")
            # Try OLS API as fallback
            return self._get_term_label_from_ols(term_id, ontology)

    def _get_term_label_from_ols(self, term_id: str, ontology: str) -> Optional[str]:
        """Fallback method to get term label from OLS API.

        Args:
            term_id: The ontology term ID
            ontology: The ontology name

        Returns:
            The correct label, or None if not found
        """
        try:
            url = f"{self.ols_base_url}/search"
            params = {"q": term_id, "ontology": ontology, "exact": "true", "rows": 1}

            response = requests.get(url, params=params, timeout=self.timeout)  # type: ignore[arg-type]

            if response.status_code == 200:
                data = response.json()
                if data.get("response") and data["response"].get("docs"):
                    docs = data["response"]["docs"]
                    if docs:
                        return docs[0].get("label")
        except Exception:
            pass

        return None

    def _is_term_obsolete(self, term_id: str) -> bool:
        """Check if a term is obsolete.

        Args:
            term_id: The ontology term ID

        Returns:
            True if the term is obsolete, False otherwise
        """
        # Check cache first
        if self.use_cache and term_id in self._obsolete_cache:
            return self._obsolete_cache[term_id]

        # Determine the ontology from the prefix
        prefix = term_id.split(":")[0] if ":" in term_id else term_id.split("_")[0]
        ontology = self.ONTOLOGY_MAP.get(prefix)

        if not ontology:
            return False

        try:
            # Use OAK adapter for the ontology
            adapter = self._get_ontology_adapter(ontology)
            if adapter is None:
                return False  # Can't check, assume not obsolete

            # Check if term is obsolete using OAK
            is_obsolete = adapter.is_obsolete(term_id)

            # Cache the result
            if self.use_cache:
                self._obsolete_cache[term_id] = is_obsolete

            return is_obsolete

        except Exception:
            pass

        return False

    def _get_go_adapter(self):
        """Get or create the GO adapter for oaklib.

        Returns:
            The GO adapter using SQLite backend
        """
        if self._go_adapter is None:
            # Use SQLite backend for GO
            self._go_adapter = get_adapter("sqlite:obo:go")
        return self._go_adapter

    def _get_ontology_adapter(self, ontology_name: str):
        """Get or create an adapter for a specific ontology.

        Args:
            ontology_name: The ontology name (e.g., 'go', 'hp', 'mondo')

        Returns:
            The ontology adapter using SQLite backend
        """
        if ontology_name not in self._ontology_adapters:
            try:
                # Use SQLite backend for the ontology
                self._ontology_adapters[ontology_name] = get_adapter(
                    f"sqlite:obo:{ontology_name}"
                )
            except Exception as e:
                print(f"Could not load adapter for {ontology_name}: {e}")
                return None
        return self._ontology_adapters[ontology_name]

    def _get_cache_dir(self) -> Path:
        """Get the cache directory for validation data.

        Returns:
            Path to cache directory
        """
        cache_dir = Path.home() / ".cache" / "ai-gene-review"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def _load_go_branch_cache(self) -> Dict[str, Set[str]]:
        """Load cached GO branch descendants from disk.

        Returns:
            Dictionary of branch root to descendant sets
        """
        cache_file = self._get_cache_dir() / "go_branch_descendants.pkl"
        if cache_file.exists():
            try:
                with open(cache_file, "rb") as f:
                    return pickle.load(f)
            except Exception as e:
                print(f"Could not load GO branch cache: {e}")
        return {}

    def _save_go_branch_cache(self, cache: Dict[str, Set[str]]) -> None:
        """Save GO branch descendants cache to disk.

        Args:
            cache: Dictionary of branch root to descendant sets
        """
        cache_file = self._get_cache_dir() / "go_branch_descendants.pkl"
        try:
            with open(cache_file, "wb") as f:
                pickle.dump(cache, f)
        except Exception as e:
            print(f"Could not save GO branch cache: {e}")

    def _get_go_branch_descendants(self, root: str) -> Set[str]:
        """Get all descendants of a GO branch root.

        Args:
            root: The GO root term (e.g., GO:0003674 for molecular_function)

        Returns:
            Set of all descendant GO terms
        """
        # Check in-memory cache first
        if root in self._go_branch_cache:
            return self._go_branch_cache[root]

        # Try to load from disk cache
        if not self._go_branch_cache:
            self._go_branch_cache = self._load_go_branch_cache()
            if root in self._go_branch_cache:
                return self._go_branch_cache[root]

        # If not cached, compute descendants
        try:
            print(
                f"Computing GO descendants for {root} (this may take a moment, but will be cached)..."
            )
            adapter = self._get_go_adapter()
            # Get all descendants using IS_A relationship
            descendants = set(adapter.descendants(root, predicates=[IS_A]))
            # Include the root itself
            descendants.add(root)

            # Cache the result in memory
            self._go_branch_cache[root] = descendants

            # Save to disk for future runs
            self._save_go_branch_cache(self._go_branch_cache)

            print(f"Cached {len(descendants)} descendants for {root}")
            return descendants
        except Exception as e:
            print(f"Error getting GO descendants for {root}: {e}")
            # Return empty set if there's an error
            return set()

    def _check_go_branch(self, term_id: str, expected_root: str) -> Optional[str]:
        """Check if a GO term is in the expected branch.

        Args:
            term_id: The GO term to check
            expected_root: The expected branch root (e.g., GO:0003674 for MF)

        Returns:
            Error message if term is in wrong branch, None if correct
        """
        if not term_id.startswith("GO:"):
            return None  # Not a GO term, skip check

        # Get descendants of the expected root
        expected_descendants = self._get_go_branch_descendants(expected_root)

        if term_id not in expected_descendants:
            # Determine which branch it's actually in
            actual_branch = "unknown"
            if term_id in self._get_go_branch_descendants(self.GO_MF_ROOT):
                actual_branch = "molecular_function"
            elif term_id in self._get_go_branch_descendants(self.GO_BP_ROOT):
                actual_branch = "biological_process"
            elif term_id in self._get_go_branch_descendants(self.GO_CC_ROOT):
                actual_branch = "cellular_component"

            expected_branch = {
                self.GO_MF_ROOT: "molecular_function",
                self.GO_BP_ROOT: "biological_process",
                self.GO_CC_ROOT: "cellular_component",
            }.get(expected_root, "unknown")

            return f"GO term {term_id} is in the {actual_branch} branch but should be in the {expected_branch} branch"

        return None

    def validate_terms_in_core_functions(
        self, data: Any, path: str = ""
    ) -> List[TermValidationResult]:
        """Special validation for terms in core_functions with GO branch checking.

        Args:
            data: The data structure to validate
            path: Current path in the data structure

        Returns:
            List of validation results
        """
        results = []

        if isinstance(data, dict) and "core_functions" in data:
            for i, func in enumerate(data.get("core_functions", [])):
                if not isinstance(func, dict):
                    continue

                func_path = (
                    f"{path}.core_functions[{i}]" if path else f"core_functions[{i}]"
                )

                # Check molecular_function - should be in MF branch
                if "molecular_function" in func:
                    mf = func["molecular_function"]
                    if isinstance(mf, dict) and "id" in mf:
                        mf_id = mf["id"]
                        mf_path = f"{func_path}.molecular_function"

                        # First do normal term validation
                        result = self.validate_term(mf_id, mf.get("label"), mf_path)
                        if result.is_valid:
                            # Then check branch
                            branch_error = self._check_go_branch(mf_id, self.GO_MF_ROOT)
                            if branch_error:
                                result.is_valid = False
                                result.error_message = branch_error
                        results.append(result)

                # Check directly_involved_in - should be in BP branch
                if "directly_involved_in" in func:
                    for j, bp_term in enumerate(func.get("directly_involved_in", [])):
                        if isinstance(bp_term, dict) and "id" in bp_term:
                            bp_id = bp_term["id"]
                            bp_path = f"{func_path}.directly_involved_in[{j}]"

                            # First do normal term validation
                            result = self.validate_term(
                                bp_id, bp_term.get("label"), bp_path
                            )
                            if result.is_valid:
                                # Then check branch
                                branch_error = self._check_go_branch(
                                    bp_id, self.GO_BP_ROOT
                                )
                                if branch_error:
                                    result.is_valid = False
                                    result.error_message = branch_error
                            results.append(result)

                # Check locations - should be in CC branch
                if "locations" in func:
                    for j, cc_term in enumerate(func.get("locations", [])):
                        if isinstance(cc_term, dict) and "id" in cc_term:
                            cc_id = cc_term["id"]
                            cc_path = f"{func_path}.locations[{j}]"

                            # First do normal term validation
                            result = self.validate_term(
                                cc_id, cc_term.get("label"), cc_path
                            )
                            if result.is_valid:
                                # Then check branch
                                branch_error = self._check_go_branch(
                                    cc_id, self.GO_CC_ROOT
                                )
                                if branch_error:
                                    result.is_valid = False
                                    result.error_message = branch_error
                            results.append(result)

        return results


def validate_yaml_file(yaml_file: Path) -> Tuple[bool, List[TermValidationResult]]:
    """Validate all ontology terms in a YAML file.

    Args:
        yaml_file: Path to the YAML file to validate

    Returns:
        Tuple of (all_valid, list_of_results)

    Example:
        >>> # Would require a real file
        >>> # all_valid, results = validate_yaml_file(Path("test.yaml"))
        >>> # for result in results:
        >>> #     if not result.is_valid:
        >>> #         print(result.error_message)
    """
    with open(yaml_file) as f:
        data = yaml.safe_load(f)

    validator = TermValidator()
    results = validator.validate_terms_in_data(data)

    all_valid = all(r.is_valid for r in results)

    return all_valid, results
