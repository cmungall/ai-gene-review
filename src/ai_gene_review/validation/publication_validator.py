"""Validator for publication references to prevent hallucination."""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

from ai_gene_review.etl.publication import extract_pmid, fetch_pubmed_data


@dataclass
class PublicationValidationResult:
    """Result of validating a single publication reference."""

    pmid: str
    provided_title: Optional[str]
    correct_title: Optional[str]
    is_valid: bool
    error_message: Optional[str] = None
    path: Optional[str] = None
    from_cache: bool = False


@dataclass
class PublicationValidator:
    """Validates publication PMID/title pairs against actual publication data.

    This validator helps prevent hallucination by ensuring that when a
    publication is referenced by PMID, its title matches the actual
    publication title.

    It first checks the local cache in publications/ directory, then
    falls back to fetching from PubMed if not cached.

    Example:
        >>> validator = PublicationValidator()
        >>> result = validator.validate_publication(
        ...     "PMID:29727692",
        ...     "C11orf70 Mutations Disrupting the Intraflagellar Transport-Dependent Assembly of Multiple Axonemal Dyneins Cause Primary Ciliary Dyskinesia."
        ... )
        >>> result.is_valid
        True
    """

    # Cache directory for publications
    cache_dir: Path = field(default_factory=lambda: Path("publications"))

    # In-memory cache for publication data
    _title_cache: Dict[str, Optional[str]] = field(default_factory=dict)

    # Whether to use cache
    use_cache: bool = True

    # Whether to auto-fetch missing publications
    auto_fetch: bool = True

    def validate_publication(
        self,
        pmid_str: str,
        provided_title: Optional[str] = None,
        path: Optional[str] = None,
    ) -> PublicationValidationResult:
        """Validate a single publication reference.

        Args:
            pmid_str: The PMID (e.g., "PMID:12345" or "12345")
            provided_title: The title provided by the user (optional)
            path: Path in the document where this reference appears

        Returns:
            PublicationValidationResult with validation details
        """
        # Extract clean PMID
        pmid = extract_pmid(pmid_str)

        if not pmid or not pmid.isdigit():
            return PublicationValidationResult(
                pmid=pmid_str,
                provided_title=provided_title,
                correct_title=None,
                is_valid=False,
                error_message=f"Invalid PMID format: {pmid_str}",
                path=path,
            )

        # Get the correct title
        correct_title, from_cache = self._get_publication_title(pmid)

        if correct_title is None:
            return PublicationValidationResult(
                pmid=pmid,
                provided_title=provided_title,
                correct_title=None,
                is_valid=False,
                error_message=f"Could not verify PMID:{pmid} (publication not found)",
                path=path,
                from_cache=from_cache,
            )

        # If no title was provided, we just validate the PMID exists
        if provided_title is None:
            return PublicationValidationResult(
                pmid=pmid,
                provided_title=None,
                correct_title=correct_title,
                is_valid=True,
                path=path,
                from_cache=from_cache,
            )

        # Compare titles (normalize for comparison)
        is_valid = self._titles_match(provided_title, correct_title)
        error_message = None

        if not is_valid:
            # Check if it's a substring match (sometimes titles are truncated)
            if self._is_substring_match(provided_title, correct_title):
                # Warning rather than error for partial matches
                error_message = f"Partial title match for PMID:{pmid}: provided title appears truncated"
                is_valid = True  # Consider partial matches as valid
            else:
                error_message = f"Title mismatch for PMID:{pmid}"
                if len(provided_title) < 50 and len(correct_title) > 100:
                    error_message += " (provided title may be truncated)"

        return PublicationValidationResult(
            pmid=pmid,
            provided_title=provided_title,
            correct_title=correct_title,
            is_valid=is_valid,
            error_message=error_message,
            path=path,
            from_cache=from_cache,
        )

    def validate_publications_in_data(
        self, data: Any, path: str = "", invalid_pmids: Optional[set] = None
    ) -> List[PublicationValidationResult]:
        """Recursively validate all publication references in a data structure.

        Args:
            data: The data structure to validate
            path: Current path in the data structure
            invalid_pmids: Set of PMIDs marked as invalid (collected on first pass)

        Returns:
            List of validation results for all publications found
        """
        # First pass - collect all invalid PMIDs from references section
        if invalid_pmids is None:
            invalid_pmids = set()
            if isinstance(data, dict) and "references" in data:
                for ref in data.get("references", []):
                    if isinstance(ref, dict) and ref.get("is_invalid", False):
                        ref_id = ref.get("id", "")
                        if ref_id:
                            # Add both with and without prefix
                            invalid_pmids.add(ref_id)
                            if ref_id.startswith("PMID:"):
                                invalid_pmids.add(ref_id[5:])
                            else:
                                invalid_pmids.add(f"PMID:{ref_id}")

        results = []

        if isinstance(data, dict):
            # Check if this is a reference entry
            if "id" in data and isinstance(data["id"], str):
                # Check if it's a PMID
                if data["id"].startswith("PMID") or data["id"].isdigit():
                    # Skip validation if marked as invalid
                    if data.get("is_invalid", False):
                        # Don't validate invalid references
                        pass
                    else:
                        result = self.validate_publication(
                            data["id"], data.get("title"), path
                        )
                        results.append(result)

            # Also check for original_reference_id fields
            if "original_reference_id" in data:
                ref_id = data["original_reference_id"]
                if ref_id and (
                    ref_id.startswith("PMID")
                    or (isinstance(ref_id, str) and ref_id.isdigit())
                ):
                    # Skip validation if this PMID is marked as invalid in references
                    if ref_id in invalid_pmids:
                        # Don't validate PMIDs marked as invalid
                        pass
                    else:
                        # For original_reference_id, we usually don't have a title to validate
                        # Just check that the PMID exists
                        result = self.validate_publication(
                            ref_id,
                            None,  # No title to validate
                            f"{path}.original_reference_id"
                            if path
                            else "original_reference_id",
                        )
                        results.append(result)

            # Recurse into dict values
            for key, value in data.items():
                new_path = f"{path}.{key}" if path else key
                results.extend(
                    self.validate_publications_in_data(value, new_path, invalid_pmids)
                )

        elif isinstance(data, list):
            # Recurse into list items
            for i, item in enumerate(data):
                new_path = f"{path}[{i}]"
                results.extend(
                    self.validate_publications_in_data(item, new_path, invalid_pmids)
                )

        return results

    def _get_publication_title(self, pmid: str) -> tuple[Optional[str], bool]:
        """Get the correct title for a publication.

        Args:
            pmid: Clean PMID (without prefix)

        Returns:
            Tuple of (title, from_cache) or (None, False) if not found
        """
        # Check in-memory cache first
        if self.use_cache and pmid in self._title_cache:
            return self._title_cache[pmid], True

        # Check local file cache
        cached_file = self.cache_dir / f"PMID_{pmid}.md"
        if cached_file.exists():
            try:
                content = cached_file.read_text()
                # Extract title from frontmatter
                if content.startswith("---"):
                    # Find the end of frontmatter
                    end_marker = content.find("---", 3)
                    if end_marker > 0:
                        frontmatter_str = content[3:end_marker]
                        frontmatter = yaml.safe_load(frontmatter_str)
                        if frontmatter and "title" in frontmatter:
                            title = frontmatter["title"]
                            # Cache in memory
                            if self.use_cache:
                                self._title_cache[pmid] = title
                            return title, True
            except Exception as e:
                print(f"Error reading cached publication {pmid}: {e}")

        # If not in cache and auto_fetch is enabled, fetch from PubMed
        if self.auto_fetch:
            try:
                pub = fetch_pubmed_data(pmid)
                if pub:
                    title = pub.title
                    # Cache in memory
                    if self.use_cache:
                        self._title_cache[pmid] = title
                    # Optionally cache to file
                    # (not doing this automatically to avoid unexpected side effects)
                    return title, False
            except Exception as e:
                print(f"Error fetching publication {pmid}: {e}")

        return None, False

    def _titles_match(self, title1: str, title2: str) -> bool:
        """Check if two titles match (with normalization).

        Args:
            title1: First title
            title2: Second title

        Returns:
            True if titles match (ignoring case, punctuation, and whitespace)
        """

        def normalize(s: str) -> str:
            # Remove punctuation and extra whitespace
            s = re.sub(r"[^\w\s]", " ", s.lower())
            s = re.sub(r"\s+", " ", s)
            return s.strip()

        return normalize(title1) == normalize(title2)

    def _is_substring_match(self, short_title: str, full_title: str) -> bool:
        """Check if short_title is a substring of full_title.

        Args:
            short_title: Potentially truncated title
            full_title: Full title

        Returns:
            True if short_title is a significant substring of full_title
        """

        def normalize(s: str) -> str:
            s = re.sub(r"[^\w\s]", " ", s.lower())
            s = re.sub(r"\s+", " ", s)
            return s.strip()

        short_norm = normalize(short_title)
        full_norm = normalize(full_title)

        # Check if it's a prefix (common for truncated titles)
        if full_norm.startswith(short_norm) and len(short_norm) >= 20:
            return True

        # Check if it's contained within (less common)
        if short_norm in full_norm and len(short_norm) >= 30:
            return True

        return False


def mark_invalid_pmids(
    yaml_file: Path, pmids: List[str], output_file: Optional[Path] = None
) -> int:
    """Mark specific PMIDs as invalid in a YAML file.

    This function adds is_invalid: true to references that can't be retrieved,
    so they won't be validated in future runs.

    Args:
        yaml_file: Path to the YAML file to update
        pmids: List of PMID strings (e.g., ["PMID:12345", "34521819"])
        output_file: Optional output path (defaults to input file)

    Returns:
        Number of references marked as invalid
    """
    import yaml

    # Normalize PMIDs
    normalized_pmids = set()
    for pmid in pmids:
        if pmid.startswith("PMID:"):
            normalized_pmids.add(pmid)
            normalized_pmids.add(pmid[5:])  # Also add without prefix
        else:
            normalized_pmids.add(f"PMID:{pmid}")
            normalized_pmids.add(pmid)

    # Load YAML file
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)

    marked_count = 0

    # Update references section
    if "references" in data and isinstance(data["references"], list):
        for ref in data["references"]:
            if isinstance(ref, dict) and "id" in ref:
                ref_id = ref["id"]
                # Check if this PMID should be marked invalid
                if ref_id in normalized_pmids or (
                    ref_id.startswith("PMID:") and ref_id[5:] in normalized_pmids
                ):
                    if not ref.get("is_invalid", False):
                        ref["is_invalid"] = True
                        marked_count += 1
                        # Add a note about why it's invalid
                        if "title" in ref and ref["title"] == "TODO: Fetch title":
                            ref["title"] = (
                                "INVALID: Cannot retrieve from PubMed (possibly replaced or retracted)"
                            )

    # Save the updated YAML
    if output_file is None:
        output_file = yaml_file

    with open(output_file, "w") as f:
        yaml.dump(
            data, f, default_flow_style=False, sort_keys=False, allow_unicode=True
        )

    return marked_count


def validate_yaml_file_publications(
    yaml_file: Path,
) -> tuple[bool, List[PublicationValidationResult]]:
    """Validate all publication references in a YAML file.

    Args:
        yaml_file: Path to the YAML file to validate

    Returns:
        Tuple of (all_valid, list_of_results)
    """
    with open(yaml_file) as f:
        data = yaml.safe_load(f)

    validator = PublicationValidator()
    results = validator.validate_publications_in_data(data)

    all_valid = all(r.is_valid for r in results)

    return all_valid, results
