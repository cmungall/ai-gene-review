"""Gene data ETL module for fetching UniProt and GOA data.

This module provides functionality to fetch gene data from UniProt and GOA APIs
and save them in a structured directory format.

Example:
    >>> from ai_gene_review.etl.gene import fetch_gene_data
    >>> fetch_gene_data(("human", "CFAP300"), uniprot_id="Q9BRQ4")  # doctest: +SKIP

    This creates:
    genes/
      human/
        CFAP300/
          CFAP300-uniprot.txt
          CFAP300-goa.csv
"""

from pathlib import Path
from typing import Tuple, Optional, List, Dict, Any
import requests
import yaml


def fetch_gene_data(
    gene_info: Tuple[str, str],
    uniprot_id: Optional[str] = None,
    base_path: Optional[Path] = None,
    seed_annotations: bool = True,
    fetch_titles: bool = True,
    alias: Optional[str] = None,
) -> Dict[str, Any]:
    """Fetch gene data from UniProt and GOA APIs and save to files.

    Creates/updates an ai-review.yaml file with GOA annotations if it doesn't exist
    or is missing annotations.

    Args:
        gene_info: Tuple of (organism, gene_name) e.g. ("human", "CFAP300")
        uniprot_id: Optional UniProt accession ID. If not provided, will attempt to resolve.
        base_path: Base directory for output files. Defaults to current directory.
        seed_annotations: If True, creates/seeds ai-review.yaml with GOA annotations.
        fetch_titles: If True, fetch actual titles from PubMed when seeding (default: True).
        alias: Optional alias to use for directory name and file prefixes instead of gene_name.

    Returns:
        Dictionary with status information:
            - yaml_created: bool - True if ai-review.yaml was created
            - yaml_existed: bool - True if ai-review.yaml already existed
            - annotations_added: int - Number of annotations added
            - references_added: int - Number of references added

    Raises:
        ValueError: If UniProt ID cannot be resolved or data cannot be fetched.

    Example:
        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as tmpdir:
        ...     base = Path(tmpdir)
        ...     result = fetch_gene_data(("human", "TP53"), uniprot_id="P04637", base_path=base)  # doctest: +SKIP
        ...     assert (base / "genes" / "human" / "TP53").exists()  # doctest: +SKIP
    """
    organism, gene_name = gene_info
    
    # Use alias for directory name and file prefixes if provided
    dir_name = alias if alias else gene_name
    file_prefix = alias if alias else gene_name

    if base_path is None:
        base_path = Path.cwd()

    # Create directory structure
    gene_dir = base_path / "genes" / organism / dir_name
    gene_dir.mkdir(parents=True, exist_ok=True)

    # Resolve UniProt ID if not provided
    if uniprot_id is None:
        uniprot_id = resolve_gene_to_uniprot(gene_name, organism)

    # Fetch and save UniProt data
    uniprot_data = fetch_uniprot_data(uniprot_id)
    uniprot_file = gene_dir / f"{file_prefix}-uniprot.txt"
    uniprot_file.write_text(uniprot_data)

    # Fetch and save GOA data
    goa_data = fetch_goa_data(uniprot_id)
    goa_file = gene_dir / f"{file_prefix}-goa.tsv"
    goa_file.write_text(goa_data)

    # Initialize result status
    result = {
        "yaml_created": False,
        "yaml_existed": False,
        "annotations_added": 0,
        "references_added": 0,
    }

    # Seed ai-review.yaml with GOA annotations if requested
    if seed_annotations:
        yaml_file = gene_dir / f"{file_prefix}-ai-review.yaml"

        # Import here to avoid circular dependency
        from ai_gene_review.validation.goa_validator import GOAValidator

        # Check if file already exists
        yaml_existed = yaml_file.exists()
        result["yaml_existed"] = yaml_existed

        # Get taxon ID and proper label from organism
        organism_to_taxon = {
            "human": ("NCBITaxon:9606", "Homo sapiens"),
            "mouse": ("NCBITaxon:10090", "Mus musculus"),
            "rat": ("NCBITaxon:10116", "Rattus norvegicus"),
            "yeast": ("NCBITaxon:559292", "Saccharomyces cerevisiae"),
            "fly": ("NCBITaxon:7227", "Drosophila melanogaster"),
            "worm": ("NCBITaxon:6239", "Caenorhabditis elegans"),
            "zebrafish": ("NCBITaxon:7955", "Danio rerio"),
        }

        # Check if we have a predefined mapping
        if organism.lower() in organism_to_taxon:
            taxon_info = organism_to_taxon[organism.lower()]
        # Check if it's a UniProt organism code
        elif organism.isupper() and len(organism) <= 5:
            taxon_id = resolve_organism_code_to_taxon(organism)
            if taxon_id:
                # Get organism name from UniProt data if available
                organism_name = get_organism_name_from_uniprot(uniprot_id) or organism
                taxon_info = (f"NCBITaxon:{taxon_id}", organism_name)
            else:
                # Default fallback
                taxon_info = (f"NCBITaxon:{organism}", organism.capitalize())
        else:
            # Default for unknown organisms
            taxon_info = (f"NCBITaxon:{organism}", organism.capitalize())
        taxon_id, taxon_label = taxon_info

        # Create minimal YAML structure if file doesn't exist
        if not yaml_existed:
            yaml_data = {
                "id": uniprot_id,
                "gene_symbol": gene_name,
                "taxon": {"id": taxon_id, "label": taxon_label},
                "description": f"TODO: Add description for {gene_name}",
            }

            # Write initial YAML
            with open(yaml_file, "w") as f:
                yaml.dump(
                    yaml_data,
                    f,
                    default_flow_style=False,
                    sort_keys=False,
                    allow_unicode=True,
                )
            result["yaml_created"] = True

        # Now seed missing GOA annotations
        validator = GOAValidator()
        added_count, _, refs_added = validator.seed_missing_annotations(
            yaml_file, goa_file, fetch_titles=fetch_titles
        )
        result["annotations_added"] = added_count
        result["references_added"] = refs_added

        if added_count > 0:
            print(
                f"  ✓ Seeded {added_count} GOA annotations in {file_prefix}-ai-review.yaml"
            )
        else:
            if yaml_existed:
                print(
                    f"  - {file_prefix}-ai-review.yaml already contains all GOA annotations"
                )

    return result


def expand_organism_name(organism: str) -> str:
    """Expand organism code or common name to full scientific name.

    Args:
        organism: Organism code (e.g., "PSEPK"), common name (e.g., "human"),
                  or taxon ID

    Returns:
        Full scientific name of the organism, or original value if not found

    Examples:
        >>> expand_organism_name("human")
        'Homo sapiens'
        >>> expand_organism_name("PSEPK")  # doctest: +SKIP
        'Pseudomonas putida (strain ATCC 47054 / DSM 6125 / CFBP 8728 / NCIMB 11950 / KT2440)'
        >>> expand_organism_name("Homo sapiens")
        'Homo sapiens'
    """
    # Common name mappings
    common_to_scientific = {
        "human": "Homo sapiens",
        "mouse": "Mus musculus",
        "rat": "Rattus norvegicus",
        "yeast": "Saccharomyces cerevisiae",
        "fly": "Drosophila melanogaster",
        "worm": "Caenorhabditis elegans",
        "zebrafish": "Danio rerio",
        "chicken": "Gallus gallus",
        "pig": "Sus scrofa",
        "dog": "Canis lupus familiaris",
        "cow": "Bos taurus",
        "ecoli": "Escherichia coli",
        "e.coli": "Escherichia coli",
        "e coli": "Escherichia coli",
    }

    # Check common names first
    if organism.lower() in common_to_scientific:
        return common_to_scientific[organism.lower()]

    # Check if it's a UniProt organism code (uppercase, <= 5 chars)
    if organism.isupper() and len(organism) <= 5:
        try:
            # Query UniProt taxonomy API
            url = "https://rest.uniprot.org/taxonomy/stream"
            params = {"query": f"mnemonic:{organism}", "format": "json", "size": "1"}
            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()

            results = response.json().get("results", [])
            if results:
                return results[0].get("scientificName", organism)
        except Exception:
            pass

    # Check if it's a taxon ID (all digits)
    if organism.isdigit():
        try:
            # Query UniProt taxonomy API with taxon ID
            url = "https://rest.uniprot.org/taxonomy/stream"
            params = {"query": f"id:{organism}", "format": "json", "size": "1"}
            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()

            results = response.json().get("results", [])
            if results:
                return results[0].get("scientificName", organism)
        except Exception:
            pass

    # Return original if already looks like scientific name or can't expand
    return organism


def get_organism_name_from_uniprot(uniprot_id: str) -> Optional[str]:
    """Get the organism scientific name from a UniProt ID.

    Args:
        uniprot_id: UniProt accession ID

    Returns:
        Scientific name of the organism, or None if not found
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?format=json&fields=organism_name"
        response = requests.get(url, timeout=10)
        response.raise_for_status()

        data = response.json()
        return data.get("organism", {}).get("scientificName")
    except Exception:
        return None


def resolve_organism_code_to_taxon(organism_code: str) -> Optional[str]:
    """Resolve a UniProt organism code to an NCBI taxonomy ID.

    Args:
        organism_code: UniProt organism code (e.g., "PSEPK", "ECOLI")

    Returns:
        NCBI taxonomy ID as string, or None if not found

    Example:
        >>> taxon_id = resolve_organism_code_to_taxon("ECOLI")
        >>> taxon_id  # doctest: +SKIP
        '511145'
    """
    try:
        # Query UniProt taxonomy API with the mnemonic
        url = "https://rest.uniprot.org/taxonomy/stream"
        params = {"query": f"mnemonic:{organism_code}", "format": "json", "size": "1"}
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()

        results = response.json().get("results", [])
        if results:
            return str(results[0].get("taxonId"))

    except Exception:
        # Silently fail and return None
        pass

    return None


def swissprot_by_gene_taxon(
    gene_symbol: str, organism: str, limit: int = 10
) -> List[Dict[str, str]]:
    """Look up reviewed (Swiss-Prot) UniProtKB entries by gene symbol and organism.

    Wrapper for backwards compatibility.

    Args:
        gene_symbol: Gene symbol to search for
        organism: Organism name (e.g., "human", "mouse", "yeast") or NCBI taxon ID
        limit: Maximum number of results to return

    Returns:
        List of dictionaries with UniProt entry information

    Raises:
        ValueError: If API request fails
    """
    return uniprot_by_gene_taxon(gene_symbol, organism, limit, reviewed_only=True)


def uniprot_by_gene_taxon(
    gene_symbol: str, organism: str, limit: int = 10, reviewed_only: bool = True
) -> List[Dict[str, str]]:
    """Look up UniProtKB entries by gene symbol and organism.

    Args:
        gene_symbol: Gene symbol to search for
        organism: Organism name (e.g., "human", "mouse"), UniProt organism code (e.g., "PSEPK"), or NCBI taxon ID
        limit: Maximum number of results to return
        reviewed_only: If True, only return Swiss-Prot (reviewed) entries. If False, include TrEMBL.

    Returns:
        List of dictionaries with UniProt entry information

    Raises:
        ValueError: If API request fails
    """
    # Map common organism names to NCBI taxonomy IDs
    organism_to_taxon = {
        "human": "9606",
        "mouse": "10090",
        "rat": "10116",
        "yeast": "559292",
        "fly": "7227",
        "worm": "6239",
        "zebrafish": "7955",
        "chicken": "9031",
        "pig": "9823",
        "dog": "9615",
        "cow": "9913",
    }

    # Determine if we have a taxon ID or need to look it up
    taxon = organism_to_taxon.get(organism.lower(), organism)

    # Check if it's already a taxon ID (all digits)
    if taxon.isdigit():
        org_field = f"organism_id:{taxon}"
    # Check if it's a UniProt organism code (uppercase letters possibly with digits)
    elif organism.isupper() and len(organism) <= 5:
        # Look up taxon ID from organism code
        taxon_id = resolve_organism_code_to_taxon(organism)
        if taxon_id:
            org_field = f"organism_id:{taxon_id}"
        else:
            # Fall back to organism name search
            org_field = f'organism_name:"{organism}"'
    else:
        # Try as organism name
        org_field = f'organism_name:"{organism}"'

    # Build query
    reviewed_filter = " AND (reviewed:true)" if reviewed_only else ""
    query = f"(gene_exact:{gene_symbol}) AND ({org_field}){reviewed_filter}"

    # API parameters
    params = {
        "query": query,
        "format": "json",
        "fields": "accession,id,gene_primary,organism_name,reviewed",
        "size": str(limit),
        "compressed": "false",
        "download": "true",
    }

    # Make request
    url = "https://rest.uniprot.org/uniprotkb/stream"
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()

    # Parse results
    hits = []
    for rec in response.json().get("results", []):
        hits.append(
            {
                "accession": rec["primaryAccession"],
                "entry_name": rec.get("uniProtkbId", ""),
                "gene": next(
                    (
                        g.get("geneName", {}).get("value")
                        for g in rec.get("genes", [])
                        if g.get("geneName")
                    ),
                    None,
                ),
                "organism": rec.get("organism", {}).get("scientificName"),
            }
        )

    return hits


def resolve_gene_to_uniprot(gene_name: str, organism: str) -> str:
    """Resolve a gene name to a UniProt accession ID.

    For human genes, only Swiss-Prot (reviewed) entries are accepted.
    For other organisms, falls back to TrEMBL if no Swiss-Prot entry is found.

    Special case: If gene_name looks like a UniProt accession (6+ alphanumeric chars
    starting with letter), verify it exists directly rather than treating as gene name.

    Args:
        gene_name: Name of the gene or UniProt accession
        organism: Organism name (e.g., "human", "mouse", "yeast")

    Returns:
        UniProt accession ID (Swiss-Prot preferred, TrEMBL as fallback)

    Raises:
        ValueError: If gene cannot be resolved to UniProt ID
    """
    # Check if gene_name looks like a UniProt accession
    # UniProt accessions: 6+ chars, start with letter, alphanumeric + optional underscore/dash
    import re
    if re.match(r'^[A-Z][A-Z0-9_-]{5,}$', gene_name.upper()):
        # Looks like an accession, try to validate it exists
        try:
            # Try to fetch the UniProt data to verify it exists
            # This will also work for TrEMBL entries
            fetch_uniprot_data(gene_name)
            return gene_name
        except ValueError:
            # If fetch fails, fall through to treat as gene name
            pass

    # First try Swiss-Prot entries
    hits = uniprot_by_gene_taxon(gene_name, organism, limit=1, reviewed_only=True)

    if hits:
        return hits[0]["accession"]

    # For human, Swiss-Prot is required - no fallback to TrEMBL
    if organism.lower() == "human":
        raise ValueError(
            f"Could not find Swiss-Prot UniProt ID for gene {gene_name} in {organism}. "
            "All human genes should have Swiss-Prot entries."
        )

    # For other organisms, try TrEMBL as fallback
    hits = uniprot_by_gene_taxon(gene_name, organism, limit=1, reviewed_only=False)

    if not hits:
        raise ValueError(
            f"Could not find any UniProt ID (Swiss-Prot or TrEMBL) for gene {gene_name} in {organism}"
        )

    # Return the first hit (could be Swiss-Prot or TrEMBL)
    return hits[0]["accession"]


def fetch_uniprot_data(uniprot_id: str) -> str:
    """Fetch UniProt data in text format.

    Args:
        uniprot_id: UniProt accession ID (e.g., "Q9BRQ4")

    Returns:
        UniProt entry in text format

    Raises:
        ValueError: If UniProt data cannot be fetched
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    response = requests.get(url)

    if response.status_code != 200:
        raise ValueError(
            f"Failed to fetch UniProt data for {uniprot_id}: HTTP {response.status_code}"
        )

    return response.text


def fetch_goa_data(uniprot_id: str) -> str:
    """Fetch Gene Ontology Annotation (GOA) data from QuickGO API.

    Args:
        uniprot_id: UniProt accession ID

    Returns:
        GOA data in TSV format

    Raises:
        ValueError: If GOA data cannot be fetched
    """
    # Use QuickGO search endpoint with JSON format
    url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"

    # Build URL with properly formatted parameters
    # QuickGO expects multiple includeFields parameters, not an array
    params = [
        ("geneProductId", uniprot_id),
        ("includeFields", "goName"),
        ("includeFields", "taxonName"),
        ("includeFields", "name"),
        ("limit", "100"),  # API limit
    ]

    # Request JSON format
    headers = {"Accept": "application/json"}

    try:
        response = requests.get(url, params=params, headers=headers, timeout=30)
        response.raise_for_status()

        # Parse JSON response
        data = response.json()
        results = data.get("results", [])

        if not results:
            # Return header only if no data
            return "GENE PRODUCT DB\tGENE PRODUCT ID\tSYMBOL\tQUALIFIER\tGO TERM\tGO NAME\tGO ASPECT\tECO ID\tGO EVIDENCE CODE\tREFERENCE\tWITH/FROM\tTAXON ID\tTAXON NAME\tASSIGNED BY\tGENE NAME\tDATE\n"

        # Sort results: IBA first, then by most recent date (descending), then by GO ID
        sorted_results = sorted(
            results,
            key=lambda x: (
                0 if x.get("goEvidence", "") == "IBA" else 1,  # IBA annotations first
                -(
                    int(x.get("date", "0")) if x.get("date", "").isdigit() else 0
                ),  # Negative for descending date
                x.get("goId", ""),  # Then by GO ID
            ),
        )

        # Convert JSON to TSV format
        tsv_lines = [
            "GENE PRODUCT DB\tGENE PRODUCT ID\tSYMBOL\tQUALIFIER\tGO TERM\tGO NAME\tGO ASPECT\tECO ID\tGO EVIDENCE CODE\tREFERENCE\tWITH/FROM\tTAXON ID\tTAXON NAME\tASSIGNED BY\tGENE NAME\tDATE"
        ]

        for result in sorted_results:
            # Extract database and ID from geneProductId
            gene_product_id = result.get("geneProductId", "")
            if ":" in gene_product_id:
                db, prod_id = gene_product_id.split(":", 1)
            else:
                db, prod_id = "", gene_product_id

            # Extract WITH/FROM information
            with_from_list = []
            with_from_data = result.get("withFrom")
            if with_from_data:
                for wf in with_from_data:
                    for xref in wf.get("connectedXrefs", []):
                        with_from_list.append(
                            f"{xref.get('db', '')}:{xref.get('id', '')}"
                        )
            with_from = "|".join(with_from_list) if with_from_list else ""

            # Build TSV row
            row = [
                db,  # GENE PRODUCT DB
                prod_id,  # GENE PRODUCT ID
                result.get("symbol", ""),  # SYMBOL
                result.get("qualifier", ""),  # QUALIFIER
                result.get("goId", ""),  # GO TERM
                result.get("goName", ""),  # GO NAME
                result.get("goAspect", ""),  # GO ASPECT
                result.get("evidenceCode", ""),  # ECO ID
                result.get("goEvidence", ""),  # GO EVIDENCE CODE
                result.get("reference", ""),  # REFERENCE
                with_from,  # WITH/FROM
                str(result.get("taxonId", "")),  # TAXON ID
                result.get("taxonName", ""),  # TAXON NAME
                result.get("assignedBy", ""),  # ASSIGNED BY
                result.get("name", ""),  # GENE NAME
                result.get("date", ""),  # DATE
            ]

            tsv_lines.append("\t".join(row))

        return "\n".join(tsv_lines) + "\n"

    except requests.exceptions.RequestException:
        # Return header only if request fails
        return "GENE PRODUCT DB\tGENE PRODUCT ID\tSYMBOL\tQUALIFIER\tGO TERM\tGO NAME\tGO ASPECT\tECO ID\tGO EVIDENCE CODE\tREFERENCE\tWITH/FROM\tTAXON ID\tTAXON NAME\tASSIGNED BY\tGENE NAME\tDATE\n"
