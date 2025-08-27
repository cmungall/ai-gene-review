"""ETL for fetching and caching PubMed/PMC publications."""

import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml
from Bio import Entrez  # type: ignore[import-untyped]

# Set email for NCBI (required for Entrez)
Entrez.email = "ai-gene-review@example.com"  # type: ignore


@dataclass
class Publication:
    """Represents a publication with metadata and content.
    
    Example:
        >>> pub = Publication(
        ...     pmid="12345",
        ...     title="Test Article",
        ...     authors=["Smith J", "Doe J"],
        ...     journal="Test Journal",
        ...     year="2024",
        ...     abstract="This is a test abstract."
        ... )
        >>> frontmatter = pub.to_frontmatter_dict()
        >>> frontmatter['pmid']
        '12345'
        >>> len(pub.to_markdown()) > 0
        True
    """
    pmid: str
    title: str
    authors: List[str]
    journal: str
    year: str
    abstract: str
    full_text: Optional[str] = None
    pmcid: Optional[str] = None
    doi: Optional[str] = None
    keywords: Optional[List[str]] = None
    
    def to_frontmatter_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML frontmatter."""
        data = {
            'pmid': self.pmid,
            'title': self.title,
            'authors': self.authors,
            'journal': self.journal,
            'year': self.year,
        }
        if self.pmcid:
            data['pmcid'] = self.pmcid
        if self.doi:
            data['doi'] = self.doi
        if self.keywords:
            data['keywords'] = self.keywords
        return data
    
    def to_markdown(self) -> str:
        """Generate markdown content with frontmatter."""
        # Create frontmatter
        frontmatter = yaml.dump(
            self.to_frontmatter_dict(), 
            default_flow_style=False, 
            allow_unicode=True, 
            sort_keys=False
        )
        
        # Create markdown body
        body_parts = [f"# {self.title}\n"]
        
        if self.authors:
            body_parts.append(f"**Authors:** {', '.join(self.authors)}\n")
        
        body_parts.append(f"**Journal:** {self.journal} ({self.year})\n")
        
        if self.doi:
            body_parts.append(f"**DOI:** [{self.doi}](https://doi.org/{self.doi})\n")
        
        if self.pmcid:
            body_parts.append(f"**PMC:** [{self.pmcid}](https://www.ncbi.nlm.nih.gov/pmc/articles/{self.pmcid}/)\n")
        
        body_parts.append(f"\n## Abstract\n\n{self.abstract}\n")
        
        if self.full_text:
            body_parts.append(f"\n## Full Text\n\n{self.full_text}\n")
        
        return f"---\n{frontmatter}---\n\n{''.join(body_parts)}"


def extract_pmid(pmid_str: str) -> str:
    """Extract PMID from various formats.
    
    Examples:
        >>> extract_pmid('12345678')
        '12345678'
        >>> extract_pmid('PMID:12345678')
        '12345678'
        >>> extract_pmid('pmid: 12345678')
        '12345678'
        >>> extract_pmid('PMID12345678')
        '12345678'
        >>> extract_pmid(' PMID: 12345678 ')
        '12345678'
    """
    # Remove PMID prefix and clean up
    pmid = re.sub(r'^(PMID|pmid)[:\s]*', '', pmid_str.strip())
    return pmid.strip()


def get_cached_title(pmid: str, cache_dir: Path = Path("publications")) -> Optional[str]:
    """Get just the title from cached publication if available.
    
    Args:
        pmid: PubMed ID (without PMID prefix)
        cache_dir: Directory containing cached publications
        
    Returns:
        Title string if cached, None otherwise
    """
    cache_file = cache_dir / f"PMID_{pmid}.md"
    if cache_file.exists():
        # Parse the markdown file to get title
        content = cache_file.read_text()
        
        # Extract title from first # header
        for line in content.split('\n'):
            if line.startswith('# '):
                return line[2:].strip()
    
    return None


def get_cached_publication(pmid: str, cache_dir: Path = Path("publications")) -> Optional[Publication]:
    """Get publication from cache if available.
    
    Args:
        pmid: PubMed ID (without PMID prefix)
        cache_dir: Directory containing cached publications
        
    Returns:
        Publication object if cached, None otherwise
    """
    # For now, return None to always fetch fresh for full data
    # Could implement full markdown parsing if needed
    return None


def fetch_pubmed_data(pmid: str, use_cache: bool = True, cache_dir: Path = Path("publications")) -> Optional[Publication]:
    """Fetch publication data from PubMed, using cache if available.
    
    Args:
        pmid: PubMed ID (without PMID prefix)
        use_cache: Whether to use cached data if available
        cache_dir: Directory for caching publications
        
    Returns:
        Publication object with abstract, or None if not found
        
    Example:
        >>> # This would require network access
        >>> # pub = fetch_pubmed_data("29727692")
        >>> # if pub:
        >>> #     assert pub.pmid == "29727692"
    """
    # Check cache first
    if use_cache:
        cached = get_cached_publication(pmid, cache_dir)
        if cached:
            return cached
    
    try:
        # Fetch summary from PubMed
        handle = Entrez.esummary(db="pubmed", id=pmid, retmode="xml")
        summary_records = Entrez.read(handle)
        handle.close()
        
        if not summary_records or 'error' in summary_records[0]:
            return None
        
        record = summary_records[0]
        
        # Extract basic metadata (convert Bio.Entrez objects to strings)
        title = str(record.get('Title', 'No title'))
        authors = [str(author) for author in record.get('AuthorList', [])]
        journal = str(record.get('Source', 'Unknown journal'))
        year = str(record.get('PubDate', 'Unknown'))[:4] if 'PubDate' in record else 'Unknown'
        doi = str(record.get('DOI', '')) if record.get('DOI') else None
        
        # Check for PMC ID
        pmcid = None
        
        # Method 1: Check ArticleIds
        if 'ArticleIds' in record:
            article_ids = record['ArticleIds']
            if isinstance(article_ids, dict):
                for id_type, id_value in article_ids.items():
                    if 'pmc' in str(id_type).lower():
                        pmcid = str(id_value)
                        break
        
        # Method 2: Use Entrez elink to find PMC ID
        if not pmcid:
            try:
                handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
                link_records = Entrez.read(handle)
                handle.close()
                
                if link_records and link_records[0].get('LinkSetDb'):
                    for linkset in link_records[0]['LinkSetDb']:
                        if linkset.get('DbTo') == 'pmc' and linkset.get('Link'):
                            pmcid = 'PMC' + str(linkset['Link'][0]['Id'])
                            break
            except Exception:
                pass
        
        # Fetch abstract from PubMed
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract_text = handle.read()
        handle.close()
        
        # Clean up abstract
        abstract = abstract_text.strip()
        if not abstract or abstract == "":
            abstract = "No abstract available."
        
        publication = Publication(
            pmid=pmid,
            title=title,
            authors=authors,
            journal=journal,
            year=year,
            abstract=abstract,
            pmcid=pmcid,
            doi=doi
        )
        
        # Try to fetch full text from PMC if available
        if pmcid:
            full_text = fetch_pmc_fulltext(pmcid)
            if full_text:
                publication.full_text = full_text
        
        # Cache the publication if we fetched it
        if use_cache:
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_file = cache_dir / f"PMID_{pmid}.md"
            try:
                cache_file.write_text(publication.to_markdown())
            except Exception:
                # Silently fail on cache write errors
                pass
        
        return publication
        
    except Exception as e:
        print(f"Error fetching PMID {pmid}: {e}")
        return None


def fetch_pmc_fulltext(pmcid: str) -> Optional[str]:
    """Fetch full text from PMC.
    
    Args:
        pmcid: PMC ID (with or without PMC prefix)
        
    Returns:
        Full text as string, or None if not available
    """
    try:
        # Remove PMC prefix if present
        pmcid = pmcid.replace('PMC', '')
        
        # Fetch full text XML from PMC
        handle = Entrez.efetch(db="pmc", id=pmcid, rettype="full", retmode="xml")
        xml_data = handle.read()
        handle.close()
        
        # Parse XML to extract text
        root = ET.fromstring(xml_data)
        
        # Extract body text
        body_texts = []
        for elem in root.iter():
            if elem.tag in ['p', 'sec']:
                text = ''.join(elem.itertext()).strip()
                if text:
                    body_texts.append(text)
        
        if body_texts:
            return '\n\n'.join(body_texts)
        
    except Exception as e:
        print(f"Could not fetch full text for PMC{pmcid}: {e}")
    
    return None


def cache_publication(
    pmid: str, 
    output_dir: Path = Path("publications"),
    force: bool = False
) -> bool:
    """Cache a single publication to the filesystem.
    
    Args:
        pmid: PubMed ID (with or without PMID prefix)
        output_dir: Directory to save cached publication
        force: Force re-download even if already cached
        
    Returns:
        True if successfully cached, False otherwise
        
    Example:
        >>> # This would require network access
        >>> # success = cache_publication("PMID:29727692", Path("test_publications"))
        >>> # if success:
        >>> #     assert (Path("test_publications") / "PMID_29727692.md").exists()
    """
    # Clean PMID
    pmid = extract_pmid(pmid)
    
    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if already cached
    output_file = output_dir / f"PMID_{pmid}.md"
    if output_file.exists() and not force:
        print(f"PMID {pmid} already cached (use force=True to re-download)")
        return True
    
    # Fetch publication data
    print(f"Fetching PMID {pmid}...")
    publication = fetch_pubmed_data(pmid)
    
    if not publication:
        print(f"Failed to fetch PMID {pmid}")
        return False
    
    # Write to file
    output_file.write_text(publication.to_markdown())
    
    if publication.full_text:
        print(f"Cached PMID {pmid} with full text from PMC")
    else:
        print(f"Cached PMID {pmid} with abstract only")
    
    return True


def cache_publications(
    pmids: List[str],
    output_dir: Path = Path("publications"),
    force: bool = False,
    delay: float = 0.5
) -> int:
    """Cache multiple publications.
    
    Args:
        pmids: List of PubMed IDs
        output_dir: Directory to save cached publications
        force: Force re-download even if already cached
        delay: Delay between requests in seconds (be polite to NCBI)
        
    Returns:
        Number of successfully cached publications
        
    Example:
        >>> # This would require network access
        >>> # count = cache_publications(["29727692", "29727693"])
        >>> # assert count >= 0
    """
    success_count = 0
    
    for i, pmid in enumerate(pmids):
        if i > 0:
            time.sleep(delay)  # Be polite to NCBI servers
        
        if cache_publication(pmid, output_dir, force):
            success_count += 1
    
    print(f"Cached {success_count}/{len(pmids)} publications")
    return success_count


def extract_pmids_from_yaml(yaml_file: Path) -> List[str]:
    """Extract PMIDs from a gene review YAML file.
    
    Args:
        yaml_file: Path to gene review YAML file
        
    Returns:
        List of PMIDs found in the file
        
    Example:
        >>> import tempfile
        >>> data = {
        ...     "references": [
        ...         {"id": "PMID:12345"},
        ...         {"id": "PMID:67890"}
        ...     ],
        ...     "existing_annotations": [
        ...         {"original_reference_id": "PMID:11111"}
        ...     ]
        ... }
        >>> with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        ...     yaml.dump(data, f)
        ...     temp_file = Path(f.name)
        >>> pmids = extract_pmids_from_yaml(temp_file)
        >>> sorted(pmids)
        ['11111', '12345', '67890']
        >>> temp_file.unlink()
    """
    pmids = set()
    
    with open(yaml_file) as f:
        data = yaml.safe_load(f)
    
    if not data:
        return []
    
    # Extract from references
    if 'references' in data and data['references']:
        for ref in data['references']:
            if isinstance(ref, dict) and 'id' in ref:
                ref_id = ref['id']
                if ref_id and ref_id.startswith('PMID'):
                    pmids.add(extract_pmid(ref_id))
    
    # Extract from existing_annotations
    if 'existing_annotations' in data and data['existing_annotations']:
        for annotation in data['existing_annotations']:
            if isinstance(annotation, dict) and 'original_reference_id' in annotation:
                ref_id = annotation['original_reference_id']
                if ref_id and ref_id.startswith('PMID'):
                    pmids.add(extract_pmid(ref_id))
    
    return list(pmids)