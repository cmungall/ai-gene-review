"""ETL for fetching and caching Reactome pathway information."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests
import yaml


@dataclass
class ReactomePathway:
    """Represents a Reactome pathway with metadata.

    Example:
        >>> pathway = ReactomePathway(
        ...     stable_id="R-HSA-6785807",
        ...     display_name="Interleukin-4 and Interleukin-13 signaling",
        ...     species="Homo sapiens"
        ... )
        >>> len(pathway.to_markdown()) > 0
        True
    """

    stable_id: str
    display_name: str
    species: Optional[str] = None
    summary: Optional[str] = None

    def to_frontmatter_dict(self) -> dict:
        """Convert to dictionary for YAML frontmatter."""
        data = {
            "stable_id": self.stable_id,
            "display_name": self.display_name,
        }
        if self.species:
            data["species"] = self.species
        if self.summary:
            data["summary"] = self.summary
        return data

    def to_markdown(self) -> str:
        """Generate markdown content with frontmatter."""
        # Create frontmatter
        frontmatter = yaml.dump(
            self.to_frontmatter_dict(),
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )

        # Create markdown body
        body_parts = [f"# {self.display_name}\n"]

        body_parts.append(
            f"**Reactome ID:** [{self.stable_id}](https://reactome.org/content/detail/{self.stable_id})\n"
        )

        if self.species:
            body_parts.append(f"**Species:** {self.species}\n")

        if self.summary:
            body_parts.append(f"\n## Summary\n\n{self.summary}\n")

        return f"---\n{frontmatter}---\n\n{''.join(body_parts)}"


def extract_reactome_id(reactome_str: str) -> str:
    """Extract Reactome ID from various formats.

    Examples:
        >>> extract_reactome_id('R-HSA-6785807')
        'R-HSA-6785807'
        >>> extract_reactome_id('Reactome:R-HSA-6785807')
        'R-HSA-6785807'
        >>> extract_reactome_id('reactome: R-HSA-6785807')
        'R-HSA-6785807'
        >>> extract_reactome_id(' Reactome:R-HSA-6785807 ')
        'R-HSA-6785807'
    """
    # Remove Reactome prefix and clean up
    reactome_id = reactome_str.strip()
    if ":" in reactome_id:
        reactome_id = reactome_id.split(":", 1)[1].strip()
    return reactome_id


def get_cached_title(
    reactome_id: str, cache_dir: Path = Path("reactome")
) -> Optional[str]:
    """Get just the title from cached Reactome pathway if available.

    Args:
        reactome_id: Reactome stable ID (e.g., R-HSA-6785807)
        cache_dir: Directory containing cached Reactome pathways

    Returns:
        Title string if cached, None otherwise
    """
    cache_file = cache_dir / f"{reactome_id}.md"
    if cache_file.exists():
        # Parse the markdown file to get title
        content = cache_file.read_text()

        # Extract title from first # header
        for line in content.split("\n"):
            if line.startswith("# "):
                return line[2:].strip()

    return None


def fetch_reactome_data(
    reactome_id: str, use_cache: bool = True, cache_dir: Path = Path("reactome")
) -> Optional[ReactomePathway]:
    """Fetch pathway data from Reactome API, using cache if available.

    Args:
        reactome_id: Reactome stable ID (e.g., R-HSA-6785807)
        use_cache: Whether to use cached data if available
        cache_dir: Directory for caching pathways

    Returns:
        ReactomePathway object, or None if not found

    Example:
        >>> # This would require network access
        >>> # pathway = fetch_reactome_data("R-HSA-6785807")
        >>> # if pathway:
        >>> #     assert pathway.stable_id == "R-HSA-6785807"
    """
    # Clean the ID
    reactome_id = extract_reactome_id(reactome_id)

    # Check cache first
    if use_cache:
        cached_title = get_cached_title(reactome_id, cache_dir)
        if cached_title:
            # Return minimal pathway object from cache
            return ReactomePathway(stable_id=reactome_id, display_name=cached_title)

    try:
        # Fetch from Reactome API
        url = f"https://reactome.org/ContentService/data/query/{reactome_id}"
        response = requests.get(url, timeout=10)

        if response.status_code != 200:
            return None

        data = response.json()

        # Extract pathway information
        display_name = data.get("displayName", "Unknown pathway")
        species_name = None
        summary_text = None

        # Get species if available
        if "species" in data and isinstance(data["species"], list) and data["species"]:
            species = data["species"][0]
            if isinstance(species, dict):
                species_name = species.get("displayName")

        # Get summary if available
        if (
            "summation" in data
            and isinstance(data["summation"], list)
            and data["summation"]
        ):
            summation = data["summation"][0]
            if isinstance(summation, dict):
                summary_text = summation.get("text")

        pathway = ReactomePathway(
            stable_id=reactome_id,
            display_name=display_name,
            species=species_name,
            summary=summary_text,
        )

        # Cache the pathway if we fetched it
        if use_cache:
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_file = cache_dir / f"{reactome_id}.md"
            try:
                cache_file.write_text(pathway.to_markdown())
            except Exception:
                # Silently fail on cache write errors
                pass

        return pathway

    except Exception as e:
        print(f"Error fetching Reactome {reactome_id}: {e}")
        return None


def cache_reactome_pathway(
    reactome_id: str, output_dir: Path = Path("reactome"), force: bool = False
) -> bool:
    """Cache a single Reactome pathway to the filesystem.

    Args:
        reactome_id: Reactome stable ID (with or without Reactome: prefix)
        output_dir: Directory to save cached pathway
        force: Force re-download even if already cached

    Returns:
        True if successfully cached, False otherwise

    Example:
        >>> # This would require network access
        >>> # success = cache_reactome_pathway("Reactome:R-HSA-6785807")
        >>> # if success:
        >>> #     assert Path("reactome/R-HSA-6785807.md").exists()
    """
    # Clean ID
    reactome_id = extract_reactome_id(reactome_id)

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if already cached
    output_file = output_dir / f"{reactome_id}.md"
    if output_file.exists() and not force:
        print(f"Reactome {reactome_id} already cached (use force=True to re-download)")
        return True

    # Fetch pathway data
    print(f"Fetching Reactome {reactome_id}...")
    pathway = fetch_reactome_data(reactome_id, use_cache=False)

    if not pathway:
        print(f"Failed to fetch Reactome {reactome_id}")
        return False

    # Write to file
    output_file.write_text(pathway.to_markdown())
    print(f"Cached Reactome {reactome_id}: {pathway.display_name}")

    return True
