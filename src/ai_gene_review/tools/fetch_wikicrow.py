#!/usr/bin/env python3
"""
DEPRECATED: This script is deprecated and non-functional.

Wikicrow requires JavaScript rendering which cannot be handled by simple HTTP requests.
Proper implementation would require browser automation tools like Selenium or Playwright.
"""

import argparse
import sys
from pathlib import Path
import subprocess
import json
from typing import Optional


def fetch_wikicrow_with_browser(gene_symbol: str) -> Optional[str]:
    """
    Fetch gene information from Wikicrow using browser automation.
    
    This uses a headless browser approach to handle JavaScript-rendered content.
    
    Args:
        gene_symbol: Human gene symbol (e.g., TP53, BRCA1)
        
    Returns:
        Markdown formatted content from Wikicrow or None if fetch fails
    """
    url = f"https://wikicrow.ai/{gene_symbol}"
    
    # Try using curl with JavaScript evaluation hints
    # Note: This is a fallback - proper implementation would use Selenium or Playwright
    try:
        # First try with curl to see if we can get anything
        result = subprocess.run(
            ['curl', '-s', '-L', 
             '-H', 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36',
             '-H', 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
             url],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0:
            content = result.stdout
            
            # Check if we got JavaScript app shell (common with React/Next.js)
            if '__NEXT_DATA__' in content or 'window.__INITIAL_STATE__' in content:
                # Try to extract JSON data
                import re
                
                # Look for Next.js data
                next_data_match = re.search(r'<script id="__NEXT_DATA__"[^>]*>(.*?)</script>', content, re.DOTALL)
                if next_data_match:
                    try:
                        json_data = json.loads(next_data_match.group(1))
                        # Extract relevant content from JSON structure
                        # This would need to be adapted based on actual structure
                        return f"# Wikicrow: {gene_symbol}\n\nSource: {url}\n\nNote: Content requires JavaScript rendering. Please visit the URL directly in a browser to view the full content."
                    except json.JSONDecodeError:
                        pass
            
            # If we can't parse JavaScript data, return a placeholder
            return f"""# Wikicrow: {gene_symbol}

Source: {url}

Note: Wikicrow requires JavaScript rendering for full content display.
To view the complete gene information, please visit the URL directly in your browser.

The Wikicrow page likely contains:
- Gene overview and function
- Associated diseases and phenotypes
- Protein structure and domains
- Expression patterns
- Interactions and pathways
- References to scientific literature

For automated fetching, consider using browser automation tools like Selenium or Playwright."""
            
    except subprocess.TimeoutExpired:
        print(f"Timeout while fetching {url}", file=sys.stderr)
    except Exception as e:
        print(f"Error fetching from Wikicrow: {e}", file=sys.stderr)
    
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Fetch gene information from Wikicrow',
        epilog='Note: Wikicrow uses JavaScript rendering. Full content may require browser access.'
    )
    parser.add_argument('gene', help='Human gene symbol (e.g., TP53)')
    parser.add_argument('--output-dir', default='genes/human', 
                       help='Output directory for the gene folder')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing file')
    
    args = parser.parse_args()
    
    # Ensure it's a human gene
    gene_symbol = args.gene.upper()
    
    # Create output path
    output_dir = Path(args.output_dir) / gene_symbol
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{gene_symbol}-wikicrow.md"
    
    # Check if file exists
    if output_file.exists() and not args.force:
        print(f"File {output_file} already exists. Use --force to overwrite.")
        return 1
    
    print(f"Attempting to fetch Wikicrow content for {gene_symbol}...")
    print("Note: Full content requires JavaScript rendering.")
    
    content = fetch_wikicrow_with_browser(gene_symbol)
    
    if content:
        output_file.write_text(content)
        print(f"Saved placeholder/partial content to {output_file}")
        print(f"Visit https://wikicrow.ai/{gene_symbol} in your browser for full content.")
        return 0
    else:
        print(f"Failed to fetch content for {gene_symbol}")
        return 1


if __name__ == '__main__':
    sys.exit(main())