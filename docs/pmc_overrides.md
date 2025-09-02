# PMC ID Override System

## Overview

The PMC ID override system addresses incorrect linkages in NCBI's database where PubMed IDs (PMIDs) are incorrectly associated with PMC (PubMed Central) IDs. This can happen when:

- NCBI's elink service incorrectly cross-references papers
- Multiple papers about the same topic get confused in the linkage database
- Database entries are corrupted or misaligned

## How It Works

1. **Override Table**: A TSV file (`src/ai_gene_review/etl/pmc_overrides.tsv`) contains manual corrections
2. **Format**: Simple two-column format: `PMID<tab>PMCID`
3. **Loading**: The override table is loaded and cached on first use
4. **Priority**: Overrides take precedence over NCBI's elink results

## File Format

The override file is a tab-separated values (TSV) file with the following format:

```tsv
# Comments start with #
# PMID	PMCID	Optional comment
2001740		# No PMC version exists (NCBI incorrectly links to PMC11824087)
12345	PMC67890	# Correct PMC ID for this paper
```

- **PMID**: The PubMed ID (without "PMID:" prefix)
- **PMCID**: The correct PMC ID, or empty if no PMC version exists
- **Comments**: Lines starting with # are ignored
- **Blank PMCID**: Indicates the paper has no PMC version

## Usage

### Adding an Override

1. Edit `src/ai_gene_review/etl/pmc_overrides.tsv`
2. Add a new line with the PMID and correct PMCID (or leave blank)
3. Add a comment explaining the issue
4. Report the issue to NCBI at: https://www.ncbi.nlm.nih.gov/home/about/contact/

### Example Cases

#### Case 1: No PMC Version
```tsv
2001740		# 1991 FEBS Letters paper, no PMC version available
```

#### Case 2: Wrong PMC Linked
```tsv
12345	PMC67890	# NCBI links to PMC99999 which is a different paper
```

#### Case 3: PMC Exists but Not Linked
```tsv
54321	PMC11111	# NCBI doesn't show PMC link but it exists
```

## Implementation Details

The override system is implemented in `src/ai_gene_review/etl/publication.py`:

- `load_pmc_overrides()`: Loads and caches the override table
- `fetch_pubmed_data()`: Checks overrides before querying NCBI's elink

## Testing

Run tests with:
```bash
uv run pytest tests/test_pmc_overrides.py
```

## Reporting Issues to NCBI

When you find an incorrect PMC linkage:

1. Add it to the override table (immediate fix)
2. Report to NCBI:
   - URL: https://www.ncbi.nlm.nih.gov/home/about/contact/
   - Include both PMIDs and PMC IDs involved
   - Provide paper titles and years to help identify the issue
   
## Example Report to NCBI

```
Subject: Incorrect PMC ID linkage for PMID 2001740

The PubMed entry for PMID 2001740 (Yamaguchi & Sakurai, FEBS Lett 1991) 
is incorrectly linked to PMC11824087.

PMC11824087 actually corresponds to PMID 25230901 (Yamaguchi, J Cancer 
Res Clin Oncol 2014), which is a different paper about the same protein.

The 1991 paper appears to have no PMC version available, as expected for 
older subscription journal articles.

Please correct this linkage in the NCBI database.
```