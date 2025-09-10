
# AI Gene Review

AI-assisted tool for reviewing and curating gene annotations. This project provides a structured workflow for validating existing Gene Ontology (GO) annotations using AI-driven analysis combined with literature research and bioinformatics evidence.

## Overview

The AI Gene Review tool helps researchers and curators:
- **Review existing GO annotations** using strict, defined criteria
- **Synthesize high-quality annotations** from multiple evidence sources
- **Fetch and organize** gene data from UniProt and GOA databases
- **Validate annotation files** against LinkML schemas
- **Manage references** and supporting literature

## Quick Start

### Installation

1. Install [uv](https://docs.astral.sh/uv/) for dependency management
2. Clone the repository and install dependencies:
   ```bash
   git clone https://github.com/cmungall/ai-gene-review.git
   cd ai-gene-review
   uv sync --group dev
   ```

### Basic Usage

**Fetch gene data:**
```bash
uv run ai-gene-review fetch-gene human TP53
```

**Validate a gene review file:**
```bash
uv run ai-gene-review validate genes/human/TP53/TP53-ai-review.yaml
```

**Fetch publications for a gene:**
```bash
uv run ai-gene-review fetch-gene-pmids genes/human/TP53/TP53-ai-review.yaml
```

**Generate statistics report:**
```bash
just stats                # Generate HTML report
just stats-open           # Generate and open in browser
```

## Workflow Overview

1. **Fetch Gene Data**: Download UniProt records and GO annotations
2. **Literature Research**: Gather supporting publications and evidence
3. **Create Review**: Structure annotations using the YAML schema
4. **Validate**: Check against LinkML schema and best practices
5. **Iterate**: Refine annotations based on validation results

## Key Features

- 🧬 **Multi-organism support**: Human, mouse, worm, and other model organisms
- 📚 **Literature integration**: Automatic PubMed citation fetching and caching
- ✅ **Schema validation**: LinkML-based validation for consistency
- 🛡️ **Anti-hallucination validation**: ID/label tuple checksums prevent AI fabrication of terms
- 🔄 **Batch processing**: Handle multiple genes efficiently
- 📊 **Structured reviews**: YAML-based gene annotation reviews
- 🔍 **Evidence tracking**: Detailed provenance and supporting text


## Resources

### Documentation & Visualization

- **Documentation Website**: [https://monarch-initiative.github.io/ai-gene-review](https://monarch-initiative.github.io/ai-gene-review)
- **Interactive Web App**: [https://ai4curation.github.io/ai-gene-review/app/index.html](https://ai4curation.github.io/ai-gene-review/app/index.html) - Browse and explore gene annotation reviews
- **Statistics Dashboard**: [https://ai4curation.github.io/ai-gene-review/docs/stats_report.html](https://ai4curation.github.io/ai-gene-review/docs/stats_report.html) - Summary Stats

## Gene Review Structure

Each gene review follows a structured YAML format containing:

- **Gene metadata**: UniProt ID, gene symbol, taxon information
- **Description**: Comprehensive summary of gene function
- **References**: Literature and bioinformatics sources
- **Existing annotations**: Review of current GO annotations with actions (ACCEPT, MODIFY, REMOVE, etc.)
- **Core functions**: Curated essential gene functions

Example structure:
```yaml
id: Q9BRQ4
gene_symbol: CFAP300
taxon:
  id: NCBITaxon:9606
  label: Homo sapiens
description: >-
  CFAP300 is a cilium- and flagellum-specific protein...
existing_annotations:
  - term:
      id: GO:0005515
      label: protein binding
    action: MODIFY
    reason: "While evidence is strong, 'protein binding' is uninformative..."
```

## Example Data

The repository includes example gene reviews for:
- **Human**: BRCA1, CFAP300, RBFOX3, TP53
- **Mouse**: Various examples
- **Worm**: lrx-1

Browse the `genes/` directory to see complete examples.

## Case Studies

### PedH (Pseudomonas putida KT2440) - Lanthanide-Dependent Alcohol Dehydrogenase

The review of **pedH** revealed several important curation insights:

#### Key Discoveries

1. **Lanthanide vs Calcium Dependency**: PedH was incorrectly annotated with "calcium ion binding" (GO:0005509) when it actually requires lanthanide ions (La³⁺, Ce³⁺, Pr³⁺, Nd³⁺, Sm³⁺) for activity. This highlights the importance of reviewing automated annotations based on sequence similarity.

2. **Cellular Localization Precision**: Bioinformatics analysis confirmed PedH is a **soluble periplasmic enzyme**, not membrane-associated:
   - Signal peptide (aa 1-25) directs export, then is cleaved
   - No transmembrane regions in mature protein
   - Functions throughout periplasmic space, not just at membrane boundaries
   - Led to choosing GO:0042597 (periplasmic space) over GO:0030288 (outer membrane-bounded periplasmic space)

3. **Dual Functional Roles**: PedH serves both as:
   - **Metabolic enzyme**: Oxidizes alcohols in 2-phenylethanol degradation pathway
   - **Regulatory sensor**: Part of lanthanide-sensing system controlling gene expression via PedS2/PedR2 two-component system

4. **Missing GO Terms Identified**: The review revealed gaps in GO:
   - No term for "lanthanide ion binding" (distinct from transition metal binding)
   - No term for "lanthanide-dependent alcohol dehydrogenase activity"

#### Lessons for Curation

- **Verify metal cofactors carefully** - Don't assume calcium when other metals are possible
- **Consider protein mobility** - Soluble vs membrane-associated matters for localization terms
- **Look for regulatory functions** - Enzymes may have sensory/regulatory roles beyond catalysis
- **Use bioinformatics to validate** - Signal peptide and TM predictions can clarify localization

## Anti-Hallucination Validation

The AI Gene Review system implements a robust **anti-hallucination validation mechanism** using **ID/label tuple checksums** to prevent AI systems from fabricating or misusing ontological terms.

### How It Works

Every ontology term in the system requires both an `id` (semantic identifier) and `label` (human-readable name):

```yaml
term:
  id: GO:0005515      # Ontology identifier
  label: protein binding  # Canonical label
```

### Validation Process

The `TermValidator` performs multi-layer validation:

1. **Format Validation**: Ensures IDs follow proper CURIE patterns (`PREFIX:NUMBER`)
2. **Existence Validation**: Verifies terms exist in authoritative ontologies via OAK/OLS APIs
3. **Label Matching**: Cross-references provided labels against canonical ontology labels
4. **Branch Validation**: Ensures GO terms are in correct ontological branches (MF/BP/CC)
5. **Obsolescence Checking**: Flags outdated terms

### Why This Prevents Hallucination

✅ **Dual Verification**: Both ID and label must be correct and consistent
✅ **External Truth Source**: Validates against authoritative ontologies (GO, HP, MONDO, etc.)
✅ **Real-time Checking**: Uses live API calls to catch fabricated terms
✅ **Semantic Consistency**: Ensures terms make sense in their context

### Examples

```yaml
# ❌ This would be caught as invalid
term:
  id: GO:0005515
  label: "DNA binding"  # Wrong label for GO:0005515

# ✅ This passes validation
term:
  id: GO:0005515  
  label: "protein binding"  # Correct canonical label

# ❌ This would be flagged as fabricated
term:
  id: GO:9999999
  label: "made up function"  # Non-existent term
```

### Supported Ontologies

The validator supports 10+ major ontologies:
- **GO**: Gene Ontology (molecular functions, biological processes, cellular components)
- **HP**: Human Phenotype Ontology
- **MONDO**: Mondo Disease Ontology  
- **CL**: Cell Ontology
- **UBERON**: Uberon Anatomy Ontology
- **CHEBI**: Chemical Entities of Biological Interest
- **PR**: Protein Ontology
- **SO**: Sequence Ontology
- **PATO**: Phenotype And Trait Ontology
- **NCBITaxon**: NCBI Taxonomy

This validation system represents a **novel approach to preventing ontological hallucination** in AI curation workflows and could serve as a model for other AI applications working with structured biological knowledge.

## Repository Structure

* **[genes/](genes/)** - Gene review data organized by organism
  * `human/`, `mouse/`, `worm/` - Species-specific gene directories
  * Each gene folder contains: YAML review, UniProt data, GO annotations, notes
* **[docs/](docs/)** - MkDocs-managed documentation
* **[src/ai_gene_review/](src/ai_gene_review/)** - Core Python package
  * `cli.py` - Command-line interface
  * `schema/` - LinkML schema definitions
  * `etl/` - Data extraction and loading modules
* **[tests/](tests/)** - Python tests and example data
* **[publications/](publications/)** - Cached PubMed articles

## Developer Tools

This project uses [just](https://github.com/casey/just/) command runner for development tasks.

**Available commands:**
```bash
just --list           # Show all available commands
just test             # Run tests, type checking, and linting
just format           # Run code formatting checks
just install          # Install project dependencies
```

**CLI Commands:**
```bash
uv run ai-gene-review --help                    # Show CLI help
uv run ai-gene-review fetch-gene human BRCA1    # Fetch gene data
uv run ai-gene-review validate <yaml-file>      # Validate review file
uv run ai-gene-review batch-fetch <input-file>  # Process multiple genes
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed contribution guidelines including:
- Code of conduct and best practices
- Understanding LinkML schemas
- Pull request workflow
- Development setup

## Credits

This project uses the template [monarch-project-copier](https://github.com/monarch-initiative/monarch-project-copier)
