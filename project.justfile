## Add your own just recipes here. This is imported by the main justfile.

all: validate-all test

# Fetch gene data from UniProt and GOA
fetch-gene organism gene *args="":
    uv run ai-gene-review fetch-gene {{organism}} {{gene}} --output-dir . {{args}}

# Conduct deep research on a gene using OpenAI Deep Research API
deep-research organism gene *args="":
    uv run python src/ai_gene_review/tools/deep_research.py {{gene}} --organism {{organism}} --output-dir genes/{{organism}}/{{gene}} {{args}}

# Fetch a specific PMID
fetch-pmid pmid output_dir="publications":
    uv run ai-gene-review fetch-pmid {{pmid}} --output-dir {{output_dir}}

# Fetch all PMIDs referenced in a gene's review file
fetch-gene-pmids organism gene output_dir="publications":
    uv run ai-gene-review fetch-gene-pmids {{organism}} {{gene}} --output-dir {{output_dir}}

# Fetch PMIDs from a file
fetch-pmids-from-file file output_dir="publications":
    uv run ai-gene-review fetch-pmids-from-file {{file}} --output-dir {{output_dir}}

# Validate gene review YAML files against LinkML schema (for multiple files or patterns)
validate-files files:
    uv run ai-gene-review validate {{files}}

# Validate a specific gene's review file (short alias)
validate organism gene:
    uv run ai-gene-review validate --verbose genes/{{organism}}/{{gene}}/{{gene}}-ai-review.yaml

# Validate a specific gene's review file (long name for clarity)
validate-gene organism gene:
    uv run ai-gene-review validate genes/{{organism}}/{{gene}}/{{gene}}-ai-review.yaml

# Validate with verbose output
validate-gene-verbose organism gene:
    uv run ai-gene-review validate genes/{{organism}}/{{gene}}/{{gene}}-ai-review.yaml --verbose

# Validate all genes for an organism (shows detailed errors by default)
validate-organism organism:
    uv run ai-gene-review validate --verbose "genes/{{organism}}/*/*-ai-review.yaml"

# Validate all gene review files (shows detailed errors by default)
validate-all:
    uv run ai-gene-review validate --verbose "genes/*/*/*-ai-review.yaml"

# Validate all gene review files (summary only, no details)
validate-all-summary:
    uv run ai-gene-review validate "genes/*/*/*-ai-review.yaml"

# Validate all gene review files (strict mode - warnings become errors)
validate-all-strict:
    uv run ai-gene-review validate --verbose --strict "genes/*/*/*-ai-review.yaml"

# Seed missing GOA annotations in a gene's review file
seed-goa organism gene:
    uv run ai-gene-review seed-goa genes/{{organism}}/{{gene}}/{{gene}}-ai-review.yaml

# Validate GOA annotations for a specific gene
validate-goa organism gene:
    uv run ai-gene-review validate-goa genes/{{organism}}/{{gene}}/{{gene}}-ai-review.yaml

# Export existing_annotations to CSV format
export-annotations output_file="exports/exported_annotations.csv":
    @mkdir -p exports
    uv run python -c "from ai_gene_review.export import TabularExporter; from pathlib import Path; exporter = TabularExporter(); files = list(Path('genes').glob('**/*-ai-review.yaml')); print(f'Found {len(files)} files'); exporter.export_to_csv(files, '{{output_file}}'); print(f'Exported to {{output_file}}')"

# Export existing_annotations to TSV format  
export-annotations-tsv output_file="exports/exported_annotations.tsv":
    @mkdir -p exports
    uv run python -c "from ai_gene_review.export import TabularExporter; from pathlib import Path; exporter = TabularExporter(); files = list(Path('genes').glob('**/*-ai-review.yaml')); print(f'Found {len(files)} files'); exporter.export_to_tsv(files, '{{output_file}}'); print(f'Exported to {{output_file}}')"

# Export annotations for a specific organism
export-organism-annotations organism output_file="exports/exported_annotations.csv":
    @mkdir -p exports
    uv run python -c "from ai_gene_review.export import TabularExporter; from pathlib import Path; exporter = TabularExporter(); files = list(Path('genes/{{organism}}').glob('**/*-ai-review.yaml')); print(f'Found {len(files)} files for {{organism}}'); exporter.export_to_csv(files, '{{output_file}}'); print(f'Exported to {{output_file}}')"


# Batch fetch genes from a file
batch-fetch input_file output_dir=".":
    uv run ai-gene-review batch-fetch {{input_file}} --output-dir {{output_dir}}

# Example: Fetch common human genes
fetch-examples:
    uv run ai-gene-review fetch-gene human TP53
    uv run ai-gene-review fetch-gene human BRCA1
    uv run ai-gene-review fetch-gene human EGFR
    uv run ai-gene-review fetch-gene human CFAP300

# Build static site from gene markdown files
build-site:
    uv run build-site

# Build and serve the site
build-and-serve: build-site
    uv run mkdocs serve

# Generate project artifacts from LinkML schema
gen-project:
    uv run gen-project src/ai_gene_review/schema/gene_review.yaml -d assets

pydantic:
    uv run gen-pydantic src/ai_gene_review/schema/gene_review.yaml > src/ai_gene_review/datamodel/gene_review_model.py.tmp && mv src/ai_gene_review/datamodel/gene_review_model.py.tmp src/ai_gene_review/datamodel/gene_review_model.py

gen-all: gen-project pydantic