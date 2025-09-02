#!/usr/bin/env python3
"""CLI wrapper for OpenAI Deep Research API for gene research."""

import sys
from pathlib import Path
from typing import Optional

import click
from openai import OpenAI

from ai_gene_review.etl.gene import expand_organism_name


@click.command()
@click.argument("gene_symbol", type=str)
@click.option("--organism", default="human", help="Organism (default: human)")
@click.option(
    "--model", default="o3-deep-research-2025-06-26", help="Deep research model to use"
)
@click.option(
    "--output-dir",
    type=click.Path(),
    help="Output directory for results (defaults to genes/<organism>/<gene>/)",
)
@click.option(
    "--api-key",
    envvar="OPENAI_API_KEY",
    help="OpenAI API key (or set OPENAI_API_KEY env var)",
)
@click.option(
    "--save-reasoning",
    is_flag=True,
    default=False,
    help="Save reasoning steps to debug file (default: off)",
)
@click.option(
    "--system-prompt-file",
    type=click.Path(exists=True),
    help="File containing custom system prompt (optional)",
)
def research_gene(
    gene_symbol: str,
    organism: str,
    model: str,
    output_dir: Optional[str],
    api_key: Optional[str],
    save_reasoning: bool,
    system_prompt_file: Optional[str],
):
    """Research a gene using OpenAI's Deep Research API.

    This tool performs comprehensive research on a given gene symbol,
    gathering information about its function, structure, disease associations,
    and relevant literature.

    Examples:
        deep-research CFAP300
        deep-research TP53 --organism human --output-dir custom/path
    """
    if not api_key:
        click.echo(
            "Error: OpenAI API key required. Set OPENAI_API_KEY environment variable or use --api-key option.",
            err=True,
        )
        sys.exit(1)

    # Initialize OpenAI client
    client = OpenAI(api_key=api_key)

    # Determine output directory
    if output_dir:
        output_path = Path(output_dir)
    else:
        output_path = Path("genes") / organism / gene_symbol

    output_path.mkdir(parents=True, exist_ok=True)

    # Load system prompt (custom or default)
    if system_prompt_file:
        with open(system_prompt_file, "r", encoding="utf-8") as f:
            system_message = f.read()
    else:
        system_message = """You are a molecular biologist and gene annotation expert conducting comprehensive research to support GO annotation curation.

Provide detailed, well-cited information focusing on:
1. Gene function and molecular mechanisms
2. Cellular localization and subcellular components  
3. Biological processes involvement
4. Disease associations and phenotypes
5. Protein domains and structural features
6. Expression patterns and regulation
7. Evolutionary conservation
8. Key experimental evidence and literature

Format as a comprehensive research report with citations suitable for Gene Ontology annotation curation."""

    # Expand organism name for better research results
    organism_full = expand_organism_name(organism)

    user_query = f"""Research the {organism_full} gene {gene_symbol}. Provide a comprehensive report covering function, localization, processes, domains, disease associations, expression, conservation, and relevant GO terms."""

    click.echo(f"Starting deep research on {gene_symbol} ({organism_full})...")
    click.echo(f"Output will be saved to: {output_path}")

    try:
        # Make the Deep Research API call
        click.echo("⏳ This may take several minutes to complete...")
        response = client.responses.create(
            model=model,
            input=[
                {
                    "role": "developer",
                    "content": [{"type": "input_text", "text": system_message}],
                },
                {
                    "role": "user",
                    "content": [{"type": "input_text", "text": user_query}],
                },
            ],
            tools=[{"type": "web_search_preview"}],
        )

        # Extract the final report
        final_output = response.output[-1]
        if hasattr(final_output, "content") and final_output.content:  # type: ignore
            content = final_output.content  # type: ignore
            if isinstance(content, list) and len(content) > 0:
                first_content = content[0]
                if hasattr(first_content, "text"):
                    report_text = first_content.text  # type: ignore
                else:
                    report_text = str(first_content)
            else:
                report_text = str(content)
        else:
            report_text = str(final_output)

        # Save the research report
        report_file = output_path / f"{gene_symbol}-deep-research.md"
        with open(report_file, "w", encoding="utf-8") as f:
            f.write(f"# Deep Research Report: {gene_symbol} ({organism})\n\n")
            f.write("Generated using OpenAI Deep Research API\n\n")
            f.write("---\n\n")
            f.write(report_text)

        click.echo(f"✅ Research completed and saved to: {report_file}")

        # Extract and save citations if available
        if hasattr(final_output, "content") and final_output.content:  # type: ignore
            content = final_output.content  # type: ignore
            if isinstance(content, list) and len(content) > 0:
                first_content = content[0]
                if hasattr(first_content, "annotations"):
                    annotations = first_content.annotations  # type: ignore
                    if annotations:
                        citations_file = output_path / f"{gene_symbol}-citations.md"
                        with open(citations_file, "w", encoding="utf-8") as f:
                            f.write(f"# Citations for {gene_symbol} Research\n\n")
                            for i, annotation in enumerate(annotations, 1):
                                f.write(f"{i}. {annotation}\n")
                        click.echo(f"📚 Citations saved to: {citations_file}")

        # Save reasoning steps if available and requested
        if save_reasoning:
            reasoning_steps = [
                item
                for item in response.output
                if hasattr(item, "type") and item.type == "reasoning"
            ]  # type: ignore
            if reasoning_steps:
                reasoning_file = output_path / f"{gene_symbol}-reasoning.md"
                with open(reasoning_file, "w", encoding="utf-8") as f:
                    f.write(f"# Reasoning Steps for {gene_symbol} Research\n\n")
                    for i, step in enumerate(reasoning_steps, 1):
                        f.write(f"## Step {i}\n\n")
                        if hasattr(step, "content"):  # type: ignore
                            step_content = step.content  # type: ignore
                            if isinstance(step_content, list) and len(step_content) > 0:
                                if hasattr(step_content[0], "text"):
                                    f.write(f"{step_content[0].text}\n\n")  # type: ignore
                                else:
                                    f.write(f"{step_content[0]}\n\n")
                            else:
                                f.write(f"{step_content}\n\n")
                        else:
                            f.write(f"{step}\n\n")
                click.echo(f"🧠 Reasoning steps saved to: {reasoning_file}")
            else:
                click.echo("ℹ️  No reasoning steps available to save")

    except Exception as e:
        click.echo(f"❌ Error during research: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    research_gene()
