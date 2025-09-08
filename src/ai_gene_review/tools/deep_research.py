#!/usr/bin/env python3
"""CLI wrapper for OpenAI Deep Research API for gene research."""

import sys
import time
import traceback
from pathlib import Path
from typing import Optional

import click
import httpx
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
@click.option(
    "--timeout",
    type=int,
    default=600,
    help="Request timeout in seconds (default: 600, i.e., 10 minutes)",
)
@click.option(
    "--max-retries",
    type=int,
    default=3,
    help="Maximum number of retries on timeout (default: 3)",
)
def research_gene(
    gene_symbol: str,
    organism: str,
    model: str,
    output_dir: Optional[str],
    api_key: Optional[str],
    save_reasoning: bool,
    system_prompt_file: Optional[str],
    timeout: int,
    max_retries: int,
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

    # Initialize OpenAI client with custom timeout
    http_client = httpx.Client(
        timeout=httpx.Timeout(
            connect=30.0,  # Connection timeout
            read=timeout,  # Read timeout
            write=30.0,  # Write timeout
            pool=30.0,  # Pool timeout
        )
    )
    client = OpenAI(api_key=api_key, http_client=http_client)

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
    click.echo(f"Timeout set to: {timeout} seconds")
    click.echo(f"Max retries: {max_retries}")

    response = None
    last_error = None
    
    for attempt in range(max_retries):
        try:
            # Make the Deep Research API call
            if attempt > 0:
                click.echo(f"\n🔄 Retry attempt {attempt + 1}/{max_retries}...")
            else:
                click.echo("\n⏳ This may take several minutes to complete...")
            
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
            break  # Success, exit retry loop
            
        except httpx.ReadTimeout as e:
            last_error = e
            click.echo(f"\n⚠️  Timeout error (attempt {attempt + 1}/{max_retries}): {e}", err=True)
            if attempt < max_retries - 1:
                click.echo("Will retry in 5 seconds...", err=True)
                import time
                time.sleep(5)
            
        except httpx.HTTPStatusError as e:
            last_error = e
            click.echo(f"\n⚠️  HTTP error (attempt {attempt + 1}/{max_retries}): {e}", err=True)
            click.echo(f"Status code: {e.response.status_code}", err=True)
            if attempt < max_retries - 1:
                click.echo("Will retry in 5 seconds...", err=True)
                import time
                time.sleep(5)
                
        except Exception as e:
            last_error = e
            click.echo(f"\n⚠️  Unexpected error (attempt {attempt + 1}/{max_retries}): {type(e).__name__}: {e}", err=True)
            if attempt < max_retries - 1:
                click.echo("Will retry in 5 seconds...", err=True)
                import time
                time.sleep(5)
    
    if response is None:
        click.echo(f"\n❌ Failed after {max_retries} attempts", err=True)
        click.echo("\n📋 Full stacktrace:", err=True)
        click.echo("-" * 60, err=True)
        traceback.print_exc()
        click.echo("-" * 60, err=True)
        
        if isinstance(last_error, httpx.ReadTimeout):
            click.echo("\n💡 Suggestions:", err=True)
            click.echo("  - Try increasing the timeout with --timeout option (e.g., --timeout 1200 for 20 minutes)", err=True)
            click.echo("  - The deep research models can take 10-20 minutes for complex queries", err=True)
            click.echo("  - Check OpenAI service status at https://status.openai.com/", err=True)
        
        sys.exit(1)
    
    try:

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
        click.echo(f"\n❌ Error processing response: {type(e).__name__}: {e}", err=True)
        click.echo("\n📋 Full stacktrace:", err=True)
        click.echo("-" * 60, err=True)
        traceback.print_exc()
        click.echo("-" * 60, err=True)
        sys.exit(1)
    
    finally:
        # Clean up HTTP client
        if 'http_client' in locals():
            http_client.close()


if __name__ == "__main__":
    research_gene()
