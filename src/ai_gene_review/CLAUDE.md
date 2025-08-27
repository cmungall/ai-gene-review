## IMPORTANT INSTRUCTIONS FOR CODE

### Best Practices

- Use test driven development, write tests first before implementing a feature
- Use doctests liberally - these serve as both explanatory examples for humans and as unit tests
- For longer examples, write pytest tests
- Always write pytest functional style rather than unittest OO style
- Use modern pytest idioms, including `@pytest.mark.parametrize` to test for combinations of inputs
- NEVER write mock tests unless requested. I need to rely on tests to know if something breaks
- For tests that have external dependencies, use `@pytest.mark.integration`
- Do not "fix" issues by changing or weakening test conditions. Try harder, or ask questions if a test fails
- Always run tests
- Use docstrings

### Code Style

- Avoid try/except blocks, these can mask bugs (except when interfacing with external systems)
- Fail fast is a good principle
- Follow the DRY principle
- Avoid repeating chunks of code, but also avoid premature over-abstraction
- Pydantic or LinkML is favored for data objects
- For state in engine-style OO classes, dataclasses is favored
- Declarative principles are favored
- Always use type hints, always document methods and classes
- **NEVER write scripts that just tell users to go to a website** - if you can't programmatically call a tool, be honest about it and write scripts that parse output instead
- Don't write misleading code that pretends to do something it doesn't actually do

## Essential Commands


### Building and Testing
- `just --list` - See all commands
- `just test` - Performs unit tests, doctests, ruff/linting
- `just test-full` - As above plus integration tests
- `just pytest` - Run Python tests only
- `just mypy` - Run type checking
- `just format` - Run ruff linting/formatting checks
- `uv run pytest tests/test_simple.py::test_simple` - Run a specific test

You can run the underlying commands (with `uv run ...`) but in general justfile targets should be favored.

### Running the CLI
- `uv run ai-gene-review --help` - Run the CLI tool with options

### Documentation
- `just _serve` - Run local documentation server with mkdocs

## Project Architecture

### Layout
- **src/ai_gene_review/** - Code goes here
  - `cli.py` - Typer-based CLI interface, entry point for the application
- **tests/** - Test suite using pytest with parametrized tests
  - `input/` - Example files for testing
- **docs/** - MkDocs documentation
- **mkdocs.yml** - Index of docs

### Technology Stack
- **Python 3.10+** with `uv` for dependency management
- **LinkML** for data modeling (linkml-runtime)
- **Typer** for CLI interface
- **pytest** for testing
- **mypy** for type checking
- **ruff** for linting and formatting
- **MkDocs Material** for documentation

### Key Configuration Files
- `pyproject.toml` - Python project configuration, dependencies, and tool settings
- `justfile` - Command runner recipes for common development tasks
- `mkdocs.yml` - Documentation configuration
- `uv.lock` - Locked dependency versions

## Development Workflow

1. Dependencies are managed via `uv` - use `uv add` for new dependencies
2. All commands are run through `just` or `uv run`
3. The project uses dynamic versioning from git tags
4. Documentation is auto-deployed to GitHub Pages at https://monarch-initiative.github.io/my-awesome-tool
