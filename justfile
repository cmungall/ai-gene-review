# ============ Hint for for Windows Users ============

# On Windows the "sh" shell that comes with Git for Windows should be used.
# If it is not on path, provide the path to the executable in the following line.
#set windows-shell := ["C:/Program Files/Git/usr/bin/sh", "-cu"]

# ============ Variables used in recipes ============

# Set shebang line for cross-platform Python recipes (assumes presence of launcher on Windows)
shebang := if os() == 'windows' {
  'py'
} else {
  '/usr/bin/env python3'
}


# ============== Project recipes ==============

# List all commands as default command. The prefix "_" hides the command.
_default: _status
    @just --list

# Initialize a new project (use this for projects not yet under version control)
[group('project management')]
setup: _git-init _ai-instructions install _git-add
  git commit -m "Initialise git with minimal project" -a

  
# Install project dependencies
[group('project management')]
install:
  uv sync --group dev


# Run all tests
[group('model development')]
test: pytest doctest mypy format

test-full: test pytest-integration

pytest:
  uv run pytest tests

# include integration tests
pytest-integration:
	$(RUN) pytest -m ""

doctest:
  uv run pytest  --doctest-modules src

mypy:
  uv run mypy src tests

format:
	uv run ruff check .

# ============== Hidden internal recipes ==============

_status:
  @echo "OK"

# Update project template
_update-template:
  copier update --trust --skip-answered


# Run documentation server
_serve:
  uv run mkdocs serve

# Initialize git repository
_git-init:
  git init

# Add files to git
_git-add:
  git add .

# Commit files to git
_git-commit:
  git commit -m 'chore: just setup was run' -a

# Show git status
_git-status:
  git status

goosehints:
  [ -f .goosehints ] || ln -s CLAUDE.md .goosehints

copilot-instructions:
  [ -f .github/copilot-instructions.md ] || cd .github && ln -s ../CLAUDE.md copilot-instructions.md

_ai-instructions: goosehints copilot-instructions

gh-add-topics:
  gh repo edit --add-topic "ai-gene-review,monarchinitiative,linkml"

gh-add-secrets:
  gh secret set PAT_FOR_PR --body "$PAT_FOR_PR"
  gh secret set ANTHROPIC_API_KEY --body "$ANTHROPIC_API_KEY"
  gh secret set OPENAI_API_KEY --body "$OPENAI_API_KEY"
  gh secret set CBORG_API_KEY --body "$CBORG_API_KEY"
  gh secret set CLAUDE_CODE_OATH_TOKEN --body "$CLAUDE_CODE_OATH_TOKEN"

gh-invite-the-dragon:
  gh api repos/monarch-initiative/ai-gene-review/collaborators/dragon-ai-agent -X PUT -f permission=push

# ============== Include project-specific recipes ==============

import "project.justfile"
