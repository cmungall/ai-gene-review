# Claude Code Hooks

This directory contains hooks that integrate with Claude Code to provide automated validation and other features.

## validate_ai_review_hook.py

This hook automatically validates any `-ai-review.yaml` files when they are written or edited using Claude Code.

### Features
- Automatically runs validation when saving ai-review.yaml files
- Blocks file modifications if validation fails
- Shows detailed validation output in the Claude Code interface
- Filters out noise from warning messages for cleaner output

### How it works
1. Intercepts Write, Edit, and MultiEdit operations on files ending with `-ai-review.yaml`
2. Runs the validation command: `uv run ai-gene-review validate <file>`
3. Displays validation results
4. Returns exit code 1 to block the operation if validation fails

### Testing
Run the test script to verify the hook is working:
```bash
python3 .claude/hooks/test_validation_hook.py
```

### Protection Hook
The `protect_files_hook.py` prevents modifications to certain protected paths including:
- publications/
- genes/ directories (except ai-review.yaml files)
- deep-research.md files
- goa.tsv files
- .uniprot.txt files