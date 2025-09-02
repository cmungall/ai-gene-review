#!/usr/bin/env python3
"""
Test script for the validation hook.
This script simulates the hook environment to test different scenarios.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path


def test_hook(tool_name: str, file_path: str, expected_exit_code: int = 0):
    """Test the validation hook with given parameters."""
    # Prepare the input data that Claude Code would send
    hook_input = {"tool_name": tool_name, "tool_input": {"file_path": file_path}}

    # Run the hook
    result = subprocess.run(
        [sys.executable, ".claude/hooks/validate_ai_review_hook.py"],
        input=json.dumps(hook_input),
        text=True,
        capture_output=True,
    )

    print(f"\nTest: {tool_name} on {file_path}")
    print(f"Expected exit code: {expected_exit_code}, Got: {result.returncode}")

    if result.stderr:
        print("Hook output:")
        print(result.stderr)

    assert result.returncode == expected_exit_code, (
        f"Expected exit code {expected_exit_code}, got {result.returncode}"
    )
    return result


def main():
    print("Testing validation hook...")

    # Test 1: Non-ai-review file (should pass through)
    print("\n" + "=" * 60)
    print("Test 1: Non-ai-review file - should be ignored")
    test_hook("Write", "/some/path/regular-file.yaml", expected_exit_code=0)

    # Test 2: Non-Write/Edit tool (should pass through)
    print("\n" + "=" * 60)
    print("Test 2: Read tool on ai-review file - should be ignored")
    test_hook("Read", "/some/path/test-ai-review.yaml", expected_exit_code=0)

    # Test 3: Valid ai-review file (should pass)
    print("\n" + "=" * 60)
    print("Test 3: Valid ai-review file - should pass validation")
    # We'll use an existing valid file
    valid_file = "genes/human/RBFOX3/RBFOX3-ai-review.yaml"
    if Path(valid_file).exists():
        test_hook("Write", valid_file, expected_exit_code=0)
    else:
        print(f"Skipping - {valid_file} not found")

    # Test 4: Invalid ai-review file (should fail)
    print("\n" + "=" * 60)
    print("Test 4: Invalid ai-review file - should block operation")

    # Create a temporary invalid ai-review file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix="-ai-review.yaml", delete=False
    ) as f:
        f.write("""
id: TEST001
gene_symbol: INVALID
organism: human
gene_id: "999999999"
gene_name: "Invalid test gene"
existing_annotations:
  - go_id: GO:9999999  # Invalid GO ID
    go_term_label: "Invalid GO term"
    evidence_type: IEA
    original_reference_id: PMID:12345678
    date: "2024-01-01"
    aspect: molecular_function
    review:
      rationale: "Test"
      is_consistent: true
      suggested_changes: []
""")
        temp_file = f.name

    try:
        # This should fail validation and return exit code 1
        test_hook("Write", temp_file, expected_exit_code=1)
    finally:
        # Clean up
        Path(temp_file).unlink(missing_ok=True)

    # Test 5: MultiEdit tool
    print("\n" + "=" * 60)
    print("Test 5: MultiEdit tool on ai-review file")
    test_hook("MultiEdit", "/some/path/test-ai-review.yaml", expected_exit_code=0)

    print("\n" + "=" * 60)
    print("✅ All hook tests passed!")


if __name__ == "__main__":
    main()
