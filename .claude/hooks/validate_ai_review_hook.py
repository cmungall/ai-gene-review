#!/usr/bin/env python3
"""
Hook to automatically validate ai-review.yaml files after they are written or edited.
This hook runs the validation tool and displays the results to provide immediate feedback.
"""

import sys
import json
import subprocess
import os
from pathlib import Path


def main():
    # Read the hook input from stdin
    data = json.load(sys.stdin)

    # Extract the file path from the tool input
    tool_name = data.get("tool_name", "")
    file_path = data.get("tool_input", {}).get("file_path", "")

    # Only process Write and Edit tool calls
    if tool_name not in ["Write", "Edit", "MultiEdit"]:
        sys.exit(0)

    # Check if this is an ai-review.yaml file
    if not file_path.endswith("-ai-review.yaml"):
        sys.exit(0)

    # Convert to Path object for easier manipulation
    file_path = Path(file_path)

    # Check if the file exists (it should after Write/Edit)
    if not file_path.exists():
        print(f"⚠️ File not found: {file_path}", file=sys.stderr)
        sys.exit(0)

    # Run the validation command
    try:
        # Build the validation command
        cmd = ["uv", "run", "ai-gene-review", "validate", str(file_path)]

        # Run the command and capture output
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=os.path.dirname(
                os.path.dirname(os.path.dirname(__file__))
            ),  # Project root
        )

        # Display the validation output
        print("\n" + "=" * 60, file=sys.stderr)
        print(f"🔍 Validation Results for {file_path.name}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        # Show stdout (the actual validation results)
        if result.stdout:
            # Filter out the uv warning lines
            lines = result.stdout.split("\n")
            filtered_lines = []
            skip_next = False
            for line in lines:
                if "warning: Failed to uninstall package" in line:
                    skip_next = True
                    continue
                elif skip_next and ("Uninstalled" in line or "Installed" in line):
                    skip_next = False
                    continue
                elif "/eutils/__init__.py" in line and "UserWarning" in line:
                    continue
                elif "pkg_resources is deprecated" in line:
                    continue
                else:
                    filtered_lines.append(line)

            output = "\n".join(filtered_lines).strip()
            if output:
                print(output, file=sys.stderr)

        # Show any errors
        if result.returncode != 0 and result.stderr:
            # Filter stderr similarly
            lines = result.stderr.split("\n")
            filtered_lines = []
            for line in lines:
                if "/eutils/__init__.py" not in line:
                    filtered_lines.append(line)

            error_output = "\n".join(filtered_lines).strip()
            if error_output:
                print("\n⚠️ Validation errors:", file=sys.stderr)
                print(error_output, file=sys.stderr)

        print("=" * 60 + "\n", file=sys.stderr)

        # Return non-zero exit code if validation failed
        if result.returncode != 0:
            print("❌ Validation failed - blocking file modification", file=sys.stderr)
            print("Fix validation errors before saving the file.", file=sys.stderr)
            sys.exit(1)  # Block the operation

    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to run validation: {e}", file=sys.stderr)
        sys.exit(1)  # Block on validation errors
    except Exception as e:
        print(f"❌ Unexpected error during validation: {e}", file=sys.stderr)
        # Allow operation to proceed on unexpected errors (don't block on hook failures)
        sys.exit(0)

    # Exit 0 if validation passed
    sys.exit(0)


if __name__ == "__main__":
    main()
