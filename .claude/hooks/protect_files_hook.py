#!/usr/bin/env python3

import sys
import json

PROTECTED_PATH_FRAGMENTS = [
    "publications",
    "deep-research.md",
    "goa.tsv",
    ".uniprot.txt",
]

data = json.load(sys.stdin)
path = data.get("tool_input", {}).get("file_path", "")

sys.exit(2 if any(p in path for p in PROTECTED_PATH_FRAGMENTS) else 0)
