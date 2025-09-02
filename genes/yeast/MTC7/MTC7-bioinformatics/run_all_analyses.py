#!/usr/bin/env python3
"""
Master script to run all MTC7 bioinformatics analyses.
"""

import subprocess
import sys
from pathlib import Path


def run_analysis(script_name, description):
    """Run an analysis script and handle errors."""
    print(f"\n{'=' * 60}")
    print(f"Running: {description}")
    print('=' * 60)
    
    try:
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        if result.stderr:
            print("Warnings:", result.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}:")
        print(e.stdout)
        print(e.stderr)
        return False


def main():
    print("MTC7 Comprehensive Bioinformatics Analysis")
    print("=" * 60)
    
    analyses = [
        ("analyze_transmembrane.py", "Transmembrane Region Analysis"),
        ("analyze_domains.py", "Domain and Motif Analysis"),
        ("analyze_conservation.py", "Conservation Analysis"),
        ("analyze_structure.py", "Structural Feature Analysis")
    ]
    
    results = []
    
    for script, description in analyses:
        success = run_analysis(script, description)
        results.append((description, success))
    
    print("\n" + "=" * 60)
    print("Analysis Summary:")
    print("=" * 60)
    
    for description, success in results:
        status = "✓ Complete" if success else "✗ Failed"
        print(f"{status}: {description}")
    
    print("\nGenerated files:")
    output_files = [
        "tm_analysis_results.json",
        "mtc7_hydropathy.png",
        "domain_analysis_results.json",
        "mtc7_charge_distribution.png",
        "conservation_analysis_results.json",
        "structure_analysis_results.json",
        "mtc7_structural_features.png"
    ]
    
    for file in output_files:
        if Path(file).exists():
            print(f"  ✓ {file}")
        else:
            print(f"  ✗ {file} (not found)")
    
    print("\nAnalysis complete! Check RESULTS.md for a comprehensive summary.")


if __name__ == "__main__":
    main()