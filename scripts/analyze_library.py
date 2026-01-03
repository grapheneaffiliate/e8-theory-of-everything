#!/usr/bin/env python3
"""
GSM Library Analyzer
====================

Analyzes the refined library from the Self-Correcting Engine.
Shows best discoveries for each physical constant target.

Usage:
    python analyze_library.py
    python analyze_library.py gsm_refined_library.json
"""

import json
import sys
import os

# Default library file
LIBRARY_FILE = "gsm_refined_library.json"


def load_library(filepath):
    """Load the library from JSON file."""
    if not os.path.exists(filepath):
        print(f"Error: Library file '{filepath}' not found.")
        print("\nTo create a library, run:")
        print("    python gsm_self_correcting_engine.py")
        return None
    
    with open(filepath, 'r', encoding='utf-8') as f:
        return json.load(f)


def analyze_library(data):
    """Analyze and display library contents."""
    
    print("=" * 70)
    print("GSM REFINED LIBRARY ANALYSIS")
    print("=" * 70)
    
    metadata = data.get('metadata', {})
    theorems = data.get('theorems', [])
    
    print(f"\nLast updated: {metadata.get('last_updated', 'Unknown')}")
    print(f"Version: {metadata.get('version', 'Unknown')}")
    print(f"Total theorems: {len(theorems)}")
    
    # Count tags
    tags = {}
    for t in theorems:
        for tag in t.get('tags', []):
            tags[tag] = tags.get(tag, 0) + 1
    
    print(f"\nTag distribution:")
    for tag, count in sorted(tags.items(), key=lambda x: -x[1]):
        print(f"  {tag}: {count}")
    
    # Group by target
    targets = {}
    for t in theorems:
        for tag in t.get('tags', []):
            if 'Target:' in tag:
                target_name = tag.replace('Target:', '').replace('_Refined', '').replace('_Rough', '')
                if target_name not in targets:
                    targets[target_name] = []
                targets[target_name].append(t)
    
    # Show best for each target
    print("\n" + "=" * 70)
    print("BEST DISCOVERIES BY TARGET")
    print("=" * 70)
    
    for target_name in sorted(targets.keys()):
        results = targets[target_name]
        # Sort by error
        results.sort(key=lambda x: x.get('error', float('inf')))
        
        best = results[0]
        
        print(f"\n{target_name.upper()}")
        print("-" * 50)
        print(f"  Expression: {best['expression']}")
        print(f"  Value:      {best['value']:.15f}")
        print(f"  Error:      {best.get('error', 0):.2e}")
        print(f"  Error ppm:  {best.get('error_ppm', 0):.4f}")
        print(f"  Complexity: {best.get('complexity', 'N/A')}")
        
        if 'base_expression' in best:
            print(f"  Base:       {best['base_expression']}")
        if 'correction_term' in best:
            print(f"  Correction: {best['correction_term']}")
        
        # Show top 3 alternatives
        if len(results) > 1:
            print(f"\n  Alternatives ({len(results)-1} more):")
            for alt in results[1:4]:
                ppm = alt.get('error_ppm', 0)
                print(f"    {ppm:.4f} ppm: {alt['expression'][:60]}...")

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    all_errors = [t.get('error_ppm', float('inf')) for t in theorems if 'error_ppm' in t]
    if all_errors:
        sub_ppm = sum(1 for e in all_errors if e < 1)
        sub_10ppm = sum(1 for e in all_errors if e < 10)
        sub_100ppm = sum(1 for e in all_errors if e < 100)
        
        print(f"\nPrecision breakdown:")
        print(f"  < 1 ppm:   {sub_ppm} theorems")
        print(f"  < 10 ppm:  {sub_10ppm} theorems")
        print(f"  < 100 ppm: {sub_100ppm} theorems")
        print(f"  Total:     {len(all_errors)} theorems with error data")


def export_best_formulas(data, output_file="best_formulas.txt"):
    """Export the best formula for each target to a text file."""
    
    theorems = data.get('theorems', [])
    
    targets = {}
    for t in theorems:
        for tag in t.get('tags', []):
            if 'Target:' in tag:
                target_name = tag.replace('Target:', '').replace('_Refined', '').replace('_Rough', '')
                if target_name not in targets:
                    targets[target_name] = []
                targets[target_name].append(t)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("GSM Self-Correcting Engine - Best Formulas\n")
        f.write("=" * 50 + "\n\n")
        
        for target_name in sorted(targets.keys()):
            results = sorted(targets[target_name], key=lambda x: x.get('error', float('inf')))
            best = results[0]
            
            f.write(f"{target_name}:\n")
            f.write(f"  {best['expression']}\n")
            f.write(f"  Value: {best['value']:.15f}\n")
            f.write(f"  Error: {best.get('error_ppm', 0):.4f} ppm\n\n")
    
    print(f"\nBest formulas exported to: {output_file}")


def main():
    # Get library file from command line or use default
    library_file = sys.argv[1] if len(sys.argv) > 1 else LIBRARY_FILE
    
    # Load library
    data = load_library(library_file)
    if data is None:
        return
    
    # Analyze
    analyze_library(data)
    
    # Export best formulas
    export_best_formulas(data)


if __name__ == "__main__":
    main()
