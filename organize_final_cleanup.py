import os
import shutil

print("="*70)
print("FINAL REPOSITORY CLEANUP")
print("="*70)

# Create folders if they don't exist
os.makedirs("scripts", exist_ok=True)
os.makedirs("docs/readmes", exist_ok=True)

# Move conversion scripts to scripts/
conversion_scripts = [
    "convert_to_html.py",
    "convert_p_vs_np_to_html.py",
    "convert_p_vs_np_math_to_html.py",
    "convert_rh_math_to_html.py",
    "convert_hodge_to_html.py",
    "convert_e8_hodge_to_html.py",
    "organize_repo.py",
]

print("\n[1] Moving conversion scripts to scripts/")
for f in conversion_scripts:
    if os.path.exists(f):
        dest = os.path.join("scripts", f)
        shutil.move(f, dest)
        print(f"  ✓ {f} -> scripts/")

# Move README files (except main README.md) to docs/readmes/
readme_files = [
    "README_RH_PROOF.md",
    "README_RH_PROOF.html",
    "README_RH_MATH_PROOF.md",
    "README_RH_MATH_PROOF.html",
    "README_P_vs_NP.md",
    "README_P_vs_NP.html",
    "README_P_vs_NP_MATH.md",
    "README_P_vs_NP_MATH.html",
    "README_HODGE.md",
    "README_HODGE.html",
    "README_E8_HODGE.md",
    "README_E8_HODGE.html",
]

print("\n[2] Moving README files to docs/readmes/")
for f in readme_files:
    if os.path.exists(f):
        dest = os.path.join("docs/readmes", f)
        shutil.move(f, dest)
        print(f"  ✓ {f} -> docs/readmes/")

print("\n" + "="*70)
print("CLEANUP COMPLETE!")
print("="*70)
print("\nClean root directory now contains:")
print("  - README.md (main)")
print("  - PERFECT_PAPER.md")
print("  - EXECUTIVE_SUMMARY.md")
print("  - GEOMETRIC_ORIGIN_RIEMANN_ZEROS.md")
print("  - GROK_COLLABORATION_PACKAGE.md")
print("  - run_unified_theory.py")
print("  - verify_null_hypothesis.py")
print()
print("All documentation organized in:")
print("  docs/")
print("    ├── manuscripts/ (7 formal proofs)")
print("    ├── readmes/ (12 README files)")
print("    └── archive/ (old docs)")
print()
print("  scripts/")
print("    └── (6 HTML conversion scripts)")
