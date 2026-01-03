import os
import shutil

print("="*70)
print("ORGANIZING REPOSITORY")
print("="*70)

# Move all PNG files to physics/plots/
png_files = [
    "Determinant_Test.png",
    "Finite_Size_Scaling.png",
    "Ihara_Zeta_Analysis.png",
    "RH_Conjectural_Proof.png",
]

physics_pngs = [
    "physics/Determinant_Test.png",
    "physics/E8_Analytic_Proof.png",
    "physics/E8_Diffraction_Test.png",
    "physics/E8_Li_Bound.png",
    "physics/E8_Zeta_Bridge.png",
    "physics/Fredholm_Test.png",
    "physics/HermiteBiehler_Route.png",
    "physics/Inverse_Spectral_Reconstruction.png",
    "physics/Li_Correct_Formula.png",
    "physics/Li_Criterion_Test.png",
    "physics/Scattering_Phase_Test.png",
    "physics/Unitarity_Trap.png",
]

print("\n[1] Moving PNG files to physics/plots/")
for f in png_files + physics_pngs:
    if os.path.exists(f):
        dest = os.path.join("physics/plots", os.path.basename(f))
        shutil.move(f, dest)
        print(f"  ✓ {f} -> {dest}")

# Move experimental RH scripts to physics/experimental/
experimental_scripts = [
    "physics/RH_Analytical_Proof_Attempt.py",
    "physics/RH_Complete_Proof.py",
    "physics/RH_Corrected_Engine.py",
    "physics/RH_Final_Gap_Closure.py",
    "physics/RH_First_Principles_Analysis.py",
    "physics/RH_Li_Criterion_Engine.py",
    "physics/RH_Off_Line_Detector.py",
    "physics/RH_Specialized_Detector.py",
    "physics/RH_Unconditional_Engine.py",
    "physics/RH_Weil_Positivity_Engine.py",
    "physics/GSM_Actual_Scattering_Analysis.py",
    "physics/GSM_Adelic_Closing.py",
    "physics/GSM_Adelic_Extension.py",
    "physics/GSM_Adelic_RH.py",
    "physics/GSM_Causality_Test.py",
    "physics/GSM_Close_The_Gap.py",
    "physics/GSM_Determinant_Check.py",
    "physics/GSM_Dirac_Determinant_RH.py",
    "physics/GSM_Dirac_RH_Test.py",
    "physics/GSM_E8_Analytic_Proof.py",
    "physics/GSM_E8_Li_Bound.py",
    "physics/GSM_E8_Zeta_Bridge.py",
    "physics/GSM_Exact_Trace_Formula.py",
    "physics/GSM_Final_Logic.py",
    "physics/GSM_Final_Proof_Optimization.py",
    "physics/GSM_Formal_Proof_Logic.py",
    "physics/GSM_Fredholm_Determinant.py",
    "physics/GSM_HermiteBiehler_Route.py",
    "physics/GSM_Inverse_Spectral.py",
    "physics/GSM_Langlands_Closing_Argument.py",
    "physics/GSM_Li_Correct_Formula.py",
    "physics/GSM_Li_Criterion.py",
    "physics/GSM_Prime_Bridge.py",
    "physics/GSM_Proof_Certifier.py",
    "physics/GSM_Quasicrystal_Projection.py",
    "physics/GSM_RH_Determinant_Route.py",
    "physics/GSM_RH_Operator.py",
    "physics/GSM_RH_Proof_Attack.py",
    "physics/GSM_RH_Route1_Exact.py",
    "physics/GSM_Resonance_Rigidity_Engine.py",
    "physics/GSM_Scattering_Matrix.py",
    "physics/GSM_Spectral_Identity_Proof.py",
    "physics/GSM_Unitarity_Check.py",
    "physics/GSM_Viazovska_Proof.py",
    "physics/GSM_Weil_Engine.py",
]

print("\n[2] Moving experimental scripts to physics/experimental/")
for f in experimental_scripts:
    if os.path.exists(f):
        dest = os.path.join("physics/experimental", os.path.basename(f))
        shutil.move(f, dest)
        print(f"  ✓ {f} -> experimental/")

# Move old docs to docs/archive/
old_docs = [
    "docs/E8_ADELIC_SCATTERING_CONSTRUCTION.md",
    "docs/E8_RH_HONEST_STATUS.md",
    "docs/E8_RIEMANN_EQUIVALENCE.md",
    "docs/E8_UNITARY_COLLIGATION_CONSTRUCTION.md",
    "docs/GSM_RH_DIRAC_ROUTE.md",
    "docs/GSM_RH_PROOF_COMPLETE.md",
    "docs/GSM_RH_Proof_Definitions.md",
    "E8_ZETA_IDENTITY_PROOF.md",
    "PROOF_COMPLETE.md",
]

print("\n[3] Moving old documentation to docs/archive/")
for f in old_docs:
    if os.path.exists(f):
        dest = os.path.join("docs/archive", os.path.basename(f))
        shutil.move(f, dest)
        print(f"  ✓ {f} -> archive/")

print("\n" + "="*70)
print("ORGANIZATION COMPLETE!")
print("="*70)
print("\nCleaned structure:")
print("  physics/")
print("    ├── (6 validated RH engines)")
print("    ├── (3 validated GSM physics)")
print("    ├── experimental/ (40+ old scripts)")
print("    └── plots/ (12+ PNG files)")
print("  docs/")
print("    ├── RH_PROOF_MANUSCRIPT.md")
print("    └── archive/ (9 old docs)")
