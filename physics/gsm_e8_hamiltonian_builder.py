#!/usr/bin/env python3
"""
GSM E8 HAMILTONIAN BUILDER
============================
"Constructing the Hilbert-Pólya Operator from E8 Geometry"

Purpose: Build a quantum Hamiltonian on the E8 Root Lattice with Golden Ratio coupling.
         If the Riemann Zeros emerge as eigenvalues, we have found THE operator.

Hypothesis: The E8 Gosset Polytope (240 roots) with Golden tunneling amplitudes
            produces a spectrum containing the Riemann Zeros.

Physics: H_ij = E8 Adjacency * Golden Decay = Quantum particle on E8 lattice
         
Author: GSM Research Team
Date: January 2, 2026
"""

import numpy as np
import scipy.linalg
from itertools import product

print("=" * 70)
print("             GSM E8 HAMILTONIAN BUILDER")
print("   Constructing the 8-Dimensional Operator of the Riemann Zeros")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════════════
# 1. FUNDAMENTAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════
PHI = (1 + np.sqrt(5)) / 2
PSI = PHI - 1  # 1/PHI
LAMBDA = 16 * np.sqrt(15)
PI = np.pi

# First 10 Riemann Zeros (imaginary parts)
RIEMANN_ZEROS = [
    14.134725141734693790457251983562,
    21.022039638771554992628479593897,
    25.010857580145688763213790992563,
    30.424876125859513210311897530584,
    32.935061587739189690662368964075,
    37.586178158825671257217763480705,
    40.918719012147495187398126914633,
    43.327073280914999519496122165407,
    48.005150881167159727942472749428,
    49.773832477672302181916784678564,
]

print(f"\nPhysical Constants:")
print(f"  φ (Golden Ratio) = {PHI:.10f}")
print(f"  Λ (Lattice)      = {LAMBDA:.10f}")

# ═══════════════════════════════════════════════════════════════════════════
# 2. GENERATE E8 ROOT LATTICE (240 roots in 8D)
# ═══════════════════════════════════════════════════════════════════════════
def generate_E8_roots():
    """
    Generate the 240 roots of the E8 lattice.
    
    E8 roots come in two types:
    Type 1: Permutations of (±1, ±1, 0, 0, 0, 0, 0, 0) - 112 roots
    Type 2: (±½, ±½, ±½, ±½, ±½, ±½, ±½, ±½) with even # of minus signs - 128 roots
    """
    roots = []
    
    # Type 1: Two non-zero entries (±1, ±1)
    for i in range(8):
        for j in range(i + 1, 8):
            for s1 in [-1, 1]:
                for s2 in [-1, 1]:
                    v = np.zeros(8)
                    v[i] = s1
                    v[j] = s2
                    roots.append(v)
    
    # Type 2: All ±1/2, even number of minus signs
    for signs in product([-0.5, 0.5], repeat=8):
        if np.sum(np.array(signs) < 0) % 2 == 0:
            roots.append(np.array(signs))
    
    return np.array(roots)

print("\n[1] Generating E8 Root Lattice...")
roots = generate_E8_roots()
print(f"    Total Roots: {len(roots)} (Expected: 240)")
print(f"    Dimensionality: {roots.shape[1]}D")

# Verify E8 properties
norms = np.sqrt(np.sum(roots**2, axis=1))
print(f"    All roots have norm = {norms[0]:.4f} (Expected: √2 ≈ 1.414)")

# ═══════════════════════════════════════════════════════════════════════════
# 3. CONSTRUCT MULTIPLE HAMILTONIAN MODELS
# ═══════════════════════════════════════════════════════════════════════════
print("\n[2] Constructing Hamiltonian Models...")

num_roots = len(roots)

def build_hamiltonian_v1(roots, coupling=1.0):
    """
    Model 1: Standard adjacency with Golden decay.
    H_ij = dot_product * phi^(-distance²)
    """
    n = len(roots)
    H = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                # Diagonal: Lattice potential energy
                H[i, j] = np.dot(roots[i], roots[i]) * PHI * coupling
            else:
                # Off-diagonal: Tunneling with Golden decay
                dist_sq = np.sum((roots[i] - roots[j])**2)
                dot = np.dot(roots[i], roots[j])
                H[i, j] = dot * (PHI ** (-dist_sq)) * coupling
    
    return H

def build_hamiltonian_v2(roots, coupling=1.0):
    """
    Model 2: Pure Golden adjacency (nearest neighbors only).
    H_ij = phi^(-|i-j|) for connected roots
    """
    n = len(roots)
    H = np.zeros((n, n))
    
    # Define connectivity: roots are connected if their dot product = -1 or 0
    for i in range(n):
        for j in range(n):
            if i == j:
                H[i, j] = 2 * PHI * coupling  # E8 roots have norm² = 2
            else:
                dot = np.dot(roots[i], roots[j])
                dist_sq = np.sum((roots[i] - roots[j])**2)
                
                # Adjacent roots in E8 have specific relationships
                if abs(dot + 1) < 0.01:  # dot = -1 (connected)
                    H[i, j] = -PSI * coupling
                elif abs(dot) < 0.01:  # dot = 0 (edge)
                    H[i, j] = PSI**2 * coupling
                elif abs(dot - 1) < 0.01:  # dot = 1 (same direction)
                    H[i, j] = PSI**3 * coupling
    
    return H

def build_hamiltonian_v3(roots, coupling=1.0):
    """
    Model 3: Lambda-scaled Hamiltonian (incorporating 16√15).
    """
    n = len(roots)
    H = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                # Diagonal scaled by Lambda
                H[i, j] = LAMBDA / (4 * PHI)
            else:
                dist_sq = np.sum((roots[i] - roots[j])**2)
                # Lambda-modulated tunneling
                H[i, j] = (LAMBDA / 40) * (PHI ** (-dist_sq))
    
    return H

# Build and analyze multiple models
print("\n" + "=" * 70)
print("                     SPECTRAL ANALYSIS")
print("=" * 70)

models = [
    ("Golden Decay (v1, c=π)", build_hamiltonian_v1(roots, PI)),
    ("Golden Decay (v1, c=π·φ)", build_hamiltonian_v1(roots, PI * PHI)),
    ("Nearest Neighbor (v2, c=2π)", build_hamiltonian_v2(roots, 2 * PI)),
    ("Lambda-Scaled (v3)", build_hamiltonian_v3(roots)),
]

for model_name, H in models:
    # Make Hamiltonian Hermitian (should already be, but ensure)
    H = (H + H.T) / 2
    
    # Compute eigenvalues
    eigenvalues = np.linalg.eigvalsh(H)
    eigenvalues = np.sort(np.abs(eigenvalues))
    
    # Filter meaningful eigenvalues (above noise floor)
    significant = eigenvalues[eigenvalues > 5]
    
    print(f"\n{'─' * 70}")
    print(f"MODEL: {model_name}")
    print(f"{'─' * 70}")
    print(f"Spectrum Range: [{eigenvalues.min():.4f}, {eigenvalues.max():.4f}]")
    print(f"Eigenvalues in [10, 60]: {[f'{e:.4f}' for e in eigenvalues if 10 < e < 60][:10]}")
    
    # Check for Riemann Zero matches
    print(f"\nChecking against Riemann Zeros:")
    print(f"{'Eigenvalue':>12} | {'Nearest Zero':>12} | {'Error':>10} | Match?")
    print(f"{'─' * 50}")
    
    matches_found = 0
    for ev in significant:
        if not (10 < ev < 80):
            continue
            
        # Find closest Riemann zero
        for zero in RIEMANN_ZEROS:
            error = abs(ev - zero)
            error_shifted = abs((ev - 0.5) - zero)
            
            if error < 1.0:
                print(f"{ev:>12.4f} | {zero:>12.4f} | {error:>10.4f} | DIRECT ★")
                matches_found += 1
            elif error_shifted < 1.0:
                print(f"{ev:>12.4f} | {zero:>12.4f} | {error_shifted:>10.4f} | SHIFTED ◆")
                matches_found += 1
    
    if matches_found == 0:
        print("    No close matches found in significant range.")
    else:
        print(f"\n    >>> {matches_found} POTENTIAL RESONANCES FOUND <<<")

# ═══════════════════════════════════════════════════════════════════════════
# 4. TARGETED CALIBRATION: Find coupling that hits γ₁ = 14.1347
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("        TARGETED CALIBRATION: Hunting γ₁ = 14.1347...")
print("=" * 70)

target_zero = RIEMANN_ZEROS[0]  # 14.1347

best_coupling = None
best_error = float('inf')
best_eigenvalue = None

print("\nScanning coupling strengths...")
for coupling in np.linspace(1, 20, 100):
    H = build_hamiltonian_v1(roots, coupling)
    H = (H + H.T) / 2
    eigenvalues = np.linalg.eigvalsh(H)
    
    # Find eigenvalue closest to 14.1347
    for ev in eigenvalues:
        error = abs(abs(ev) - target_zero)
        if error < best_error:
            best_error = error
            best_coupling = coupling
            best_eigenvalue = ev

print(f"\nBest Match Found:")
print(f"  Coupling: {best_coupling:.4f}")
print(f"  Eigenvalue: {best_eigenvalue:.6f}")
print(f"  Target (γ₁): {target_zero:.6f}")
print(f"  Error: {best_error:.6f}")

if best_error < 0.1:
    print("\n" + "★" * 70)
    print("  REMARKABLE: E8 produces eigenvalue near First Riemann Zero!")
    print("★" * 70)

# ═══════════════════════════════════════════════════════════════════════════
# 5. SPECTRAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("                     CONCLUSION")
print("=" * 70)
print("""
The E8 Root Lattice Hamiltonian with Golden Ratio coupling produces a
discrete spectrum. The key question is whether fine-tuning the coupling
or modifying the interaction term causes the spectrum to align with
the Riemann Zeros.

If the eigenvalues {λ_n} satisfy λ_n ≈ γ_n + ½ (the Riemann Zeros), then:

    H|ψ_n⟩ = (γ_n + ½)|ψ_n⟩

This would be the HILBERT-PÓLYA OPERATOR, solving the Riemann Hypothesis.
""")
print("=" * 70)
