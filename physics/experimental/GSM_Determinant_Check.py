#!/usr/bin/env python3
"""
GSM STEP 6: THE DETERMINANT GAP
================================
Testing Identity: det(H - E*I) ~ ξ(1/2 + iE)

If the dips (zeros) align, we have empirically closed the gap!
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import product

print("="*70)
print("GSM STEP 6: THE DETERMINANT GAP")
print("Testing Identity: det(H - E*I) ~ ζ(1/2 + iE)")
print("="*70)

# 1. CONSTANTS & E8 SETUP
PHI = (1 + np.sqrt(5)) / 2
LAMBDA_SCALING = 14.1347 / 7.52  # Match first zero

def get_E8_roots():
    roots = []
    # Type I: ±e_i ± e_j (i<j)
    for i in range(8):
        for j in range(i + 1, 8):
            for s1 in [-1, 1]:
                for s2 in [-1, 1]:
                    v = np.zeros(8)
                    v[i] = s1
                    v[j] = s2
                    roots.append(v)
    # Type II: (±1/2, ..., ±1/2) with even number of minus signs
    for signs in product([-0.5, 0.5], repeat=8):
        if np.sum(np.array(signs) < 0) % 2 == 0:
            roots.append(np.array(signs))
    return np.array(roots)

# 2. BUILD HAMILTONIAN (Golden Tunneling Model)
print("\n[1] Building E8 Hamiltonian with Golden Tunneling...")
roots = get_E8_roots()
n_roots = len(roots)
print(f"    E8 roots: {n_roots}")

H = np.zeros((n_roots, n_roots))
for i in range(n_roots):
    for j in range(n_roots):
        if i == j:
            H[i,j] = np.dot(roots[i], roots[i]) * PHI
        else:
            dist_sq = np.sum((roots[i] - roots[j])**2)
            dot = np.dot(roots[i], roots[j])
            H[i,j] = dot * (PHI ** (-dist_sq))

# Scale H to physical units
H_phys = H * LAMBDA_SCALING

# 3. DIAGONALIZE ONCE
print("[2] Diagonalizing H...")
eigenvalues = np.linalg.eigvalsh(H_phys)
evals = np.sort(eigenvalues)

print(f"    Eigenvalue range: [{evals.min():.4f}, {evals.max():.4f}]")
print(f"    First 5 eigenvalues: {evals[:5]}")

# 4. COMPUTE SPECTRAL DETERMINANT vs ZETA
print("\n[3] Computing Spectral Determinant vs Zeta comparison...")

# Try to use mpmath for zeta, fallback to scipy
try:
    from mpmath import zeta as mp_zeta, mpc
    print("    Using mpmath for Riemann zeta")
    use_mpmath = True
except ImportError:
    try:
        from scipy.special import zeta as scipy_zeta
        print("    Using scipy for Riemann zeta (less accurate at critical line)")
        use_mpmath = False
    except ImportError:
        print("    WARNING: No zeta function available, using placeholder")
        use_mpmath = None

# Riemann zeros for reference
riemann_zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151]

# Scan range covering first ~5 zeros
scan_range = np.linspace(10, 50, 500)
det_vals = []
zeta_vals = []

for E in scan_range:
    # Spectral Determinant: sum log|eigenvalue - E|
    # Dips occur when E is close to an eigenvalue
    spectral_val = np.sum(np.log(np.abs(evals - E) + 1e-50))
    det_vals.append(spectral_val)
    
    # Riemann Zeta at 1/2 + iE
    if use_mpmath:
        s = mpc(0.5, E)
        z = mp_zeta(s)
        zeta_vals.append(np.log(float(abs(z)) + 1e-50))
    elif use_mpmath is False:
        # scipy zeta doesn't work on critical line easily
        # use asymptotic approximation
        zeta_vals.append(-E * np.log(E) / (2 * np.pi) + E/2)  # rough approx
    else:
        zeta_vals.append(0)

det_vals = np.array(det_vals)
zeta_vals = np.array(zeta_vals)

# 5. NORMALIZE AND PLOT
print("[4] Normalizing and plotting...")

# Z-score normalization
det_norm = (det_vals - np.mean(det_vals)) / (np.std(det_vals) + 1e-10)
zeta_norm = (zeta_vals - np.mean(zeta_vals)) / (np.std(zeta_vals) + 1e-10)

# Create plot
plt.style.use('dark_background')
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Top: Overlay plot
ax1 = axes[0]
ax1.plot(scan_range, det_norm, color='cyan', label='Log |det(H - E)| (E8 Operator)', lw=2)
ax1.plot(scan_range, zeta_norm, color='magenta', linestyle='--', label='Log |ζ(1/2 + iE)|', lw=1.5, alpha=0.8)

# Mark Riemann zeros
for i, gamma in enumerate(riemann_zeros):
    if 10 <= gamma <= 50:
        ax1.axvline(gamma, color='gold', linestyle=':', alpha=0.5, lw=1)
        if i < 5:
            ax1.text(gamma, ax1.get_ylim()[1]*0.9, f'γ_{i+1}', color='gold', fontsize=8, ha='center')

ax1.set_title("Step 6: The Determinant Gap Test\nDo the Operator Zeros match the Zeta Zeros?", 
              fontsize=14, color='gold')
ax1.set_xlabel("Energy E (Imaginary part of s)", fontsize=12)
ax1.set_ylabel("Normalized Log Magnitude", fontsize=12)
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.2)

# Bottom: Eigenvalue distribution vs Riemann zeros
ax2 = axes[1]

# Eigenvalue histogram
ax2.hist(evals, bins=50, color='cyan', alpha=0.6, label='E8 Eigenvalues', density=True)

# Mark Riemann zeros
for i, gamma in enumerate(riemann_zeros[:5]):
    ax2.axvline(gamma, color='magenta', linestyle='--', lw=2, alpha=0.8)
    ax2.text(gamma, ax2.get_ylim()[1]*0.8 if ax2.get_ylim()[1] > 0 else 0.1, 
             f'γ_{i+1}={gamma:.2f}', color='magenta', fontsize=9, rotation=90, va='bottom')

ax2.set_title("Eigenvalue Distribution vs Riemann Zeros", fontsize=14, color='gold')
ax2.set_xlabel("Energy", fontsize=12)
ax2.set_ylabel("Density", fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('Determinant_Test.png', dpi=150, bbox_inches='tight')
print(f"\n[5] Plot saved to 'Determinant_Test.png'")

# 6. NUMERICAL ANALYSIS: Do dips align?
print("\n" + "="*70)
print("NUMERICAL ANALYSIS: Zero Alignment")
print("="*70)

# Find local minima in det_norm (the dips)
from scipy.signal import find_peaks

# Invert to find minima as peaks
det_minima_idx, _ = find_peaks(-det_norm, prominence=0.5)
det_minima_E = scan_range[det_minima_idx]

print(f"\nE8 Spectral Determinant dips (local minima) at E =")
for i, E in enumerate(det_minima_E[:10]):
    # Find closest Riemann zero
    closest_gamma = min(riemann_zeros, key=lambda g: abs(g - E))
    diff = E - closest_gamma
    print(f"    Dip {i+1}: E = {E:.4f} | Closest γ = {closest_gamma:.4f} | Δ = {diff:+.4f}")

# Find zeta minima for comparison
if use_mpmath:
    zeta_minima_idx, _ = find_peaks(-zeta_norm, prominence=0.5)
    zeta_minima_E = scan_range[zeta_minima_idx]
    print(f"\nZeta dips (should be at γ_n):")
    for i, E in enumerate(zeta_minima_E[:5]):
        closest_gamma = min(riemann_zeros, key=lambda g: abs(g - E))
        diff = E - closest_gamma
        print(f"    Dip {i+1}: E = {E:.4f} | γ = {closest_gamma:.4f} | Δ = {diff:+.4f}")

# 7. DIRECT EIGENVALUE vs ZERO COMPARISON
print("\n" + "="*70)
print("DIRECT COMPARISON: E8 Eigenvalues vs Riemann Zeros")
print("="*70)

# Scale eigenvalues to match first zero
scale = riemann_zeros[0] / evals[evals > 10][0] if any(evals > 10) else 1.0
print(f"\nScaling factor: {scale:.6f}")

print("\nScaled E8 Eigenvalues in range [10, 50]:")
evals_scaled = evals * scale
evals_in_range = evals_scaled[(evals_scaled > 10) & (evals_scaled < 50)]

for i, e in enumerate(evals_in_range[:10]):
    closest_gamma = min(riemann_zeros, key=lambda g: abs(g - e))
    diff = e - closest_gamma
    print(f"    λ_{i+1} (scaled) = {e:.4f} | Closest γ = {closest_gamma:.4f} | Δ = {diff:+.4f}")

# 8. CONCLUSION
print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

# Calculate alignment score
if len(det_minima_E) > 0:
    alignment_errors = []
    for E in det_minima_E[:5]:
        closest_gamma = min(riemann_zeros, key=lambda g: abs(g - E))
        alignment_errors.append(abs(E - closest_gamma))
    mean_error = np.mean(alignment_errors)
    print(f"\nMean alignment error (first 5 dips): {mean_error:.4f}")
    
    if mean_error < 1.0:
        print("\n✅ STRONG ALIGNMENT: Dips within 1 unit of Riemann zeros!")
        print("   This suggests det(H - E) ~ ζ(1/2 + iE) numerically.")
    elif mean_error < 3.0:
        print("\n⚠️ MODERATE ALIGNMENT: Some correspondence but not exact.")
    else:
        print("\n❌ WEAK ALIGNMENT: E8 eigenvalues don't match Riemann zeros.")
else:
    print("\nNo clear dips found in spectral determinant.")

print("\nFINAL STATUS:")
print("- If dips align: Numerical evidence for det(H-E) ~ ξ(s)")
print("- Exact identity still requires ANALYTIC proof")
print("="*70)
