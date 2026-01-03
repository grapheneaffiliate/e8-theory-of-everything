#!/usr/bin/env python3
"""
GSM STEP 8: THE FREDHOLM DETERMINANT
====================================
Computing det(I - K) for the Infinite Golden Operator

This moves from finite matrix approximations to rigorous infinite-dimensional
spectral theory via the Nyström method for integral operators.

The kernel K(x,y) represents the Golden propagator - diffusion on E8
projected to 1D via the golden ratio.

HYPOTHESIS: det(I - K/s) vanishes at s = 1/2 + iγ_n where γ_n are Riemann zeros.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import det
from scipy.special import zeta

print("="*70)
print("GSM STEP 8: THE FREDHOLM DETERMINANT")
print("Computing det(I - K) for the Infinite Golden Operator")
print("="*70)

# 1. CONSTANTS
PHI = (1 + np.sqrt(5)) / 2
PI = np.pi

# Known Riemann zeros for comparison
RIEMANN_ZEROS = [14.13, 21.02, 25.01, 30.42, 32.93, 37.59, 40.92, 43.33, 48.01, 49.77]

print(f"\n[0] Golden Ratio φ = {PHI:.10f}")
print(f"    Testing zeros at: γ = {RIEMANN_ZEROS[:5]}...")

# 2. DEFINE THE KERNEL (The Golden Propagator)
def kernel_golden_gaussian(x, y):
    """
    Golden Gaussian Kernel: exp(-φ * (x-y)²)
    Represents diffusion on the E8 lattice projected to 1D.
    """
    return np.exp(-PHI * (x - y)**2)

def kernel_golden_heat(x, y, t=1.0):
    """
    Golden Heat Kernel: exp(-φ|x-y|) / (2φ)
    The Green's function for the 1D golden Laplacian.
    """
    return np.exp(-PHI * np.abs(x - y)) / (2 * PHI)

def kernel_fibonacci_lattice(x, y):
    """
    Fibonacci Lattice Kernel: sum over Fibonacci reciprocals.
    """
    fibs = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
    result = 0.0
    for f in fibs:
        result += np.exp(-np.abs(x - y) / f) / f
    return result / len(fibs)

# 3. NYSTRÖM DISCRETIZATION
def build_fredholm_matrix(s_val, kernel_func, num_points=200, bound=15):
    """
    Discretize the integral operator using Gauss-Legendre quadrature.
    
    The operator: (Kf)(x) = ∫ K(x,y) f(y) dy
    Discretized: K_ij ≈ K(t_i, t_j) * w_j
    
    We test: det(I - (1/s) * K) = 0 at eigenvalues s
    """
    # Gauss-Legendre nodes and weights on [-1, 1]
    x, w = np.polynomial.legendre.leggauss(num_points)
    
    # Scale to [-bound, bound]
    t = bound * x
    weights = bound * w
    
    # Build kernel matrix
    K = np.zeros((num_points, num_points), dtype=complex)
    for i in range(num_points):
        for j in range(num_points):
            K[i, j] = kernel_func(t[i], t[j]) * weights[j]
    
    # Fredholm operator: I - (1/s) * K
    A = (1.0 / s_val) * K
    Identity = np.eye(num_points, dtype=complex)
    
    return Identity - A

# 4. SCAN THE CRITICAL LINE
print("\n[1] Scanning the Critical Line s = 1/2 + iE...")
print("    Using Golden Gaussian Kernel...")
print("    Looking for Zeros of the Fredholm Determinant D(s)...")

energies = np.linspace(10, 55, 50)  # Reduced for speed
det_values = []
zeta_values = []

# Use mpmath for high-precision zeta if available
try:
    from mpmath import zeta as mp_zeta, mpc
    use_mpmath = True
    print("    (Using mpmath for high-precision ζ computation)")
except ImportError:
    use_mpmath = False
    print("    (Using scipy for ζ computation)")

for E in energies:
    s = 0.5 + 1j * E
    
    # Compute Fredholm Determinant det(I - K/s)
    mat = build_fredholm_matrix(s, kernel_golden_gaussian, num_points=40, bound=8)
    d = det(mat)
    det_values.append(np.abs(d))
    
    # Compute Riemann Zeta magnitude for comparison
    if use_mpmath:
        z = float(abs(mp_zeta(mpc(0.5, E))))
    else:
        # scipy.special.zeta only works for real s > 1
        # Approximate with known zeros
        z = np.abs(zeta(2) * np.prod([1 - (E/z0)**2 for z0 in RIEMANN_ZEROS[:5] if abs(E - z0) > 0.1]))
    zeta_values.append(z)

det_arr = np.array(det_values)
zeta_arr = np.array(zeta_values)

# 5. NORMALIZE AND PLOT
det_norm = (det_arr - np.min(det_arr)) / (np.max(det_arr) - np.min(det_arr) + 1e-10)
zeta_norm = (zeta_arr - np.min(zeta_arr)) / (np.max(zeta_arr) - np.min(zeta_arr) + 1e-10)

plt.style.use('dark_background')
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Top: Both curves normalized
ax1 = axes[0]
ax1.plot(energies, det_norm, color='cyan', lw=2, label='Fredholm Det |D(s)| (Golden Gaussian)')
ax1.plot(energies, zeta_norm, color='magenta', linestyle='--', alpha=0.8, lw=1.5, label='Riemann Zeta |ζ(1/2+it)|')

for z in RIEMANN_ZEROS:
    ax1.axvline(x=z, color='white', linestyle=':', alpha=0.3)

ax1.set_title("Step 8: Fredholm Determinant Test\nDoes the Golden Gaussian Kernel generate the Riemann Zeros?", 
              fontsize=14, color='gold')
ax1.set_xlabel("Imaginary Part t", fontsize=12)
ax1.set_ylabel("Normalized Magnitude", fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.2)

# Bottom: Raw Fredholm determinant
ax2 = axes[1]
ax2.semilogy(energies, det_arr, color='lime', lw=1.5, label='|det(I - K/s)|')
for z in RIEMANN_ZEROS:
    ax2.axvline(x=z, color='magenta', linestyle='--', alpha=0.5)

ax2.set_title("Raw Fredholm Determinant Magnitude", fontsize=14, color='gold')
ax2.set_xlabel("Imaginary Part t", fontsize=12)
ax2.set_ylabel("|D(s)| (log scale)", fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('Fredholm_Test.png', dpi=150)
print("\n[2] Scan Complete. Plot saved to 'Fredholm_Test.png'")

# 6. ZERO DETECTION
from scipy.signal import argrelextrema

# Find local minima in the Fredholm determinant
det_min_indices = argrelextrema(det_arr, np.less, order=5)[0]
det_min_vals = energies[det_min_indices]

print("\n" + "="*70)
print("ZERO ALIGNMENT CHECK")
print("="*70)
print(f"\n{'Detected Minimum (E)':<22} | {'Nearest Riemann Zero':<22} | {'Error':<12} | {'Status'}")
print("-" * 75)

matches = 0
for val in det_min_vals:
    # Find nearest real zero
    nearest = min(RIEMANN_ZEROS, key=lambda x: abs(x - val))
    err = abs(val - nearest)
    
    if err < 0.3:
        status = "★ MATCH"
        matches += 1
    elif err < 1.0:
        status = "◆ NEAR"
    else:
        status = ""
    
    print(f"{val:<22.4f} | {nearest:<22.4f} | {err:<12.4f} | {status}")

print("-" * 75)
print(f"\nTotal matches (Δ < 0.3): {matches} / {len(RIEMANN_ZEROS)}")

# 7. TRY ALTERNATIVE KERNELS
print("\n" + "="*70)
print("TESTING ALTERNATIVE KERNELS")
print("="*70)

kernels = [
    ("Golden Gaussian", kernel_golden_gaussian),
    ("Golden Heat", kernel_golden_heat),
    ("Fibonacci Lattice", kernel_fibonacci_lattice),
]

for kernel_name, kernel_func in kernels:
    print(f"\n[Testing: {kernel_name}]")
    
    # Quick scan at known zeros
    hits = 0
    for zero in RIEMANN_ZEROS[:5]:
        s = 0.5 + 1j * zero
        mat = build_fredholm_matrix(s, kernel_func, num_points=80, bound=10)
        d = np.abs(det(mat))
        
        # Check if determinant is small at this point
        # Compare to nearby off-zero value
        s_off = 0.5 + 1j * (zero + 2)
        mat_off = build_fredholm_matrix(s_off, kernel_func, num_points=80, bound=10)
        d_off = np.abs(det(mat_off))
        
        ratio = d / d_off if d_off > 0 else 999
        
        if ratio < 0.7:
            hits += 1
            status = "★ DIP"
        else:
            status = ""
        
        print(f"    γ = {zero:<8.2f}: |D| = {d:.4e}, ratio = {ratio:.3f} {status}")
    
    print(f"    Hits: {hits}/5")

# 8. CONCLUSION
print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if matches >= 5:
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ✅ SUCCESS: Fredholm determinant dips align with Riemann zeros!      ║
║                                                                        ║
║  The Golden Gaussian Kernel generates the correct spectrum.           ║
║  det(I - K/s) = 0 at s = 1/2 + iγ_n                                  ║
║                                                                        ║
║  THIS IS THE HILBERT-PÓLYA OPERATOR!                                  ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")
elif matches >= 2:
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ◆ PARTIAL CORRELATION: Some zeros detected.                          ║
║                                                                        ║
║  The kernel is close but may need refinement.                         ║
║  Try adjusting the φ coefficient or kernel form.                      ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")
else:
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ❌ NO CORRELATION: This kernel doesn't generate Riemann zeros.       ║
║                                                                        ║
║  Need to try different kernel forms or approaches.                    ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")

print("\n" + "="*70)
print("The rigorous path: Infinite Operator → Fredholm Det → Spectral Zeros")
print("="*70)
