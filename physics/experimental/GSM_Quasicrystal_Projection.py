#!/usr/bin/env python3
"""
GSM STEP 7: E8 QUASICRYSTAL PROJECTION
======================================
Testing if Riemann Zeros are the diffraction pattern of a Golden Slice of E8

HYPOTHESIS: The Riemann zeros are the "spectral shadow" of E8 - 
the diffraction peaks when the 8D lattice is cut at a golden angle.

This moves from Quantum Mechanics (Hamiltonian eigenvalues) to 
Spectral Geometry (Fourier transform of point sets).
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import product

print("="*70)
print("GSM STEP 7: E8 QUASICRYSTAL PROJECTION")
print("Testing if Riemann Zeros are the diffraction pattern of a Golden Slice of E8")
print("="*70)

# 1. CONSTANTS & SETUP
PHI = (1 + np.sqrt(5)) / 2
PI = np.pi

# Known Riemann Zeros for verification (Imaginary parts)
RIEMANN_ZEROS = [
    14.1347, 21.0220, 25.0108, 30.4248, 32.9350,
    37.5861, 40.9187, 43.3270, 48.0051, 49.7738
]

# 2. GENERATE E8 ROOTS
def get_E8_roots():
    roots = []
    # Type 1: ±e_i ± e_j
    for i in range(8):
        for j in range(i + 1, 8):
            for s1 in [-1, 1]:
                for s2 in [-1, 1]:
                    v = np.zeros(8)
                    v[i] = s1
                    v[j] = s2
                    roots.append(v)
    # Type 2: (±1/2, ..., ±1/2) with even number of minus signs
    for signs in product([-0.5, 0.5], repeat=8):
        if np.sum(np.array(signs) < 0) % 2 == 0:
            roots.append(np.array(signs))
    return np.array(roots)

roots = get_E8_roots()
print(f"\n[1] Generated {len(roots)} E8 Roots in 8D.")

# 3. DEFINE THE "GOLDEN SLICE" PROJECTION VECTOR
# We project onto a line defined by powers of Phi to ensure irrationality
# V_proj = (1, phi, phi^2, ... phi^7) normalized
v_proj = np.array([PHI**i for i in range(8)])
v_proj = v_proj / np.linalg.norm(v_proj)

print(f"[2] Defined Golden Projection Vector:")
print(f"    V = (1, φ, φ², ..., φ⁷) normalized")
print(f"    First 3 components: {v_proj[:3]}")

# 4. PERFORM THE CUT-AND-PROJECT
# Map 8D roots to 1D points: x_i = root_i . v_proj
points_1d = np.array([np.dot(r, v_proj) for r in roots])

# The points tend to cluster. Let's define the quasicrystal locations.
# We look at the absolute values of the projected points.
qc_points = np.abs(points_1d)
# Filter out zeros and sort
qc_points = np.sort(qc_points[qc_points > 1e-5])

print(f"[3] Projected to 1D Quasicrystal. Found {len(qc_points)} distinct sites.")
print(f"    Point range: [{qc_points.min():.4f}, {qc_points.max():.4f}]")

# 5. CALCULATE STRUCTURE FACTOR (DIFFRACTION PATTERN)
# S(k) = | Sum_j exp(i * k * x_j) |^2
# This simulates firing X-rays at the 1D points and looking for bright spots.
print("\n[4] Calculating Diffraction Pattern (Structure Factor)...")

# Scan range in k-space corresponding to Riemann Zeros
k_scan = np.linspace(10, 55, 2000)
structure_factor = np.zeros_like(k_scan)

for i, k in enumerate(k_scan):
    # Sum of phases from all points
    phase_sum = np.sum(np.exp(1j * k * qc_points))
    structure_factor[i] = np.abs(phase_sum)**2

# Normalize factor for plotting
structure_factor = structure_factor / np.max(structure_factor)

# 6. ANALYSIS & PLOTTING
print("[5] Analyzing Peaks against Riemann Zeros...")

plt.style.use('dark_background')
plt.figure(figsize=(14, 7))

# Plot Diffraction Pattern
plt.plot(k_scan, structure_factor, color='cyan', lw=1.5, label='E8 Quasicrystal Diffraction S(k)')

# Overlay known Riemann Zeros
for i, z in enumerate(RIEMANN_ZEROS):
    label = 'Riemann Zeros γₙ' if i == 0 else None
    plt.axvline(x=z, color='magenta', linestyle='--', alpha=0.7, lw=1.5, label=label)

plt.title("The Shadow of E8: Quasicrystal Diffraction vs Riemann Zeros\n"
          "Do the Cyan Peaks align with the Magenta Lines?", fontsize=14, color='gold')
plt.xlabel("Wavevector k (Imaginary part of s)", fontsize=12)
plt.ylabel("Normalized Intensity |S(k)|²", fontsize=12)
plt.legend(loc='upper right')
plt.grid(True, alpha=0.2)
plt.xlim(10, 55)

plt.tight_layout()
plt.savefig('E8_Diffraction_Test.png', dpi=150)
print("    -> Plot saved to 'E8_Diffraction_Test.png'")

# 7. PEAK DETECTION AND COMPARISON
from scipy.signal import find_peaks

peaks, properties = find_peaks(structure_factor, height=0.1, prominence=0.05)
peak_k_vals = k_scan[peaks]
peak_heights = structure_factor[peaks]

# Sort by height (brightest first)
sort_idx = np.argsort(peak_heights)[::-1]
peak_k_vals_sorted = peak_k_vals[sort_idx]

print(f"\n[6] Found {len(peaks)} diffraction peaks above threshold.")
print(f"    Brightest 10 peaks at k =")
for i, k in enumerate(peak_k_vals_sorted[:10]):
    print(f"        {i+1}. k = {k:.4f}")

print("\n" + "="*70)
print("PEAK ALIGNMENT CHECK")
print("="*70)
print(f"\n{'Riemann Zero (γ)':<18} | {'Closest Peak (k)':<18} | {'Difference':<12} | {'Status'}")
print("-" * 70)

matches = 0
near_matches = 0
alignment_errors = []

for zero in RIEMANN_ZEROS:
    if len(peak_k_vals) > 0:
        # Find nearest peak
        idx = (np.abs(peak_k_vals - zero)).argmin()
        closest_peak = peak_k_vals[idx]
        diff = abs(closest_peak - zero)
        alignment_errors.append(diff)
        
        if diff < 0.5:
            mark = "★ MATCH"
            matches += 1
        elif diff < 1.5:
            mark = "◆ NEAR"
            near_matches += 1
        else:
            mark = ""
            
        print(f"{zero:<18.4f} | {closest_peak:<18.4f} | {diff:<12.4f} | {mark}")
    else:
        print(f"{zero:<18.4f} | {'No peaks found':<18} | {'N/A':<12} |")

print("-" * 70)
print(f"\nSTATISTICS:")
print(f"    Strong matches (Δ < 0.5): {matches}/{len(RIEMANN_ZEROS)}")
print(f"    Near matches (Δ < 1.5):   {near_matches}/{len(RIEMANN_ZEROS)}")
if alignment_errors:
    print(f"    Mean alignment error:     {np.mean(alignment_errors):.4f}")
    print(f"    Min alignment error:      {np.min(alignment_errors):.4f}")
    print(f"    Max alignment error:      {np.max(alignment_errors):.4f}")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if matches >= 7:
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ✅ STRONG SUCCESS: Diffraction pattern strongly correlates with       ║
║     Riemann Zeros!                                                     ║
║                                                                        ║
║  INTERPRETATION:                                                       ║
║  The Riemann zeros ARE the diffraction peaks of the E8 lattice        ║
║  projected onto the critical line via the Golden Ratio.               ║
║                                                                        ║
║  This supports the hypothesis:                                         ║
║  ξ(s) ~ Structure Factor of E8 Quasicrystal                           ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")
elif matches >= 5:
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ⚠️ MODERATE SUCCESS: Partial correlation detected.                   ║
║                                                                        ║
║  The projection vector or quasicrystal construction may need tuning.  ║
║  Try different golden-based projection vectors or include more        ║
║  lattice points from E8 shells.                                       ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")
elif matches >= 2:
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ◆ WEAK CORRELATION: Some matches found.                              ║
║                                                                        ║
║  The basic structure is there but needs refinement:                   ║
║  - Try including multiple E8 shells (norms 2, 4, 6, ...)             ║
║  - Use a different projection onto H4 first, then to 1D              ║
║  - Consider the DUAL lattice E8*                                     ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")
else:
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ❌ NO CORRELATION: The simple golden projection doesn't work.        ║
║                                                                        ║
║  This doesn't disprove the hypothesis - it means:                     ║
║  - The projection needs to be more sophisticated                       ║
║  - Try: E8 → H4 → Penrose 1D (two-step projection)                   ║
║  - Or: Use the VORONOI cells of E8 instead of root points            ║
║  - Or: Include the imaginary part via analytic continuation          ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")

# 8. ADDITIONAL TEST: Try multiple projection vectors
print("\n" + "="*70)
print("BONUS: Testing Alternative Projection Vectors")
print("="*70)

def test_projection(v_proj, name):
    """Test a projection vector and return match count."""
    v_norm = v_proj / np.linalg.norm(v_proj)
    pts = np.abs(np.array([np.dot(r, v_norm) for r in roots]))
    pts = np.sort(pts[pts > 1e-5])
    
    # Structure factor
    sf = np.zeros_like(k_scan)
    for i, k in enumerate(k_scan):
        sf[i] = np.abs(np.sum(np.exp(1j * k * pts)))**2
    sf = sf / np.max(sf)
    
    # Find peaks
    pks, _ = find_peaks(sf, height=0.1)
    pk_vals = k_scan[pks]
    
    # Count matches
    m = 0
    for z in RIEMANN_ZEROS:
        if len(pk_vals) > 0:
            if np.min(np.abs(pk_vals - z)) < 0.5:
                m += 1
    return m

# Try various golden-based projections
projections = [
    ("φ powers (1, φ, φ², ...)", np.array([PHI**i for i in range(8)])),
    ("φ⁻¹ powers", np.array([PHI**(-i) for i in range(8)])),
    ("Alternating φ", np.array([PHI**((-1)**i * i) for i in range(8)])),
    ("Fibonacci weights", np.array([1, 1, 2, 3, 5, 8, 13, 21])),
    ("φ-mixed", np.array([PHI, 1, PHI, 1, PHI, 1, PHI, 1])),
    ("Random (control)", np.random.randn(8)),
]

print(f"\n{'Projection Vector':<25} | {'Matches':<10}")
print("-" * 40)
for name, v in projections:
    m = test_projection(v, name)
    marker = "★" if m >= 5 else ("◆" if m >= 3 else "")
    print(f"{name:<25} | {m:<10} {marker}")

print("\n" + "="*70)
print("TEST COMPLETE")
print("="*70)
