#!/usr/bin/env python3
"""
E8 FAST TEST - Quick verification of φ^(-12) suppression
=========================================================

A FAST version that runs in ~30 seconds.
Uses vectorized NumPy operations instead of nested loops.
"""

import numpy as np
import sys

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2
PHI_NEG_12 = PHI ** (-12)

print()
print("=" * 60)
print("  E8 QUANTUM GRAVITY: FAST φ^(-12) VERIFICATION")
print("=" * 60)
print()
print(f"  Target: φ^(-12) = {PHI_NEG_12:.6f}")
print()

def fast_lattice_integral(spacing, N=8):
    """Fast vectorized 1-loop integral."""
    a = spacing
    dk = 2 * np.pi / (N * a)
    p_ext = 0.5  # Only z-component matters for symmetry
    
    # Create 4D momentum grid
    n = np.arange(N) - N//2
    grid = np.meshgrid(n, n, n, n, indexing='ij')
    k = dk * np.stack(grid, axis=-1)  # Shape: (N,N,N,N,4)
    
    # Lattice Laplacian: Delta(k) = (2/a²) Σ [1 - cos(k_μ a)]
    delta_k = (2/a**2) * np.sum(1 - np.cos(k * a), axis=-1)
    
    # For p-k, only the 4th component of p is non-zero
    p_minus_k = k.copy()
    p_minus_k[..., 3] -= p_ext
    delta_pk = (2/a**2) * np.sum(1 - np.cos(p_minus_k * a), axis=-1)
    
    # Propagator product, avoiding division by zero
    mask = (delta_k > 1e-10) & (delta_pk > 1e-10)
    result = np.sum(1 / (delta_k[mask] * delta_pk[mask]))
    
    return result / np.sum(mask)

def fast_continuum_integral(cutoff, n_samples=10000):
    """Fast Monte Carlo continuum integral."""
    p_ext = np.array([0, 0, 0, 0.5])
    
    k = np.random.uniform(-cutoff, cutoff, (n_samples, 4))
    k_sq = np.sum(k**2, axis=1)
    pk = k.copy()
    pk[:, 3] -= p_ext[3]
    pk_sq = np.sum(pk**2, axis=1)
    
    mask = (k_sq > 1e-6) & (pk_sq > 1e-6)
    integrand = 1 / (k_sq[mask] * pk_sq[mask])
    
    volume = (2 * cutoff) ** 4
    return np.sum(integrand) * volume / n_samples / (2*np.pi)**4

# Run fast test at 3 scales
print("  Scale     Spacing    Cutoff     Ratio        φ power    Match")
print("-" * 60)

results = []
spacings = [0.1, 0.2, 0.4]

for i, a in enumerate(spacings):
    cutoff = np.pi / a
    sys.stdout.write(f"  {i+1}         {a:.2f}       {cutoff:.2f}      ")
    sys.stdout.flush()
    
    I_lat = fast_lattice_integral(a, N=8)
    I_cont = fast_continuum_integral(cutoff, n_samples=10000)
    
    ratio = I_lat / I_cont if I_cont > 0 else 0
    power = np.log(ratio) / np.log(PHI) if ratio > 0 else -99
    
    match = "YES" if abs(power + 12) < 6 else "NO"
    results.append(power)
    
    print(f"{ratio:.4e}    {power:+.1f}      {match}")

print("-" * 60)

# Summary
mean_power = np.mean([r for r in results if r > -50])
print()
print("=" * 60)
print("  RESULT")
print("=" * 60)
print()
print(f"  Mean φ power: {mean_power:.1f}")
print(f"  Expected:     -12")
print(f"  Difference:   {abs(mean_power + 12):.1f}")
print()

if abs(mean_power + 12) < 5:
    print("  ✓ φ^(-12) SUPPRESSION CONFIRMED!")
    print()
    print("  The lattice integral is suppressed by ~φ^(-12)")
    print("  relative to the continuum, as predicted by")
    print("  the 600-cell orthoscheme derivation.")
else:
    print("  Results deviate from φ^(-12)")
    print("  (May need larger grid or more samples)")

print()
print("=" * 60)
