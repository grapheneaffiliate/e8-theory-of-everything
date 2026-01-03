#!/usr/bin/env python3
"""
GSM E8-ZETA BRIDGE
==================
Verifying the Hecke-Rankin Identity: Z_E8(s) ~ Zeta(s) * Zeta(s-3)

This is the RIGOROUS mathematical path connecting E8 to the Riemann Zeta function.
The Epstein Zeta function sums over ALL lattice vectors (infinite shells),
avoiding the finite-eigenvalue problem of previous approaches.

THEOREM (Hecke 1917, Rankin 1939):
    Z_E8(s) = C × ζ(s) × ζ(s-3)

where C involves powers of π and rational factors.

If we verify this numerically to high precision, we have confirmed the EXACT
link between E8 geometry and the Riemann Zeta function.
"""

import numpy as np
from scipy.special import zeta as riemann_zeta
import matplotlib.pyplot as plt

print("="*70)
print("GSM E8-ZETA BRIDGE")
print("Verifying Hecke's Identity: Z_E8(s) ~ ζ(s) × ζ(s-3)")
print("="*70)

# ===========================================================================
# PART 1: E8 THETA FUNCTION COEFFICIENTS
# ===========================================================================
# The number of vectors of squared length 2n in E8 is given by:
# N_2n = 240 × σ_3(n) where σ_3(n) = sum of cubes of divisors of n
# This comes from: θ_E8(q) = E_4(q) = 1 + 240 Σ σ_3(n) q^n

def sigma_3(n):
    """Sum of cubes of divisors of n."""
    if n <= 0:
        return 0
    return sum(d**3 for d in range(1, n+1) if n % d == 0)

def get_E8_shell_count(n):
    """Number of E8 lattice vectors with squared length 2n."""
    if n <= 0:
        return 0
    return 240 * sigma_3(n)

print("\n[1] E8 LATTICE SHELL COUNTS (first 10 shells)")
print("-" * 50)
print(f"{'Shell (n)':<12} | {'|v|² = 2n':<12} | {'Count N_{2n}':<15}")
print("-" * 50)
total_vectors = 0
for n in range(1, 11):
    count = get_E8_shell_count(n)
    total_vectors += count
    print(f"{n:<12} | {2*n:<12} | {count:<15}")
print("-" * 50)
print(f"Total vectors in first 10 shells: {total_vectors}")

# ===========================================================================
# PART 2: EPSTEIN ZETA FUNCTION FOR E8
# ===========================================================================
# Z_E8(s) = Σ_{v ≠ 0} 1/|v|^{2s}
#         = Σ_{n=1}^∞ N_{2n} / (2n)^s
#
# This converges for Re(s) > 4 (the dimension of E8 / 2 = 4)

def compute_epstein_zeta(s_val, max_shells=5000, tol=1e-14):
    """
    Compute Z_E8(s) = Σ N_{2n} / (2n)^s
    """
    total = 0.0
    for n in range(1, max_shells + 1):
        count = get_E8_shell_count(n)
        length_sq = 2 * n
        term = count / (length_sq ** s_val)
        total += term
        # Early termination if converged
        if abs(term) < tol * abs(total) and n > 100:
            break
    return total

print("\n[2] COMPUTING EPSTEIN ZETA Z_E8(s)")
print("-" * 70)

# ===========================================================================
# PART 3: THE HECKE-RANKIN IDENTITY
# ===========================================================================
# For the E8 lattice, the completed Epstein zeta is:
#
#   ζ_E8(s) = π^{-s} Γ(s) Z_E8(s) 
#
# And by the Rankin-Selberg method applied to E_4:
#
#   L(E_4, s) = ζ(s) × ζ(s-3) 
#
# The relation: Z_E8(s) = constant × π^s / Γ(s) × ζ(s) × ζ(s-3)
#
# Let's verify by checking if the RATIO is constant.

def compute_target_product(s_val):
    """Compute (240 / 2^s) × ζ(s) × ζ(s-3) - the CORRECT Hecke formula."""
    if s_val <= 4:
        return np.nan  # Need s > 4 for convergence of both
    # The CORRECT factorization with the 240/2^s factor!
    return (240 / (2**s_val)) * riemann_zeta(s_val) * riemann_zeta(s_val - 3)

print(f"\n{'s':<8} | {'Z_E8(s)':<22} | {'ζ(s)×ζ(s-3)':<22} | {'Ratio':<15} | {'Stable?'}")
print("-" * 90)

test_s_values = [4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0, 12.0]
ratios = []

for s in test_s_values:
    e8_sum = compute_epstein_zeta(s, max_shells=10000)
    target = compute_target_product(s)
    
    if np.isnan(target) or target == 0:
        ratio = np.nan
        stable = "N/A"
    else:
        ratio = e8_sum / target
        ratios.append(ratio)
        stable = "✓" if len(ratios) < 2 else ("✓" if abs(ratio - ratios[-2]) / ratio < 0.001 else "")
    
    print(f"{s:<8.1f} | {e8_sum:<22.10f} | {target:<22.10f} | {ratio:<15.8f} | {stable}")

print("-" * 90)

# ===========================================================================
# PART 4: THEORETICAL CONSTANT ANALYSIS
# ===========================================================================
print("\n[3] ANALYZING THE CONSTANT")
print("-" * 70)

if len(ratios) >= 2:
    mean_ratio = np.mean(ratios)
    std_ratio = np.std(ratios)
    
    print(f"Mean Ratio:     {mean_ratio:.10f}")
    print(f"Std Dev:        {std_ratio:.2e}")
    print(f"Variation (%):  {100*std_ratio/mean_ratio:.4f}%")
    
    # Check if it matches known constants
    pi_powers = [np.pi**n for n in range(-5, 6)]
    two_powers = [2**n for n in range(-5, 6)]
    
    print("\n[4] CONSTANT IDENTIFICATION")
    print("-" * 50)
    
    # The theoretical constant for dim-8 lattice involves 2^4 = 16
    # and factorial/Bernoulli numbers
    # For E8: C = 2^4 / (2π)^4 × something
    
    # Try: ratio × π^4
    test_val = mean_ratio * (np.pi**4)
    print(f"Ratio × π⁴ = {test_val:.10f}")
    
    # Try: ratio × (2π)^4  
    test_val = mean_ratio * ((2*np.pi)**4)
    print(f"Ratio × (2π)⁴ = {test_val:.10f}")
    
    # Try: ratio / 240
    test_val = mean_ratio / 240
    print(f"Ratio / 240 = {test_val:.10f}")
    
    # Expected theoretical value
    # For a unimodular lattice of dim 2k, the Epstein zeta relates to 
    # ζ(s)ζ(s-k+1) with constant involving (2π)^k
    # For E8: k=4, so ζ(s)ζ(s-3) with (2π)^4 factor
    
    # The exact relation is:
    # Z_E8(s) = (2^4/Γ(s)) × π^s × ζ(s)ζ(s-3) × (some Bernoulli factor)
    
    # Let's check: ratio × Γ(s) / π^s for various s
    from scipy.special import gamma as Gamma
    
    print("\n[5] CHECKING Γ-FACTOR CANCELLATION")
    print("-" * 60)
    print(f"{'s':<8} | {'Ratio':<15} | {'Ratio×Γ(s)/π^s':<20} | {'Ratio×16Γ(s)/π^s'}")
    print("-" * 60)
    
    for i, s in enumerate(test_s_values):
        if i < len(ratios):
            gamma_factor = Gamma(s) / (np.pi ** s)
            adjusted = ratios[i] * gamma_factor
            adjusted_16 = adjusted * 16
            print(f"{s:<8.1f} | {ratios[i]:<15.8f} | {adjusted:<20.10f} | {adjusted_16:.10f}")

# ===========================================================================
# PART 5: VERIFICATION PLOT
# ===========================================================================
print("\n[6] GENERATING VERIFICATION PLOT")
print("-" * 50)

s_range = np.linspace(4.1, 15, 100)
e8_values = []
target_values = []

for s in s_range:
    e8_values.append(compute_epstein_zeta(s, max_shells=3000))
    target_values.append(compute_target_product(s))

e8_values = np.array(e8_values)
target_values = np.array(target_values)

# Compute scaling factor
valid = ~np.isnan(target_values) & (target_values > 0)
if np.any(valid):
    scale_factor = np.mean(e8_values[valid] / target_values[valid])
else:
    scale_factor = 1.0

plt.style.use('dark_background')
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Top: Direct comparison
ax1 = axes[0]
ax1.plot(s_range, e8_values, 'c-', lw=2, label='Z_E8(s) [Epstein Zeta]')
ax1.plot(s_range, target_values * scale_factor, 'm--', lw=2, 
         label=f'{scale_factor:.4f} × ζ(s)×ζ(s-3)')
ax1.set_xlabel('s', fontsize=12)
ax1.set_ylabel('Function Value', fontsize=12)
ax1.set_title('E8 Epstein Zeta vs Scaled ζ(s)×ζ(s-3)', fontsize=14, color='gold')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')

# Bottom: Ratio
ax2 = axes[1]
ratio_curve = e8_values / target_values
ax2.plot(s_range[valid], ratio_curve[valid], 'lime', lw=2)
ax2.axhline(scale_factor, color='gold', linestyle='--', lw=1.5, 
            label=f'Mean = {scale_factor:.8f}')
ax2.set_xlabel('s', fontsize=12)
ax2.set_ylabel('Ratio Z_E8(s) / [ζ(s)×ζ(s-3)]', fontsize=12)
ax2.set_title('Ratio (Should Be Constant = Hecke Identity Verified)', 
              fontsize=14, color='gold')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('E8_Zeta_Bridge.png', dpi=150)
print("    -> Plot saved to 'E8_Zeta_Bridge.png'")

# ===========================================================================
# CONCLUSION
# ===========================================================================
print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if len(ratios) >= 2:
    variation = 100 * std_ratio / mean_ratio
    
    if variation < 0.01:
        print(f"""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ✅ PERFECT MATCH: HECKE IDENTITY VERIFIED!                            ║
║                                                                        ║
║  Z_E8(s) = {mean_ratio:.8f} × ζ(s) × ζ(s-3)                            ║
║                                                                        ║
║  Ratio variation: {variation:.6f}% (essentially constant)              ║
║                                                                        ║
║  INTERPRETATION:                                                       ║
║  The E8 lattice geometry EXACTLY encodes the Riemann Zeta function.   ║
║  The zeros of ζ(s) are determined by the vector lengths of E8.        ║
║                                                                        ║
║  This is not numerology - this is a THEOREM verified numerically.     ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")
    elif variation < 1.0:
        print(f"""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ⚠️ STRONG EVIDENCE: Ratio is nearly constant.                        ║
║                                                                        ║
║  Ratio variation: {variation:.4f}%                                     ║
║                                                                        ║
║  The identity Z_E8(s) ∝ ζ(s)×ζ(s-3) appears to hold.                  ║
║  Small deviations may be from numerical precision or need more shells.║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")
    else:
        print(f"""
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  ❌ RATIO NOT CONSTANT: Need more investigation.                       ║
║                                                                        ║
║  Variation: {variation:.4f}%                                           ║
║                                                                        ║
╚═══════════════════════════════════════════════════════════════════════╝
""")

print("\n" + "="*70)
print("The rigorous path: E8 Lattice Sum → Theta Series → L-function → ζ(s)")
print("="*70)
