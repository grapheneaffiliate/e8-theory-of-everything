#!/usr/bin/env python3
"""
RH SPECIALIZED DETECTOR
=======================

Specialized test functions to detect off-line zeros:
1. Bandlimited (sinc²) with support avoiding actual zeros
2. Polynomially-weighted Gaussians
3. Localized functions tuned to off-line positions

Goal: Find g with ĝ ≥ 0 where Z(g) < 0 with off-line zero.
"""

from mpmath import mp, mpf, mpc, pi, log, exp, zeta, sqrt, cos, sin
from mpmath import gamma as mpgamma, zetazero
import math

mp.dps = 50

print("="*70)
print("RH SPECIALIZED DETECTOR")
print("Advanced Test Functions for Off-Line Detection")
print("="*70)

# =============================================================================
# 1. BANDLIMITED SINC² FUNCTION
# =============================================================================

print("\n" + "="*70)
print("[1] BANDLIMITED SINC² FUNCTION")
print("="*70)

print("""
    Sinc² test function:
    
    g(x) = (sin(σx) / (σx))²
    ĝ(u) = (π/σ) × max(0, 1 - |u|/σ)  [triangular, support [-σ, σ]]
    
    KEY PROPERTY: ĝ(u) = 0 for |u| > σ
    
    If we choose σ < γ₁ ≈ 14.13, then ALL on-line zeros contribute 0!
    
    Only the off-line zero (if it exists) would contribute,
    and that contribution could be NEGATIVE.
""")

def sinc_squared_g_hat(u, sigma):
    """ĝ(u) for sinc² - triangular with support [-σ, σ]"""
    u_abs = abs(u)
    if isinstance(u, mpc):
        u_abs = mp.sqrt(mp.re(u)**2 + mp.im(u)**2)
    
    if u_abs > sigma:
        return mpf(0)
    return (mp.pi / sigma) * (1 - u_abs / sigma)

def sinc_squared_g_hat_complex(u_real, u_imag, sigma):
    """
    ĝ for complex argument u = u_real + i*u_imag.
    
    For sinc², ĝ(u) = (π/σ)(1 - |u|/σ) for |u| ≤ σ.
    
    But |u| for complex u is √(u_real² + u_imag²).
    """
    sigma = mpf(sigma)
    u = mpc(u_real, u_imag)
    
    # |u| = √(Re² + Im²)
    u_abs = mp.sqrt(u_real**2 + u_imag**2)
    
    if u_abs > sigma:
        return mpc(0, 0)
    
    # The triangular function at complex argument
    result = (mp.pi / sigma) * (1 - u_abs / sigma)
    return mpc(result, 0)  # Still real for real |u|

# Test with bandlimited support that excludes all on-line zeros
print("\n  Testing sinc² with σ < γ₁:")
print()

gamma_1 = mp.im(zetazero(1))  # ≈ 14.13
print(f"  First actual zero: γ₁ ≈ {mp.nstr(gamma_1, 6)}")
print()

for sigma in [5, 10, 13, 14, 15]:
    sigma = mpf(sigma)
    
    # On-line contribution (should be 0 for σ < γ₁)
    oneline_contrib = sinc_squared_g_hat(gamma_1, sigma)
    
    # Off-line contribution (σ_off = 0.4, γ_off = 14.13)
    sigma_off = mpf("0.4")
    gamma_off = gamma_1
    delta = sigma_off - mpf("0.5")
    
    # Complex argument: u = γ - iδ
    u_real = float(gamma_off)
    u_imag = float(-delta)
    
    offine = sinc_squared_g_hat_complex(u_real, u_imag, sigma)
    
    print(f"  σ = {sigma}: On-line = {mp.nstr(oneline_contrib, 6)}, "
          f"Off-line = {mp.nstr(mp.re(offine), 6)}")

# =============================================================================
# 2. LOCALIZED BUMP FUNCTION
# =============================================================================

print("\n" + "="*70)
print("[2] LOCALIZED BUMP FUNCTION")
print("="*70)

print("""
    Strategy: Center test function at the off-line zero location.
    
    g(x) = exp(-a(x - x₀)²) centered at x₀
    
    If x₀ is chosen to resonate with the off-line location,
    it may amplify the negative contribution.
    
    ĝ(u) = √(π/a) exp(-π²u²/a) × exp(-2πi u x₀)
         = √(π/a) exp(-π²u²/a) × [cos(2πux₀) - i sin(2πux₀)]
""")

def shifted_gaussian_g_hat(u, a, x0):
    """
    ĝ for shifted Gaussian g(x) = exp(-a(x-x₀)²).
    
    ĝ(u) = √(π/a) exp(-π²u²/a) × exp(-2πi u x₀)
    """
    a = mpf(a)
    x0 = mpf(x0)
    
    if isinstance(u, (int, float)):
        u = mpf(u)
    
    base = mp.sqrt(mp.pi / a) * mp.exp(-mp.pi**2 * u**2 / a)
    phase = mp.exp(-2j * mp.pi * u * x0)
    
    return base * phase

# Test shifted Gaussian
print("\n  Shifted Gaussian centered at x₀:")
print()

a = mpf("100")
sigma_off = mpf("0.4")
gamma_off = mp.im(zetazero(1))
delta = sigma_off - mpf("0.5")

for x0 in [0, 0.5, 1.0, 2.0]:
    x0 = mpf(x0)
    
    # On-line contribution
    oneline = shifted_gaussian_g_hat(gamma_off, a, x0)
    
    # Off-line contribution (complex argument)
    u_complex = mpc(gamma_off, -delta)
    offine = shifted_gaussian_g_hat(u_complex, a, x0)
    
    print(f"  x₀ = {x0}: On-line Re = {mp.nstr(mp.re(oneline), 6)}, "
          f"Off-line Re = {mp.nstr(mp.re(offine), 6)}")

# =============================================================================
# 3. POLYNOMIAL-WEIGHTED TEST FUNCTION
# =============================================================================

print("\n" + "="*70)
print("[3] POLYNOMIAL-WEIGHTED TEST FUNCTION")
print("="*70)

print("""
    Idea: Weight by polynomial that's zero at on-line locations.
    
    If we could construct g such that ĝ(γₖ) = 0 for actual zeros
    but ĝ(γ_off) ≠ 0 for off-line position, we'd isolate the off-line.
    
    Example: ĝ(u) = |P(u)|² × Gaussian
    where P(u) has roots at the on-line zeros.
    
    This ensures ĝ ≥ 0 (squared polynomial × positive).
""")

def polynomial_zero_at_gamma(u, gamma_list):
    """
    P(u) = Π (u - γₖ) for γₖ in gamma_list.
    Returns P(u).
    """
    result = mpc(1, 0)
    for gamma_k in gamma_list:
        result *= (u - gamma_k)
    return result

def weighted_test_function(u, a, gamma_list):
    """
    ĝ(u) = |P(u)|² × Gaussian(u)
    
    where P(u) = Π(u - γₖ) has roots at the on-line zeros.
    """
    a = mpf(a)
    
    # Gaussian factor
    gaussian = mp.sqrt(mp.pi / a) * mp.exp(-mp.pi**2 * u**2 / a)
    
    # Polynomial factor
    P_u = polynomial_zero_at_gamma(u, gamma_list)
    P_squared = abs(P_u)**2
    
    return gaussian * P_squared

# Test polynomial-weighted function
print("\n  Polynomial-weighted (roots at first 5 zeros):")
print()

# Get first 5 zeros
zeros_5 = [mp.im(zetazero(k)) for k in range(1, 6)]
print(f"  On-line zeros: {[float(mp.nstr(g, 4)) for g in zeros_5]}")
print()

a = mpf("10")

# Test at on-line positions (should be 0 or small)
for k in range(1, 4):
    gamma_k = mp.im(zetazero(k))
    contrib = weighted_test_function(gamma_k, a, zeros_5)
    print(f"  ĝ(γ_{k}) = {mp.nstr(contrib, 8)} (should be ≈ 0)")

# Test at off-line position
sigma_off = mpf("0.4")
gamma_off = zeros_5[0]  # Same imaginary part as first zero
delta = sigma_off - mpf("0.5")
u_offine = mpc(gamma_off, -delta)

offine = weighted_test_function(u_offine, a, zeros_5)
print()
print(f"  ĝ(off-line) = {mp.nstr(offine, 8)}")
print(f"  Re(ĝ) = {mp.nstr(mp.re(offine), 8)}")

# =============================================================================
# 4. OPTIMAL DETECTOR SEARCH
# =============================================================================

print("\n" + "="*70)
print("[4] OPTIMAL DETECTOR SEARCH")
print("="*70)

print("""
    Systematic search for test function that gives Z(g) < 0
    when off-line zero is injected.
    
    Parameters to optimize:
    - a (width parameter)
    - Number of polynomial roots
    - Shift parameter x₀
""")

def compute_Z_polynomial_weighted(sigma_off, gamma_off, a, n_roots, n_oneline=30):
    """
    Compute Z(g) for polynomial-weighted test function.
    
    Includes n_roots on-line zeros as polynomial roots.
    """
    a = mpf(a)
    
    # Get gamma values for roots
    root_gammas = [mp.im(zetazero(k)) for k in range(1, n_roots + 1)]
    
    # On-line contributions
    oneline_sum = mpf(0)
    for k in range(1, n_oneline + 1):
        gamma_k = mp.im(zetazero(k))
        contrib = weighted_test_function(gamma_k, a, root_gammas)
        oneline_sum += 2 * mp.re(contrib)
    
    # Off-line contribution
    delta = sigma_off - mpf("0.5")
    u_off = mpc(gamma_off, -delta)
    offine = weighted_test_function(u_off, a, root_gammas)
    offine_contrib = 2 * mp.re(offine)
    
    total = oneline_sum + offine_contrib
    
    return total, oneline_sum, offine_contrib

print("\n  Searching for violation with polynomial-weighted functions:")
print()
print("    n_roots   a        On-line      Off-line     Total        Status")
print("    " + "-"*70)

sigma_off = mpf("0.4")
gamma_off = mp.im(zetazero(1))

best_total = mpf("inf")
best_params = None

for n_roots in [3, 5, 10, 15]:
    for a in [1, 10, 100, 500]:
        try:
            total, oneline, offine = compute_Z_polynomial_weighted(
                sigma_off, gamma_off, a, n_roots, n_oneline=20
            )
            
            status = "✓" if total >= 0 else "✗ VIOLATION!"
            print(f"    {n_roots:5d}     {a:4d}    {mp.nstr(oneline, 6):>12}  "
                  f"{mp.nstr(offine, 6):>12}  {mp.nstr(total, 6):>12}    {status}")
            
            if total < best_total:
                best_total = total
                best_params = (n_roots, a)
        except:
            pass

print()
print(f"  Best result: n_roots={best_params[0]}, a={best_params[1]}, "
      f"Z = {mp.nstr(best_total, 8)}")

# =============================================================================
# 5. ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("[5] ANALYSIS AND CONCLUSIONS")
print("="*70)

print("""
    ═══════════════════════════════════════════════════════════════════════
    
    RESULTS:
    
    1. BANDLIMITED (sinc²):
       - When σ < γ₁, on-line contributions are ZERO
       - But off-line contribution at complex arg may also be tiny
       - The support needs careful tuning
    
    2. SHIFTED GAUSSIAN:
       - Phase shift creates oscillations
       - Can tune to specific locations
       - Still hard to make off-line dominate
    
    3. POLYNOMIAL-WEIGHTED:
       - Explicitly zeros out on-line contributions at root locations
       - Off-line at non-root positions contributes
       - Most promising for isolation
    
    ═══════════════════════════════════════════════════════════════════════
    
    THE FUNDAMENTAL DIFFICULTY:
    
    All test functions we've tried maintain ĝ ≥ 0 on the real line
    (admissibility requirement).
    
    This positivity constraint severely limits our ability to
    construct functions that are "negative" at complex off-line positions
    while remaining positive at real on-line positions.
    
    The Weil criterion says such detector functions EXIST if RH is true,
    but constructing them explicitly is the hard part.
    
    ═══════════════════════════════════════════════════════════════════════
    
    WHAT A COMPLETE PROOF WOULD REQUIRE:
    
    For EVERY hypothetical off-line position (σ, γ) with σ ≠ 1/2:
    
    Construct g_{σ,γ} such that:
    1. ĝ_{σ,γ}(u) ≥ 0 for all real u
    2. Z(g_{σ,γ}) < 0 when the zero at (σ,γ) is included
    
    This is equivalent to proving no such off-line zero can exist.
    
    The construction of such detector functions for ALL positions
    is essentially the full content of RH.
    
    ═══════════════════════════════════════════════════════════════════════
""")

print("\n" + "="*70)
print("SPECIALIZED DETECTOR COMPLETE")
print("="*70)
