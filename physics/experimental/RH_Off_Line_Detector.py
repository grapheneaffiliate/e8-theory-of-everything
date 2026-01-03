#!/usr/bin/env python3
"""
RH OFF-LINE ZERO DETECTOR
=========================

THE REMAINING CHALLENGE:
If ANY zero is off-line, show Z(g) < 0 for SOME admissible g.

STRATEGY: Construct test functions that "resonate" with off-line zeros
and produce negative contributions that dominate.

KEY INSIGHT: For off-line zeros at ρ = σ ± iγ (σ ≠ 1/2):
- The contribution (ρ - 1/2)/i is COMPLEX, not real
- This creates OSCILLATING contributions that can go negative
"""

from mpmath import mp, mpf, mpc, pi, log, exp, zeta, sqrt, cos, sin
from mpmath import gamma as mpgamma, zetazero
import math

mp.dps = 50

print("="*70)
print("RH OFF-LINE ZERO DETECTOR")
print("Testing Positivity Violation for Off-Line Zeros")
print("="*70)

# =============================================================================
# 1. THE CRITICAL INSIGHT
# =============================================================================

print("\n" + "="*70)
print("[1] THE CRITICAL INSIGHT")
print("="*70)

print("""
    For zeros on the critical line: ρ = 1/2 + iγ
    
    The zero-sum contribution:
    ĝ((ρ - 1/2)/i) = ĝ(γ)  [REAL argument]
    
    For ĝ ≥ 0 (like Gaussian), this is automatically ≥ 0.
    
    ════════════════════════════════════════════════════════
    
    For HYPOTHETICAL off-line zeros: ρ = σ + iγ with σ ≠ 1/2
    
    The argument becomes:
    (ρ - 1/2)/i = (σ - 1/2)/i + γ = γ - i(σ - 1/2)  [COMPLEX!]
    
    So: ĝ((ρ - 1/2)/i) = ĝ(γ - i(σ - 1/2))
    
    For Gaussian: ĝ(u) = √(π/a) exp(-π²u²/a)
    
    At complex u = γ - iδ (where δ = σ - 1/2):
    
    u² = γ² - 2iγδ - δ²
    
    exp(-π²u²/a) = exp(-π²(γ² - δ²)/a) × exp(2iπ²γδ/a)
    
    This has OSCILLATING PHASE!
    
    The real part can become NEGATIVE when cos(2π²γδ/a) < 0.
""")

# =============================================================================
# 2. INJECTING A HYPOTHETICAL OFF-LINE ZERO
# =============================================================================

print("\n" + "="*70)
print("[2] HYPOTHETICAL OFF-LINE ZERO ANALYSIS")
print("="*70)

def analyze_offline_contribution(sigma, gamma, a=1.0):
    """
    Compute the contribution from a hypothetical off-line zero
    at ρ = sigma + i*gamma to Z(g) for Gaussian test function.
    """
    sigma = mpf(sigma)
    gamma = mpf(gamma)
    a = mpf(a)
    
    # The argument: (ρ - 1/2)/i = γ - i(σ - 1/2)
    delta = sigma - mpf("0.5")
    u_real = gamma
    u_imag = -delta  # Note: (σ-1/2)/i = -i(σ-1/2)
    
    u = mpc(u_real, u_imag)
    
    # ĝ(u) = √(π/a) exp(-π²u²/a)
    u_squared = u * u  # = γ² - 2iγδ - δ²
    exponent = -mp.pi**2 * u_squared / a
    
    g_hat_u = mp.sqrt(mp.pi / a) * mp.exp(exponent)
    
    return g_hat_u, u, u_squared

# Test with hypothetical off-line zero
print("\n  Hypothetical off-line zero at σ = 0.4, γ = 14.13:")
print()

sigma_off = 0.4
gamma_off = 14.134725  # Similar to first actual zero

for a in [0.1, 1.0, 10.0, 100.0, 1000.0]:
    contrib, u, u_sq = analyze_offline_contribution(sigma_off, gamma_off, a)
    re_part = mp.re(contrib)
    im_part = mp.im(contrib)
    
    print(f"    a = {a:6.1f}: ĝ(u) = {mp.nstr(re_part, 8)} + {mp.nstr(im_part, 8)}i")
    print(f"             |ĝ(u)| = {mp.nstr(abs(contrib), 8)}, Re(ĝ) = {mp.nstr(re_part, 8)}")
    
    if re_part < 0:
        print(f"             *** NEGATIVE CONTRIBUTION! ***")
    print()

# =============================================================================
# 3. PHASE ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("[3] PHASE ANALYSIS FOR OFF-LINE ZEROS")
print("="*70)

print("""
    For ρ = σ + iγ with σ ≠ 1/2, let δ = σ - 1/2.
    
    The complex argument: u = γ - iδ
    
    u² = γ² - 2iγδ - δ² = (γ² - δ²) - 2iγδ
    
    exp(-π²u²/a) = exp(-π²(γ² - δ²)/a) × exp(2iπ²γδ/a)
                 = R × e^{iθ}
    
    where:
    R = exp(-π²(γ² - δ²)/a)    [amplitude]
    θ = 2π²γδ/a                [phase]
    
    The real part of ĝ(u):
    Re(ĝ(u)) = √(π/a) × R × cos(θ)
    
    This is NEGATIVE when cos(θ) < 0, i.e., when:
    π/2 + nπ < θ < 3π/2 + nπ  (for integer n)
    
    Or: π/2 + nπ < 2π²γδ/a < 3π/2 + nπ
""")

def find_negative_phase(sigma, gamma):
    """
    Find values of a where the off-line contribution is negative.
    """
    delta = abs(sigma - 0.5)
    
    # θ = 2π²γδ/a
    # We want cos(θ) < 0, i.e., θ in (π/2, 3π/2) + 2πn
    
    theta_target = mp.pi/2 + mp.pi/4  # Center of first negative region
    
    # θ = 2π²γδ/a → a = 2π²γδ/θ
    a_for_negative = 2 * mp.pi**2 * gamma * delta / theta_target
    
    return a_for_negative

print("\n  For hypothetical zero at σ = 0.4, γ = 14.13:")
print()

sigma_off = 0.4
gamma_off = 14.134725
delta = abs(sigma_off - 0.5)

# Scan a values to find negative phases
print("    Scanning phase behavior:")
print()
print("    a          θ (phase)      cos(θ)        Re(ĝ) sign")
print("    " + "-"*60)

for a in [10, 50, 100, 200, 500, 1000, 2000, 5000]:
    theta = 2 * mp.pi**2 * gamma_off * delta / a
    cos_theta = mp.cos(theta)
    contrib, _, _ = analyze_offline_contribution(sigma_off, gamma_off, a)
    re_sign = "POSITIVE" if mp.re(contrib) >= 0 else "NEGATIVE ←"
    
    print(f"    {a:5d}      {mp.nstr(theta, 6):>10}     {mp.nstr(cos_theta, 6):>10}    {re_sign}")

# =============================================================================
# 4. CONSTRUCT DETECTOR TEST FUNCTION
# =============================================================================

print("\n" + "="*70)
print("[4] DETECTOR TEST FUNCTION")
print("="*70)

print("""
    To detect an off-line zero, we need a test function where:
    
    1. ĝ ≥ 0 on the real line (admissibility)
    2. The off-line contribution Re(ĝ(γ - iδ)) < 0 (negative)
    3. The on-line contributions don't overwhelm it
    
    For Gaussian with large a, the off-line phase can produce cos(θ) < 0.
    
    But we need the TOTAL sum to go negative, not just one term.
""")

def compute_total_Z_with_offline(sigma_off, gamma_off, a, N_oneline=100):
    """
    Compute Z(g) assuming:
    - N_oneline zeros on the critical line (actual zeros)
    - One extra hypothetical off-line zero at (sigma_off, gamma_off)
    """
    a = mpf(a)
    total = mpf(0)
    
    # On-line contributions (all positive)
    oneline_sum = mpf(0)
    for k in range(1, N_oneline + 1):
        rho_k = zetazero(k)
        gamma_k = mp.im(rho_k)
        
        # ĝ(γ) for on-line zero
        term = mp.sqrt(mp.pi / a) * mp.exp(-mp.pi**2 * gamma_k**2 / a)
        oneline_sum += 2 * term  # Factor of 2 for conjugate pair
    
    # Off-line contribution (can be negative!)
    offine_contrib, _, _ = analyze_offline_contribution(sigma_off, gamma_off, a)
    offine_real = mp.re(offine_contrib)
    
    # Total: on-line + off-line pair (factor of 2)
    total = oneline_sum + 2 * offine_real
    
    return total, oneline_sum, 2 * offine_real

print("\n  Computing Z(g) with hypothetical off-line zero:")
print("  (σ = 0.4, γ = 14.13)")
print()
print("    a         On-line sum     Off-line       Total Z(g)       Status")
print("    " + "-"*70)

sigma_off = 0.4
gamma_off = 14.134725

for a in [100, 500, 1000, 2000, 5000]:
    total, oneline, offine = compute_total_Z_with_offline(sigma_off, gamma_off, a, N_oneline=50)
    status = "≥ 0 ✓" if total >= 0 else "< 0 ✗ VIOLATION!"
    
    print(f"    {a:5d}     {mp.nstr(oneline, 6):>12}   {mp.nstr(offine, 6):>12}   {mp.nstr(total, 6):>12}   {status}")

# =============================================================================
# 5. FINDING THE VIOLATION
# =============================================================================

print("\n" + "="*70)
print("[5] SEARCHING FOR POSITIVITY VIOLATION")
print("="*70)

print("""
    Strategy: Vary parameters to find Z(g) < 0 with off-line zero.
    
    Key factors:
    - Large a makes off-line phase significant
    - But large a also suppresses on-line contributions (good!)
    - The off-line zero at low γ has more impact
""")

print("\n  Fine-grained search for violation:")
print()

# Try different off-line positions
violations_found = []

for sigma_off in [0.3, 0.4, 0.45]:
    for gamma_off in [14.13, 21.02, 25.01]:  # Near actual zeros
        delta = abs(sigma_off - 0.5)
        
        # Scan a values
        for a in range(100, 10001, 100):
            total, oneline, offine = compute_total_Z_with_offline(
                sigma_off, gamma_off, a, N_oneline=30
            )
            
            if total < 0:
                violations_found.append((sigma_off, gamma_off, a, total))

if violations_found:
    print("    VIOLATIONS FOUND:")
    for v in violations_found[:5]:
        print(f"      σ={v[0]}, γ={v[1]:.2f}, a={v[2]}: Z = {mp.nstr(v[3], 8)}")
else:
    print("    No violations found in search range.")
    print("    (Off-line contribution may be too small compared to on-line)")

# =============================================================================
# 6. ANALYSIS OF RESULTS
# =============================================================================

print("\n" + "="*70)
print("[6] ANALYSIS")
print("="*70)

print("""
    ═══════════════════════════════════════════════════════════════════════
    
    KEY FINDING:
    
    The off-line contribution CAN be negative (due to complex phase).
    
    However, making the TOTAL Z(g) negative requires:
    - Off-line contribution large enough to overwhelm on-line sum
    - This is hard because on-line zeros near γ₁ ≈ 14.13 contribute positively
    
    THE CHALLENGE:
    
    For Gaussian test functions, the on-line contributions scale similarly
    to the off-line magnitude, making it hard to get Z < 0.
    
    This suggests we need MORE SPECIALIZED test functions:
    - Functions that suppress on-line contributions
    - While amplifying the off-line negative phase
    
    POSSIBLE APPROACHES:
    
    1. Bandlimited functions (sinc²) with support avoiding actual zeros
    2. Polynomially-weighted Gaussians
    3. Functions tuned to specific off-line positions
    
    ═══════════════════════════════════════════════════════════════════════
    
    PARTIAL RESULT:
    
    We've shown OFF-LINE zeros produce COMPLEX contributions
    with potentially NEGATIVE real parts.
    
    Full proof requires showing this negativity dominates for SOME test function.
    
    This is the core difficulty of the Weil approach to RH.
    
    ═══════════════════════════════════════════════════════════════════════
""")

print("\n" + "="*70)
print("OFF-LINE DETECTOR COMPLETE")
print("="*70)
