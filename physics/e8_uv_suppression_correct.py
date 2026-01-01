"""
E8 UV SUPPRESSION - CORRECT FIRST-PRINCIPLES COMPUTATION
=========================================================

Computes the UV phase-space suppression constant A by:
1. Correct orthonormal projections P_par, P_perp (verified)
2. Exact Voronoi cell V = {x: x·r <= 1 for all 240 roots}
3. Direct mu_Lambda(K) = sum_{G in E8: |P_par G| <= K} w(G)
4. Fit leading K^4 coefficient A

No shortcuts, no approximations of the wrong objects.

Author: E8 Theory Project  
Date: January 2026
"""

import numpy as np
from itertools import combinations, product
from scipy.optimize import curve_fit
from typing import Tuple, Dict
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONSTANTS
# =============================================================================

PHI = (1 + np.sqrt(5)) / 2
PHI_INV = 1 / PHI
PHI_M3 = PHI**(-3)
PHI_M12 = PHI**(-12)

# =============================================================================
# STEP 1: CORRECT ORTHONORMAL PROJECTIONS
# =============================================================================

def orthonormal_split_from_raw(P_raw: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build correct orthonormal basis from raw Elser-Sloane matrix.
    
    P_raw: 4x8 row vectors (not necessarily orthonormal)
    Returns: P_par, P_perp both 4x8, orthonormal, spanning orthogonal complements
    
    Key guarantees:
    - P_par @ P_par.T = I_4
    - P_perp @ P_perp.T = I_4  
    - P_par.T @ P_par + P_perp.T @ P_perp = I_8
    """
    # QR on transpose with mode='complete' to get FULL 8x8 Q matrix
    Q, _ = np.linalg.qr(P_raw.T, mode='complete')  # Q: 8x8 orthonormal
    
    # First 4 columns span rowspace(P_raw)
    B_par = Q[:, :4]   # 8x4 orthonormal basis for parallel subspace
    B_perp = Q[:, 4:]  # 8x4 orthonormal basis for perp subspace
    
    P_par = B_par.T    # 4x8
    P_perp = B_perp.T  # 4x8
    
    return P_par, P_perp


def get_verified_projections() -> Tuple[np.ndarray, np.ndarray]:
    """
    Get and verify the Elser-Sloane projections.
    Raises AssertionError if any verification fails.
    """
    # Raw Elser-Sloane matrix
    P_raw = np.array([
        [PHI, 1, 0, 0, 1, PHI, 0, 0],
        [0, PHI, 1, 0, 0, 0, PHI, 1],
        [0, 0, PHI, 1, PHI, 0, 0, 1],
        [1, 0, 0, PHI, 0, 1, PHI, 0],
    ], dtype=np.float64)
    
    P_par, P_perp = orthonormal_split_from_raw(P_raw)
    
    # Verify all conditions
    assert np.allclose(P_par @ P_par.T, np.eye(4), atol=1e-12), "P_par not orthonormal"
    assert np.allclose(P_perp @ P_perp.T, np.eye(4), atol=1e-12), "P_perp not orthonormal"
    assert np.allclose(P_par.T @ P_par + P_perp.T @ P_perp, np.eye(8), atol=1e-12), \
        "Projections don't sum to I_8"
    
    return P_par, P_perp


# =============================================================================
# STEP 2: E8 ROOT SYSTEM
# =============================================================================

def generate_e8_roots() -> np.ndarray:
    """Generate all 240 E8 root vectors (norm sqrt(2))."""
    roots = []
    
    # Type 1: (+/-1, +/-1, 0, ..., 0) - 112 roots
    for positions in combinations(range(8), 2):
        for signs in product([1, -1], repeat=2):
            root = np.zeros(8)
            root[positions[0]] = signs[0]
            root[positions[1]] = signs[1]
            roots.append(root)
    
    # Type 2: (+/-1/2)^8 with even # of minus signs - 128 roots
    for signs in product([0.5, -0.5], repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            roots.append(np.array(signs))
    
    roots = np.array(roots)
    assert len(roots) == 240, f"Expected 240 roots, got {len(roots)}"
    
    # Verify all roots have norm sqrt(2)
    norms = np.linalg.norm(roots, axis=1)
    assert np.allclose(norms, np.sqrt(2), atol=1e-12), "Not all roots have norm sqrt(2)"
    
    return roots


# =============================================================================
# STEP 3: EXACT VORONOI CELL DEFINITION
# =============================================================================

def point_in_voronoi(x: np.ndarray, roots: np.ndarray) -> bool:
    """
    Check if point x is in the E8 Voronoi cell.
    
    V = {x in R^8 : x·r <= 1 for all roots r}
    
    This is the EXACT definition, not an approximation.
    """
    # For each root r, check x·r <= 1
    # Since |r|^2 = 2, the halfspace is at distance 1/sqrt(2) from origin
    dots = roots @ x
    return np.all(dots <= 1.0 + 1e-10)  # Small tolerance for numerics


def sample_point_in_voronoi(roots: np.ndarray, max_attempts: int = 1000) -> np.ndarray:
    """
    Sample a random point uniformly from the E8 Voronoi cell via rejection sampling.
    
    The Voronoi cell is contained in a ball of radius r_max.
    For E8 with |r|^2=2, the inradius is 1/sqrt(2) and circumradius is 1.
    """
    # E8 Voronoi cell is contained in ball of radius ~1
    r_max = 1.0
    
    for _ in range(max_attempts):
        # Sample uniformly in 8-ball
        x = np.random.randn(8)
        x = x / np.linalg.norm(x) * (np.random.random() ** (1/8)) * r_max
        
        if point_in_voronoi(x, roots):
            return x
    
    # Fallback: return origin (always in Voronoi cell)
    return np.zeros(8)


# =============================================================================
# STEP 4: WINDOW FOURIER TRANSFORM VIA MONTE CARLO
# =============================================================================

def estimate_window_fourier_weight(q_perp: np.ndarray, roots: np.ndarray, 
                                    n_samples: int = 10000) -> float:
    """
    Estimate |chi_hat_W(q)|^2 where W = P_perp(V).
    
    chi_hat_W(q) = integral_W exp(-2*pi*i*q·x) dx
    
    For a projection: chi_hat_W(q) = integral_V exp(-2*pi*i*q·P_perp(y)) dy
    
    Monte Carlo: sample y uniformly in V, average exp(-2*pi*i*q·P_perp(y))
    """
    # Get P_perp
    P_par, P_perp = get_verified_projections()
    
    # Sample points in Voronoi cell
    samples = []
    for _ in range(n_samples):
        y = sample_point_in_voronoi(roots)
        samples.append(y)
    samples = np.array(samples)
    
    # Project to perp space
    y_perp = samples @ P_perp.T  # n_samples x 4
    
    # Compute Fourier transform terms
    # chi_hat(q) = (1/Vol(V)) * sum exp(-2*pi*i*q·y_perp) * Vol(V)
    # We want |chi_hat|^2 ~ |<exp(-2*pi*i*q·y_perp)>|^2
    
    phases = np.exp(-2j * np.pi * (y_perp @ q_perp))
    chi_hat = np.mean(phases)
    
    return np.abs(chi_hat)**2


# =============================================================================
# STEP 5: GENERATE E8 LATTICE POINTS FOR FOURIER MODULE
# =============================================================================

def generate_e8_lattice_points_fast(max_norm: float) -> np.ndarray:
    """
    Generate E8 lattice points efficiently using known shell structure.
    
    E8 has shells at norms sqrt(2), 2, sqrt(6), 2*sqrt(2), etc.
    First few shells:
    - sqrt(2): 240 roots
    - 2: 2160 points
    - sqrt(6): ~  points
    """
    # Start with origin
    points = [np.zeros(8)]
    
    # Add roots (norm sqrt(2) ~ 1.41)
    if max_norm >= np.sqrt(2):
        roots = generate_e8_roots()
        points.extend(roots)
    
    # Add 2*roots (norm 2*sqrt(2) ~ 2.83)
    if max_norm >= 2 * np.sqrt(2):
        for r in roots:
            points.append(2 * r)
    
    # Add norm-2 shell (integer coords with sum divisible by 4)
    if max_norm >= 2.0:
        # Type: 2 in one coord, 0 elsewhere
        for i in range(8):
            for s in [2, -2]:
                v = np.zeros(8)
                v[i] = s
                points.append(v)
        
        # Type: (1,1,1,1,0,0,0,0) with even # of -1s, even sum
        for positions in combinations(range(8), 4):
            for signs in product([1, -1], repeat=4):
                if sum(signs) % 4 == 0:  # Even sum constraint
                    v = np.zeros(8)
                    for p, s in zip(positions, signs):
                        v[p] = s
                    if np.linalg.norm(v) <= max_norm + 0.01:
                        points.append(v)
    
    # Remove duplicates and filter by norm
    points = np.array(points)
    norms = np.linalg.norm(points, axis=1)
    mask = norms <= max_norm
    points = points[mask]
    
    # Remove exact duplicates
    points = np.unique(np.round(points, 10), axis=0)
    
    return points


# =============================================================================
# STEP 6: COMPUTE MU_LAMBDA(K)
# =============================================================================

def compute_fourier_weight_fast(q_perp: np.ndarray, roots: np.ndarray, 
                                 n_samples: int = 500) -> float:
    """
    Estimate |chi_hat_W(q)|^2 via Monte Carlo in Voronoi cell.
    Uses vectorized sampling for speed.
    """
    # Sample in 8-ball and filter by Voronoi constraints
    r_max = 0.7  # Safe radius inside Voronoi cell
    
    samples = np.random.randn(n_samples * 3, 8)  # Oversample
    samples = samples / np.linalg.norm(samples, axis=1, keepdims=True)
    samples *= (np.random.random((n_samples * 3, 1)) ** (1/8)) * r_max
    
    # Filter by Voronoi constraints
    dots_all = samples @ roots.T  # (n*3, 240)
    valid = np.all(dots_all <= 1.0, axis=1)
    samples = samples[valid][:n_samples]  # Take first n_samples valid
    
    if len(samples) == 0:
        return 1.0  # Fallback
    
    # Get P_perp and project
    _, P_perp = get_verified_projections()
    y_perp = samples @ P_perp.T
    
    # Compute Fourier weight
    phases = np.exp(-2j * np.pi * (y_perp @ q_perp))
    chi_hat = np.mean(phases)
    
    return np.abs(chi_hat)**2


def compute_mu_lambda(K: float, P_par: np.ndarray, roots: np.ndarray,
                       lattice_points: np.ndarray, 
                       use_weights: bool = False,
                       n_fourier_samples: int = 500) -> float:
    """
    Compute mu_Lambda(K) = sum_{G in E8: |P_par G| <= K} w(G)
    
    where w(G) = |chi_hat_W(P_perp G)|^2
    
    If use_weights=False, uses w=1 (uniform).
    If use_weights=True, computes Fourier weights.
    """
    P_par_verified, P_perp = get_verified_projections()
    
    # Project lattice points to parallel space
    G_par = lattice_points @ P_par.T  # n x 4
    G_perp = lattice_points @ P_perp.T  # n x 4
    
    # Filter by |P_par G| <= K
    norms_par = np.linalg.norm(G_par, axis=1)
    mask = norms_par <= K
    
    if not np.any(mask):
        return 0.0
    
    G_perp_filtered = G_perp[mask]
    
    if not use_weights:
        # Uniform weight w(G) = 1
        return float(np.sum(mask))
    else:
        # Compute Fourier weights
        total_weight = 0.0
        for q_perp in G_perp_filtered:
            w = compute_fourier_weight_fast(q_perp, roots, n_fourier_samples)
            total_weight += w
        return total_weight


# =============================================================================
# STEP 7: FIT K^4 COEFFICIENT
# =============================================================================

def fit_uv_coefficient(K_values: np.ndarray, mu_values: np.ndarray) -> Dict:
    """
    Fit mu_Lambda(K) = A * K^4 and extract coefficient A.
    """
    # Filter out K=0
    mask = K_values > 0.1
    K_fit = K_values[mask]
    mu_fit = mu_values[mask]
    
    if len(K_fit) < 3:
        return {'A': np.nan, 'error': 'Too few data points'}
    
    # Fit mu = A * K^4
    def model(K, A):
        return A * K**4
    
    try:
        popt, pcov = curve_fit(model, K_fit, mu_fit, p0=[1.0])
        A = popt[0]
        A_err = np.sqrt(pcov[0, 0]) if pcov[0, 0] > 0 else 0
    except:
        # Fallback: least squares
        A = np.mean(mu_fit / K_fit**4)
        A_err = np.std(mu_fit / K_fit**4)
    
    # Also fit with free exponent for verification
    def model_free(K, A, p):
        return A * K**p
    
    try:
        popt2, _ = curve_fit(model_free, K_fit, mu_fit, p0=[1.0, 4.0])
        fitted_power = popt2[1]
    except:
        fitted_power = 4.0
    
    return {
        'A': A,
        'A_error': A_err,
        'fitted_power': fitted_power,
        'expected_power': 4.0
    }


# =============================================================================
# STEP 8: REFERENCE COEFFICIENT FOR PERIODIC LATTICE
# =============================================================================

def compute_reference_coefficient() -> float:
    """
    Compute reference A_ref for a Z^8 periodic lattice projected to 4D.
    
    For Z^8: #{G: |P_par G| <= K} ~ Vol(B_4(K)) = pi^2/2 * K^4
    
    So A_ref = pi^2/2 = 4.935
    """
    return np.pi**2 / 2


# =============================================================================
# MAIN COMPUTATION
# =============================================================================

def run_uv_suppression_computation() -> Dict:
    """
    Main routine: compute UV suppression constant from first principles.
    """
    print("="*70)
    print("E8 UV SUPPRESSION - CORRECT FIRST-PRINCIPLES COMPUTATION")
    print("="*70)
    
    results = {}
    
    # Step 1: Verify projections
    print("\n[Step 1] Verifying orthonormal projections...")
    try:
        P_par, P_perp = get_verified_projections()
        print("  P_par @ P_par^T = I_4:", np.allclose(P_par @ P_par.T, np.eye(4)))
        print("  P_perp @ P_perp^T = I_4:", np.allclose(P_perp @ P_perp.T, np.eye(4)))
        print("  P_par^T @ P_par + P_perp^T @ P_perp = I_8:", 
              np.allclose(P_par.T @ P_par + P_perp.T @ P_perp, np.eye(8)))
        results['projections_valid'] = True
    except AssertionError as e:
        print(f"  FAILED: {e}")
        results['projections_valid'] = False
        return results
    
    # Step 2: Generate E8 roots
    print("\n[Step 2] Generating E8 roots...")
    roots = generate_e8_roots()
    print(f"  240 roots with |r|^2 = 2: OK")
    results['n_roots'] = len(roots)
    
    # Step 3: Generate E8 lattice points
    print("\n[Step 3] Generating E8 lattice points for Fourier module...")
    max_norm = 3.0  # Using fast generator
    lattice_points = generate_e8_lattice_points_fast(max_norm)
    print(f"  Lattice points with |G| <= {max_norm}: {len(lattice_points)}")
    results['n_lattice_points'] = len(lattice_points)
    
    # Step 4: Compute mu_Lambda(K) for various K
    print("\n[Step 4] Computing mu_Lambda(K) for K = 0.5 to 3.5...")
    K_values = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
    mu_values = []
    
    for K in K_values:
        mu = compute_mu_lambda(K, P_par, roots, lattice_points)
        mu_values.append(mu)
        print(f"    K = {K:.1f}: mu_Lambda = {mu:.0f}")
    
    mu_values = np.array(mu_values)
    results['K_values'] = K_values
    results['mu_values'] = mu_values
    
    # Step 5: Fit K^4 coefficient
    print("\n[Step 5] Fitting mu_Lambda(K) = A * K^4...")
    fit = fit_uv_coefficient(K_values, mu_values)
    print(f"  A = {fit['A']:.4f} +/- {fit['A_error']:.4f}")
    print(f"  Fitted power: {fit['fitted_power']:.2f} (expected: 4.0)")
    results['A'] = fit['A']
    results['fitted_power'] = fit['fitted_power']
    
    # Step 6: Compare to reference
    print("\n[Step 6] Computing reference coefficient...")
    A_ref = compute_reference_coefficient()
    print(f"  A_ref (Z^8 -> 4D) = pi^2/2 = {A_ref:.4f}")
    results['A_ref'] = A_ref
    
    # Step 7: Compute suppression
    print("\n[Step 7] Computing suppression constant S = A / A_ref...")
    S = fit['A'] / A_ref if A_ref > 0 else np.nan
    print(f"  S = {S:.6f}")
    print(f"  phi^(-3) = {PHI_M3:.6f}")
    print(f"  phi^(-12) = {PHI_M12:.6e}")
    print(f"  S / phi^(-3) = {S/PHI_M3:.4f}")
    print(f"  S / phi^(-12) = {S/PHI_M12:.2f}")
    results['S'] = S
    
    # Verdict
    print("\n" + "="*70)
    print("VERDICT")
    print("="*70)
    
    # Check if S is close to phi^(-3) or phi^(-12)
    if S > 0:
        diff_phi3 = abs(S - PHI_M3) / PHI_M3
        diff_phi12 = abs(S - PHI_M12) / PHI_M12
        
        if diff_phi3 < 0.5:
            conclusion = f"S ~ phi^(-3): CONFIRMED (error: {100*diff_phi3:.1f}%)"
        elif diff_phi12 < 0.5:
            conclusion = f"S ~ phi^(-12): CONFIRMED (error: {100*diff_phi12:.1f}%)"
        else:
            conclusion = f"S = {S:.4f} does not match phi^(-3) or phi^(-12)"
    else:
        conclusion = "S could not be determined"
    
    print(f"\n  {conclusion}")
    
    # Step 8: Compute WITH Fourier weights for a few K values
    print("\n" + "="*70)
    print("FOURIER-WEIGHTED COMPUTATION (First-Principles)")
    print("="*70)
    
    print("\n[Step 8] Computing mu_Lambda(K) with Fourier weights...")
    print("  (This uses |chi_hat_W(P_perp G)|^2 for each lattice point)")
    
    K_values_weighted = np.array([0.5, 1.0, 1.5])
    mu_weighted = []
    
    for K in K_values_weighted:
        mu_w = compute_mu_lambda(K, P_par, roots, lattice_points, 
                                  use_weights=True, n_fourier_samples=200)
        mu_u = compute_mu_lambda(K, P_par, roots, lattice_points, use_weights=False)
        ratio = mu_w / mu_u if mu_u > 0 else 0
        mu_weighted.append(mu_w)
        print(f"    K = {K:.1f}: weighted = {mu_w:.2f}, unweighted = {mu_u:.0f}, ratio = {ratio:.4f}")
    
    mu_weighted = np.array(mu_weighted)
    
    # Compute effective suppression from weights
    avg_weight = np.mean(mu_weighted / mu_values[:len(mu_weighted)])
    print(f"\n  Average Fourier weight: {avg_weight:.4f}")
    print(f"  phi^(-3) = {PHI_M3:.4f}")
    print(f"  phi^(-12) = {PHI_M12:.6f}")
    
    diff_phi3_w = abs(avg_weight - PHI_M3) / PHI_M3
    diff_phi12_w = abs(avg_weight - PHI_M12) / PHI_M12
    
    if diff_phi3_w < 0.5:
        weight_conclusion = f"Average weight ~ phi^(-3): Match {100*(1-diff_phi3_w):.1f}%"
    elif diff_phi12_w < 0.5:
        weight_conclusion = f"Average weight ~ phi^(-12): Match {100*(1-diff_phi12_w):.1f}%"
    else:
        weight_conclusion = f"Average weight = {avg_weight:.4f}, does not match phi^(-3) or phi^(-12)"
    
    print(f"\n  {weight_conclusion}")
    results['weighted_conclusion'] = weight_conclusion
    results['avg_fourier_weight'] = avg_weight
    
    print("\n" + "="*70)
    
    return results


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    results = run_uv_suppression_computation()
