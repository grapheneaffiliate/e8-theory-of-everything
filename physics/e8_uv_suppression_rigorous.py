"""
E8 UV SUPPRESSION - RIGOROUS FIRST-PRINCIPLES COMPUTATION
==========================================================

Correct approach:
1. Sphere enumeration for E8 and Z^8 (not just shells)
2. Increasing R cutoffs with convergence check
3. Ratio S = A(E8) / A(Z^8) (normalization-independent)
4. Verify fitted power ≈ 4.0 to confirm asymptotic regime

NO WEIGHTS - just clean counting ratio.

Author: E8 Theory Project  
Date: January 2026
"""

import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple, List, Dict
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONSTANTS
# =============================================================================

PHI = (1 + np.sqrt(5)) / 2
PHI_INV = 1 / PHI

# =============================================================================
# PROJECTION MATRICES (VERIFIED ORTHONORMAL)
# =============================================================================

def get_orthonormal_projections() -> Tuple[np.ndarray, np.ndarray]:
    """
    Get verified orthonormal projections P_par, P_perp.
    Both 4x8, satisfying P_par P_par^T = I_4, P_perp P_perp^T = I_4,
    P_par^T P_par + P_perp^T P_perp = I_8.
    """
    P_raw = np.array([
        [PHI, 1, 0, 0, 1, PHI, 0, 0],
        [0, PHI, 1, 0, 0, 0, PHI, 1],
        [0, 0, PHI, 1, PHI, 0, 0, 1],
        [1, 0, 0, PHI, 0, 1, PHI, 0],
    ], dtype=np.float64)
    
    Q, _ = np.linalg.qr(P_raw.T, mode='complete')
    P_par = Q[:, :4].T   # 4x8
    P_perp = Q[:, 4:].T  # 4x8
    
    return P_par, P_perp


# =============================================================================
# SPHERE ENUMERATION
# =============================================================================

def enumerate_z8_in_ball(R: float) -> np.ndarray:
    """
    Enumerate all points in Z^8 with |x| <= R.
    Uses bounded integer enumeration.
    """
    R_int = int(np.ceil(R))
    points = []
    
    # 8D bounded enumeration
    # For efficiency, we iterate and filter
    def recurse(partial: List[int], remaining_dims: int, remaining_r2: float):
        if remaining_dims == 0:
            points.append(np.array(partial, dtype=float))
            return
        
        # Max absolute value for this coordinate
        max_val = int(np.floor(np.sqrt(remaining_r2) + 1e-9))
        
        for v in range(-max_val, max_val + 1):
            new_r2 = remaining_r2 - v**2
            if new_r2 >= -1e-9:
                recurse(partial + [v], remaining_dims - 1, new_r2)
    
    recurse([], 8, R**2)
    
    if len(points) == 0:
        return np.zeros((0, 8))
    
    return np.array(points)


def enumerate_e8_in_ball(R: float) -> np.ndarray:
    """
    Enumerate all points in E8 lattice with |x| <= R.
    
    E8 = {x in Z^8 : sum(x) ≡ 0 (mod 2)}
         ∪ {x in (Z+1/2)^8 : sum(x) ≡ 0 (mod 2)}
    """
    R_int = int(np.ceil(R))
    points = []
    
    # Type 1: Integer coordinates with even sum
    def recurse_int(partial: List[int], remaining_dims: int, remaining_r2: float, partial_sum: int):
        if remaining_dims == 0:
            if partial_sum % 2 == 0:
                points.append(np.array(partial, dtype=float))
            return
        
        max_val = int(np.floor(np.sqrt(remaining_r2) + 1e-9))
        
        for v in range(-max_val, max_val + 1):
            new_r2 = remaining_r2 - v**2
            if new_r2 >= -1e-9:
                recurse_int(partial + [v], remaining_dims - 1, new_r2, partial_sum + v)
    
    recurse_int([], 8, R**2, 0)
    n_type1 = len(points)
    
    # Type 2: Half-integer coordinates with even sum of integer parts
    def recurse_half(partial: List[int], remaining_dims: int, remaining_r2: float, partial_sum: int):
        if remaining_dims == 0:
            if partial_sum % 2 == 0:
                # Convert integer parts to half-integers
                half_vec = np.array(partial, dtype=float) + 0.5
                points.append(half_vec)
            return
        
        # We enumerate integer parts k where x = k + 0.5
        # |x| = |k + 0.5|, so need to check k and k+1 ranges
        max_k = int(np.floor(np.sqrt(remaining_r2) + 0.5))
        
        for k in range(-max_k - 1, max_k + 1):
            x_val = k + 0.5
            new_r2 = remaining_r2 - x_val**2
            if new_r2 >= -1e-9:
                recurse_half(partial + [k], remaining_dims - 1, new_r2, partial_sum + k)
    
    recurse_half([], 8, R**2, 0)
    
    if len(points) == 0:
        return np.zeros((0, 8))
    
    return np.array(points)


# =============================================================================
# COUNTING FUNCTION μ(K)
# =============================================================================

def compute_mu(points: np.ndarray, P_par: np.ndarray, K_values: np.ndarray) -> np.ndarray:
    """
    Compute μ_R(K) = #{x in L : |P_par x| <= K}.
    """
    if len(points) == 0:
        return np.zeros(len(K_values))
    
    # Project to parallel space
    k_norms = np.linalg.norm(points @ P_par.T, axis=1)
    
    # Count for each K
    mu = np.array([np.sum(k_norms <= K) for K in K_values])
    
    return mu


# =============================================================================
# K^4 FITTING WITH CONVERGENCE CHECK
# =============================================================================

def fit_k4_coefficient(K_values: np.ndarray, mu_values: np.ndarray, 
                       fit_range: Tuple[float, float] = None) -> Dict:
    """
    Fit μ(K) = A * K^4 and verify the power is close to 4.
    
    Returns A, fitted_power, and fit quality metrics.
    """
    # Filter to fit range
    if fit_range is not None:
        mask = (K_values >= fit_range[0]) & (K_values <= fit_range[1])
        K_fit = K_values[mask]
        mu_fit = mu_values[mask]
    else:
        K_fit = K_values[K_values > 0.1]
        mu_fit = mu_values[K_values > 0.1]
    
    if len(K_fit) < 3 or np.all(mu_fit == 0):
        return {'A': np.nan, 'power': np.nan, 'R2': 0, 'valid': False}
    
    # Fit μ = A * K^p (free power)
    def model_free(K, A, p):
        return A * K**p
    
    try:
        # Initial guess
        p0_A = mu_fit[-1] / K_fit[-1]**4 if K_fit[-1] > 0 else 1.0
        popt, pcov = curve_fit(model_free, K_fit, mu_fit, p0=[p0_A, 4.0], 
                              bounds=([0, 0], [np.inf, 10]))
        A = popt[0]
        power = popt[1]
        
        # Compute R^2
        mu_pred = model_free(K_fit, *popt)
        ss_res = np.sum((mu_fit - mu_pred)**2)
        ss_tot = np.sum((mu_fit - np.mean(mu_fit))**2)
        R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        
        # Valid if power is close to 4
        valid = abs(power - 4.0) < 1.0 and R2 > 0.95
        
    except Exception as e:
        return {'A': np.nan, 'power': np.nan, 'R2': 0, 'valid': False, 'error': str(e)}
    
    return {
        'A': A,
        'power': power,
        'R2': R2,
        'valid': valid,
        'K_range': (K_fit[0], K_fit[-1]),
        'n_points': len(K_fit)
    }


# =============================================================================
# MAIN COMPUTATION
# =============================================================================

def run_rigorous_computation(R_values: List[float] = None, 
                             K_max_fraction: float = 0.6) -> Dict:
    """
    Run the rigorous suppression ratio computation.
    
    For each R:
    1. Enumerate E8 and Z^8 points with |x| <= R
    2. Compute μ_R(K) for K in (0, K_max) where K_max = R * K_max_fraction
    3. Fit A * K^4 in the stable regime
    4. Check convergence across R values
    5. Compute ratio S = A(E8) / A(Z^8)
    """
    if R_values is None:
        R_values = [3.0, 4.0, 5.0, 6.0]  # Can increase for better convergence
    
    print("="*70)
    print("E8 UV SUPPRESSION - RIGOROUS FIRST-PRINCIPLES COMPUTATION")
    print("="*70)
    
    # Get projections
    P_par, P_perp = get_orthonormal_projections()
    print(f"\n[Projections verified]")
    print(f"  P_par @ P_par^T = I_4: {np.allclose(P_par @ P_par.T, np.eye(4))}")
    print(f"  P_perp @ P_perp^T = I_4: {np.allclose(P_perp @ P_perp.T, np.eye(4))}")
    
    results = {
        'R_values': R_values,
        'E8': {},
        'Z8': {},
    }
    
    for R in R_values:
        print(f"\n{'='*70}")
        print(f"R = {R}")
        print(f"{'='*70}")
        
        # K grid: up to K_max_fraction * R to avoid truncation effects
        K_max = R * K_max_fraction
        K_values = np.linspace(0.1, K_max, 50)
        
        # E8 enumeration
        print(f"\n[E8] Enumerating points with |x| <= {R}...")
        e8_points = enumerate_e8_in_ball(R)
        print(f"  E8 points: {len(e8_points)}")
        
        mu_e8 = compute_mu(e8_points, P_par, K_values)
        fit_e8 = fit_k4_coefficient(K_values, mu_e8, fit_range=(0.5, K_max * 0.8))
        
        print(f"  μ(K_max): {mu_e8[-1]:.0f}")
        print(f"  Fit: A = {fit_e8['A']:.4f}, power = {fit_e8['power']:.3f}, R² = {fit_e8['R2']:.4f}")
        print(f"  Valid (power ≈ 4): {fit_e8['valid']}")
        
        results['E8'][R] = {
            'n_points': len(e8_points),
            'K_values': K_values,
            'mu': mu_e8,
            'fit': fit_e8
        }
        
        # Z^8 enumeration
        print(f"\n[Z^8] Enumerating points with |x| <= {R}...")
        z8_points = enumerate_z8_in_ball(R)
        print(f"  Z^8 points: {len(z8_points)}")
        
        mu_z8 = compute_mu(z8_points, P_par, K_values)
        fit_z8 = fit_k4_coefficient(K_values, mu_z8, fit_range=(0.5, K_max * 0.8))
        
        print(f"  μ(K_max): {mu_z8[-1]:.0f}")
        print(f"  Fit: A = {fit_z8['A']:.4f}, power = {fit_z8['power']:.3f}, R² = {fit_z8['R2']:.4f}")
        print(f"  Valid (power ≈ 4): {fit_z8['valid']}")
        
        results['Z8'][R] = {
            'n_points': len(z8_points),
            'K_values': K_values,
            'mu': mu_z8,
            'fit': fit_z8
        }
        
        # Compute ratio at this R
        if fit_e8['valid'] and fit_z8['valid'] and fit_z8['A'] > 0:
            ratio = fit_e8['A'] / fit_z8['A']
            print(f"\n  >>> S(R={R}) = A(E8) / A(Z8) = {ratio:.6f}")
            results['E8'][R]['ratio'] = ratio
        else:
            print(f"\n  >>> S(R={R}): Cannot compute (fit not valid)")
            results['E8'][R]['ratio'] = np.nan
    
    # Convergence analysis
    print(f"\n{'='*70}")
    print("CONVERGENCE ANALYSIS")
    print(f"{'='*70}")
    
    ratios = []
    for R in R_values:
        r = results['E8'][R].get('ratio', np.nan)
        if not np.isnan(r):
            ratios.append((R, r))
            print(f"  R = {R}: S = {r:.6f}")
    
    if len(ratios) >= 2:
        R_last, S_last = ratios[-1]
        R_prev, S_prev = ratios[-2]
        rel_change = abs(S_last - S_prev) / S_prev if S_prev > 0 else np.inf
        print(f"\n  Relative change (R={R_prev} → {R_last}): {100*rel_change:.2f}%")
        
        if rel_change < 0.1:  # Less than 10% change
            print(f"  [OK] Ratio appears to be converging")
            results['converged'] = True
            results['S_final'] = S_last
        else:
            print(f"  [!] Ratio NOT YET CONVERGED - need larger R")
            results['converged'] = False
    
    # Final verdict
    print(f"\n{'='*70}")
    print("FINAL VERDICT")
    print(f"{'='*70}")
    
    if results.get('converged', False):
        S = results['S_final']
        print(f"\n  Converged suppression ratio:")
        print(f"  S = A(E8) / A(Z^8) = {S:.6f}")
        print(f"\n  Comparison to φ powers:")
        print(f"    φ^(-1) = {PHI**(-1):.6f}")
        print(f"    φ^(-2) = {PHI**(-2):.6f}")
        print(f"    φ^(-3) = {PHI**(-3):.6f}")
        print(f"    φ^(-4) = {PHI**(-4):.6f}")
        print(f"    φ^(-6) = {PHI**(-6):.6f}")
        print(f"    φ^(-12) = {PHI**(-12):.6e}")
        
        # Find closest φ power
        for n in range(1, 13):
            phi_n = PHI**(-n)
            if abs(S - phi_n) / phi_n < 0.2:
                print(f"\n  >>> S ≈ φ^(-{n}) = {phi_n:.6f} (match within 20%)")
                results['phi_power'] = n
                break
    else:
        print(f"\n  Ratio not converged. Increase R_values for definitive result.")
    
    print(f"\n{'='*70}")
    
    return results


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # Start with smaller R values for speed, can increase
    results = run_rigorous_computation(R_values=[2.5, 3.0, 3.5, 4.0])
