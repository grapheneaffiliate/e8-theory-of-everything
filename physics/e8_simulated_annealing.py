"""
DEFINITIVE TEST: Simulated Annealing on E8 Stiefel Manifold
============================================================

Tests whether golden ratio geometry emerges SPONTANEOUSLY from E8
without any forcing term (μ = 0).

Action: S[P] = λ Σᵣ ||P·r||⁴  (pure quartic, NO locking)

Question: At the minimum, what is ⟨cos θ⟩ for nearest neighbors?

If ⟨cos θ⟩_min ≈ 0.447 (1/√5) or 0.618 (1/φ):
   → Golden ratio geometry is DYNAMICALLY FAVORED by E8
   
If ⟨cos θ⟩_min ≈ 0 or random:
   → H4 structure requires locking as INPUT

Author: Timothy McGirl
Date: December 31, 2025
"""

import numpy as np
from scipy.linalg import qr, expm
from itertools import product

# Constants
PHI = (1 + np.sqrt(5)) / 2
TARGET_1_SQRT5 = 1.0 / np.sqrt(5)  # 0.447214
TARGET_1_PHI = 1.0 / PHI  # 0.618034
TARGET_PHI_2 = PHI / 2  # 0.809017

# Generate E8 roots
def generate_e8_roots():
    """Generate the 240 roots of E8."""
    roots = []
    for i in range(8):
        for j in range(i + 1, 8):
            for s1, s2 in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                r = np.zeros(8)
                r[i], r[j] = s1, s2
                roots.append(r)
    for signs in product([0.5, -0.5], repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            roots.append(np.array(signs))
    return np.array(roots)

E8_ROOTS = generate_e8_roots()


def random_stiefel():
    """Random point on V₄(ℝ⁸)."""
    A = np.random.randn(4, 8)
    P, _ = qr(A.T)
    return P.T[:4, :]


def pure_quartic_action(P, roots, lambd=1.0):
    """
    PURE quartic action - NO locking term.
    
    S[P] = λ Σᵣ ||P·r||⁴
    
    This tests whether golden ratio emerges from geometry alone.
    """
    projected = P @ roots.T  # 4 x 240
    norms_sq = np.sum(projected**2, axis=0)
    return lambd * np.sum(norms_sq**2)


def geodesic_step(P, step_size=0.1):
    """Move along geodesic on Stiefel manifold."""
    A = np.random.randn(8, 8) * step_size
    A = (A - A.T) / 2  # Skew-symmetric
    exp_A = expm(A)
    P_new = P @ exp_A
    P_new, _ = qr(P_new.T)
    return P_new.T[:4, :]


def measure_cos_theta(P, roots, n_neighbors=6):
    """Measure mean cos θ for nearest neighbor pairs."""
    projected = P @ roots.T
    n_roots = len(roots)
    
    cos_thetas = []
    for i in range(min(120, n_roots)):
        distances = np.linalg.norm(roots - roots[i], axis=1)
        nearest_idx = np.argsort(distances)[1:n_neighbors+1]
        
        for j in nearest_idx:
            if i < j:
                v_i = projected[:, i]
                v_j = projected[:, j]
                norm_i = np.linalg.norm(v_i) + 1e-10
                norm_j = np.linalg.norm(v_j) + 1e-10
                cos_theta = np.dot(v_i, v_j) / (norm_i * norm_j)
                cos_thetas.append(cos_theta)
    
    return np.mean(cos_thetas) if cos_thetas else 0.0


def simulated_annealing(P_init, roots, T_start=10.0, T_end=0.01, 
                        n_steps=5000, step_size=0.1, lambd=1.0):
    """
    Simulated annealing to find minimum of PURE quartic action.
    
    Temperature schedule: exponential cooling
    """
    P = P_init.copy()
    S = pure_quartic_action(P, roots, lambd)
    
    P_best = P.copy()
    S_best = S
    
    n_accept = 0
    cooling_rate = (T_end / T_start) ** (1.0 / n_steps)
    T = T_start
    
    for step in range(n_steps):
        # Propose move
        P_new = geodesic_step(P, step_size)
        S_new = pure_quartic_action(P_new, roots, lambd)
        
        # Metropolis criterion
        delta_S = S_new - S
        if delta_S < 0 or np.random.rand() < np.exp(-delta_S / T):
            P = P_new
            S = S_new
            n_accept += 1
            
            if S < S_best:
                P_best = P.copy()
                S_best = S
        
        # Cool down
        T *= cooling_rate
    
    return P_best, S_best, n_accept / n_steps


def run_test(n_runs=10):
    """
    Run simulated annealing from multiple random starts.
    
    If all converge to same ⟨cos θ⟩ ≈ golden ratio: STRONG EVIDENCE
    If they diverge or give ≈ 0: THEORY NEEDS LOCKING TERM
    """
    print("=" * 70)
    print("DEFINITIVE TEST: Does Golden Ratio Emerge from Pure E8 Geometry?")
    print("=" * 70)
    print()
    print("Action: S[P] = λ Σᵣ ||P·r||⁴  (NO locking term)")
    print("Question: At minimum, is ⟨cos θ⟩ ≈ golden ratio value?")
    print()
    print("=" * 70)
    print(f"Running {n_runs} independent optimizations from random starts...")
    print("=" * 70)
    print()
    
    roots = E8_ROOTS
    results = []
    
    for run in range(n_runs):
        print(f"Run {run + 1}/{n_runs}:", end=" ", flush=True)
        
        P_init = random_stiefel()
        P_min, S_min, accept_rate = simulated_annealing(
            P_init, roots, 
            T_start=10.0, T_end=0.001,
            n_steps=10000, step_size=0.08, lambd=1.0
        )
        
        cos_theta = measure_cos_theta(P_min, roots)
        results.append({
            'cos_theta': cos_theta,
            'action': S_min,
            'accept_rate': accept_rate
        })
        
        print(f"⟨cos θ⟩ = {cos_theta:.4f}, S = {S_min:.2f}, accept = {accept_rate:.1%}")
    
    # Analyze convergence
    cos_values = [r['cos_theta'] for r in results]
    mean_cos = np.mean(cos_values)
    std_cos = np.std(cos_values)
    
    print()
    print("=" * 70)
    print("RESULTS:")
    print("=" * 70)
    print(f"Mean ⟨cos θ⟩: {mean_cos:.4f} ± {std_cos:.4f}")
    print()
    print("Comparison to golden ratio values:")
    print(f"  1/√5 = {TARGET_1_SQRT5:.4f}  →  error = {abs(mean_cos - TARGET_1_SQRT5):.4f} ({abs(mean_cos - TARGET_1_SQRT5)/TARGET_1_SQRT5*100:.1f}%)")
    print(f"  1/φ  = {TARGET_1_PHI:.4f}  →  error = {abs(mean_cos - TARGET_1_PHI):.4f} ({abs(mean_cos - TARGET_1_PHI)/TARGET_1_PHI*100:.1f}%)")
    print(f"  φ/2  = {TARGET_PHI_2:.4f}  →  error = {abs(mean_cos - TARGET_PHI_2):.4f} ({abs(mean_cos - TARGET_PHI_2)/TARGET_PHI_2*100:.1f}%)")
    print()
    
    # Determine verdict
    convergent = std_cos < 0.05  # All runs converge to same value
    
    errors = {
        '1/√5': abs(mean_cos - TARGET_1_SQRT5) / TARGET_1_SQRT5,
        '1/φ': abs(mean_cos - TARGET_1_PHI) / TARGET_1_PHI,
        'φ/2': abs(mean_cos - TARGET_PHI_2) / TARGET_PHI_2
    }
    best_match = min(errors, key=errors.get)
    best_error = errors[best_match]
    
    print("=" * 70)
    print("VERDICT:")
    print("=" * 70)
    
    if not convergent:
        print("✗ INCONCLUSIVE: Runs did not converge to same value")
        print(f"  Standard deviation: {std_cos:.4f} (need < 0.05)")
        print("  → Theory neither confirmed nor refuted")
        verdict = "INCONCLUSIVE"
    elif best_error < 0.05:
        print(f"✓ STRONG EVIDENCE: Golden ratio ({best_match}) emerges spontaneously!")
        print(f"  All {n_runs} runs converged to ⟨cos θ⟩ ≈ {best_match}")
        print(f"  Error: {best_error*100:.1f}% (< 5% threshold)")
        print("  → E8 geometry DYNAMICALLY FAVORS icosahedral structure")
        verdict = "CONFIRMED"
    elif best_error < 0.10:
        print(f"⚠ SUGGESTIVE: Converges near {best_match} but not within 5%")
        print(f"  Error: {best_error*100:.1f}%")
        print("  → Worth investigating further")
        verdict = "SUGGESTIVE"
    else:
        print("✗ NEGATIVE: No convergence to golden ratio values")
        print(f"  Closest match: {best_match} with {best_error*100:.1f}% error")
        print("  → H4 structure requires locking term as INPUT")
        verdict = "NEGATIVE"
    
    print("=" * 70)
    
    return results, verdict


if __name__ == "__main__":
    np.random.seed(42)  # For reproducibility
    results, verdict = run_test(n_runs=10)
