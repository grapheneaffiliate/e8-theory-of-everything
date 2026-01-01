#!/usr/bin/env python3
"""
E8 INHERENT GOLDEN RATIO TEST
=============================

HYPOTHESIS: The icosahedral angle 1/sqrt5 = 0.447 is INHERENT to E8 geometry,
not imposed by the Elser-Sloane projection.

TEST: Take 1000 completely RANDOM 4x8 projections and compute <cos theta>
for E8 nearest-neighbor pairs. If this consistently gives ~0.447,
the golden ratio is TOPOLOGICAL, not dynamical.

Author: Timothy McGirl
Date: December 31, 2025
"""

import numpy as np
from scipy.linalg import qr
from itertools import product
from typing import Tuple, List, Dict

# Constants
PHI = (1 + np.sqrt(5)) / 2
TARGET_COS = 1.0 / np.sqrt(5)  # 0.447214
ICOSAHEDRAL_DIHEDRAL = np.degrees(np.arccos(TARGET_COS))  # 63.43 deg


def generate_e8_roots() -> np.ndarray:
    """Generate all 240 E8 root vectors."""
    roots = []
    
    # Type 1: 112 integer roots (+/-1, +/-1, 0...0)
    for i in range(8):
        for j in range(i + 1, 8):
            for s1, s2 in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                r = np.zeros(8)
                r[i], r[j] = s1, s2
                roots.append(r)
    
    # Type 2: 128 half-integer roots (+/-1/2)^8 with even minus signs
    for signs in product([0.5, -0.5], repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            roots.append(np.array(signs))
    
    return np.array(roots)


def find_e8_nearest_neighbors(roots: np.ndarray) -> List[Tuple[int, int]]:
    """
    Find nearest-neighbor pairs in E8 (distance = sqrt2).
    
    In E8, nearest neighbors have distance sqrt2 between them.
    This is a TOPOLOGICAL property of the root lattice.
    """
    n = len(roots)
    pairs = []
    target_dist = np.sqrt(2)
    
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(roots[i] - roots[j])
            if abs(d - target_dist) < 0.01:
                pairs.append((i, j))
    
    return pairs


def random_projection_4x8() -> np.ndarray:
    """Generate a random orthonormal 4x8 projection matrix."""
    # Generate random 8x4 matrix and take QR decomposition
    A = np.random.randn(8, 4)
    Q, R = qr(A, mode='economic')  # Q is 8x4, R is 4x4
    # Fix sign ambiguity
    signs = np.sign(np.diag(R))
    signs[signs == 0] = 1
    Q = Q * signs  # Broadcast correctly
    return Q.T  # Shape (4, 8)


def compute_cos_theta(P: np.ndarray, roots: np.ndarray, 
                       pairs: List[Tuple[int, int]]) -> float:
    """
    Compute mean cos theta between projected nearest-neighbor pairs.
    
    This is the KEY observable: what angle do E8 neighbors make
    when projected to 4D?
    """
    projected = P @ roots.T  # Shape (4, 240)
    
    cos_values = []
    for i, j in pairs:
        v_i = projected[:, i]
        v_j = projected[:, j]
        
        norm_i = np.linalg.norm(v_i)
        norm_j = np.linalg.norm(v_j)
        
        if norm_i > 1e-10 and norm_j > 1e-10:
            cos_theta = np.dot(v_i, v_j) / (norm_i * norm_j)
            cos_values.append(cos_theta)
    
    return np.mean(cos_values) if cos_values else 0.0


def run_inherent_test(n_trials: int = 1000) -> Dict:
    """
    Run the inherent golden ratio test.
    
    Hypothesis: <cos theta> ~ 1/sqrt5 = 0.447 for ANY random projection.
    """
    print("=" * 70)
    print("E8 INHERENT GOLDEN RATIO TEST")
    print("=" * 70)
    print()
    print("Hypothesis: The icosahedral angle 1/sqrt5 ~ 0.447 emerges from")
    print("            ANY random 4x8 projection of E8 nearest neighbors.")
    print()
    print("If true -> Golden ratio is TOPOLOGICAL (inherent to E8)")
    print("If false -> Golden ratio requires specific Elser-Sloane projection")
    print()
    
    # Generate E8 roots and find neighbors
    roots = generate_e8_roots()
    print(f"Generated {len(roots)} E8 roots")
    
    pairs = find_e8_nearest_neighbors(roots)
    print(f"Found {len(pairs)} nearest-neighbor pairs (distance = sqrt2)")
    print()
    
    # Run many random projections
    print(f"Testing {n_trials} random projections...")
    print()
    
    results = []
    checkpoints = [10, 50, 100, 500, 1000]
    
    for trial in range(n_trials):
        P = random_projection_4x8()
        cos_theta = compute_cos_theta(P, roots, pairs)
        results.append(cos_theta)
        
        if (trial + 1) in checkpoints and (trial + 1) <= n_trials:
            current_mean = np.mean(results)
            current_std = np.std(results)
            error = abs(current_mean - TARGET_COS) / TARGET_COS * 100
            print(f"  After {trial + 1:4d} trials: <cos theta> = {current_mean:.4f} +/- {current_std:.4f} "
                  f"(error from 1/sqrt5: {error:.2f}%)")
    
    # Final statistics
    results = np.array(results)
    mean_cos = np.mean(results)
    std_cos = np.std(results)
    mean_angle = np.degrees(np.arccos(np.clip(mean_cos, -1, 1)))
    error_pct = abs(mean_cos - TARGET_COS) / TARGET_COS * 100
    
    print()
    print("=" * 70)
    print("FINAL RESULTS")
    print("=" * 70)
    print()
    print(f"Number of random projections tested: {n_trials}")
    print(f"Number of E8 neighbor pairs: {len(pairs)}")
    print()
    print(f"Mean <cos theta> = {mean_cos:.6f} +/- {std_cos:.6f}")
    print(f"Target 1/sqrt5  = {TARGET_COS:.6f}")
    print(f"Error        = {error_pct:.2f}%")
    print()
    print(f"Mean angle   = {mean_angle:.2f} deg")
    print(f"Icosahedral dihedral = {ICOSAHEDRAL_DIHEDRAL:.2f} deg")
    print()
    
    # Verdict
    print("=" * 70)
    print("VERDICT")
    print("=" * 70)
    
    if error_pct < 5.0:
        verdict = "CONFIRMED"
        print()
        print("[OK] THE GOLDEN RATIO IS INHERENT TO E8")
        print()
        print("The icosahedral angle cos theta = 1/sqrt5 ~ 0.447 emerges from")
        print("E8 nearest-neighbor pairs under ANY random 4x8 projection.")
        print()
        print("This means:")
        print("  1. The H4 locking term is NOT arbitrary")
        print("  2. The golden ratio is TOPOLOGICAL, not dynamical")
        print("  3. The Elser-Sloane projection REVEALS structure already in E8")
        print()
        print("The icosahedral/quasicrystal geometry is a PROPERTY OF E8 ITSELF.")
    elif error_pct < 10.0:
        verdict = "SUGGESTIVE"
        print()
        print("[!] SUGGESTIVE: Close but not conclusive")
        print(f"   Error {error_pct:.1f}% is between 5% and 10%")
    else:
        verdict = "NOT_CONFIRMED"
        print()
        print("✗ Golden ratio does NOT emerge from random projections")
        print("   The Elser-Sloane projection is special/required")
    
    print("=" * 70)
    
    return {
        'n_trials': n_trials,
        'n_pairs': len(pairs),
        'mean_cos': mean_cos,
        'std_cos': std_cos,
        'target_cos': TARGET_COS,
        'error_pct': error_pct,
        'mean_angle_deg': mean_angle,
        'target_angle_deg': ICOSAHEDRAL_DIHEDRAL,
        'verdict': verdict,
        'all_results': results
    }


if __name__ == "__main__":
    np.random.seed(42)  # For reproducibility
    results = run_inherent_test(n_trials=1000)
    
    # Additional analysis: histogram of results
    print()
    print("Distribution of <cos theta> values:")
    hist, edges = np.histogram(results['all_results'], bins=20)
    max_count = max(hist)
    for i, count in enumerate(hist):
        bar = "█" * int(40 * count / max_count)
        print(f"  {edges[i]:.3f}-{edges[i+1]:.3f}: {bar} ({count})")
