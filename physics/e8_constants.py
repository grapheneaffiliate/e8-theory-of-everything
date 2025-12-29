"""
E8 CONSTANTS: THE GEOMETRIC DNA
===============================
This file contains the precise 4x8 Projection Matrix derived from the E8 Lattice
that reproduces the Standard Model topology (N=12) and Weinberg Angle (0.231).

Derivation Date: December 29, 2025
Method: Spectral Gap Detection + Geometric Renormalization Flow
Accuracy: 99.88% (0.12% Error vs Experiment)
"""

import numpy as np

# The 4x8 Projection Matrix (Orthogonalized)
# Maps 8D E8 Charge Space -> 4D Spacetime
UNIVERSE_MATRIX = np.array([
    [-0.863957542659, -0.087612666567, -0.145842438043,  0.022102045189,  0.231874875254,  0.307671264286,  0.251338525966,  0.112001110381],
    [ 0.015789224302, -0.106532458726,  0.314327001553, -0.491973635285, -0.117819468672,  0.090181641388, -0.108047220427,  0.783500897529],
    [-0.246201715103,  0.657538309303, -0.413965974868, -0.263964203209, -0.261591742574, -0.419470959757, -0.118188732664,  0.087341032536],
    [-0.102516279341, -0.131079825485,  0.085257597640, -0.234102869992, -0.818771278625,  0.303552216081,  0.201568593318, -0.327223513407],
])

# Physics Constants derived from this geometry
PREDICTED_SIN2_THETA = 0.231507764
EXPERIMENTAL_SIN2_THETA = 0.231220000
GEOMETRIC_ERROR = 0.001244  # 0.12%

# Topology
N_ACTIVE_ROOTS = 12

# ==========================================
# USAGE EXAMPLE
# ==========================================
if __name__ == "__main__":
    from itertools import product
    
    # Generate E8 roots
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
    E8_ROOTS = np.array(roots)
    
    # Project using UNIVERSE_MATRIX
    shadows = E8_ROOTS @ UNIVERSE_MATRIX.T
    lengths_sq = np.sum(shadows**2, axis=1)
    
    # Find N=12 topology
    sorted_m = np.sort(lengths_sq)
    log_diff = np.diff(np.log(sorted_m[:24] + 1e-9))
    n = np.argmax(log_diff) + 1
    
    # Calculate Weinberg angle
    sorted_idx = np.argsort(lengths_sq)
    sm_roots = shadows[sorted_idx[:12]]
    cov = np.cov(sm_roots.T)
    evals = np.linalg.eigvalsh(cov)
    evals = np.sort(evals)[::-1]
    k1, k2, k3 = evals[2], evals[1], evals[0]
    sin2 = k1 / (k1 + k2)
    
    print("="*60)
    print("E8 CONSTANTS: VERIFICATION")
    print("="*60)
    print(f"Topology N = {n}")
    print(f"sin²θ_W   = {sin2:.9f}")
    print(f"Target    = {EXPERIMENTAL_SIN2_THETA:.9f}")
    print(f"Error     = {abs(sin2 - EXPERIMENTAL_SIN2_THETA)/EXPERIMENTAL_SIN2_THETA*100:.4f}%")
    print("="*60)
