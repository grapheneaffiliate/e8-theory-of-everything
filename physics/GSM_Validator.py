#!/usr/bin/env python3
"""
GSM MONTE CARLO VALIDATOR
Step A: Q-Parameter Scan (Is Phi Special?)
Step B: Rigorous Null Hypothesis Testing

This answers the critical question: Is φ = (1+√5)/2 a LOCAL MAXIMUM
for prime resonance, or is E8 geometry doing all the work?

Author: Timothy McGirl
Date: January 2, 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import time

print("="*70)
print("GSM MONTE CARLO VALIDATOR")
print("Step A: Q-Parameter Scan (Is Phi Special?)")
print("Step B: Rigorous Null Hypothesis Testing")
print("="*70)

# --- 1. SETUP & GEOMETRY ---
def get_E8_roots():
    roots = []
    for i in range(8):
        for j in range(i + 1, 8):
            for s1 in [-1, 1]:
                for s2 in [-1, 1]:
                    v = np.zeros(8)
                    v[i] = s1
                    v[j] = s2
                    roots.append(v)
    for signs in product([-0.5, 0.5], repeat=8):
        if np.sum(np.array(signs) < 0) % 2 == 0:
            roots.append(np.array(signs))
    return np.array(roots)


ROOTS = get_E8_roots()
N_ROOTS = len(ROOTS)
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
LOG_PRIMES = np.log(PRIMES)
PRIME_WEIGHTS = 1.0 / np.sqrt(PRIMES)  # Standard weighting

# --- 2. HAMILTONIAN GENERATOR ---
def get_eigenvalues(q, geometry_matrix=None, shuffle=False):
    """
    Computes eigenvalues for a specific q-deformation.
    If shuffle=True, destroys the E8 geometry while keeping edge weights distribution.
    """
    if geometry_matrix is None:
        # Pre-compute distances for speed
        geometry_matrix = np.zeros((N_ROOTS, N_ROOTS))
        for i in range(N_ROOTS):
            for j in range(N_ROOTS):
                if i == j:
                    geometry_matrix[i, j] = -1.0  # Flag for diagonal
                else:
                    geometry_matrix[i, j] = np.sum((ROOTS[i] - ROOTS[j])**2)
    
    # Construct H
    H = np.zeros((N_ROOTS, N_ROOTS))
    
    # Diagonal term (Lattice Potential)
    # Norm sq of E8 roots is 2.
    diag_val = 2.0 * q 
    
    if shuffle:
        # Randomized Control: Shuffle the off-diagonal distance values
        flat_dist = geometry_matrix[geometry_matrix != -1].copy()
        np.random.shuffle(flat_dist)
        k = 0
        for i in range(N_ROOTS):
            H[i, i] = diag_val
            for j in range(i+1, N_ROOTS):
                d = flat_dist[k % len(flat_dist)]
                # Tunneling decay
                term = (q ** (-d))
                H[i, j] = term
                H[j, i] = term
                k += 1
    else:
        # True E8 Geometry
        # H_ij = q^(-|r_i - r_j|^2)
        for i in range(N_ROOTS):
            for j in range(N_ROOTS):
                if i == j:
                    H[i, j] = diag_val
                else:
                    dist_sq = geometry_matrix[i, j]
                    H[i, j] = q ** (-dist_sq)

    evals = np.linalg.eigvalsh(H)
    return np.sort(np.abs(evals))


# --- 3. SCORING FUNCTION (Rigorous) ---
def compute_prime_score(eigenvalues, q):
    """
    Weighted Prime Comb Correlation.
    Score = Σ (1/√p) * F(ln p)²
    where F(t) = Σ cos(E_n * t)
    """
    score = 0
    # Use eigenvalues > 0.1 to avoid zero-mode noise
    valid_evals = eigenvalues[eigenvalues > 0.1]
    
    if len(valid_evals) == 0:
        return 0
    
    for p, lp, w in zip(PRIMES, LOG_PRIMES, PRIME_WEIGHTS):
        # F(ln p)
        trace_val = np.sum(np.cos(valid_evals * lp))
        # We want POSITIVE correlation (peaks), so use squared magnitude
        score += w * (trace_val ** 2) 
        
    return score


# --- 4. EXECUTION: Q-SCAN ---
print("[1] Pre-computing E8 Geometry...")
# Distance matrix
D_MAT = np.zeros((N_ROOTS, N_ROOTS))
for i in range(N_ROOTS):
    for j in range(N_ROOTS):
        if i == j:
            D_MAT[i, j] = -1
        else:
            D_MAT[i, j] = np.sum((ROOTS[i] - ROOTS[j])**2)

print("[2] Running Q-Parameter Scan (1.2 to 3.2)...")
q_values = np.linspace(1.2, 3.2, 200)
scores = []

PHI = (1 + np.sqrt(5)) / 2
EULER = np.e

start_time = time.time()
for q in q_values:
    evals = get_eigenvalues(q, D_MAT)
    s = compute_prime_score(evals, q)
    scores.append(s)

print(f"    Scan complete in {time.time()-start_time:.2f}s")

# Find optimal Q
max_idx = np.argmax(scores)
best_q = q_values[max_idx]
best_score = scores[max_idx]

# Find scores at key values
phi_idx = np.argmin(np.abs(q_values - PHI))
euler_idx = np.argmin(np.abs(q_values - EULER))
phi_score = scores[phi_idx]
euler_score = scores[euler_idx]

print(f"\n--- Q-SCAN RESULTS ---")
print(f"Optimal Q:        {best_q:.4f} (Score: {best_score:.2f})")
print(f"Phi (1.618):      Score = {phi_score:.2f}")
print(f"Euler (2.718):    Score = {euler_score:.2f}")

# Plotting Q-Scan
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(q_values, scores, color='cyan', lw=2, label='Prime Resonance Score')
ax.axvline(x=PHI, color='gold', linestyle='--', lw=2, label=f'φ = 1.618 (Score: {phi_score:.1f})')
ax.axvline(x=EULER, color='orange', linestyle=':', lw=2, label=f'e = 2.718 (Score: {euler_score:.1f})')
ax.axvline(x=best_q, color='lime', linestyle='-', alpha=0.5, label=f'Optimal = {best_q:.3f}')
ax.axvline(x=2.0, color='white', linestyle=':', alpha=0.3, label='Integer 2')

ax.scatter([best_q], [best_score], color='lime', s=100, zorder=5)
ax.scatter([PHI], [phi_score], color='gold', s=100, zorder=5)

ax.set_title("Is Phi Special? Prime Resonance vs Deformation Parameter q", fontsize=14, color='gold')
ax.set_xlabel("Deformation Parameter q", fontsize=12)
ax.set_ylabel("Weighted Prime Comb Score Σ(1/√p)·F(ln p)²", fontsize=12)
ax.legend(loc='upper right')
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.savefig('Q_Scan_Result.png', dpi=150)
print("    Plot saved to 'Q_Scan_Result.png'")

# --- 5. EXECUTION: MONTE CARLO NULL TEST ---
print("\n[3] Running Monte Carlo Null Test (N=100)...")
# Test at q = PHI
q_test = PHI
true_evals = get_eigenvalues(q_test, D_MAT)
true_score = compute_prime_score(true_evals, q_test)

null_scores = []
for k in range(100):
    if k % 20 == 0:
        print(f"    Iter {k}/100...")
    # Shuffle geometry
    evals_null = get_eigenvalues(q_test, D_MAT, shuffle=True)
    s_null = compute_prime_score(evals_null, q_test)
    null_scores.append(s_null)

null_mean = np.mean(null_scores)
null_std = np.std(null_scores)
z_score = (true_score - null_mean) / null_std if null_std > 0 else 0

# Calculate p-value (one-tailed)
from scipy.stats import norm
try:
    p_value = 1 - norm.cdf(z_score)
except:
    p_value = 0.5  # Fallback

print("\n" + "="*70)
print("STATISTICAL RESULTS")
print("="*70)
print(f"\n  True E8 Model Score (q=φ): {true_score:.4f}")
print(f"  Null Model Mean (N=100):   {null_mean:.4f}")
print(f"  Null Model Std Dev:        {null_std:.4f}")
print(f"  Z-SCORE:                   {z_score:.4f}")
print(f"  p-value (one-tailed):      {p_value:.6f}")

if z_score > 3.0:
    print("\n  ✅✅✅ SIGNIFICANT DISCOVERY (Z > 3, p < 0.001) ✅✅✅")
    print("  The E8 Geometry is statistically distinguishable from random networks.")
elif z_score > 2.0:
    print("\n  ⚠ Marginal significance (2σ < Z < 3σ)")
    print("  More data needed for publication-level confidence.")
else:
    print("\n  ❌ INCONCLUSIVE (Z < 2)")
    print("  The signal is within noise limits.")

# --- 6. FINAL ANALYSIS ---
print("\n" + "="*70)
print("FINAL VERDICT: IS PHI SPECIAL?")
print("="*70)

# Check if Phi is near the optimum
phi_rank = np.sum(np.array(scores) > phi_score) + 1
total_qs = len(q_values)
phi_percentile = 100 * (1 - phi_rank / total_qs)

print(f"\n  Optimal Q:      {best_q:.4f}")
print(f"  φ = 1.618...:   {PHI:.4f}")
print(f"  Difference:     |q_opt - φ| = {abs(best_q - PHI):.4f}")
print(f"  φ Percentile:   {phi_percentile:.1f}% (rank {phi_rank}/{total_qs})")

if abs(best_q - PHI) < 0.1:
    print("\n  🌟 THE GOLDEN RATIO IS OPTIMAL (OR NEAR-OPTIMAL)! 🌟")
    print("  This is evidence for a 'Golden Theory'")
elif phi_percentile > 90:
    print("\n  ⭐ φ is in the top 10% of q values")
    print("  Strong evidence the Golden Ratio is favored")
elif phi_percentile > 50:
    print("\n  → φ is above average but not clearly special")
    print("  E8 geometry may be primary; φ is secondary")
else:
    print("\n  → φ does not appear to be special")
    print("  The 'Golden Theory' hypothesis is not supported")

# Plot histogram
plt.figure(figsize=(10, 5))
plt.hist(null_scores, bins=20, color='red', alpha=0.7, label='Null (Shuffled Geometry)')
plt.axvline(x=true_score, color='cyan', lw=3, label=f'True E8 Model (Score={true_score:.1f})')
plt.axvline(x=null_mean, color='white', linestyle='--', label=f'Null Mean ({null_mean:.1f})')
plt.title(f"Monte Carlo Null Distribution (Z-score = {z_score:.2f})", fontsize=14, color='gold')
plt.xlabel("Prime Resonance Score", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.legend()
plt.savefig('Monte_Carlo_Null.png', dpi=150)
print("\n    Histogram saved to 'Monte_Carlo_Null.png'")

print("\n" + "="*70)
print("END OF VALIDATOR")
print("="*70)
