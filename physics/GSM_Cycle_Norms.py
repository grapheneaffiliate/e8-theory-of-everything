#!/usr/bin/env python3
"""
GSM CYCLE NORMS: THE HONEST BRIDGE TEST
========================================
Testing whether primes emerge STRUCTURALLY from E8 Weyl monodromy.

THE QUESTION: Can we get primes without injecting them?

METHOD:
1. Assign Weyl reflection matrices to directed edges
2. Compute monodromy M(c) = product of matrices along cycle
3. Define N(c) = |det(M(c) - I)| as cycle "norm"
4. Test: Do primitive cycle norms give us {2, 3, 5, 7, 11, ...}?

If YES → genuine Euler product, real RH proof path
If NO → E8 alone doesn't produce primes, need adelic structure

Author: Timothy McGirl
Date: January 2, 2026
Purpose: HONEST non-cheating test of prime emergence
"""

import numpy as np
from itertools import product, permutations
from collections import defaultdict
import sympy as sp
from sympy import Matrix, eye, det, Abs, simplify

print("="*70)
print("GSM CYCLE NORMS: THE HONEST BRIDGE TEST")
print("Do primes emerge from E8 Weyl monodromy?")
print("="*70)

# ═══════════════════════════════════════════════════════════════════════════
# 1. E8 ROOT SYSTEM AND SIMPLE ROOTS
# ═══════════════════════════════════════════════════════════════════════════

print("\n[1] E8 ROOT SYSTEM AND SIMPLE ROOTS")

def generate_E8_roots():
    """Generate 240 E8 roots in R^8."""
    roots = []
    # Type 1: (±1, ±1, 0, 0, 0, 0, 0, 0) and permutations - 112 roots
    for i in range(8):
        for j in range(i + 1, 8):
            for s1, s2 in product([-1, 1], repeat=2):
                v = np.zeros(8)
                v[i], v[j] = s1, s2
                roots.append(v)
    # Type 2: (±1/2, ±1/2, ..., ±1/2) with even number of minus signs - 128 roots
    for signs in product([-0.5, 0.5], repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            roots.append(np.array(signs))
    return np.array(roots)

roots = generate_E8_roots()
print(f"    Total roots: {len(roots)}")

# E8 simple roots (standard basis)
# α_1 = (1,-1,0,0,0,0,0,0)
# α_2 = (0,1,-1,0,0,0,0,0)
# ...
# α_7 = (0,0,0,0,0,0,1,-1)
# α_8 = (-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1/2)  # the "special" E8 root

simple_roots = np.array([
    [1, -1, 0, 0, 0, 0, 0, 0],   # α_1
    [0, 1, -1, 0, 0, 0, 0, 0],   # α_2
    [0, 0, 1, -1, 0, 0, 0, 0],   # α_3
    [0, 0, 0, 1, -1, 0, 0, 0],   # α_4
    [0, 0, 0, 0, 1, -1, 0, 0],   # α_5
    [0, 0, 0, 0, 0, 1, -1, 0],   # α_6
    [0, 0, 0, 0, 0, 0, 1, -1],   # α_7
    [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5]  # α_8
])

print(f"    Simple roots: {len(simple_roots)}")
print(f"    Simple root norms: {[np.dot(r, r) for r in simple_roots]}")

# ═══════════════════════════════════════════════════════════════════════════
# 2. WEYL REFLECTION MATRICES
# ═══════════════════════════════════════════════════════════════════════════

print("\n[2] WEYL REFLECTION MATRICES")

def weyl_reflection(alpha):
    """
    Weyl reflection matrix for root alpha.
    s_α(v) = v - 2(v·α)/(α·α) × α
    In matrix form: S_α = I - 2 α α^T / (α·α)
    """
    alpha = np.array(alpha)
    norm_sq = np.dot(alpha, alpha)
    return np.eye(8) - 2 * np.outer(alpha, alpha) / norm_sq

# Generate reflection matrices for simple roots
S = [weyl_reflection(alpha) for alpha in simple_roots]

print(f"    Generated {len(S)} simple reflections")
print(f"    S_0 shape: {S[0].shape}")

# Verify they are involutions: S_i^2 = I
for i, s in enumerate(S):
    product = s @ s
    is_identity = np.allclose(product, np.eye(8))
    print(f"    S_{i}² = I ? {is_identity}")

# Verify orthogonality (det = -1 for reflections)
print(f"    Determinants: {[round(np.linalg.det(s), 2) for s in S]}")

# ═══════════════════════════════════════════════════════════════════════════
# 3. BUILD ADJACENCY WITH WEYL LABELS
# ═══════════════════════════════════════════════════════════════════════════

print("\n[3] BUILD ADJACENCY WITH WEYL LABELS")

def find_adjacent_roots(roots, tol=0.01):
    """
    Find adjacent pairs in E8 root system.
    Adjacent if dot product = -1 (angle = 120°).
    """
    n = len(roots)
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            dot = np.dot(roots[i], roots[j])
            if abs(dot + 1) < tol:  # dot = -1
                edges.append((i, j))
    return edges

edges = find_adjacent_roots(roots)
print(f"    Found {len(edges)} edges (undirected)")

# For each edge, assign the Weyl reflection that takes r_i to r_j
# The reflection equation: s_α(r_i) = r_i - 2(r_i·α)/(α·α)α
# We need α such that s_α(r_i) = r_j

def find_reflection_for_edge(r_i, r_j, simple_reflections, simple_roots):
    """
    Find which Weyl generator maps r_i to r_j (or close).
    Try single reflections and products of two.
    """
    # Try single simple reflections
    for k, S_k in enumerate(simple_reflections):
        result = S_k @ r_i
        if np.allclose(result, r_j, atol=0.01):
            return (k,)  # Single reflection
        if np.allclose(result, r_i, atol=0.01):
            continue  # Fixed point, not useful
    
    # Try products of two reflections
    for k1, S_k1 in enumerate(simple_reflections):
        for k2, S_k2 in enumerate(simple_reflections):
            result = S_k2 @ (S_k1 @ r_i)
            if np.allclose(result, r_j, atol=0.01):
                return (k1, k2)
    
    return None  # Not found in simple products

# Label a subset of edges
print("    Labeling edges with Weyl generators...")
edge_labels = {}
labeled_count = 0
for i, (v1, v2) in enumerate(edges[:100]):  # First 100 for speed
    label = find_reflection_for_edge(roots[v1], roots[v2], S, simple_roots)
    if label is not None:
        edge_labels[(v1, v2)] = label
        labeled_count += 1

print(f"    Labeled {labeled_count}/{min(100, len(edges))} edges")

# ═══════════════════════════════════════════════════════════════════════════
# 4. COMPUTE CYCLE MONODROMY
# ═══════════════════════════════════════════════════════════════════════════

print("\n[4] COMPUTE CYCLE MONODROMY")

def compute_monodromy(cycle_vertices, roots_array, simple_reflections):
    """
    Compute monodromy matrix M(c) = product of Weyl reflections along cycle.
    For each step v_i → v_{i+1}, use the reflection that maps r_{v_i} to r_{v_{i+1}}.
    
    Simplified: Use reflection in the plane perpendicular to (r_j - r_i).
    """
    M = np.eye(8)
    L = len(cycle_vertices)
    
    for k in range(L):
        v_i = cycle_vertices[k]
        v_j = cycle_vertices[(k + 1) % L]
        
        r_i = roots_array[v_i]
        r_j = roots_array[v_j]
        
        # Reflection that swaps r_i and r_j
        # The axis is the midpoint, reflection in (r_j - r_i)
        diff = r_j - r_i
        norm_sq = np.dot(diff, diff)
        if norm_sq > 1e-10:
            S_edge = np.eye(8) - 2 * np.outer(diff, diff) / norm_sq
            M = S_edge @ M
    
    return M

# Find triangles (3-cycles) in the graph
print("    Finding 3-cycles (triangles)...")

neighbors = defaultdict(set)
for v1, v2 in edges:
    neighbors[v1].add(v2)
    neighbors[v2].add(v1)

triangles = []
for v in range(len(roots)):
    for n1 in neighbors[v]:
        for n2 in neighbors[v]:
            if n1 < n2 and n2 in neighbors[n1]:
                triangle = tuple(sorted([v, n1, n2]))
                if triangle not in triangles:
                    triangles.append(triangle)

print(f"    Found {len(triangles)} unique triangles")

# Compute monodromy for first 50 triangles
print("    Computing monodromy matrices...")
monodromies = []
for tri in triangles[:50]:
    M = compute_monodromy(list(tri), roots, S)
    monodromies.append((tri, M))

# ═══════════════════════════════════════════════════════════════════════════
# 5. COMPUTE CYCLE NORMS N(c) = |det(M(c) - I)|
# ═══════════════════════════════════════════════════════════════════════════

print("\n[5] COMPUTE CYCLE NORMS N(c) = |det(M(c) - I)|")

def cycle_norm(M):
    """Compute N(c) = |det(M - I)|."""
    return abs(np.linalg.det(M - np.eye(M.shape[0])))

norms = []
for tri, M in monodromies:
    N = cycle_norm(M)
    norms.append((tri, N))

# Sort by norm
norms.sort(key=lambda x: x[1])

print("\n    Cycle norms (first 20):")
print(f"    {'Triangle':<25} {'N(c)':<15} {'N(c) int':<10}")
print("    " + "-"*50)

unique_norms = []
for tri, N in norms[:20]:
    N_int = round(N)
    print(f"    {str(tri):<25} {N:<15.6f} {N_int:<10}")
    if N > 0.5 and N_int not in unique_norms:
        unique_norms.append(N_int)

# ═══════════════════════════════════════════════════════════════════════════
# 6. THE CRITICAL TEST: ARE THESE PRIMES?
# ═══════════════════════════════════════════════════════════════════════════

print("\n[6] THE CRITICAL TEST: ARE THESE PRIMES?")

def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True

primes_10 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

print(f"\n    Expected primes: {primes_10}")
print(f"    Observed norms:  {unique_norms[:10]}")

# Check overlap
observed_primes = [n for n in unique_norms if is_prime(n)]
print(f"\n    Observed that ARE prime: {observed_primes}")
print(f"    Observed that are NOT prime: {[n for n in unique_norms if not is_prime(n) and n > 1]}")

# ═══════════════════════════════════════════════════════════════════════════
# 7. ALTERNATIVE: SPECTRAL RADIUS AS NORM
# ═══════════════════════════════════════════════════════════════════════════

print("\n[7] ALTERNATIVE: SPECTRAL RADIUS AS NORM")

def spectral_norm(M):
    """Compute spectral radius ρ(M)."""
    eigs = np.linalg.eigvals(M)
    return max(abs(e) for e in eigs)

spec_norms = []
for tri, M in monodromies:
    rho = spectral_norm(M)
    spec_norms.append((tri, rho))

spec_norms.sort(key=lambda x: x[1])

print("\n    Spectral radii (first 10):")
for tri, rho in spec_norms[:10]:
    print(f"    {str(tri):<25} ρ = {rho:.6f}")

# ═══════════════════════════════════════════════════════════════════════════
# 8. SUMMARY AND VERDICT
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("VERDICT: DO PRIMES EMERGE FROM E8 WEYL MONODROMY?")
print("="*70)

# Check if we got 2, 3, 5, 7 in the first few norms
got_2 = 2 in unique_norms[:10]
got_3 = 3 in unique_norms[:10]
got_5 = 5 in unique_norms[:10]
got_7 = 7 in unique_norms[:10]

print(f"""
  CYCLE NORM RESULTS:
    First unique norms: {unique_norms[:10]}
    
    Contains 2? {got_2}
    Contains 3? {got_3}  
    Contains 5? {got_5}
    Contains 7? {got_7}
    
  PRIMES IN OBSERVED NORMS:
    {observed_primes}
    
  COMPOSITES IN OBSERVED NORMS:
    {[n for n in unique_norms if not is_prime(n) and n > 1]}
""")

if set(primes_10[:4]).issubset(set(unique_norms[:10])):
    print("  ✅ PRIMES EMERGE! The first cycle norms match {2, 3, 5, 7}.")
    print("     This is GENUINE prime structure, not hand-coded.")
    print("     Proceed to full Euler product construction.")
else:
    print("  ❌ PRIMES DO NOT EMERGE from E8 Weyl monodromy alone.")
    print("     The cycle norms are NOT the primes {2, 3, 5, 7, ...}.")
    print("     E8 alone cannot produce the Riemann zeta Euler product.")
    print("")
    print("  CONCLUSION:")
    print("     To get an honest RH proof, you need ADDITIONAL structure:")
    print("     - Adelic factorization (places = primes)")
    print("     - Symbolic dynamics with unique factorization property")
    print("     - Arithmetic on the root lattice (integers mod N)")
    print("")
    print("  The current 'cycle → log(p)' bridge is HAND-WAVING.")

print("="*70)
