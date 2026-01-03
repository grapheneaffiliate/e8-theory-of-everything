#!/usr/bin/env python3
"""
GSM SYMBOLIC DYNAMICS: THE BRIDGE
==================================
Converting integer cycle lengths to ln(p) structurally.

The critical missing piece: We have cycles of length {3,4,5,...}
but Riemann needs {ln(2), ln(3), ln(5), ln(7),...}

APPROACH: Define a NEW metric on E8 where cycle lengths ARE ln(p).

Author: Timothy McGirl  
Date: January 2, 2026
"""

import numpy as np
from itertools import product
from collections import defaultdict
import matplotlib.pyplot as plt

print("="*70)
print("GSM SYMBOLIC DYNAMICS: THE BRIDGE")
print("Converting Cycle Lengths to ln(p)")
print("="*70)

PHI = (1 + np.sqrt(5)) / 2

# ═══════════════════════════════════════════════════════════════════════════
# 1. PRIME NUMBER GENERATION
# ═══════════════════════════════════════════════════════════════════════════

def sieve_primes(n):
    """Generate primes up to n using Sieve of Eratosthenes."""
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n + 1) if is_prime[i]]

primes = sieve_primes(100)
ln_primes = [np.log(p) for p in primes]

print("\n[1] Prime Numbers and Their Logarithms")
print(f"    First 15 primes: {primes[:15]}")
print(f"    Their logs:      {[f'{x:.4f}' for x in ln_primes[:15]]}")

# ═══════════════════════════════════════════════════════════════════════════
# 2. E8 ROOT GENERATION (same as Ihara script)
# ═══════════════════════════════════════════════════════════════════════════

def generate_E8_roots():
    """Generate the 240 roots of the E8 lattice."""
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

print("\n[2] Generating E8 Lattice...")
roots = generate_E8_roots()
N = len(roots)
print(f"    {N} roots generated")

# ═══════════════════════════════════════════════════════════════════════════
# 3. THE KEY IDEA: PRIME-WEIGHTED METRIC
# ═══════════════════════════════════════════════════════════════════════════

print("\n[3] CONSTRUCTING PRIME-WEIGHTED METRIC")
print("""
    IDEA: Replace Euclidean distance with a metric where
    edge lengths ARE logarithms of primes.
    
    Standard E8: all edges have length √2 (roots at angle 120°)
    New metric:  assign edge (i,j) → ln(p_k) where k = f(i,j)
    
    The function f must be:
    - Symmetric: f(i,j) = f(j,i)
    - Consistent: cycles close with defined length
    - Structural: derives from E8 geometry (not arbitrary)
""")

# Build adjacency matrix
def build_adjacency(roots):
    n = len(roots)
    adj = []
    for i in range(n):
        for j in range(i+1, n):
            dot = np.dot(roots[i], roots[j])
            if abs(dot + 1) < 0.01:  # Adjacent if dot = -1
                adj.append((i, j))
    return adj

edges = build_adjacency(roots)
print(f"    Total edges: {len(edges)}")

# ═══════════════════════════════════════════════════════════════════════════
# 4. PRIME ORBIT LABELING SCHEME
# ═══════════════════════════════════════════════════════════════════════════

print("\n[4] PRIME ORBIT LABELING")
print("""
    CONSTRUCTION: The E8 graph has automorphism group W(E8) of order 696,729,600.
    This group acts transitively on:
    - Vertices (240 roots)
    - Edges (6720 edges)
    
    HYPOTHESIS: There exists a subset of primitive orbits that can be
    bijected to prime powers p^k such that:
    
    orbit_length(cycle_p) = k * ln(p)
    
    where cycle_p corresponds to prime p, and k is the multiplicity.
""")

# Count edge classes under rotational symmetry
# Simplified: group edges by their "type" based on vertex properties

def classify_edge(i, j, roots):
    """Classify edge by structural properties."""
    r_i, r_j = roots[i], roots[j]
    # Type based on whether vertices are "Type 1" (±1, ±1, 0...) or "Type 2" (±1/2...)
    type_i = 1 if np.max(np.abs(r_i)) == 1 else 2
    type_j = 1 if np.max(np.abs(r_j)) == 1 else 2
    return (min(type_i, type_j), max(type_i, type_j))

edge_classes = defaultdict(list)
for i, j in edges:
    cls = classify_edge(i, j, roots)
    edge_classes[cls].append((i, j))

print(f"    Edge classes found: {len(edge_classes)}")
for cls, edges_in_class in edge_classes.items():
    print(f"    Class {cls}: {len(edges_in_class)} edges")

# ═══════════════════════════════════════════════════════════════════════════
# 5. GOLDEN SYMBOLIC SHIFT SPACE
# ═══════════════════════════════════════════════════════════════════════════

print("\n[5] GOLDEN SYMBOLIC SHIFT SPACE")
print("""
    Define: Σ_E8 = {sequences (x_0, x_1, ...) in {0,1,...,239}^Z | consecutive symbols adjacent}
    
    The shift map σ: Σ_E8 → Σ_E8 given by σ(x)_n = x_{n+1}
    
    Periodic orbits of σ = primitive cycles in E8 graph
    
    KEY STEP: Define a CEILING FUNCTION that rounds cycle lengths:
    
    ceiling_p(ℓ) = ln(p) where p is the smallest prime with ln(p) ≥ ℓ/c
    
    for some constant c determined by E8 geometry.
""")

# Find cycle length distribution
def find_short_cycles(adj_list, max_len=6):
    """Find cycles up to given length."""
    n = 240
    # Build adjacency list
    neighbors = defaultdict(list)
    for i, j in adj_list:
        neighbors[i].append(j)
        neighbors[j].append(i)
    
    cycles_by_len = defaultdict(int)
    
    # DFS for cycles
    for start in range(n):
        stack = [(start, [start], {start})]
        while stack:
            node, path, visited = stack.pop()
            if len(path) > max_len:
                continue
            for neighbor in neighbors[node]:
                if neighbor == start and len(path) >= 3:
                    cycles_by_len[len(path)] += 1
                elif neighbor not in visited:
                    stack.append((neighbor, path + [neighbor], visited | {neighbor}))
    
    # Each cycle counted 2*len times (once per vertex, forward and backward)
    for length in cycles_by_len:
        cycles_by_len[length] //= (2 * length)
    
    return dict(cycles_by_len)

print("    Computing cycle distribution...")
cycle_dist = find_short_cycles(edges, max_len=6)
print(f"    Cycles found: {cycle_dist}")

# ═══════════════════════════════════════════════════════════════════════════
# 6. THE PRIME MAPPING FUNCTION
# ═══════════════════════════════════════════════════════════════════════════

print("\n[6] THE PRIME MAPPING FUNCTION")
print("""
    Given integer cycle length ℓ, map to ln(p) via:
    
    Option A (Linear): L(ℓ) = c * ℓ where c = ln(2)/3 ≈ 0.231
                       Then L(3) = ln(2), L(4) ≈ ln(2.52), L(5) ≈ ln(3.17)
    
    Option B (Exponential): L(ℓ) = ln(prime(ℓ-2))
                            Then L(3) = ln(2), L(4) = ln(3), L(5) = ln(5)
    
    Option B is exact but arbitrary.
    Option A is natural but doesn't hit primes exactly.
    
    THE INSIGHT: E8's 240 = 2^4 * 3 * 5 suggests Option B has structural roots!
    240 = 16 * 15 = F(12) * L(5) (Fibonacci-Lucas product nearby: F_12=144, L_5=11)
""")

# Option A: Linear scaling
c_linear = np.log(2) / 3

print("\n    Option A (Linear scaling c = ln(2)/3):")
for L in range(3, 10):
    mapped = c_linear * L
    closest_prime_idx = np.argmin([abs(np.log(p) - mapped) for p in primes[:20]])
    closest_prime = primes[closest_prime_idx]
    print(f"    L={L} → {mapped:.4f} ≈ ln({closest_prime}) = {np.log(closest_prime):.4f}")

# Option B: Direct prime mapping
print("\n    Option B (Direct prime mapping):")
for L in range(3, 10):
    if L - 2 < len(primes):
        p = primes[L - 3]  # L=3 → prime[0]=2, L=4 → prime[1]=3, etc.
        print(f"    L={L} → ln({p}) = {np.log(p):.4f}")

# ═══════════════════════════════════════════════════════════════════════════
# 7. SPECTRAL COMPARISON
# ═══════════════════════════════════════════════════════════════════════════

print("\n[7] SPECTRAL COMPARISON")
print("    Comparing trace formula with both mappings...")

# Build Hamiltonian and get eigenvalues
def build_hamiltonian(roots):
    n = len(roots)
    H = np.zeros((n, n))
    for i in range(n):
        H[i, i] = 2 * PHI
        for j in range(i+1, n):
            dot = np.dot(roots[i], roots[j])
            dist_sq = np.sum((roots[i] - roots[j])**2)
            val = dot * PHI**(-dist_sq)
            H[i, j] = val
            H[j, i] = val
    return H

H = build_hamiltonian(roots)
eigenvalues = np.linalg.eigvalsh(H)
eigenvalues = np.sort(eigenvalues[eigenvalues > 0.1])

# Compute spectral trace
t_vals = np.linspace(0.5, 3.5, 1000)
spectral_trace = np.zeros_like(t_vals)
for lam in eigenvalues:
    spectral_trace += np.cos(lam * t_vals)

# Compute geometric trace with prime mapping
prime_trace = np.zeros_like(t_vals)
for L, count in cycle_dist.items():
    if L >= 3:
        ln_p = np.log(primes[L - 3]) if L - 3 < len(primes) else 0
        amp = PHI**(-(2 * L))  # Simplified amplitude
        for i, t in enumerate(t_vals):
            prime_trace[i] += count * amp * np.cos(t * ln_p)

# Mark ln(p) positions
ln2, ln3, ln5, ln7 = np.log(2), np.log(3), np.log(5), np.log(7)

print(f"    ln(2) = {ln2:.4f}")
print(f"    ln(3) = {ln3:.4f}")
print(f"    ln(5) = {ln5:.4f}")
print(f"    ln(7) = {ln7:.4f}")

# ═══════════════════════════════════════════════════════════════════════════
# 8. PLOT
# ═══════════════════════════════════════════════════════════════════════════

print("\n[8] Generating Comparison Plot...")

plt.style.use('dark_background')
fig, axes = plt.subplots(2, 1, figsize=(12, 8))

# Spectral trace
ax1 = axes[0]
ax1.plot(t_vals, spectral_trace, color='cyan', lw=1.5, label='Spectral: Σcos(λt)')
for ln_p, p in [(ln2, 2), (ln3, 3), (ln5, 5), (ln7, 7)]:
    ax1.axvline(x=ln_p, color='magenta', linestyle='--', alpha=0.7, label=f'ln({p})' if p==2 else '')
ax1.set_title("Spectral Trace vs ln(p) Positions", color='gold', fontsize=14)
ax1.set_xlabel("t")
ax1.set_ylabel("Trace")
ax1.legend()
ax1.grid(True, alpha=0.2)

# Compare mappings
ax2 = axes[1]
for L in range(3, 8):
    # Linear mapping
    linear_pos = c_linear * L
    ax2.scatter([linear_pos], [L], color='orange', s=100, marker='o', label='Linear' if L==3 else '')
    # Prime mapping
    if L - 3 < len(primes):
        prime_pos = np.log(primes[L - 3])
        ax2.scatter([prime_pos], [L], color='lime', s=100, marker='^', label='Prime' if L==3 else '')
        ax2.annotate(f'L={L}→ln({primes[L-3]})', (prime_pos, L), xytext=(5, 5), 
                    textcoords='offset points', color='white', fontsize=9)

# Add ln(p) reference lines
for ln_p, p in [(ln2, 2), (ln3, 3), (ln5, 5), (ln7, 7), (np.log(11), 11)]:
    ax2.axvline(x=ln_p, color='cyan', linestyle=':', alpha=0.5)
    ax2.text(ln_p, 8.5, f'ln({p})', ha='center', color='cyan', fontsize=8)

ax2.set_title("Cycle Length → ln(p) Mappings", color='gold', fontsize=14)
ax2.set_xlabel("Mapped Length (log scale)")
ax2.set_ylabel("Original Cycle Length L")
ax2.set_xlim([0.5, 2.8])
ax2.set_ylim([2.5, 9])
ax2.legend()
ax2.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('Symbolic_Dynamics_Bridge.png', dpi=150)
print("    Plot saved to 'Symbolic_Dynamics_Bridge.png'")

# ═══════════════════════════════════════════════════════════════════════════
# 9. SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("SYMBOLIC DYNAMICS SUMMARY")
print("="*70)

print("""
THE BRIDGE PROBLEM:
  Current: Cycle lengths = {3, 4, 5, 6, ...}
  Needed:  Cycle lengths = {ln(2), ln(3), ln(5), ln(7), ...}

PROPOSED SOLUTIONS:

1. PRIME INDEXING (Option B):
   - Define L(ℓ) = ln(prime(ℓ-2))
   - L(3) = ln(2), L(4) = ln(3), L(5) = ln(5), ...
   - EXACT but requires justifying the index shift

2. LINEAR SCALING (Option A):
   - Define L(ℓ) = c * ℓ with c = ln(2)/3
   - Approximates primes but not exact
   - More "natural" continuous deformation

3. SYMBOLIC DYNAMICS (True Solution):
   - Redefine the metric on E8 graph
   - Edge weights become ln(p_e) for some assignment p_e
   - Cycle length = Σ(edge weights) must equal ln(p_cycle)
   - THIS REQUIRES MODULAR ARITHMETIC / NUMBER THEORY

THE HARD THEOREM:
  Prove there exists a consistent edge labeling of E8 such that
  primitive cycle lengths = {ln(p) : p prime}

This is the UNSOLVED PART of the Hilbert-Pólya approach.
""")

print("="*70)
print("STATUS: Bridge construction is THEORETICAL - needs proof")
print("="*70)
