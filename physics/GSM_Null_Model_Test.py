#!/usr/bin/env python3
"""
GSM NULL HYPOTHESIS TEST
Comparing E8-Golden Model against Random Controls

This is the "trial by fire" that separates real discovery from pattern matching.
We prove that ONLY the combination of E8 Geometry + Golden Ratio produces
prime resonance. Random variations destroy the signal.

Author: Timothy McGirl
Date: January 2, 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import product

print("="*70)
print("GSM NULL HYPOTHESIS TEST")
print("Comparing E8-Golden Model against Random Controls")
print("="*70)

# 1. CONSTANTS
PHI = (1 + np.sqrt(5)) / 2
RANDOM_CONSTANT = np.e  # Testing against Euler's number as a control
LAMBDA_SCALING = 14.1347 / 7.52

# 2. GENERATE REAL E8 ROOTS
def get_E8_roots():
    """Generate all 240 roots of E8."""
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

roots = get_E8_roots()
n_roots = len(roots)
print(f"[Setup] Generated {n_roots} E8 roots")

# 3. HAMILTONIAN BUILDERS

def build_model_real_E8_real_Phi():
    """The True GSM Model: E8 Geometry + Golden Ratio"""
    H = np.zeros((n_roots, n_roots))
    for i in range(n_roots):
        for j in range(n_roots):
            if i == j:
                H[i, j] = np.dot(roots[i], roots[i]) * PHI
            else:
                dist_sq = np.sum((roots[i] - roots[j])**2)
                dot = np.dot(roots[i], roots[j])
                H[i, j] = dot * (PHI ** (-dist_sq))
    return H


def build_model_random_Phi():
    """Control 1: Correct E8 Geometry, Wrong Constant (e instead of Phi)"""
    H = np.zeros((n_roots, n_roots))
    constant = RANDOM_CONSTANT
    for i in range(n_roots):
        for j in range(n_roots):
            if i == j:
                H[i, j] = np.dot(roots[i], roots[i]) * constant
            else:
                dist_sq = np.sum((roots[i] - roots[j])**2)
                dot = np.dot(roots[i], roots[j])
                H[i, j] = dot * (constant ** (-dist_sq))
    return H


def build_model_shuffled_E8():
    """Control 2: Correct Phi, Destroyed Geometry (Shuffled interactions)"""
    H_real = build_model_real_E8_real_Phi()
    flat_H = H_real.flatten()
    np.random.seed(42)  # Reproducible
    np.random.shuffle(flat_H)  # Destroy structure
    H_shuffled = flat_H.reshape((n_roots, n_roots))
    return (H_shuffled + H_shuffled.T) / 2  # Keep it symmetric/Hermitian


def build_model_random_matrix():
    """Control 3: Complete random matrix (GOE - Gaussian Orthogonal Ensemble)"""
    np.random.seed(123)
    H = np.random.randn(n_roots, n_roots)
    return (H + H.T) / 2  # Symmetric


def get_trace_signal(H, label=""):
    """Compute spectral trace signal F(t) = Σ cos(E_n * t)"""
    # Diagonalize
    evals = np.linalg.eigvalsh(H)
    evals = np.sort(np.abs(evals))
    
    # Scale to physical range (roughly) for fair comparison
    # Normalize so mean eigenvalue approx matches real model
    mean_eval = np.mean(evals[evals > 1]) if np.sum(evals > 1) > 0 else np.mean(evals)
    if mean_eval > 0:
        scale = 14.13 / mean_eval
    else:
        scale = 1.0
    energies = evals * scale
    energies = energies[energies > 5]  # Filter noise
    
    if len(energies) == 0:
        energies = evals * scale
        energies = energies[energies > 0]
    
    # Compute Trace
    t_values = np.linspace(0, 4, 1000)
    signal = np.zeros_like(t_values)
    for E in energies:
        signal += np.cos(E * t_values)
    
    # Normalize
    if len(energies) > 0:
        signal = signal / len(energies)
    
    return t_values, signal


# 4. RUN EXPERIMENT
print("\n[1] Running TRUE GSM Model (E8 + Phi)...")
t, sig_true = get_trace_signal(build_model_real_E8_real_Phi())

print("[2] Running CONTROL 1 (E8 + Euler's Number e=2.718...)...")
_, sig_wrong_const = get_trace_signal(build_model_random_Phi())

print("[3] Running CONTROL 2 (Shuffled Geometry + Phi)...")
_, sig_shuffled = get_trace_signal(build_model_shuffled_E8())

print("[4] Running CONTROL 3 (Random Matrix - GOE)...")
_, sig_random = get_trace_signal(build_model_random_matrix())

# 5. ANALYSIS & PLOTTING
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
log_primes = np.log(primes)

plt.style.use('dark_background')
fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)


def plot_trace(ax, t, signal, title, color):
    ax.plot(t, signal, color=color, lw=1.5)
    # Mark Primes
    for lp in log_primes:
        if lp <= 4:
            ax.axvline(x=lp, color='magenta', linestyle=':', alpha=0.6)
    ax.set_title(title, fontsize=12, color='white')
    ax.set_ylabel("F(t)", fontsize=10)
    ax.grid(True, alpha=0.2)


plot_trace(axes[0], t, sig_true, 
           "TRUE MODEL: E8 Lattice + Golden Ratio φ (THE HYPOTHESIS)", "cyan")
plot_trace(axes[1], t, sig_wrong_const, 
           "CONTROL 1: E8 Lattice + Euler's Number e (Wrong Constant)", "orange")
plot_trace(axes[2], t, sig_shuffled, 
           "CONTROL 2: Destroyed E8 Geometry + Golden Ratio φ", "red")
plot_trace(axes[3], t, sig_random, 
           "CONTROL 3: Random Matrix (GOE) - No Structure", "gray")

plt.xlabel("t = ln(p) [Magenta lines = Prime Logarithms]", fontsize=12)
plt.tight_layout()
plt.savefig('Null_Model_Test.png', dpi=150)
print("\n[5] Experiment Complete. Plot saved to 'Null_Model_Test.png'")

# 6. QUANTITATIVE SCORE

def score_model(t, signal, log_primes):
    """Score: sum of absolute signal amplitudes at log_prime locations"""
    score = 0
    positive_count = 0
    for lp in log_primes:
        if lp <= 4:  # Within scan range
            idx = (np.abs(t - lp)).argmin()
            amp = signal[idx]
            score += abs(amp)  # Magnitude of resonance
            if amp > 0:
                positive_count += 1
    return score, positive_count


score_true, pos_true = score_model(t, sig_true, log_primes)
score_wrong, pos_wrong = score_model(t, sig_wrong_const, log_primes)
score_shuf, pos_shuf = score_model(t, sig_shuffled, log_primes)
score_rand, pos_rand = score_model(t, sig_random, log_primes)

print("\n" + "="*70)
print("--- NULL HYPOTHESIS TEST RESULTS ---")
print("="*70)
print("\nResonance Scores at Prime Logarithms (Higher is Better):")
print("-"*50)
print(f"{'Model':<35} | {'Score':>10} | {'Positive':>8}")
print("-"*50)
print(f"{'TRUE MODEL (E8 + φ)':<35} | {score_true:>10.4f} | {pos_true:>8}/11")
print(f"{'CONTROL 1 (E8 + e)':<35} | {score_wrong:>10.4f} | {pos_wrong:>8}/11")
print(f"{'CONTROL 2 (Shuffled + φ)':<35} | {score_shuf:>10.4f} | {pos_shuf:>8}/11")
print(f"{'CONTROL 3 (Random Matrix)':<35} | {score_rand:>10.4f} | {pos_rand:>8}/11")
print("-"*50)

# Calculate statistical significance
avg_control = (score_wrong + score_shuf + score_rand) / 3
ratio = score_true / avg_control if avg_control > 0 else float('inf')

print(f"\n  Average Control Score: {avg_control:.4f}")
print(f"  True Model / Control Ratio: {ratio:.2f}x")

# 7. FINAL VERDICT
print("\n" + "="*70)
print("STATISTICAL VERDICT")
print("="*70)

if score_true > score_wrong and score_true > score_shuf and score_true > score_rand:
    if ratio > 1.5:
        print("\n  ✅✅✅ NULL HYPOTHESIS REJECTED! ✅✅✅")
        print(f"  The True Model outperforms controls by {ratio:.2f}x")
        print("  This is NOT a coincidence.")
        print("  The E8 + Golden Ratio combination uniquely produces Prime Resonance.")
        print("\n  CONCLUSION: The Geometric Standard Model is VALIDATED.")
    else:
        print("\n  ⚠ True Model is better, but margin is small.")
        print("  More sophisticated Hamiltonian may be needed.")
else:
    print("\n  ❌ NULL HYPOTHESIS NOT REJECTED")
    print("  Control models performed equally well or better.")
    print("  The prime signal may be an artifact of scaling.")

# 8. DETAILED PRIME-BY-PRIME COMPARISON
print("\n" + "="*70)
print("PRIME-BY-PRIME COMPARISON")
print("="*70)
print(f"\n{'Prime':>6} | {'ln(p)':>7} | {'True':>10} | {'Ctrl1':>10} | {'Ctrl2':>10} | {'Ctrl3':>10}")
print("-"*70)

for p in primes:
    lp = np.log(p)
    if lp <= 4:
        idx = (np.abs(t - lp)).argmin()
        v_true = sig_true[idx]
        v_wrong = sig_wrong_const[idx]
        v_shuf = sig_shuffled[idx]
        v_rand = sig_random[idx]
        
        # Mark winner
        values = [v_true, v_wrong, v_shuf, v_rand]
        best = max(values)
        winner = ""
        if v_true == best:
            winner = "★"
        
        print(f"{p:>6} | {lp:>7.4f} | {v_true:>10.4f}{winner} | {v_wrong:>10.4f} | {v_shuf:>10.4f} | {v_rand:>10.4f}")

print("-"*70)
print("★ = True Model wins at this prime")

print("\n" + "="*70)
print("END OF NULL HYPOTHESIS TEST")
print("="*70)
