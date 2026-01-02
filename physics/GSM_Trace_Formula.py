#!/usr/bin/env python3
"""
GSM TRACE FORMULA: The Spectral Fourier Transform
Testing if E8 Eigenvalues vibrate at Prime Frequencies

This is the RIGOROUS "Spectral Test" - not numerology, but spectral geometry.
If E8 eigenvalues encode Riemann zeros, their Fourier transform shows
sharp spikes at ln(prime) locations.

Author: Timothy McGirl
Date: January 2, 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import product

print("="*70)
print("GSM TRACE FORMULA: The Spectral Fourier Transform")
print("Testing if E8 Eigenvalues vibrate at Prime Frequencies")
print("="*70)

# 1. PHYSICAL CONSTANTS
PHI = (1 + np.sqrt(5)) / 2
LAMBDA_SCALING = 14.1347 / 7.52  # Calibration from previous "Gold" run

# 2. GENERATE E8 ROOTS
def get_E8_roots():
    """Generate all 240 roots of the E8 lattice."""
    roots = []
    # Type 1: Permutations of (+/-1, +/-1, 0,0,0,0,0,0)
    for i in range(8):
        for j in range(i + 1, 8):
            for s1 in [-1, 1]:
                for s2 in [-1, 1]:
                    v = np.zeros(8)
                    v[i] = s1
                    v[j] = s2
                    roots.append(v)
    # Type 2: (+/- 0.5... ) even sign changes
    for signs in product([-0.5, 0.5], repeat=8):
        if np.sum(np.array(signs) < 0) % 2 == 0:
            roots.append(np.array(signs))
    return np.array(roots)

print("[1] Generating E8 Lattice...")
roots = get_E8_roots()
n_roots = len(roots)
print(f"    -> {n_roots} roots generated (expected: 240)")

# 3. BUILD HAMILTONIAN (The "Golden Tunneling" Model)
print("[2] Constructing Hamiltonian (Golden Tunneling)...")
H = np.zeros((n_roots, n_roots))

# Interaction: Dot Product * Phi Decay
# This logic was the most promising in previous runs
for i in range(n_roots):
    for j in range(n_roots):
        if i == j:
            H[i, j] = np.dot(roots[i], roots[i]) * PHI  # Diagonal
        else:
            dist_sq = np.sum((roots[i] - roots[j])**2)
            dot = np.dot(roots[i], roots[j])
            H[i, j] = dot * (PHI ** (-dist_sq))

print(f"    -> Hamiltonian matrix: {n_roots}x{n_roots}")

# 4. COMPUTE EIGENVALUES
print("[3] Diagonalizing Matrix...")
eigenvalues = np.linalg.eigvalsh(H)
# Scale them to match the physical region of the first zero
# We apply the calibration found in the previous step
eigenvalues = eigenvalues * LAMBDA_SCALING

# Filter low noise
energies = np.sort(np.abs(eigenvalues))
energies = energies[energies > 10]
print(f"    -> Extracted {len(energies)} resonant modes (E > 10)")
print(f"    -> Energy range: [{energies.min():.4f}, {energies.max():.4f}]")

# 5. THE TRACE FORMULA (The "Prime Hunter")
# We scan "time" t. If E8 encodes primes, we see spikes at t = ln(Prime).
print("\n[4] Computing Spectral Trace Function F(t) = Sum cos(E_n * t)...")

t_values = np.linspace(0, 4, 2000)  # Scan range corresponding to small primes
trace_signal = np.zeros_like(t_values)

for E in energies:
    # The Trace Formula: Oscillations of eigenvalues sum to delta functions at primes
    trace_signal += np.cos(E * t_values)

# Normalize
trace_signal = trace_signal / len(energies)

# 6. PLOTTING & ANALYSIS
print("[5] Analyzing Peaks against Prime Logarithms...")

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
log_primes = np.log(primes)

# Find local maxima in signal
from scipy.signal import find_peaks
try:
    peaks, properties = find_peaks(trace_signal, height=np.max(trace_signal)*0.3, distance=20)
    peak_locations = t_values[peaks]
    peak_amplitudes = trace_signal[peaks]
except ImportError:
    # Fallback without scipy
    peaks = []
    peak_locations = []
    peak_amplitudes = []
    for i in range(1, len(trace_signal)-1):
        if trace_signal[i] > trace_signal[i-1] and trace_signal[i] > trace_signal[i+1]:
            if trace_signal[i] > np.max(trace_signal)*0.3:
                peaks.append(i)
                peak_locations.append(t_values[i])
                peak_amplitudes.append(trace_signal[i])
    peak_locations = np.array(peak_locations)
    peak_amplitudes = np.array(peak_amplitudes)

plt.style.use('dark_background')
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

# Top plot: Full trace
ax1.plot(t_values, trace_signal, color='cyan', lw=1.5, label='E8 Spectral Trace F(t)')
if len(peak_locations) > 0:
    ax1.scatter(peak_locations, peak_amplitudes, color='yellow', s=50, zorder=5, label='Signal Peaks')

# Mark the locations of Primes
for p, lp in zip(primes, log_primes):
    if lp <= 4:  # Only plot within range
        ax1.axvline(x=lp, color='magenta', linestyle='--', alpha=0.6)
        ax1.text(lp, ax1.get_ylim()[1]*0.95, f"ln({p})", color='magenta', ha='center', fontsize=9)

ax1.set_title("The Sound of the E8 Lattice: Spectral Trace Formula\nPeaks should align with ln(Prime) [Berry-Keating Test]", 
              fontsize=14, color='gold')
ax1.set_xlabel("t (Time / Log Scale)", fontsize=12)
ax1.set_ylabel("F(t) = Σ cos(Eₙt) / N", fontsize=12)
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.2)

# Bottom plot: Correlation analysis
# Measure how well signal correlates with prime locations
ax2.bar(range(len(primes)), [trace_signal[(np.abs(t_values - lp)).argmin()] for lp in log_primes], 
        color='magenta', alpha=0.7, label='Signal at ln(p)')
ax2.set_xticks(range(len(primes)))
ax2.set_xticklabels([f'ln({p})=' + f'{np.log(p):.2f}' for p in primes], rotation=45, fontsize=9)
ax2.set_ylabel('Trace Amplitude', fontsize=12)
ax2.set_title('Signal Amplitude at Prime Logarithm Locations', fontsize=12, color='gold')
ax2.axhline(y=0, color='white', linestyle='-', alpha=0.3)
ax2.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('E8_Trace_Formula.png', dpi=150)
print("    -> Plot saved to 'E8_Trace_Formula.png'")

# 7. TEXT ANALYSIS OF PEAKS
print("\n" + "="*70)
print("--- PRIME CORRELATION ANALYSIS ---")
print("="*70)
print("\nChecking if signal peaks align with ln(primes)...\n")

print(f"{'Prime':>6} | {'ln(p)':>8} | {'Amplitude':>12} | {'Status':>12}")
print("-"*50)

positive_count = 0
total_amp = 0
for p, lp in zip(primes, log_primes):
    if lp > 4:
        continue
    # Find signal value at this exact log-prime location
    idx = (np.abs(t_values - lp)).argmin()
    amp = trace_signal[idx]
    total_amp += amp
    status = "✓ POSITIVE" if amp > 0 else "✗ negative"
    if amp > 0:
        positive_count += 1
    print(f"{p:>6} | {lp:>8.4f} | {amp:>12.6f} | {status}")

print("-"*50)
print(f"\nSummary:")
print(f"  Positive at {positive_count}/{len([p for p in primes if np.log(p) <= 4])} prime locations")
print(f"  Mean amplitude at primes: {total_amp/len([p for p in primes if np.log(p) <= 4]):.6f}")

# Check peak-prime alignment
print("\n--- PEAK-PRIME ALIGNMENT ---")
if len(peak_locations) > 0:
    print(f"\nDetected {len(peak_locations)} peaks in signal")
    for i, (loc, amp) in enumerate(zip(peak_locations[:10], peak_amplitudes[:10])):
        # Find nearest prime
        distances = np.abs(log_primes - loc)
        nearest_idx = np.argmin(distances)
        nearest_prime = primes[nearest_idx]
        dist = distances[nearest_idx]
        match = "MATCH!" if dist < 0.1 else ""
        print(f"  Peak {i+1}: t={loc:.4f}, amp={amp:.4f} | Nearest: ln({nearest_prime})={log_primes[nearest_idx]:.4f} (Δ={dist:.4f}) {match}")
else:
    print("No significant peaks detected.")

# 8. FINAL VERDICT
print("\n" + "="*70)
print("SPECTRAL TEST CONCLUSION")
print("="*70)

mean_prime_amp = np.mean([trace_signal[(np.abs(t_values - lp)).argmin()] for lp in log_primes if lp <= 4])
mean_nonprime_amp = np.mean(trace_signal)

print(f"\n  Mean amplitude at ln(primes): {mean_prime_amp:.6f}")
print(f"  Mean amplitude (background):  {mean_nonprime_amp:.6f}")
print(f"  Ratio (prime/background):     {mean_prime_amp/mean_nonprime_amp:.2f}x")

if mean_prime_amp > mean_nonprime_amp and positive_count > len(primes)//2:
    print("\n  ✓✓✓ E8 EIGENVALUES SHOW PRIME RESONANCE! ✓✓✓")
    print("  The Hamiltonian passes the Berry-Keating spectral test.")
    print("  This supports the geometric origin of Riemann zeros.")
else:
    print("\n  ⚠ Prime resonance is weak or absent.")
    print("  Consider alternative Hamiltonian interactions (symplectic, etc.)")

print("\n" + "="*70)
print("END OF SPECTRAL ANALYSIS")
print("="*70)
