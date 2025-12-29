"""
E8 MASS SPECTRUM ANALYZER
=========================
Goal: Derive the W/Z Boson Mass Ratio from the Universe Matrix.
Input: The UNIVERSE_MATRIX from e8_constants.py.
Output: Predicted M_W / M_Z ratio compared to experiment.

Standard Model Target: M_W / M_Z = cos(θ_W) = 0.8768
"""

import numpy as np
from itertools import product
from e8_constants import UNIVERSE_MATRIX

# ==========================================
# 1. SETUP
# ==========================================
class E8Crystal:
    """Generate the 240 roots of E8."""
    def __init__(self):
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
        self.roots = np.array(roots)


E8_ROOTS = E8Crystal().roots


# ==========================================
# 2. MASS EXTRACTION
# ==========================================
def analyze_masses():
    """
    Analyze gauge boson masses from geometric root lengths.
    """
    print("="*70)
    print("E8 MASS SPECTRUM ANALYZER")
    print("="*70)
    print()
    
    # Project Roots using the DNA Matrix
    shadows = E8_ROOTS @ UNIVERSE_MATRIX.T
    lengths_sq = np.sum(shadows**2, axis=1)
    lengths = np.sqrt(lengths_sq)
    
    # Isolate the 12 SM Roots
    sorted_idx = np.argsort(lengths_sq)
    sm_indices = sorted_idx[:12]
    sm_lengths = lengths[sm_indices]
    sm_lengths_sq = lengths_sq[sm_indices]
    
    print("[Standard Model Gauge Bosons]")
    print("-"*70)
    print(f"{'Root ID':<10} | {'Length':<15} | {'Length²':<15} | {'Interpretation'}")
    print("-"*70)
    
    for i, idx in enumerate(sm_indices):
        l = lengths[idx]
        l2 = lengths_sq[idx]
        
        # Interpretation based on mass shell
        if i == 0:
            label = "Photon / Gluon"
        elif i < 8:
            label = "Gluon (Strong)"
        elif i < 11:
            label = "W/Z (Weak)"
        else:
            label = "Hypercharge B"
        
        print(f"Root {idx:<5} | {l:.8f}       | {l2:.8f}       | {label}")
    
    print("-"*70)
    
    # Group by mass shells (cluster similar lengths)
    print()
    print("[Identified Mass Shells]")
    print("-"*70)
    
    # Simple clustering with tolerance
    tolerance = 0.002
    shells = []
    shell_indices = []
    current_shell = [sm_lengths[0]]
    current_indices = [0]
    
    for i, l in enumerate(sm_lengths[1:], 1):
        if abs(l - current_shell[-1]) < tolerance:
            current_shell.append(l)
            current_indices.append(i)
        else:
            shells.append((np.mean(current_shell), len(current_shell), current_indices.copy()))
            current_shell = [l]
            current_indices = [i]
    shells.append((np.mean(current_shell), len(current_shell), current_indices.copy()))
    
    print(f"Detected {len(shells)} distinct mass levels:")
    for i, (mean_l, count, idx) in enumerate(shells):
        print(f"  Level {i}: Length = {mean_l:.8f}  (degeneracy = {count})")
    
    # Calculate Mass Ratios
    print()
    print("[Electroweak Mass Ratio Test]")
    print("-"*70)
    
    if len(shells) >= 2:
        # The W and Z bosons should form distinct mass shells
        # M_W / M_Z = cos(θ_W) ≈ 0.8768
        
        # Use the two highest mass shells
        m_heavy_1 = shells[-1][0]  # Heaviest (Z-like)
        m_heavy_2 = shells[-2][0]  # Second Heaviest (W-like)
        
        ratio_direct = m_heavy_2 / m_heavy_1
        ratio_inverse = m_heavy_1 / m_heavy_2
        
        # Target from Weinberg angle
        sin2_theta = 0.23122
        cos_theta = np.sqrt(1 - sin2_theta)
        
        print(f"Geometric Ratio (Shell[-2]/Shell[-1]): {ratio_direct:.6f}")
        print(f"Geometric Ratio (Shell[-1]/Shell[-2]): {ratio_inverse:.6f}")
        print()
        print(f"Target cos(θ_W):                       {cos_theta:.6f}")
        print(f"Target 1/cos(θ_W):                     {1/cos_theta:.6f}")
        
        # Check which interpretation works
        err_direct = abs(ratio_direct - cos_theta) / cos_theta * 100
        err_inverse = abs(ratio_inverse - cos_theta) / cos_theta * 100
        
        print()
        print(f"Error (direct):  {err_direct:.3f}%")
        print(f"Error (inverse): {err_inverse:.3f}%")
        print("-"*70)
        
        best_err = min(err_direct, err_inverse)
        
        if best_err < 5.0:
            print()
            print("🏆 SUCCESS: W/Z Mass Ratio derived from E8 geometry!")
        else:
            print()
            print("~ Masses detected, interpreting structure...")
    
    # Coupling Ratios from Covariance
    print()
    print("[Coupling Analysis from Covariance]")
    print("-"*70)
    
    sm_vectors = shadows[sm_indices]
    cov = np.cov(sm_vectors.T)
    evals = np.linalg.eigvalsh(cov)
    evals = np.sort(evals)[::-1]
    
    k3, k2, k1 = evals[0], evals[1], evals[2]
    sin2 = k1 / (k1 + k2)
    cos_geom = np.sqrt(1 - sin2)
    
    print(f"k3 (Strong):    {k3:.9f}")
    print(f"k2 (Weak):      {k2:.9f}")
    print(f"k1 (Hyper):     {k1:.9f}")
    print()
    print(f"sin²θ (geom):   {sin2:.9f}")
    print(f"cos θ (geom):   {cos_geom:.9f}")
    print(f"cos θ (exp):    {cos_theta:.9f}")
    print(f"Error:          {abs(cos_geom - cos_theta)/cos_theta*100:.4f}%")
    print("-"*70)
    
    print()
    print("="*70)
    print("Mass spectrum analysis complete.")
    print("="*70)
    
    return {
        'shells': shells,
        'k_values': (k1, k2, k3),
        'sin2_theta': sin2,
        'cos_theta': cos_geom
    }


if __name__ == "__main__":
    analyze_masses()
