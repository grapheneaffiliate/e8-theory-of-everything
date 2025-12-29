"""
E8 FERMION HUNTER
=================
Goal: Identify Standard Model Matter (Fermions) within the E8 'Dark Sector'.
Method:
1. Load the N=12 Gauge Bosons (Forces).
2. Scan the remaining 228 roots.
3. Filter for roots that INTERACT with the Weak Force (W/Z Bosons).
   Interaction = Non-zero geometric dot product.
4. Check if these candidates cluster into 3 Generations (Mass Shells).
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
# 2. FERMION HUNTING LOGIC
# ==========================================
def hunt_fermions():
    """
    Hunt for fermions in the 228 "dark" roots by testing
    geometric interaction with the Weak Force bosons.
    """
    print("="*70)
    print("E8 FERMION HUNTER: SEARCHING FOR MATTER")
    print("="*70)
    print()
    
    # 1. Project Universe
    shadows = E8_ROOTS @ UNIVERSE_MATRIX.T
    lengths_sq = np.sum(shadows**2, axis=1)
    lengths = np.sqrt(lengths_sq)
    
    # 2. Identify Sectors
    # The 12 lightest are Bosons (from Phase 1)
    sorted_idx = np.argsort(lengths_sq)
    boson_indices = sorted_idx[:12]
    candidate_indices = sorted_idx[12:]
    
    print(f"Light Sector (Bosons): {len(boson_indices)}")
    print(f"Dark Sector (Candidates): {len(candidate_indices)}")
    print("-"*70)
    print()
    
    # 3. Identify The Weak Force Bosons
    # The 3 heaviest bosons are the W⁺, W⁻, Z (SU(2) Weak)
    weak_boson_indices = boson_indices[-3:]
    weak_vectors = shadows[weak_boson_indices]
    
    print("[Weak Force Identification]")
    print("-"*70)
    print(f"W/Z Bosons (heaviest 3 of 12):")
    for i, idx in enumerate(weak_boson_indices):
        print(f"  Boson {idx}: Mass = {lengths[idx]:.6f}")
    print()
    
    # 4. Scan Dark Sector for Weak Interaction
    print("[Scanning Dark Sector for Matter...]")
    print("-"*70)
    
    fermion_candidates = []
    
    for idx in candidate_indices:
        root_vec = shadows[idx]
        length = lengths[idx]
        
        # Calculate Geometric Interaction (Dot Product with Weak Bosons)
        # Interaction Strength = Sum(|Root · W|)
        interaction = np.sum(np.abs(root_vec @ weak_vectors.T))
        
        # Filter: Fermions MUST interact with the Weak Force
        if interaction > 0.1:  # Threshold for significant interaction
            fermion_candidates.append({
                'id': idx,
                'mass': length,
                'interaction': interaction
            })
    
    print(f"Found {len(fermion_candidates)} roots interacting with Weak Force.")
    print()
    
    if not fermion_candidates:
        print("No candidates found. Adjusting threshold...")
        return None
    
    # 5. Analyze Generations (Mass Clustering)
    print("[GENERATION ANALYSIS]")
    print("-"*70)
    print("Looking for mass clustering in Fermion candidates...")
    print()
    
    # Extract and sort masses
    masses = sorted([f['mass'] for f in fermion_candidates])
    
    # Cluster masses with tolerance
    tolerance = 0.005
    shells = []
    current_shell = [masses[0]]
    
    for m in masses[1:]:
        if abs(m - current_shell[-1]) < tolerance:
            current_shell.append(m)
        else:
            shells.append(current_shell.copy())
            current_shell = [m]
    shells.append(current_shell)
    
    # 6. Report Findings
    print(f"Detected {len(shells)} distinct Fermion Mass Shells:")
    print()
    print(f"{'Gen':<5} | {'Mass (Length)':<15} | {'Count (States)':<15} | {'Interpretation'}")
    print("-"*70)
    
    generation_count = 0
    for i, shell in enumerate(shells):
        avg_mass = np.mean(shell)
        count = len(shell)
        
        # Interpret based on count
        if count == 16:
            interp = "Full Generation (16 Weyl)"
        elif count == 48:
            interp = "3 Generations × 16"
        elif count % 3 == 0 and count >= 6:
            interp = f"Color multiplet ({count//3} × 3)"
        elif count >= 8:
            interp = "Large multiplet"
        else:
            interp = f"{count} states"
        
        print(f"{i+1:<5} | {avg_mass:.6f}        | {count:<15} | {interp}")
        
        if count >= 8:  # Significant population
            generation_count += 1
    
    print("-"*70)
    print()
    
    # 7. Verdict
    print("="*70)
    print("FERMION HUNTING RESULTS")
    print("="*70)
    print()
    
    if generation_count == 3:
        print("🏆 DISCOVERY: EXACTLY 3 GENERATIONS FOUND!")
        print("   The E8 geometry naturally organizes matter into 3 mass tiers.")
        print("   This matches the Standard Model: e/μ/τ, u/c/t, d/s/b")
    elif generation_count == 4:
        print("? INTERESTING: 4 Generations found!")
        print("  Possible heavy 4th generation at GUT scale?")
    elif generation_count == 2:
        print("~ PARTIAL: Only 2 Generations detected.")
        print("  May need to adjust interaction threshold.")
    else:
        print(f"~ Result: {generation_count} mass clusters found.")
        print("  Structure is complex; further analysis needed.")
    
    print()
    print("="*70)
    
    return {
        'n_candidates': len(fermion_candidates),
        'n_generations': generation_count,
        'shells': shells
    }


if __name__ == "__main__":
    hunt_fermions()
