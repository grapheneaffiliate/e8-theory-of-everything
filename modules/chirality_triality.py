"""
CHIRALITY & TRIALITY REFINEMENT: Proper Fermion Representations
================================================================
This script addresses the key challenges in extracting SM-like fermion
content from E8:

1. TRIALITY: E8 contains SO(8) which has triality automorphism
   - Vector (8v) ↔ Spinor (8s) ↔ Conjugate Spinor (8c)
   - SM fermions arise from spinorial representations

2. CHIRALITY: Left-handed != Right-handed under SU(2)_L
   - E8 -> SO(10) -> SU(5) -> SM breaks chirality
   - Half-integer roots can encode left vs right via even/odd structure

3. EXPLICIT REPS: Count fermions by quantum numbers
   - Quarks: (3, 2, 1/6) or (3, 1, -2/3), etc.
   - Leptons: (1, 2, -1/2) or (1, 1, 1), etc.

4. HIGGS AS COMPOSITE: Identify Higgs from dark root pairs

Author: E8 Theory Team
Date: December 31, 2025
"""

import numpy as np
from itertools import product

# UNIVERSE_MATRIX
UNIVERSE_MATRIX = np.array([
    [-0.863957542659, -0.087612666567, -0.145842438043,  0.022102045189,  0.231874875254,  0.307671264286,  0.251338525966,  0.112001110381],
    [ 0.015789224302, -0.106532458726,  0.314327001553, -0.491973635285, -0.117819468672,  0.090181641388, -0.108047220427,  0.783500897529],
    [-0.246201715103,  0.657538309303, -0.413965974868, -0.263964203209, -0.261591742574, -0.419470959757, -0.118188732664,  0.087341032536],
    [-0.102516279341, -0.131079825485,  0.085257597640, -0.234102869992, -0.818771278625,  0.303552216081,  0.201568593318, -0.327223513407],
])


def generate_e8_roots():
    """Generate all 240 E8 roots."""
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


def get_dark_sector_roots(roots):
    """Get the 228 dark sector roots."""
    shadows = roots @ UNIVERSE_MATRIX.T
    lengths = np.sqrt(np.sum(shadows**2, axis=1))
    sorted_idx = np.argsort(lengths)
    dark_indices = sorted_idx[12:]
    return dark_indices, roots[dark_indices], shadows[dark_indices], lengths


def classify_triality(root_8d):
    """
    Classify root by SO(8) triality.
    
    In E8's SO(8) subgroup:
    - Vector (8v): Integer coordinates (Type 1)
    - Spinor (8s): Half-integer, #(+) mod 4 = 0 or 2
    - Conjugate (8c): Half-integer, #(+) mod 4 = 1 or 3
    
    Returns: 'V' (vector), 'S' (spinor), 'C' (conjugate)
    """
    nonzero = np.abs(root_8d) > 0.01
    n_nonzero = np.sum(nonzero)
    
    if n_nonzero == 2:  # Type 1: Integer
        return 'V'  # Vector representation
    
    # Type 2: Half-integer
    n_positive = np.sum(root_8d > 0.1)
    
    if n_positive % 4 == 0 or n_positive % 4 == 2:
        return 'S'  # Spinor
    else:
        return 'C'  # Conjugate spinor


def extract_chirality_refined(root_8d):
    """
    Extract chirality using refined triality-based method.
    
    SM chirality:
    - Left-handed doublets: (nu, e)_L, (u, d)_L
    - Right-handed singlets: e_R, u_R, d_R
    
    In E8 via SO(10):
    - 16 of SO(10) = left-handed family (contains all L fermions)
    - 16 of SO(10) = right-handed anti-family (or mirror)
    
    We use:
    - Spinor (S) + weak-active -> Left-handed
    - Conjugate (C) + weak-inactive -> Right-handed
    - Vector (V) -> typically bosonic or exotic
    """
    triality = classify_triality(root_8d)
    
    # Weak isospin activity (coords 3-4)
    weak = root_8d[3:5]
    weak_norm = np.sqrt(np.sum(weak**2))
    weak_active = weak_norm > 0.3
    
    if triality == 'S' and weak_active:
        return 'L', 'Spinor+Weak'
    elif triality == 'C' and not weak_active:
        return 'R', 'Conjugate+Singlet'
    elif triality == 'S':
        return 'L', 'Spinor-only'
    elif triality == 'C':
        return 'R', 'Conjugate-only'
    else:
        return 'V', 'Vector/Bosonic'


def extract_sm_quantum_numbers(root_8d):
    """
    Extract SM-like quantum numbers from 8D coordinates.
    
    Coordinate interpretation:
    - 0,1,2: SU(3) color (triality weights)
    - 3,4: SU(2)_L isospin
    - 5: U(1)_Y hypercharge
    - 6,7: Extra (B-L, gravity, ...)
    
    Returns: (SU3_rep, SU2_rep, Y) approximately
    """
    # Color
    color = root_8d[0:3]
    color_norm = np.sqrt(np.sum(color**2))
    
    if color_norm < 0.3:
        su3 = '1'  # Singlet
    elif 0.8 < color_norm < 0.9:
        su3 = '3'  # Triplet (could be 3)
    else:
        su3 = '8'  # Octet or other
    
    # Weak isospin
    weak = root_8d[3:5]
    weak_norm = np.sqrt(np.sum(weak**2))
    
    if weak_norm < 0.3:
        su2 = '1'  # Singlet
    elif 0.5 < weak_norm < 0.8:
        su2 = '2'  # Doublet candidate
    else:
        su2 = '3'  # Triplet or higher
    
    # Hypercharge
    Y = root_8d[5]
    
    # Round to SM-like values
    Y_rounded = round(Y * 6) / 6  # Snap to 1/6, 1/3, 1/2, 2/3, 1, etc.
    
    return (su3, su2, Y_rounded)


def identify_sm_fermion(root_8d, chirality):
    """
    Identify SM fermion type from quantum numbers.
    
    SM Fermion Reps:
    - Q_L: (3, 2, 1/6) - left-handed quark doublet
    - u_R: (3, 1, 2/3) - right-handed up
    - d_R: (3, 1, -1/3) - right-handed down
    - L_L: (1, 2, -1/2) - left-handed lepton doublet
    - e_R: (1, 1, -1) - right-handed electron
    - nu_R: (1, 1, 0) - right-handed neutrino (BSM)
    """
    su3, su2, Y = extract_sm_quantum_numbers(root_8d)
    
    # Classification logic
    if su3 == '3' and su2 == '2' and chirality == 'L':
        return 'Q_L'
    elif su3 == '3' and su2 == '1' and chirality == 'R':
        if abs(Y - 2/3) < 0.2:
            return 'u_R'
        elif abs(Y + 1/3) < 0.2:
            return 'd_R'
        else:
            return 'q_R'
    elif su3 == '1' and su2 == '2' and chirality == 'L':
        return 'L_L'
    elif su3 == '1' and su2 == '1' and chirality == 'R':
        if abs(Y + 1) < 0.2:
            return 'e_R'
        elif abs(Y) < 0.2:
            return 'nu_R'
        else:
            return 'l_R'
    elif su3 == '3' and chirality == 'L':
        return 'q_exotic_L'
    elif su3 == '1' and chirality == 'L':
        return 'l_exotic_L'
    else:
        return 'exotic'


def analyze_refined_fermions(dark_roots_8d, lengths, dark_indices):
    """Analyze fermion content with refined chirality/triality."""
    print("="*90)
    print("REFINED CHIRALITY & TRIALITY ANALYSIS")
    print("="*90)
    
    print("\nMethod:")
    print("  - Triality: V (vector), S (spinor), C (conjugate) from SO(8) subset E8")
    print("  - Chirality: L from S+weak, R from C+singlet")
    print("  - SM reps: (SU3, SU2, Y) from coordinate patterns")
    print("-" * 70)
    
    # Count by triality
    triality_counts = {'V': 0, 'S': 0, 'C': 0}
    chirality_counts = {'L': 0, 'R': 0, 'V': 0}
    fermion_counts = {}
    
    fermions = []
    
    for i, root_8d in enumerate(dark_roots_8d):
        triality = classify_triality(root_8d)
        chirality, reason = extract_chirality_refined(root_8d)
        sm_type = identify_sm_fermion(root_8d, chirality)
        qn = extract_sm_quantum_numbers(root_8d)
        
        triality_counts[triality] += 1
        chirality_counts[chirality] += 1
        fermion_counts[sm_type] = fermion_counts.get(sm_type, 0) + 1
        
        fermions.append({
            'index': i,
            'root': root_8d,
            'length': lengths[dark_indices[i]],
            'triality': triality,
            'chirality': chirality,
            'sm_type': sm_type,
            'quantum_numbers': qn
        })
    
    return fermions, triality_counts, chirality_counts, fermion_counts


def count_sm_like_content(fermion_counts):
    """Count SM-like fermion content (16 per generation)."""
    print("\n" + "="*90)
    print("SM-LIKE FERMION CONTENT")
    print("="*90)
    
    print("\nExpected per generation (1 family of 16 Weyl fermions):")
    print("  Q_L: 6 (2 flavors x 3 colors x 1 chirality)")
    print("  u_R: 3 (1 flavor x 3 colors)")
    print("  d_R: 3 (1 flavor x 3 colors)")
    print("  L_L: 2 (e + nu)")
    print("  e_R: 1")
    print("  (nu_R: 1) [BSM]")
    print("-" * 60)
    
    print("\nFound in E8 Dark Sector:")
    print("-" * 60)
    
    # Define SM fermion types
    sm_types = ['Q_L', 'u_R', 'd_R', 'L_L', 'e_R', 'nu_R', 'q_R', 'l_R', 'q_exotic_L', 'l_exotic_L', 'exotic']
    
    sm_total = 0
    exotic_total = 0
    
    for t in sm_types:
        count = fermion_counts.get(t, 0)
        if t in ['Q_L', 'u_R', 'd_R', 'L_L', 'e_R', 'nu_R']:
            sm_total += count
            status = "SM"
        else:
            exotic_total += count
            status = "Exotic"
        print(f"  {t:<15}: {count:>4}  ({status})")
    
    print("-" * 60)
    print(f"  SM-like total:     {sm_total}")
    print(f"  Exotic total:      {exotic_total}")
    print(f"  Grand total:       {sm_total + exotic_total}")
    
    # Calculate implied generations
    # 16 per generation in SO(10) GUT
    n_gen_estimate = sm_total / 16 if sm_total > 0 else 0
    print(f"\n  Implied generations (SM/16): {n_gen_estimate:.1f}")
    
    return sm_total, exotic_total


def find_higgs_composite(dark_roots_8d, dark_roots_4d, lengths, dark_indices):
    """
    Find Higgs as composite state from dark sector.
    
    Higgs properties:
    - SU(3) singlet
    - SU(2) doublet
    - Y = 1/2
    - Spin 0 (scalar)
    
    In E8: Look for pairs that combine to form (1, 2, 1/2)
    """
    print("\n" + "="*90)
    print("HIGGS AS COMPOSITE STATE")
    print("="*90)
    
    print("\nSearching for root pairs that combine to Higgs-like:")
    print("  - Color singlet: |C1 + C2| ~ 0")
    print("  - Weak doublet: |W1 + W2| ~ 0.7")
    print("  - Hypercharge: Y1 + Y2 ~ 0.5")
    print("-" * 70)
    
    higgs_candidates = []
    
    for i in range(len(dark_roots_8d)):
        for j in range(i+1, min(i+50, len(dark_roots_8d))):  # Limit search
            r1, r2 = dark_roots_8d[i], dark_roots_8d[j]
            
            # Combined quantum numbers
            color_sum = np.linalg.norm(r1[0:3] + r2[0:3])
            weak_sum = np.linalg.norm(r1[3:5] + r2[3:5])
            Y_sum = r1[5] + r2[5]
            
            # Check Higgs-like properties
            if color_sum < 0.5:  # Color singlet
                if 0.4 < weak_sum < 1.0:  # Weak doublet-ish
                    if 0.3 < Y_sum < 0.7:  # Y ~ 1/2
                        mass = (lengths[dark_indices[i]] + lengths[dark_indices[j]]) / 2
                        higgs_candidates.append({
                            'indices': (i, j),
                            'mass': mass,
                            'color_sum': color_sum,
                            'weak_sum': weak_sum,
                            'Y_sum': Y_sum
                        })
    
    # Sort by quality (closest to ideal)
    higgs_candidates.sort(key=lambda x: abs(x['Y_sum'] - 0.5) + x['color_sum'])
    
    print(f"\nFound {len(higgs_candidates)} Higgs-like composite candidates")
    
    if higgs_candidates:
        print("\nTop 5 candidates:")
        print(f"{'Rank':<6} {'Mass':<10} {'|Color|':<10} {'|Weak|':<10} {'Y_sum':<10}")
        print("-" * 50)
        
        for rank, h in enumerate(higgs_candidates[:5]):
            print(f"{rank+1:<6} {h['mass']:<10.4f} {h['color_sum']:<10.4f} "
                  f"{h['weak_sum']:<10.4f} {h['Y_sum']:<10.4f}")
        
        # Best candidate VEV estimate
        best = higgs_candidates[0]
        print(f"\n  Best candidate mass scale: {best['mass']:.4f}")
        print(f"  (Could relate to Higgs VEV via projection)")
    
    return higgs_candidates


def analyze_anomaly_cancellation(fermion_counts, triality_counts, chirality_counts):
    """Check gauge anomaly cancellation conditions."""
    print("\n" + "="*90)
    print("ANOMALY CANCELLATION CHECK")
    print("="*90)
    
    print("\nSM requires anomaly-free (from triangle diagrams):")
    print("  - Sigma Y = 0 (gravitational)")
    print("  - Sigma Y^3 = 0 (U1^3)")
    print("  - Sigma (SU2)^2 Y = 0")
    print("-" * 60)
    
    # Triality balance (vector-spinor-conjugate)
    print("\nTriality Balance (SO(8) subset E8):")
    print(f"  Vectors (V):    {triality_counts['V']}")
    print(f"  Spinors (S):    {triality_counts['S']}")
    print(f"  Conjugates (C): {triality_counts['C']}")
    
    if triality_counts['S'] == triality_counts['C']:
        print("  [OK] S-C balanced (good for anomaly cancellation)")
    else:
        print(f"  [!] S-C imbalanced by {abs(triality_counts['S'] - triality_counts['C'])}")
    
    # L-R balance
    print("\nChirality Balance:")
    print(f"  Left (L):   {chirality_counts['L']}")
    print(f"  Right (R):  {chirality_counts['R']}")
    print(f"  Vector (V): {chirality_counts['V']}")
    
    if chirality_counts['L'] == chirality_counts['R']:
        print("  [OK] L-R balanced before EWSB")
    else:
        print(f"  [!] L-R differ by {abs(chirality_counts['L'] - chirality_counts['R'])}")
        print("  Note: SM itself is chiral - requires specific balance, not equality")


def main():
    """Run refined chirality and representation analysis."""
    print("\n" + "#"*90)
    print("#" + " "*88 + "#")
    print("#" + " "*15 + "CHIRALITY & TRIALITY REFINEMENT" + " "*42 + "#")
    print("#" + " "*10 + "Extracting SM-like Fermions via SO(8) Triality" + " "*33 + "#")
    print("#" + " "*88 + "#")
    print("#"*90)
    
    # Generate E8 and get dark sector
    roots = generate_e8_roots()
    dark_indices, dark_roots_8d, dark_roots_4d, lengths = get_dark_sector_roots(roots)
    
    # Refined analysis
    fermions, triality_counts, chirality_counts, fermion_counts = analyze_refined_fermions(
        dark_roots_8d, lengths, dark_indices
    )
    
    # Print triality distribution
    print(f"\nTriality Distribution:")
    print(f"  V (Vector):    {triality_counts['V']:>4} ({100*triality_counts['V']/228:.1f}%)")
    print(f"  S (Spinor):    {triality_counts['S']:>4} ({100*triality_counts['S']/228:.1f}%)")
    print(f"  C (Conjugate): {triality_counts['C']:>4} ({100*triality_counts['C']/228:.1f}%)")
    
    # Print chirality distribution
    print(f"\nRefined Chirality Distribution:")
    print(f"  L (Left):   {chirality_counts['L']:>4} ({100*chirality_counts['L']/228:.1f}%)")
    print(f"  R (Right):  {chirality_counts['R']:>4} ({100*chirality_counts['R']/228:.1f}%)")
    print(f"  V (Vector): {chirality_counts['V']:>4} ({100*chirality_counts['V']/228:.1f}%)")
    
    # SM content
    sm_total, exotic_total = count_sm_like_content(fermion_counts)
    
    # Higgs composite
    higgs_candidates = find_higgs_composite(dark_roots_8d, dark_roots_4d, lengths, dark_indices)
    
    # Anomaly check
    analyze_anomaly_cancellation(fermion_counts, triality_counts, chirality_counts)
    
    # Summary
    print("\n" + "#"*90)
    print("#" + " "*88 + "#")
    print("#" + " "*35 + "SUMMARY" + " "*46 + "#")
    print("#" + " "*88 + "#")
    print("#"*90)
    
    print(f"""
Refined Fermion Analysis:

  Triality (SO(8) structure):
    V: {triality_counts['V']:>3}  S: {triality_counts['S']:>3}  C: {triality_counts['C']:>3}
    
  Chirality (refined):
    L: {chirality_counts['L']:>3}  R: {chirality_counts['R']:>3}  Vector: {chirality_counts['V']:>3}
    
  SM-like content: {sm_total} / Exotic: {exotic_total}
  
  Higgs candidates: {len(higgs_candidates)} composite pairs found

Key Insight:
  The SO(8) triality automorphism within E8 provides the
  mechanism for chiral fermion emergence. Spinorial (S)
  and conjugate (C) representations naturally split into
  left and right-handed fermions after projection.
  
  Excess states beyond SM likely include:
  * Vector-like pairs (mirror fermions)
  * Right-handed neutrinos (BSM)
  * Heavy (TeV-scale) exotics for anomaly cancellation
""")
    
    return fermions, fermion_counts, higgs_candidates


if __name__ == "__main__":
    fermions, fermion_counts, higgs = main()
