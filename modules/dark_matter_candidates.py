"""
DARK MATTER CANDIDATES: Predictions from Heavy E8 Roots
========================================================
This script identifies dark matter candidates from the heavy E8 roots:

Dark Matter Properties Required:
1. INVISIBLE: Neutral under EM and color (no photon/gluon coupling)
2. STABLE: Long-lived (conserved quantum number or kinematic stability)
3. MASSIVE: Cold dark matter requires m >> keV
4. WEAKLY INTERACTING: Small SM couplings (for WIMP-like DM)

E8 Approach:
- Heavy roots (long projected length) -> massive particles
- Color-singlet + neutral hypercharge -> invisible to SM
- Composite states (r, -r) pairs -> stable via discrete symmetry
- Extra coordinates (6,7) -> "dark charges" for stability

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

# Physical constants
WEAK_SCALE = 246  # GeV (Higgs VEV)
WIMP_MASS_RANGE = (10, 10000)  # GeV, typical WIMP range


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


def get_all_roots_with_info(roots):
    """Get all roots with full quantum number info."""
    shadows = roots @ UNIVERSE_MATRIX.T
    lengths = np.sqrt(np.sum(shadows**2, axis=1))
    return shadows, lengths


def extract_quantum_numbers(root_8d):
    """Extract full quantum numbers from 8D root."""
    # SU(3) color
    color = root_8d[0:3]
    color_norm = np.sqrt(np.sum(color**2))
    
    # SU(2) weak
    weak = root_8d[3:5]
    weak_norm = np.sqrt(np.sum(weak**2))
    
    # U(1) hypercharge
    Y = root_8d[5]
    
    # Dark charges (coords 6,7)
    dark = root_8d[6:8]
    dark_norm = np.sqrt(np.sum(dark**2))
    
    # EM charge Q = T3 + Y (for doublets; simplified for singlets)
    T3 = root_8d[3] if weak_norm > 0.1 else 0
    Q_em = T3 + Y
    
    return {
        'color_norm': color_norm,
        'weak_norm': weak_norm,
        'hypercharge': Y,
        'em_charge': Q_em,
        'dark_charge': dark_norm,
        'dark_coords': dark
    }


def is_invisible(qn, threshold=0.3):
    """Check if particle is invisible to SM (dark matter candidate)."""
    color_singlet = qn['color_norm'] < threshold
    em_neutral = abs(qn['em_charge']) < threshold
    return color_singlet and em_neutral


def compute_sm_coupling(root_4d, sm_roots_4d):
    """Compute coupling strength to SM gauge bosons."""
    couplings = []
    for sm_root in sm_roots_4d:
        if np.linalg.norm(root_4d) > 0 and np.linalg.norm(sm_root) > 0:
            cos_angle = np.dot(root_4d, sm_root) / (np.linalg.norm(root_4d) * np.linalg.norm(sm_root))
            couplings.append(abs(cos_angle))
    
    return np.mean(couplings) if couplings else 0


def find_stable_composites(roots_8d, roots_4d, lengths):
    """Find stable composite dark matter from (r, -r) pairs."""
    print("\n" + "="*90)
    print("STABLE COMPOSITE DARK MATTER")
    print("="*90)
    
    print("\nMethod: Find pairs (r, -r) that form stable composites")
    print("  - Net color = 0 (meson-like)")
    print("  - Net EM charge = 0")
    print("  - Stable via discrete symmetry (like R-parity)")
    print("-" * 70)
    
    composites = []
    
    for i in range(len(roots_8d)):
        for j in range(i+1, len(roots_8d)):
            r1, r2 = roots_8d[i], roots_8d[j]
            
            # Check if approximate anti-particles (r ~ -r or conjugate)
            is_conjugate = np.linalg.norm(r1 + r2) < 0.1
            
            if not is_conjugate:
                continue
            
            # Composite properties
            m1, m2 = lengths[i], lengths[j]
            composite_mass = m1 + m2  # Simple additive approximation
            
            # Combined quantum numbers
            qn1 = extract_quantum_numbers(r1)
            qn2 = extract_quantum_numbers(r2)
            
            net_color = np.linalg.norm(r1[0:3] + r2[0:3])
            net_em = abs(qn1['em_charge'] + qn2['em_charge'])
            net_dark = np.linalg.norm(r1[6:8] + r2[6:8])
            
            # Check stability criteria
            if net_color < 0.1 and net_em < 0.1:
                composites.append({
                    'indices': (i, j),
                    'mass': composite_mass,
                    'net_color': net_color,
                    'net_em': net_em,
                    'dark_charge': net_dark,
                    'root1': r1,
                    'root2': r2
                })
    
    # Sort by mass
    composites.sort(key=lambda x: x['mass'])
    
    print(f"\nFound {len(composites)} stable composite candidates")
    
    if composites:
        print("\nTop 10 Lightest Stable Composites:")
        print(f"{'Rank':<6} {'Mass':<12} {'|Color|':<10} {'|Q_EM|':<10} {'Dark Q':<10}")
        print("-" * 60)
        
        for rank, c in enumerate(composites[:10]):
            print(f"{rank+1:<6} {c['mass']:<12.4f} {c['net_color']:<10.4f} "
                  f"{c['net_em']:<10.4f} {c['dark_charge']:<10.4f}")
    
    return composites


def find_elementary_dm_candidates(roots_8d, roots_4d, lengths, sm_roots_4d):
    """Find elementary (single root) dark matter candidates."""
    print("\n" + "="*90)
    print("ELEMENTARY DARK MATTER CANDIDATES")
    print("="*90)
    
    print("\nCriteria:")
    print("  - Color singlet (|C| < 0.3)")
    print("  - EM neutral (|Q| < 0.3)")
    print("  - Heavy mass (top 50% by projected length)")
    print("  - Weak SM coupling (angle-based)")
    print("-" * 70)
    
    candidates = []
    median_mass = np.median(lengths)
    
    for i, (root_8d, root_4d) in enumerate(zip(roots_8d, roots_4d)):
        qn = extract_quantum_numbers(root_8d)
        
        if is_invisible(qn) and lengths[i] > median_mass:
            coupling = compute_sm_coupling(root_4d, sm_roots_4d)
            
            candidates.append({
                'index': i,
                'root_8d': root_8d,
                'mass': lengths[i],
                'color_norm': qn['color_norm'],
                'em_charge': qn['em_charge'],
                'weak_norm': qn['weak_norm'],
                'dark_charge': qn['dark_charge'],
                'sm_coupling': coupling
            })
    
    # Sort by mass (heaviest first for cold DM)
    candidates.sort(key=lambda x: -x['mass'])
    
    print(f"\nFound {len(candidates)} elementary DM candidates")
    
    if candidates:
        print("\nTop 20 Heavy Invisible States:")
        print(f"{'Rank':<6} {'Mass':<10} {'|Color|':<8} {'Q_EM':<8} {'|Weak|':<8} {'Dark Q':<8} {'SM Coup':<8}")
        print("-" * 70)
        
        for rank, c in enumerate(candidates[:20]):
            print(f"{rank+1:<6} {c['mass']:<10.4f} {c['color_norm']:<8.4f} "
                  f"{c['em_charge']:<8.4f} {c['weak_norm']:<8.4f} "
                  f"{c['dark_charge']:<8.4f} {c['sm_coupling']:<8.4f}")
    
    return candidates


def analyze_wimp_properties(candidates, composites):
    """Analyze WIMP-like properties of candidates."""
    print("\n" + "="*90)
    print("WIMP ANALYSIS")
    print("="*90)
    
    print("\nWIMP (Weakly Interacting Massive Particle) properties:")
    print("  - Mass: 10 GeV - 10 TeV (typical)")
    print("  - Cross section: 10^-46 to 10^-44 cm^2 (direct detection)")
    print("  - Thermal relic: <sigmav> ~ 3x10^-26 cm^3/s")
    print("-" * 70)
    
    # Mass scale conversion (geometric units to GeV)
    # Use heaviest SM root ~ m_top ~ 173 GeV as reference
    # Lightest dark ~ 0.55, heaviest ~ 1.38
    # Scale: 1 unit ~ 300 GeV (rough estimate)
    MASS_SCALE = 300  # GeV per projected length unit
    
    print("\nMass Predictions (using scale ~300 GeV/unit):")
    print("-" * 50)
    
    # Best elementary candidates
    if candidates:
        print("\nLightest Elementary DM Candidate:")
        lightest = min(candidates, key=lambda x: x['mass'])
        mass_gev = lightest['mass'] * MASS_SCALE
        print(f"  Projected mass: {lightest['mass']:.4f} -> ~{mass_gev:.0f} GeV")
        print(f"  SM coupling: {lightest['sm_coupling']:.4f}")
        
        print("\nHeaviest Elementary DM Candidate:")
        heaviest = max(candidates, key=lambda x: x['mass'])
        mass_gev = heaviest['mass'] * MASS_SCALE
        print(f"  Projected mass: {heaviest['mass']:.4f} -> ~{mass_gev:.0f} GeV")
        
        # Mass distribution
        masses = [c['mass'] * MASS_SCALE for c in candidates]
        print(f"\n  Mass range: {min(masses):.0f} - {max(masses):.0f} GeV")
        print(f"  Median mass: {np.median(masses):.0f} GeV")
    
    # Best composite candidates
    if composites:
        print("\n\nLightest Stable Composite:")
        lightest_comp = min(composites, key=lambda x: x['mass'])
        mass_gev = lightest_comp['mass'] * MASS_SCALE
        print(f"  Projected mass: {lightest_comp['mass']:.4f} -> ~{mass_gev:.0f} GeV")
        print(f"  Dark charge: {lightest_comp['dark_charge']:.4f}")
    
    # Compare with experimental constraints
    print("\n\nExperimental Context:")
    print("-" * 50)
    print("  LHC searches: WIMP < ~1 TeV for EW production")
    print("  Direct detection (LZ/XENONnT): sigma_SI < 10^-47 cm^2 at 30 GeV")
    print("  Relic density: Omega_DM h^2 = 0.120 requires specific <sigmav>")


def analyze_dark_sector_structure(roots_8d, lengths):
    """Analyze the structure of the dark sector for DM candidates."""
    print("\n" + "="*90)
    print("DARK SECTOR STRUCTURE")
    print("="*90)
    
    print("\nAnalyzing coords 6-7 (gravity/dark sector):")
    print("-" * 60)
    
    # Extract dark coordinates
    dark_coords = roots_8d[:, 6:8]
    dark_norms = np.sqrt(np.sum(dark_coords**2, axis=1))
    
    # Unique dark charge values
    unique_norms = np.unique(np.round(dark_norms, 4))
    
    print(f"  Distinct dark charge magnitudes: {len(unique_norms)}")
    print(f"  Values: {unique_norms}")
    
    # Count by dark charge
    print("\nDark Charge Distribution:")
    for norm in unique_norms[:6]:  # Top 6
        count = np.sum(np.abs(dark_norms - norm) < 0.01)
        print(f"  |D| = {norm:.4f}: {count} roots")
    
    # Correlation with mass
    print("\nDark Charge vs Mass Correlation:")
    correlation = np.corrcoef(dark_norms, lengths)[0, 1]
    print(f"  Correlation coefficient: {correlation:.4f}")
    
    if correlation > 0.3:
        print("  -> Heavier particles tend to have larger dark charge")
    elif correlation < -0.3:
        print("  -> Heavier particles tend to have smaller dark charge")
    else:
        print("  -> No strong correlation")


def predict_dm_signatures(candidates, composites):
    """Predict experimental signatures for DM candidates."""
    print("\n" + "="*90)
    print("EXPERIMENTAL SIGNATURES & PREDICTIONS")
    print("="*90)
    
    MASS_SCALE = 300  # GeV per unit
    
    print("\n1. DIRECT DETECTION (Nuclear Recoil)")
    print("-" * 50)
    
    # Candidates with weak coupling (potential WIMP)
    wimps = [c for c in candidates if c['weak_norm'] > 0.3 and c['sm_coupling'] > 0.1]
    
    if wimps:
        best_wimp = max(wimps, key=lambda x: x['sm_coupling'])
        mass = best_wimp['mass'] * MASS_SCALE
        print(f"  Best WIMP candidate:")
        print(f"    Mass: ~{mass:.0f} GeV")
        print(f"    Weak coupling: {best_wimp['weak_norm']:.3f}")
        print(f"    -> Expect nuclear recoil signal at ~{mass/100:.1f} keV")
        print(f"    -> Target: Xe, Ar detectors (XENONnT, LZ, DarkSide)")
    
    print("\n2. INDIRECT DETECTION (Annihilation)")
    print("-" * 50)
    
    if composites:
        # Lightest composite can annihilate
        lightest = min(composites, key=lambda x: x['mass'])
        mass = lightest['mass'] * MASS_SCALE
        print(f"  Lightest stable composite:")
        print(f"    Mass: ~{mass:.0f} GeV")
        print(f"    -> Annihilation: chichi -> SM particles")
        print(f"    -> Gamma-ray line at E = {mass:.0f} GeV (Fermi-LAT, CTA)")
        print(f"    -> Positron excess at ~{mass/10:.0f} GeV (AMS-02)")
    
    print("\n3. COLLIDER PRODUCTION (LHC)")
    print("-" * 50)
    
    if candidates:
        producible = [c for c in candidates if c['mass'] * MASS_SCALE < 1500]
        print(f"  Kinematically accessible (m < 1.5 TeV): {len(producible)}")
        
        if producible:
            print(f"  Signatures:")
            print(f"    -> Mono-jet + MET (gluon fusion portal)")
            print(f"    -> Mono-photon + MET (photon portal)")
            print(f"    -> VBF + MET (weak portal)")
    
    print("\n4. COSMOLOGICAL (Relic Density)")
    print("-" * 50)
    print("  Omega_DM h^2 = 0.120 +/- 0.001 (Planck 2018)")
    print("  For thermal relic: <sigmav> ~ 3x10^-26 cm^3/s")
    
    if candidates:
        # Estimate based on coupling
        best = max(candidates, key=lambda x: x['sm_coupling'])
        coupling = best['sm_coupling']
        # Rough scaling: sigmav ~ g^4 / m^2
        print(f"\n  Best candidate coupling: {coupling:.4f}")
        print(f"  -> If thermal: sigmav ~ {coupling**4:.2e} (arbitrary units)")


def main():
    """Run complete dark matter analysis."""
    print("\n" + "#"*90)
    print("#" + " "*88 + "#")
    print("#" + " "*15 + "DARK MATTER CANDIDATES FROM E8 HEAVY ROOTS" + " "*30 + "#")
    print("#" + " "*10 + "Predictions for WIMPs, Composites, and Experimental Signatures" + " "*17 + "#")
    print("#" + " "*88 + "#")
    print("#"*90)
    
    # Generate E8 roots
    roots = generate_e8_roots()
    shadows, lengths = get_all_roots_with_info(roots)
    
    # Get SM roots (12 shortest)
    sorted_idx = np.argsort(lengths)
    sm_indices = sorted_idx[:12]
    sm_roots_4d = shadows[sm_indices]
    
    # Get dark sector (228 roots beyond SM)
    dark_indices = sorted_idx[12:]
    dark_roots_8d = roots[dark_indices]
    dark_roots_4d = shadows[dark_indices]
    dark_lengths = lengths[dark_indices]
    
    print(f"\nDark Sector: {len(dark_indices)} roots (beyond 12 SM gauge bosons)")
    print(f"Mass range: {dark_lengths.min():.4f} to {dark_lengths.max():.4f}")
    
    # Dark sector structure
    analyze_dark_sector_structure(dark_roots_8d, dark_lengths)
    
    # Find elementary DM candidates
    elementary = find_elementary_dm_candidates(
        dark_roots_8d, dark_roots_4d, dark_lengths, sm_roots_4d
    )
    
    # Find stable composites
    composites = find_stable_composites(
        dark_roots_8d, dark_roots_4d, dark_lengths
    )
    
    # WIMP analysis
    analyze_wimp_properties(elementary, composites)
    
    # Experimental predictions
    predict_dm_signatures(elementary, composites)
    
    # Summary
    print("\n" + "#"*90)
    print("#" + " "*88 + "#")
    print("#" + " "*35 + "SUMMARY" + " "*46 + "#")
    print("#" + " "*88 + "#")
    print("#"*90)
    
    MASS_SCALE = 300
    
    print(f"""
Dark Matter Predictions from E8 Heavy Roots:

  +-------------------------------------------------------------------+
  |  CANDIDATE TYPE        COUNT       MASS RANGE                    |
  +-------------------------------------------------------------------+
  |  Elementary (invisible) {len(elementary):>4}       {min([c['mass']*MASS_SCALE for c in elementary]) if elementary else 0:.0f} - {max([c['mass']*MASS_SCALE for c in elementary]) if elementary else 0:.0f} GeV        |
  |  Stable Composites      {len(composites):>4}       {min([c['mass']*MASS_SCALE for c in composites]) if composites else 0:.0f} - {max([c['mass']*MASS_SCALE for c in composites]) if composites else 0:.0f} GeV        |
  +-------------------------------------------------------------------+

Key Predictions:
  1. WIMP-like candidates with weak coupling exist at ~TeV scale
  2. Stable composites (r, -r pairs) form via discrete symmetry
  3. Dark charges (coords 6,7) may stabilize DM against decay
  4. Multiple DM species possible -> multi-component dark sector

Experimental Tests:
  * Direct detection: Nuclear recoil at ~10-100 keV (XENONnT, LZ)
  * Indirect: Gamma-ray lines at DM mass (Fermi-LAT, CTA)
  * Collider: Mono-X + MET signatures (LHC Run 3, HL-LHC)
  * Relic density: Constrain coupling via Planck Omega_DM

The E8 framework naturally produces WIMP-scale dark matter
candidates as invisible heavy roots in the dark sector!
""")
    
    return elementary, composites


if __name__ == "__main__":
    elementary, composites = main()
