"""
NEUTRINO SECTOR: See-Saw Mechanism and PMNS Matrix
===================================================
Goal: Derive neutrino masses and mixing from E8 dark sector roots.

Theory:
- Light neutrino masses via Type-I see-saw: m_nu = -m_D^2 / M_R
- Heavy right-handed neutrinos from dark sector
- PMNS matrix from geometric mixing angles

EXPERIMENTAL TARGETS:
- Neutrino mass differences: Deltam^2_21 ~ 7.5x10^-^5 eV^2, Deltam^2_31 ~ 2.5x10^-^3 eV^2
- PMNS mixing angles: theta_12 ~ 33 deg, theta_23 ~ 45 deg, theta_13 ~ 8.5 deg
- CP phase: delta_CP ~ 230 deg

Author: E8 Theory Team
Date: December 31, 2025
"""

import numpy as np
from itertools import product
from typing import Dict, List, Tuple
from e8_constants import UNIVERSE_MATRIX, EXPERIMENTAL_SIN2_THETA

# Experimental neutrino data
NEUTRINO_DATA = {
    'mass_differences': {
        'delta_m21_sq': 7.53e-5,  # eV^2
        'delta_m31_sq': 2.453e-3,  # eV^2
    },
    'pmns_angles': {
        'theta_12': 33.44,  # degrees
        'theta_23': 49.0,   # degrees
        'theta_13': 8.57,   # degrees
    },
    'cp_phase': 230.0  # degrees (tentative)
}


class NeutrinoSector:
    """
    Neutrino mass generation via geometric see-saw mechanism.
    
    Type-I See-Saw Formula:
        m_nu = -m_D^T M_R^(-1) m_D
    
    where:
    - m_D: Dirac mass matrix (SM ↔ Dark sector coupling)
    - M_R: Majorana mass matrix (heavy right-handed neutrinos)
    """
    
    def __init__(self):
        """Initialize E8 geometry for neutrino analysis."""
        self.roots_8d = self._generate_e8_roots()
        self.roots_4d = self.roots_8d @ UNIVERSE_MATRIX.T
        self.lengths_4d = np.sqrt(np.sum(self.roots_4d**2, axis=1))
        
        # Classify sectors
        sorted_idx = np.argsort(self.lengths_4d)
        self.sm_indices = sorted_idx[:12]
        self.dark_indices = sorted_idx[12:]
        
        self.sm_roots = self.roots_4d[self.sm_indices]
        self.dark_roots = self.roots_4d[self.dark_indices]
        
        print("Neutrino Sector Initialized:")
        print(f"  Standard Model roots: {len(self.sm_indices)}")
        print(f"  Dark sector roots: {len(self.dark_indices)}")
    
    def _generate_e8_roots(self) -> np.ndarray:
        """Generate all 240 E8 roots."""
        roots = []
        
        # Type 1: Simple roots
        for i in range(8):
            for j in range(i + 1, 8):
                for s1, s2 in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                    r = np.zeros(8)
                    r[i], r[j] = s1, s2
                    roots.append(r)
        
        # Type 2: Half-integer roots
        for signs in product([0.5, -0.5], repeat=8):
            if sum(1 for s in signs if s < 0) % 2 == 0:
                roots.append(np.array(signs))
        
        return np.array(roots)
    
    def find_right_handed_neutrinos(self) -> List[Dict]:
        """
        Search dark sector for heavy neutral singlets.
        
        Criteria for right-handed neutrinos:
        1. Neutral under U(1)_EM (no electric charge)
        2. Singlet under SU(3)_color (no color charge)
        3. Heavy mass (dark sector)
        4. Couples weakly to SM left-handed leptons
        """
        print("\n" + "="*70)
        print("SEARCHING FOR RIGHT-HANDED NEUTRINOS")
        print("="*70)
        print("\nCriteria:")
        print("  1. Neutral (no EM coupling)")
        print("  2. Singlet (no color charge)")
        print("  3. Heavy (M_R ~ TeV scale from dark sector)")
        print("  4. Weak coupling to SM leptons")
        print("-" * 70)
        
        candidates = []
        
        # Check each dark sector root
        for i, dark_idx in enumerate(self.dark_indices[:100]):
            dark_root = self.dark_roots[i]
            
            # 1. Check neutrality (weak coupling to photon)
            # Photon is typically the first SM root
            photon_coupling = abs(np.dot(dark_root, self.sm_roots[0]))
            
            # 2. Check singlet (weak coupling to gluons)
            # Gluons are roots 1-8 typically
            gluon_couplings = [abs(np.dot(dark_root, self.sm_roots[j])) 
                              for j in range(1, min(9, len(self.sm_roots)))]
            avg_gluon = np.mean(gluon_couplings) if gluon_couplings else 0
            
            # 3. Mass (from 4D projection length)
            mass_scale = self.lengths_4d[dark_idx]
            
            # 4. Lepton coupling (should be non-zero but weak)
            # Check coupling to electroweak sector
            ew_couplings = [abs(np.dot(dark_root, self.sm_roots[j])) 
                           for j in range(9, len(self.sm_roots))]
            avg_ew = np.mean(ew_couplings) if ew_couplings else 0
            
            # Candidate selection
            is_neutral = photon_coupling < 0.05
            is_singlet = avg_gluon < 0.1
            is_heavy = mass_scale > 0.5
            has_weak_coupling = 0.01 < avg_ew < 0.5
            
            if is_neutral and is_singlet and is_heavy and has_weak_coupling:
                candidates.append({
                    'index': dark_idx,
                    'mass_scale': mass_scale,
                    'em_coupling': photon_coupling,
                    'color_coupling': avg_gluon,
                    'ew_coupling': avg_ew,
                    'root_4d': dark_root
                })
        
        # Sort by mass scale
        candidates.sort(key=lambda x: x['mass_scale'])
        
        print(f"\n[OK] Found {len(candidates)} right-handed neutrino candidates")
        
        if len(candidates) >= 3:
            print("\nTop 3 Heavy Neutrino Candidates:")
            print(f"{'#':<3} {'Index':<8} {'Mass Scale':<12} {'EM':<10} {'Color':<10} {'EW':<10}")
            print("-" * 70)
            for i, c in enumerate(candidates[:3], 1):
                print(f"{i:<3} {c['index']:<8} {c['mass_scale']:<12.6f} "
                      f"{c['em_coupling']:<10.6f} {c['color_coupling']:<10.6f} "
                      f"{c['ew_coupling']:<10.6f}")
        
        return candidates
    
    def calculate_dirac_masses(self, rh_neutrinos: List[Dict]) -> np.ndarray:
        """
        Calculate Dirac mass matrix m_D.
        
        m_D connects left-handed SM leptons to right-handed neutrinos.
        Geometric origin: overlap between lepton roots and RH neutrino roots.
        """
        if len(rh_neutrinos) < 3:
            print("[!] Need at least 3 RH neutrinos for see-saw")
            return None
        
        print("\n" + "="*70)
        print("DIRAC MASS MATRIX CALCULATION")
        print("="*70)
        
        # Assume 3 generations of left-handed leptons
        # These come from fermion shells in SM sector
        n_gen = 3
        m_D = np.zeros((n_gen, n_gen))
        
        # Use geometric overlap as mass coupling
        for i in range(n_gen):
            for j in range(n_gen):
                if j < len(rh_neutrinos):
                    # Coupling proportional to root overlap
                    rh_root = rh_neutrinos[j]['root_4d']
                    # Use SM roots as proxies for lepton doublets
                    lh_root = self.sm_roots[9 + i] if (9 + i) < len(self.sm_roots) else self.sm_roots[i]
                    
                    overlap = abs(np.dot(lh_root, rh_root))
                    m_D[i, j] = overlap * 100  # GeV scale (typical Dirac mass)
        
        print("\nDirac Mass Matrix m_D (GeV):")
        print(m_D)
        
        return m_D
    
    def calculate_majorana_masses(self, rh_neutrinos: List[Dict]) -> np.ndarray:
        """
        Calculate Majorana mass matrix M_R for heavy neutrinos.
        
        M_R ~ geometric mass scale of dark sector roots.
        """
        if len(rh_neutrinos) < 3:
            return None
        
        print("\n" + "="*70)
        print("MAJORANA MASS MATRIX CALCULATION")
        print("="*70)
        
        # Majorana masses from geometric lengths
        # Scale to TeV (typical GUT/seesaw scale)
        n_gen = 3
        M_R = np.zeros((n_gen, n_gen))
        
        for i in range(n_gen):
            if i < len(rh_neutrinos):
                M_R[i, i] = rh_neutrinos[i]['mass_scale'] * 1e12  # TeV -> eV
        
        print("\nMajorana Mass Matrix M_R (eV):")
        print(M_R)
        print(f"\nTypical scale: {np.mean(np.diag(M_R)):.2e} eV")
        
        return M_R
    
    def apply_seesaw(self, m_D: np.ndarray, M_R: np.ndarray) -> Dict:
        """
        Apply Type-I see-saw formula to get light neutrino masses.
        
        m_nu = -m_D^T M_R^(-1) m_D
        """
        if m_D is None or M_R is None:
            return None
        
        print("\n" + "="*70)
        print("TYPE-I SEE-SAW MECHANISM")
        print("="*70)
        
        # Invert Majorana matrix
        M_R_inv = np.linalg.inv(M_R)
        
        # See-saw formula
        m_nu = -m_D.T @ M_R_inv @ m_D
        
        # Diagonalize to get mass eigenvalues
        eigenvalues, eigenvectors = np.linalg.eigh(m_nu)
        masses = np.abs(eigenvalues)
        
        print("\nLight Neutrino Masses (eV):")
        for i, m in enumerate(masses, 1):
            print(f"  m_nu{i} = {m:.6e} eV")
        
        # Calculate mass differences
        if len(masses) >= 3:
            delta_m21_sq = masses[1]**2 - masses[0]**2
            delta_m31_sq = masses[2]**2 - masses[0]**2
            
            print("\nMass Differences:")
            print(f"  Deltam^2_21 = {delta_m21_sq:.6e} eV^2")
            print(f"  Deltam^2_31 = {delta_m31_sq:.6e} eV^2")
            
            print("\nExperimental Targets:")
            print(f"  Deltam^2_21 = {NEUTRINO_DATA['mass_differences']['delta_m21_sq']:.6e} eV^2")
            print(f"  Deltam^2_31 = {NEUTRINO_DATA['mass_differences']['delta_m31_sq']:.6e} eV^2")
            
            # Calculate errors
            error_21 = abs(delta_m21_sq - NEUTRINO_DATA['mass_differences']['delta_m21_sq']) / NEUTRINO_DATA['mass_differences']['delta_m21_sq'] * 100
            error_31 = abs(delta_m31_sq - NEUTRINO_DATA['mass_differences']['delta_m31_sq']) / NEUTRINO_DATA['mass_differences']['delta_m31_sq'] * 100
            
            print(f"\nErrors:")
            print(f"  Deltam^2_21 error: {error_21:.2f}%")
            print(f"  Deltam^2_31 error: {error_31:.2f}%")
        
        return {
            'masses': masses,
            'eigenvectors': eigenvectors,
            'mass_matrix': m_nu,
            'delta_m21_sq': delta_m21_sq if len(masses) >= 3 else 0,
            'delta_m31_sq': delta_m31_sq if len(masses) >= 3 else 0
        }
    
    def derive_pmns_matrix(self, eigenvectors: np.ndarray) -> Dict:
        """
        Derive PMNS (Pontecorvo-Maki-Nakagawa-Sakata) mixing matrix.
        
        The PMNS matrix U relates flavor eigenstates to mass eigenstates.
        Parameterization:
            U = U_23(theta_23) U_13(theta_13, delta) U_12(theta_12)
        """
        print("\n" + "="*70)
        print("PMNS MIXING MATRIX DERIVATION")
        print("="*70)
        
        # PMNS matrix is the eigenvector matrix (up to phases)
        U_pmns = eigenvectors
        
        # Extract mixing angles from matrix elements
        # Standard parameterization
        if U_pmns.shape[0] >= 3:
            # theta_23: atmospheric mixing
            theta_23 = np.arctan2(abs(U_pmns[2, 1]), abs(U_pmns[2, 2]))
            
            # theta_13: reactor mixing  
            theta_13 = np.arcsin(abs(U_pmns[0, 2]))
            
            # theta_12: solar mixing
            theta_12 = np.arctan2(abs(U_pmns[0, 1]), abs(U_pmns[0, 0]))
            
            # Convert to degrees
            theta_12_deg = np.degrees(theta_12)
            theta_23_deg = np.degrees(theta_23)
            theta_13_deg = np.degrees(theta_13)
            
            print("\nPMNS Mixing Angles:")
            print(f"  theta_12 (solar)      = {theta_12_deg:.2f} deg")
            print(f"  theta_23 (atmospheric) = {theta_23_deg:.2f} deg")
            print(f"  theta_13 (reactor)    = {theta_13_deg:.2f} deg")
            
            print("\nExperimental Values:")
            print(f"  theta_12 = {NEUTRINO_DATA['pmns_angles']['theta_12']:.2f} deg")
            print(f"  theta_23 = {NEUTRINO_DATA['pmns_angles']['theta_23']:.2f} deg")
            print(f"  theta_13 = {NEUTRINO_DATA['pmns_angles']['theta_13']:.2f} deg")
            
            # Calculate errors
            error_12 = abs(theta_12_deg - NEUTRINO_DATA['pmns_angles']['theta_12'])
            error_23 = abs(theta_23_deg - NEUTRINO_DATA['pmns_angles']['theta_23'])
            error_13 = abs(theta_13_deg - NEUTRINO_DATA['pmns_angles']['theta_13'])
            
            print(f"\nErrors:")
            print(f"  theta_12 error: {error_12:.2f} deg")
            print(f"  theta_23 error: {error_23:.2f} deg")
            print(f"  theta_13 error: {error_13:.2f} deg")
            
            return {
                'pmns_matrix': U_pmns,
                'theta_12': theta_12_deg,
                'theta_23': theta_23_deg,
                'theta_13': theta_13_deg,
                'error_12': error_12,
                'error_23': error_23,
                'error_13': error_13
            }
        
        return {'pmns_matrix': U_pmns}
    
    def run_complete_analysis(self) -> Dict:
        """Run complete neutrino sector analysis."""
        print("\n" + "#"*70)
        print("#" + " "*68 + "#")
        print("#" + " "*15 + "NEUTRINO SECTOR ANALYSIS" + " "*29 + "#")
        print("#" + " "*68 + "#")
        print("#"*70)
        
        results = {}
        
        # 1. Find RH neutrinos
        rh_neutrinos = self.find_right_handed_neutrinos()
        results['rh_neutrinos'] = len(rh_neutrinos)
        
        if len(rh_neutrinos) < 3:
            print("\n[!] Insufficient RH neutrinos found. Need 3 generations.")
            return results
        
        # 2. Calculate Dirac masses
        m_D = self.calculate_dirac_masses(rh_neutrinos)
        
        # 3. Calculate Majorana masses
        M_R = self.calculate_majorana_masses(rh_neutrinos)
        
        # 4. Apply see-saw
        seesaw_results = self.apply_seesaw(m_D, M_R)
        if seesaw_results:
            results.update(seesaw_results)
        
        # 5. Derive PMNS matrix
        if seesaw_results and 'eigenvectors' in seesaw_results:
            pmns_results = self.derive_pmns_matrix(seesaw_results['eigenvectors'])
            results.update(pmns_results)
        
        print("\n" + "="*70)
        print("NEUTRINO SECTOR SUMMARY")
        print("="*70)
        print(f"[OK] Right-handed neutrinos: {len(rh_neutrinos)} found")
        if 'masses' in results:
            print(f"[OK] Light neutrino masses: {len(results['masses'])} derived")
        if 'theta_12' in results:
            print(f"[OK] PMNS angles calculated")
        print("="*70)
        
        return results


def main():
    """Run neutrino sector analysis."""
    neutrino = NeutrinoSector()
    results = neutrino.run_complete_analysis()
    
    print("\n[ACHIEVEMENT] NEUTRINO ANALYSIS COMPLETE!")
    print("="*70)
    
    return results


if __name__ == "__main__":
    main()
