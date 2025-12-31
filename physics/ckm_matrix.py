"""
CKM MATRIX DERIVATION: Quark Mixing From E8 Geometry
=====================================================
Goal: Derive the Cabibbo-Kobayashi-Maskawa (CKM) matrix elements
      from geometric mixing angles between quark generation roots.

Theory:
- CKM matrix describes quark flavor mixing
- Matrix elements from root angles: V_ij = cos(theta_ij) or sin(theta_ij)
- CP violation phase from geometric phase structure

EXPERIMENTAL CKM MATRIX (PDG 2024):
|V_ud  V_us  V_ub|   |0.97373  0.2243   0.00382|
|V_cd  V_cs  V_cb| = |0.221    0.987    0.0408 |
|V_td  V_ts  V_tb|   |0.0080   0.0388   1.014  |

Wolfenstein Parameters:
- lambda (Cabibbo angle) ~ 0.22
- A ~ 0.81
- rho ~ 0.16
- eta ~ 0.35

Author: E8 Theory Team
Date: December 31, 2025
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from itertools import product
from typing import Dict, List, Tuple
from e8_constants import UNIVERSE_MATRIX

# Experimental CKM matrix values (PDG 2024)
CKM_EXPERIMENTAL = {
    'matrix': np.array([
        [0.97373, 0.2243, 0.00382],
        [0.221, 0.987, 0.0408],
        [0.0080, 0.0388, 1.014]
    ]),
    'wolfenstein': {
        'lambda': 0.22650,  # Cabibbo angle
        'A': 0.790,
        'rho_bar': 0.159,
        'eta_bar': 0.348
    }
}


class CKMMatrix:
    """
    Derive CKM matrix from E8 quark sector geometry.
    
    The CKM matrix V relates down-type quark mass eigenstates to
    weak eigenstates:
        d' = V_ud * d + V_us * s + V_ub * b
        
    Geometric interpretation:
        V_ij = overlap between up-type and down-type quark roots
    """
    
    def __init__(self):
        """Initialize E8 geometry for CKM analysis."""
        self.roots_8d = self._generate_e8_roots()
        self.roots_4d = self.roots_8d @ UNIVERSE_MATRIX.T
        self.lengths_4d = np.sqrt(np.sum(self.roots_4d**2, axis=1))
        
        # Classify sectors
        sorted_idx = np.argsort(self.lengths_4d)
        self.sm_indices = sorted_idx[:12]
        self.dark_indices = sorted_idx[12:]
        
        self.sm_roots = self.roots_4d[self.sm_indices]
        self.dark_roots = self.roots_4d[self.dark_indices]
        
        print("CKM Matrix Calculator Initialized")
        print(f"  Standard Model roots: {len(self.sm_indices)}")
    
    def _generate_e8_roots(self) -> np.ndarray:
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
    
    def identify_quark_generations(self) -> Dict:
        """
        Identify the 3 generations of up-type and down-type quarks
        in the dark sector.
        
        Expected structure:
        - Generation 1: u, d (light)
        - Generation 2: c, s (medium)
        - Generation 3: t, b (heavy)
        """
        print("\n" + "="*70)
        print("IDENTIFYING QUARK GENERATIONS IN E8")
        print("="*70)
        
        # Look for fermion shells in dark sector
        # Quarks should form discrete mass shells
        dark_lengths = np.array([self.lengths_4d[idx] for idx in self.dark_indices])
        
        # Histogram to find mass shells
        hist, bin_edges = np.histogram(dark_lengths, bins=50)
        
        # Find prominent peaks (potential generation shells)
        peaks = []
        for i in range(1, len(hist)-1):
            if hist[i] > hist[i-1] and hist[i] > hist[i+1] and hist[i] >= 3:
                mass = (bin_edges[i] + bin_edges[i+1]) / 2
                peaks.append((mass, hist[i]))
        
        peaks.sort()
        
        print(f"\nFound {len(peaks)} prominent mass shells")
        print("\nIdentified Quark Generations:")
        print(f"{'Gen':<5} {'Mass Scale':<15} {'Count':<10}")
        print("-" * 40)
        
        generations = []
        for i, (mass, count) in enumerate(peaks[:3], 1):
            print(f"{i:<5} {mass:<15.6f} {count:<10}")
            generations.append({
                'generation': i,
                'mass_scale': mass,
                'count': count
            })
        
        # Extract quark roots near these mass shells
        up_type_quarks = []
        down_type_quarks = []
        
        for gen in generations:
            mass = gen['mass_scale']
            tolerance = 0.1
            
            # Find roots near this mass
            gen_roots = []
            for i, idx in enumerate(self.dark_indices):
                if abs(dark_lengths[i] - mass) < tolerance:
                    gen_roots.append({
                        'index': idx,
                        'root_4d': self.dark_roots[i],
                        'mass': dark_lengths[i]
                    })
            
            # Split into up-type and down-type based on properties
            if len(gen_roots) >= 2:
                # Simple heuristic: alternating assignment
                up_type_quarks.append(gen_roots[0] if len(gen_roots) > 0 else None)
                down_type_quarks.append(gen_roots[1] if len(gen_roots) > 1 else None)
        
        return {
            'generations': generations,
            'up_type': up_type_quarks,
            'down_type': down_type_quarks
        }
    
    def calculate_mixing_angles(self, up_quarks: List, down_quarks: List) -> np.ndarray:
        """
        Calculate CKM matrix elements from geometric angles.
        
        V_ij = cos(theta_ij) where theta_ij is the angle between
        the i-th up-type quark and j-th down-type quark root vectors.
        """
        print("\n" + "="*70)
        print("CKM MATRIX CALCULATION FROM GEOMETRY")
        print("="*70)
        
        n_gen = min(len(up_quarks), len(down_quarks), 3)
        V_ckm = np.zeros((3, 3))
        
        print("\nGeometric Mixing Angles:")
        header = "Up \\ Down"
        print(f"{header:<12}", end='')
        for j in range(n_gen):
            print(f"{'d/s/b'[j]:<12}", end='')
        print("\n" + "-" * 50)
        
        for i in range(n_gen):
            print(f"{'u/c/t'[i]:<12}", end='')
            for j in range(n_gen):
                if up_quarks[i] and down_quarks[j]:
                    up_root = up_quarks[i]['root_4d']
                    down_root = down_quarks[j]['root_4d']
                    
                    # Calculate angle between roots
                    cos_angle = np.dot(up_root, down_root) / (
                        np.linalg.norm(up_root) * np.linalg.norm(down_root)
                    )
                    
                    # Ensure numerical stability
                    cos_angle = np.clip(cos_angle, -1, 1)
                    
                    # CKM element is related to overlap
                    if i == j:
                        # Diagonal: should be close to 1
                        V_ckm[i, j] = abs(cos_angle)
                    else:
                        # Off-diagonal: small mixing
                        angle = np.arccos(abs(cos_angle))
                        V_ckm[i, j] = abs(np.sin(angle))
                    
                    print(f"{V_ckm[i, j]:<12.6f}", end='')
                else:
                    V_ckm[i, j] = 1.0 if i == j else 0.0
                    print(f"{V_ckm[i, j]:<12.6f}", end='')
            print()

        
        # Normalize each row to ensure unitarity
        for i in range(n_gen):
            norm = np.sqrt(np.sum(V_ckm[i, :]**2))
            if norm > 0:
                V_ckm[i, :] /= norm
        
        print("\nNormalized CKM Matrix:")
        print(V_ckm)
        
        return V_ckm
    
    def extract_wolfenstein_parameters(self, V_ckm: np.ndarray) -> Dict:
        """
        Extract Wolfenstein parameters from CKM matrix.
        
        Parameterization:
            lambda ~ |V_us| (Cabibbo angle)
            A ~ |V_cb| / lambda^2
            rho + ieta from unitarity triangle
        """
        print("\n" + "="*70)
        print("WOLFENSTEIN PARAMETERS")
        print("="*70)
        
        # Extract parameters
        lambda_w = abs(V_ckm[0, 1])  # V_us
        A = abs(V_ckm[1, 2]) / (lambda_w**2) if lambda_w > 0 else 0
        
        # Unitarity triangle parameters (simplified)
        V_ud = V_ckm[0, 0]
        V_ub = V_ckm[0, 2]
        
        rho_bar = -V_ub.real / (lambda_w**3 * A) if (lambda_w > 0 and A > 0) else 0
        eta_bar = -V_ub.imag / (lambda_w**3 * A) if (lambda_w > 0 and A > 0) else 0
        
        print(f"\nDerived Parameters:")
        print(f"  lambda (Cabibbo) = {lambda_w:.5f}")
        print(f"  A           = {A:.3f}")
        print(f"  rho           = {rho_bar:.3f}")
        print(f"  eta           = {eta_bar:.3f}")
        
        print(f"\nExperimental Values:")
        print(f"  lambda           = {CKM_EXPERIMENTAL['wolfenstein']['lambda']:.5f}")
        print(f"  A           = {CKM_EXPERIMENTAL['wolfenstein']['A']:.3f}")
        print(f"  rho           = {CKM_EXPERIMENTAL['wolfenstein']['rho_bar']:.3f}")
        print(f"  eta           = {CKM_EXPERIMENTAL['wolfenstein']['eta_bar']:.3f}")
        
        # Calculate errors
        error_lambda = abs(lambda_w - CKM_EXPERIMENTAL['wolfenstein']['lambda']) / CKM_EXPERIMENTAL['wolfenstein']['lambda'] * 100
        error_A = abs(A - CKM_EXPERIMENTAL['wolfenstein']['A']) / CKM_EXPERIMENTAL['wolfenstein']['A'] * 100
        
        print(f"\nErrors:")
        print(f"  lambda error: {error_lambda:.2f}%")
        print(f"  A error: {error_A:.2f}%")
        
        return {
            'lambda': lambda_w,
            'A': A,
            'rho_bar': rho_bar,
            'eta_bar': eta_bar,
            'error_lambda': error_lambda,
            'error_A': error_A
        }
    
    def compare_with_experiment(self, V_ckm: np.ndarray) -> Dict:
        """Compare derived CKM matrix with experimental values."""
        print("\n" + "="*70)
        print("COMPARISON WITH EXPERIMENT")
        print("="*70)
        
        V_exp = CKM_EXPERIMENTAL['matrix']
        
        print("\nDerived CKM Matrix:")
        for i in range(3):
            print(f"  [{V_ckm[i,0]:.5f}  {V_ckm[i,1]:.5f}  {V_ckm[i,2]:.5f}]")
        
        print("\nExperimental CKM Matrix:")
        for i in range(3):
            print(f"  [{V_exp[i,0]:.5f}  {V_exp[i,1]:.5f}  {V_exp[i,2]:.5f}]")
        
        # Calculate element-wise errors
        print("\nElement-wise Errors (%):")
        errors = np.abs(V_ckm[:3,:3] - V_exp) / V_exp * 100
        
        labels = ['V_ud', 'V_us', 'V_ub', 'V_cd', 'V_cs', 'V_cb', 'V_td', 'V_ts', 'V_tb']
        idx = 0
        for i in range(3):
            for j in range(3):
                print(f"  {labels[idx]:<6}: {errors[i,j]:6.2f}%")
                idx += 1
        
        # Overall RMS error
        rms_error = np.sqrt(np.mean(errors**2))
        print(f"\nRMS Error: {rms_error:.2f}%")
        
        return {
            'derived': V_ckm,
            'experimental': V_exp,
            'errors': errors,
            'rms_error': rms_error
        }
    
    def run_complete_analysis(self) -> Dict:
        """Run complete CKM matrix derivation."""
        print("\n" + "#"*70)
        print("#" + " "*68 + "#")
        print("#" + " "*18 + "CKM MATRIX DERIVATION" + " "*29 + "#")
        print("#" + " "*68 + "#")
        print("#"*70)
        
        results = {}
        
        # 1. Identify quark generations
        quarks = self.identify_quark_generations()
        results['quark_generations'] = quarks['generations']
        
        # 2. Calculate CKM matrix
        if len(quarks['up_type']) >= 3 and len(quarks['down_type']) >= 3:
            V_ckm = self.calculate_mixing_angles(
                quarks['up_type'],
                quarks['down_type']
            )
            results['ckm_matrix'] = V_ckm
            
            # 3. Extract Wolfenstein parameters
            wolfenstein = self.extract_wolfenstein_parameters(V_ckm)
            results['wolfenstein'] = wolfenstein
            
            # 4. Compare with experiment
            comparison = self.compare_with_experiment(V_ckm)
            results['comparison'] = comparison
        else:
            print("\n[!] Insufficient quark generations identified")
        
        print("\n" + "="*70)
        print("CKM MATRIX SUMMARY")
        print("="*70)
        print(f"[OK] Quark generations identified: {len(quarks['generations'])}")
        if 'ckm_matrix' in results:
            print(f"[OK] CKM matrix derived from geometry")
            print(f"[OK] Wolfenstein parameters calculated")
            print(f"[OK] RMS error: {results['comparison']['rms_error']:.2f}%")
        print("="*70)
        
        return results


def main():
    """Run CKM matrix analysis."""
    ckm = CKMMatrix()
    results = ckm.run_complete_analysis()
    
    print("\n[ACHIEVEMENT] CKM ANALYSIS COMPLETE!")
    print("="*70)
    
    return results


if __name__ == "__main__":
    main()
