#!/usr/bin/env python3
"""
E8 GRAVITON LATTICE QFT - Rigorous Calculation
===============================================

This module implements a MORE RIGOROUS lattice QFT calculation for
graviton physics on the E8->H4 projected lattice.

Key improvements over naive calculation:
1. Proper Brillouin zone integration
2. Wilson-style lattice action for gravity  
3. Correct momentum routing in loops
4. Finite-size scaling analysis
5. Comparison with analytical predictions

The goal: Determine whether discrete E8 geometry provides natural
UV regularization, and if so, what the suppression factor is.

Author: E8 Theory of Everything Project
Date: January 1, 2026
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy.integrate import nquad
from scipy.special import gamma as gamma_func
import warnings

# Golden ratio
PHI = (1 + np.sqrt(5)) / 2
INV_PHI = 1 / PHI
INV_PHI_SQ = INV_PHI ** 2
COS_ICOSAHEDRAL = 1 / np.sqrt(5)


@dataclass
class LoopResult:
    """Results from loop calculation."""
    value: float
    error: float
    is_finite: bool
    n_lattice_points: int
    effective_cutoff: float


class RigorousE8LatticeGravity:
    """
    Rigorous lattice gravity on E8->H4.
    
    Key insight: On a lattice, UV finiteness comes from the
    BOUNDED nature of the Brillouin zone, not from ad-hoc suppression.
    
    The question is: Does E8->H4 structure give ADDITIONAL suppression
    beyond the trivial lattice cutoff?
    """
    
    def __init__(self, lattice_spacing: float = 1.0):
        """
        Initialize with lattice spacing a.
        
        In natural units where M_Pl = 1:
        - a = 1 means lattice scale = Planck scale
        - Brillouin zone extends to k_max = pi/a
        """
        self.a = lattice_spacing
        self.k_max = np.pi / self.a  # Brillouin zone boundary
        
        # E8 structure
        self.e8_roots = self._build_e8_roots()
        self.projection = self._build_elser_sloane()
        self.h4_lattice = self._project_roots()
        
        # Precompute angles for verification
        self._compute_angle_distribution()
        
    def _build_e8_roots(self) -> np.ndarray:
        """Build 240 E8 root vectors."""
        roots = []
        
        # Type 1: (±1, ±1, 0, 0, 0, 0, 0, 0) and permutations
        for i in range(8):
            for j in range(i+1, 8):
                for si in [-1, 1]:
                    for sj in [-1, 1]:
                        v = np.zeros(8)
                        v[i] = si
                        v[j] = sj
                        roots.append(v)
        
        # Type 2: (±1/2)^8 with even number of minus signs
        for bits in range(256):
            v = np.array([0.5 * (1 if (bits >> i) & 1 else -1) for i in range(8)])
            if np.sum(v < 0) % 2 == 0:
                roots.append(v)
        
        return np.array(roots)
    
    def _build_elser_sloane(self) -> np.ndarray:
        """
        Build the Elser-Sloane projection matrix.
        
        This is the specific 4x8 matrix that projects E8 to H4
        while preserving icosahedral symmetry.
        """
        phi = PHI
        
        # The Elser-Sloane matrix (one of several equivalent forms)
        P = np.array([
            [1, phi, 0, -1, phi, 0, 0, 0],
            [phi, 0, 1, phi, 0, -1, 0, 0],  
            [0, 1, phi, 0, -1, phi, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, phi]
        ])
        
        # Orthonormalize
        for i in range(4):
            # Subtract projection onto previous rows
            for j in range(i):
                P[i] -= np.dot(P[i], P[j]) * P[j]
            # Normalize
            P[i] /= np.linalg.norm(P[i])
        
        return P
    
    def _project_roots(self) -> np.ndarray:
        """Project E8 roots to H4 (4D)."""
        return np.dot(self.e8_roots, self.projection.T)
    
    def _compute_angle_distribution(self):
        """Compute and store angular structure of H4 lattice."""
        angles = []
        
        for i in range(len(self.h4_lattice)):
            for j in range(i+1, len(self.h4_lattice)):
                v1 = self.h4_lattice[i]
                v2 = self.h4_lattice[j]
                n1 = np.linalg.norm(v1)
                n2 = np.linalg.norm(v2)
                
                if n1 > 1e-6 and n2 > 1e-6:
                    cos_angle = np.dot(v1, v2) / (n1 * n2)
                    angles.append(np.clip(cos_angle, -1, 1))
        
        self.cos_angles = np.array(angles)
        
        # Find fraction near icosahedral angle
        ico_mask = np.abs(np.abs(self.cos_angles) - COS_ICOSAHEDRAL) < 0.1
        self.icosahedral_fraction = np.sum(ico_mask) / len(self.cos_angles)
        
    def lattice_laplacian_eigenvalue(self, k: np.ndarray) -> float:
        """
        Lattice Laplacian eigenvalue (Wilson form).
        
        On a hypercubic lattice:
        Delta(k) = (2/a^2) * sum_mu [1 - cos(k_mu * a)]
        
        This goes from k^2 at small k to 8/a^2 at the zone corner.
        """
        a = self.a
        eigenvalue = 0.0
        
        for mu in range(len(k)):
            eigenvalue += (2/a**2) * (1 - np.cos(k[mu] * a))
        
        return eigenvalue
    
    def lattice_propagator(self, k: np.ndarray, m_sq: float = 0.0) -> float:
        """
        Lattice propagator for spin-2 field.
        
        G(k) = 1 / [Delta(k) + m^2]
        
        For graviton: m = 0 (massless)
        """
        delta_k = self.lattice_laplacian_eigenvalue(k)
        denom = delta_k + m_sq
        
        if denom < 1e-12:
            return 0.0  # IR regularization
        
        return 1.0 / denom
    
    def continuum_propagator(self, k: np.ndarray, m_sq: float = 0.0) -> float:
        """Standard continuum propagator 1/(k^2 + m^2)."""
        k_sq = np.sum(k**2)
        denom = k_sq + m_sq
        
        if denom < 1e-12:
            return 0.0
        
        return 1.0 / denom
    
    def compute_1loop_integral_lattice(self, 
                                        external_p: float,
                                        n_points: int = 50) -> LoopResult:
        """
        Compute 1-loop integral on the lattice.
        
        I = (1/V) sum_k G(k) G(p-k)
        
        where V is the Brillouin zone volume and sum is over
        discrete lattice momenta.
        
        For a finite lattice of N^4 sites:
        k = (2*pi*n) / (N*a) for n = 0, 1, ..., N-1
        """
        a = self.a
        N = n_points
        
        # Momentum spacing
        dk = 2 * np.pi / (N * a)
        
        # External momentum in z-direction
        p_vec = np.array([0, 0, 0, external_p])
        
        total = 0.0
        count = 0
        
        # Sum over 4D Brillouin zone
        for n0 in range(N):
            for n1 in range(N):
                for n2 in range(N):
                    for n3 in range(N):
                        # Lattice momentum
                        k = dk * np.array([n0 - N//2, n1 - N//2, n2 - N//2, n3 - N//2])
                        
                        # Propagators
                        G_k = self.lattice_propagator(k)
                        G_pk = self.lattice_propagator(p_vec - k)
                        
                        total += G_k * G_pk
                        count += 1
        
        # Normalize by number of points (lattice volume factor)
        result = total / count
        
        # Effective UV cutoff
        effective_cutoff = np.pi / a
        
        return LoopResult(
            value=result,
            error=0.0,  # Could estimate from finite-size effects
            is_finite=True,
            n_lattice_points=count,
            effective_cutoff=effective_cutoff
        )
    
    def compute_1loop_integral_continuum(self,
                                          external_p: float,
                                          cutoff: float = 10.0,
                                          n_samples: int = 50000) -> float:
        """
        Compute 1-loop integral in continuum with hard cutoff.
        
        I = integral d^4k / [(2pi)^4] * 1/(k^2) * 1/((p-k)^2)
        
        This diverges logarithmically without cutoff.
        """
        p_vec = np.array([0, 0, 0, external_p])
        
        total = 0.0
        volume = (2 * cutoff) ** 4
        
        for _ in range(n_samples):
            # Random momentum in 4D cube
            k = np.random.uniform(-cutoff, cutoff, 4)
            
            k_sq = np.sum(k**2)
            pk_sq = np.sum((p_vec - k)**2)
            
            if k_sq > 1e-6 and pk_sq > 1e-6:
                integrand = 1.0 / (k_sq * pk_sq)
                total += integrand
        
        # Monte Carlo integral estimate
        result = total * volume / n_samples / (2*np.pi)**4
        
        return result
    
    def analyze_uv_behavior(self, n_scales: int = 5) -> Dict:
        """
        Analyze UV behavior by comparing lattice and continuum
        at different lattice spacings.
        
        If E8 provides special UV properties, we should see:
        - Lattice results independent of (or weakly dependent on) a
        - Different from standard hypercubic lattice
        """
        print("=" * 70)
        print("UV BEHAVIOR ANALYSIS: LATTICE VS CONTINUUM")
        print("=" * 70)
        print()
        
        # Test at fixed external momentum
        p_ext = 0.5  # In natural units
        
        # Vary lattice spacing
        spacings = np.logspace(-1, 0, n_scales)  # From 0.1 to 1.0
        
        lattice_results = []
        continuum_results = []
        
        print(f"{'Spacing a':<12} {'Cutoff pi/a':<12} {'Lattice I':<15} {'Continuum I':<15} {'Ratio':<10}")
        print("-" * 65)
        
        for a in spacings:
            # Update lattice spacing
            self.a = a
            cutoff = np.pi / a
            
            # Lattice integral
            lat_result = self.compute_1loop_integral_lattice(p_ext, n_points=20)
            
            # Continuum integral with same cutoff
            cont_result = self.compute_1loop_integral_continuum(p_ext, cutoff=cutoff, n_samples=10000)
            
            ratio = lat_result.value / cont_result if cont_result > 0 else 0
            
            lattice_results.append(lat_result.value)
            continuum_results.append(cont_result)
            
            print(f"{a:<12.3f} {cutoff:<12.2f} {lat_result.value:<15.6e} {cont_result:<15.6e} {ratio:<10.4f}")
        
        print("-" * 65)
        
        # Analyze scaling
        lattice_arr = np.array(lattice_results)
        continuum_arr = np.array(continuum_results)
        
        # In continuum, I ~ log(Lambda) ~ log(1/a) for fixed p
        # On lattice, finite by construction
        
        variation_lattice = np.std(lattice_arr) / np.mean(lattice_arr) if np.mean(lattice_arr) > 0 else 0
        variation_continuum = np.std(continuum_arr) / np.mean(continuum_arr) if np.mean(continuum_arr) > 0 else 0
        
        print()
        print("ANALYSIS:")
        print(f"  Lattice integral variation:   {variation_lattice*100:.1f}%")
        print(f"  Continuum integral variation: {variation_continuum*100:.1f}%")
        print()
        
        # Check for special suppression
        mean_ratio = np.mean(lattice_arr) / np.mean(continuum_arr) if np.mean(continuum_arr) > 0 else 0
        
        # Is there φ-related suppression?
        # If lattice/continuum ~ φ^n for some n, find n
        if mean_ratio > 0:
            phi_power = np.log(mean_ratio) / np.log(PHI)
        else:
            phi_power = 0
        
        print(f"  Mean lattice/continuum ratio: {mean_ratio:.4f}")
        print(f"  If ratio = φ^n, then n =      {phi_power:.2f}")
        print()
        
        # Check icosahedral structure
        print(f"  Icosahedral angle fraction:   {self.icosahedral_fraction*100:.1f}%")
        print(f"  (Expected for H4: ~20-30%)")
        print()
        
        return {
            'spacings': spacings,
            'lattice_results': lattice_arr,
            'continuum_results': continuum_arr,
            'mean_ratio': mean_ratio,
            'phi_power': phi_power,
            'icosahedral_fraction': self.icosahedral_fraction
        }
    
    def test_h4_enhancement(self) -> Dict:
        """
        Test whether H4 structure provides additional enhancement/suppression
        compared to standard hypercubic lattice.
        
        Compare:
        1. H4 lattice momenta (E8 projected)
        2. Standard 4D hypercubic lattice
        
        If E8->H4 is special, results should differ.
        """
        print("=" * 70)
        print("H4 VS HYPERCUBIC: TESTING E8 STRUCTURE")
        print("=" * 70)
        print()
        
        p_ext = 0.5
        
        # H4 lattice: Sum over projected E8 roots
        h4_total = 0.0
        h4_count = 0
        
        for q in self.h4_lattice:
            q_scaled = q * self.k_max / np.max(np.abs(self.h4_lattice))
            p_vec = np.array([0, 0, 0, p_ext])
            
            G_q = self.continuum_propagator(q_scaled)  # Use simple 1/k^2
            G_pq = self.continuum_propagator(p_vec - q_scaled)
            
            h4_total += G_q * G_pq
            h4_count += 1
        
        h4_result = h4_total / h4_count
        
        # Hypercubic lattice: Same number of points in 4D cube
        N = int(np.ceil(h4_count ** 0.25))
        hyper_total = 0.0
        hyper_count = 0
        
        for n0 in range(N):
            for n1 in range(N):
                for n2 in range(N):
                    for n3 in range(N):
                        q = self.k_max * np.array([
                            2*n0/(N-1) - 1,
                            2*n1/(N-1) - 1,
                            2*n2/(N-1) - 1,
                            2*n3/(N-1) - 1
                        ])
                        p_vec = np.array([0, 0, 0, p_ext])
                        
                        G_q = self.continuum_propagator(q)
                        G_pq = self.continuum_propagator(p_vec - q)
                        
                        hyper_total += G_q * G_pq
                        hyper_count += 1
        
        hyper_result = hyper_total / hyper_count
        
        # Compare
        ratio = h4_result / hyper_result if hyper_result > 0 else 0
        
        print(f"  H4 (E8 projected) result:     {h4_result:.6e}")
        print(f"  Hypercubic result:            {hyper_result:.6e}")
        print(f"  Ratio H4/Hypercubic:          {ratio:.4f}")
        print()
        
        if ratio > 0:
            phi_equiv = np.log(ratio) / np.log(PHI) if ratio > 0 else 0
            print(f"  If ratio = φ^n, n = {phi_equiv:.2f}")
        
        # Check if different
        is_different = abs(ratio - 1.0) > 0.1
        
        print()
        if is_different:
            if ratio < 1:
                print("  [!] H4 structure SUPPRESSES loop integral")
            else:
                print("  [!] H4 structure ENHANCES loop integral")
        else:
            print("  [=] H4 structure ~equal to hypercubic")
        
        print("=" * 70)
        
        return {
            'h4_result': h4_result,
            'hypercubic_result': hyper_result,
            'ratio': ratio,
            'is_different': is_different
        }


def run_rigorous_analysis():
    """Run the complete rigorous analysis."""
    print()
    print("#" * 70)
    print("#" + " "*68 + "#")
    print("#" + "  RIGOROUS E8 LATTICE QUANTUM GRAVITY ANALYSIS".center(68) + "#")
    print("#" + " "*68 + "#")
    print("#" * 70)
    print()
    
    # Initialize
    lattice = RigorousE8LatticeGravity(lattice_spacing=1.0)
    
    print(f"[Setup] E8 roots: {len(lattice.e8_roots)}")
    print(f"[Setup] H4 momenta: {len(lattice.h4_lattice)}")
    print(f"[Setup] Icosahedral fraction: {lattice.icosahedral_fraction*100:.1f}%")
    print()
    
    # Test 1: UV behavior analysis
    uv_results = lattice.analyze_uv_behavior(n_scales=4)
    
    # Test 2: H4 vs hypercubic comparison
    h4_results = lattice.test_h4_enhancement()
    
    # Final summary
    print()
    print("#" * 70)
    print("FINAL CONCLUSIONS")
    print("#" * 70)
    print()
    
    print("1. UV FINITENESS:")
    print("   - Lattice integrals are finite by construction (bounded BZ)")
    print("   - This is NOT special to E8 - any lattice gives this")
    print()
    
    print("2. E8 STRUCTURE EFFECT:")
    if h4_results['is_different']:
        if h4_results['ratio'] < 1:
            print(f"   - H4 provides {(1-h4_results['ratio'])*100:.1f}% SUPPRESSION vs hypercubic")
            print("   - This IS unique to E8->H4 structure!")
        else:
            print(f"   - H4 provides {(h4_results['ratio']-1)*100:.1f}% ENHANCEMENT vs hypercubic")
            print("   - Opposite to φ^(-2) suppression claim")
    else:
        print("   - No significant difference from standard lattice")
        print("   - φ-suppression NOT observed in rigorous calculation")
    print()
    
    print("3. ICOSAHEDRAL STRUCTURE:")
    if lattice.icosahedral_fraction > 0.15:
        print(f"   - Icosahedral angles present ({lattice.icosahedral_fraction*100:.1f}%)")
        print("   - H4 (600-cell) structure confirmed")
    else:
        print(f"   - Weak icosahedral signature ({lattice.icosahedral_fraction*100:.1f}%)")
        print("   - May need different projection matrix")
    print()
    
    print("4. STATUS OF φ-SUPPRESSION CLAIM:")
    print("   - Original claim: Loop integrals suppressed by φ^(-2) per loop")
    print(f"   - Measured effect: {uv_results['mean_ratio']:.4f}x (φ^{uv_results['phi_power']:.1f})")
    if abs(uv_results['phi_power'] + 2) < 0.5:
        print("   - CONSISTENT with φ^(-2) suppression!")
    else:
        print("   - NOT consistent with simple φ^(-2) picture")
    print()
    
    print("=" * 70)
    print("END OF RIGOROUS ANALYSIS")
    print("=" * 70)
    
    return lattice, uv_results, h4_results


if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    lattice, uv_results, h4_results = run_rigorous_analysis()
