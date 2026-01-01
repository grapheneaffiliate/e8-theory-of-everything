"""
E8 IDS SPECTRAL PROOF
=====================
Route A: Prove suppression via Integrated Density of States (IDS)

This script:
1. Builds a finite patch of E8->H4 quasicrystal Lambda
2. Constructs the graph Laplacian Delta_Lambda
3. Computes eigenvalues (spectrum)
4. Estimates IDS prefactor sigma via Weyl law: N(E) ~ sigma * C * E^2
5. Compares sigma to phi^(-12) numerically

The IDS prefactor IS the UV phase-space suppression in mathematically rigorous terms.

Author: E8 Theory Project
Date: January 2026
"""

import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh
from scipy.optimize import curve_fit
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONSTANTS
# =============================================================================

PHI = (1 + np.sqrt(5)) / 2
PHI_INV = 1 / PHI
PHI_M12 = PHI**(-12)


# =============================================================================
# STEP 1: E8 ROOT SYSTEM
# =============================================================================

def generate_e8_roots() -> np.ndarray:
    """
    Generate all 240 E8 root vectors.
    Type 1: all permutations of (+/-1, +/-1, 0, 0, 0, 0, 0, 0) - 112 roots
    Type 2: all (+/-1/2)^8 with even number of minus signs - 128 roots
    """
    roots = []
    
    # Type 1: (+/-1, +/-1, 0, 0, 0, 0, 0, 0) and permutations
    from itertools import combinations, product
    for positions in combinations(range(8), 2):
        for signs in product([1, -1], repeat=2):
            root = np.zeros(8)
            root[positions[0]] = signs[0]
            root[positions[1]] = signs[1]
            roots.append(root)
    
    # Type 2: (+/-1/2)^8 with even number of minus signs
    for signs in product([0.5, -0.5], repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            roots.append(np.array(signs))
    
    return np.array(roots)


# =============================================================================
# STEP 2: ELSER-SLOANE PROJECTION E8 -> H4
# =============================================================================

def elser_sloane_projection_matrix() -> np.ndarray:
    """
    The Elser-Sloane 4x8 projection matrix from E8 to H4.
    This is the mathematically canonical projection.
    """
    phi = PHI
    P = np.array([
        [phi, 1, 0, 0, 1, phi, 0, 0],
        [0, phi, 1, 0, 0, 0, phi, 1],
        [0, 0, phi, 1, phi, 0, 0, 1],
        [1, 0, 0, phi, 0, 1, phi, 0]
    ]) / np.sqrt(2 * (phi**2 + 1))
    
    return P


def project_to_4d(points_8d: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Project 8D points to 4D via P."""
    return points_8d @ P.T


# =============================================================================
# STEP 3: CUT-AND-PROJECT QUASICRYSTAL
# =============================================================================

def build_quasicrystal_patch(n_layers: int = 3, window_radius: float = 2.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build a finite patch of the E8->H4 quasicrystal using cut-and-project.
    
    The quasicrystal Lambda consists of projected points from E8 lattice
    whose perpendicular component falls within the acceptance window W.
    
    Args:
        n_layers: Number of E8 lattice shells to consider
        window_radius: Radius of acceptance window in perp space
        
    Returns:
        Lambda_par: 4D quasicrystal points (parallel space)
        Lambda_perp: 4D perpendicular components (for verification)
    """
    P = elser_sloane_projection_matrix()
    
    # Perpendicular projection (complement of P)
    # Use SVD to get orthogonal complement
    U, S, Vt = np.linalg.svd(P)
    P_perp = Vt[4:, :]  # Last 4 rows of Vt span perpendicular space
    
    # Generate E8 lattice points (roots + translations)
    roots = generate_e8_roots()
    
    # Build lattice points: sum of integer multiples of roots
    lattice_points = []
    
    # Start with origin
    lattice_points.append(np.zeros(8))
    
    # Add roots
    for root in roots:
        lattice_points.append(root)
    
    # Add sums of roots (up to n_layers)
    if n_layers >= 2:
        for i, r1 in enumerate(roots[:50]):  # Limit for efficiency
            for r2 in roots[:50]:
                if not np.allclose(r1, r2):
                    lattice_points.append(r1 + r2)
    
    if n_layers >= 3:
        for r1 in roots[:20]:
            for r2 in roots[:20]:
                for r3 in roots[:20]:
                    lattice_points.append(r1 + r2 + r3)
    
    lattice_points = np.array(lattice_points)
    
    # Remove duplicates
    lattice_points = np.unique(np.round(lattice_points, 10), axis=0)
    
    # Project to parallel and perpendicular spaces
    par = lattice_points @ P.T
    perp = lattice_points @ P_perp.T
    
    # Cut-and-project: keep only points in acceptance window
    perp_norms = np.linalg.norm(perp, axis=1)
    mask = perp_norms <= window_radius
    
    Lambda_par = par[mask]
    Lambda_perp = perp[mask]
    
    # Remove duplicates in parallel space
    unique_idx = np.unique(np.round(Lambda_par, 10), axis=0, return_index=True)[1]
    Lambda_par = Lambda_par[unique_idx]
    Lambda_perp = Lambda_perp[unique_idx]
    
    return Lambda_par, Lambda_perp


# =============================================================================
# STEP 4: GRAPH LAPLACIAN ON QUASICRYSTAL
# =============================================================================

def build_delaunay_graph_laplacian(points: np.ndarray, cutoff: float = None) -> csr_matrix:
    """
    Build the graph Laplacian for points using Delaunay triangulation
    or nearest-neighbor connectivity.
    
    L = D - A where:
        D = degree matrix (diagonal)
        A = adjacency matrix
        
    For 4D, we use distance-based connectivity (Delaunay is expensive).
    """
    n = len(points)
    
    if cutoff is None:
        # Estimate cutoff as slightly larger than typical nearest neighbor
        if n > 1:
            dists = np.linalg.norm(points[1:] - points[0], axis=1)
            cutoff = np.percentile(dists, 20) * 1.5
        else:
            cutoff = 1.0
    
    # Build adjacency matrix
    A = lil_matrix((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(points[i] - points[j])
            if d < cutoff:
                A[i, j] = 1
                A[j, i] = 1
    
    A = A.tocsr()
    
    # Degree matrix
    degrees = np.array(A.sum(axis=1)).flatten()
    D = csr_matrix(np.diag(degrees))
    
    # Laplacian
    L = D - A
    
    return L


def compute_laplacian_spectrum(L: csr_matrix, k: int = None) -> np.ndarray:
    """
    Compute eigenvalues of the Laplacian.
    
    For small matrices, use full diagonalization.
    For large matrices, use sparse eigensolver for first k eigenvalues.
    """
    n = L.shape[0]
    
    if n < 500:
        # Full diagonalization
        eigenvalues = np.linalg.eigvalsh(L.toarray())
    else:
        # Sparse: get k smallest eigenvalues
        if k is None:
            k = min(n - 2, 200)
        eigenvalues, _ = eigsh(L, k=k, which='SM')
    
    return np.sort(eigenvalues)


# =============================================================================
# STEP 5: INTEGRATED DENSITY OF STATES (IDS)
# =============================================================================

def compute_ids(eigenvalues: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the Integrated Density of States (IDS).
    
    N(E) = (1/n) * #{eigenvalues <= E}
    
    For a d-dimensional system, Weyl law predicts:
        N(E) ~ sigma * C_d * E^(d/2)
        
    where sigma is the "spectral prefactor" we want to measure.
    """
    n = len(eigenvalues)
    E_values = np.sort(eigenvalues)
    N_values = np.arange(1, n + 1) / n
    
    return E_values, N_values


def fit_weyl_law(E: np.ndarray, N: np.ndarray, d: int = 4) -> Dict:
    """
    Fit IDS to Weyl law: N(E) = sigma * C * E^(d/2)
    
    For d=4: N(E) ~ sigma * C * E^2
    
    Returns the fitted prefactor sigma.
    """
    # Filter to positive eigenvalues only
    mask = E > 0.01
    E_pos = E[mask]
    N_pos = N[mask]
    
    if len(E_pos) < 10:
        return {'sigma': np.nan, 'error': 'Too few positive eigenvalues'}
    
    # Fit N = sigma * E^2 (d=4)
    # Take log: log(N) = log(sigma) + 2*log(E)
    log_E = np.log(E_pos)
    log_N = np.log(N_pos)
    
    # Linear fit
    coeffs = np.polyfit(log_E, log_N, 1)
    power = coeffs[0]  # Should be ~d/2 = 2
    log_sigma = coeffs[1]
    sigma = np.exp(log_sigma)
    
    # Also fit with fixed power = 2
    def weyl_model(E, sigma):
        return sigma * E**2
    
    try:
        popt, pcov = curve_fit(weyl_model, E_pos, N_pos, p0=[0.1])
        sigma_fixed = popt[0]
        sigma_err = np.sqrt(pcov[0, 0])
    except:
        sigma_fixed = sigma
        sigma_err = np.nan
    
    return {
        'sigma_free': sigma,
        'fitted_power': power,
        'sigma_fixed_d4': sigma_fixed,
        'sigma_error': sigma_err,
        'expected_power': d / 2
    }


# =============================================================================
# STEP 6: REFERENCE PERIODIC LATTICE
# =============================================================================

def build_periodic_lattice_patch(n_points: int) -> np.ndarray:
    """
    Build a periodic 4D hypercubic lattice patch for reference.
    """
    # Determine grid size
    n_per_dim = int(np.ceil(n_points**(1/4)))
    
    points = []
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            for k in range(n_per_dim):
                for l in range(n_per_dim):
                    points.append([i, j, k, l])
                    if len(points) >= n_points:
                        break
                if len(points) >= n_points:
                    break
            if len(points) >= n_points:
                break
        if len(points) >= n_points:
            break
    
    return np.array(points[:n_points], dtype=float)


def compute_reference_sigma(n_points: int) -> Dict:
    """
    Compute IDS prefactor for reference periodic lattice.
    """
    # Build periodic lattice
    periodic = build_periodic_lattice_patch(n_points)
    
    # Build Laplacian
    L = build_delaunay_graph_laplacian(periodic)
    
    # Compute spectrum
    eigenvalues = compute_laplacian_spectrum(L)
    
    # Compute IDS
    E, N = compute_ids(eigenvalues)
    
    # Fit Weyl law
    fit = fit_weyl_law(E, N, d=4)
    
    return {
        'n_points': n_points,
        'sigma': fit['sigma_fixed_d4'],
        'sigma_error': fit['sigma_error'],
        'eigenvalues': eigenvalues
    }


# =============================================================================
# MAIN: IDS-BASED SUPPRESSION PROOF
# =============================================================================

def run_ids_spectral_proof(n_layers: int = 3, window_radius: float = 2.0) -> Dict:
    """
    Main routine: compute suppression factor from IDS comparison.
    
    sigma_ratio = sigma_quasicrystal / sigma_periodic
    
    If the theory is correct: sigma_ratio ~ phi^(-12)
    """
    print("="*70)
    print("E8 IDS SPECTRAL PROOF")
    print("Route A: Compute suppression from Integrated Density of States")
    print("="*70)
    
    results = {}
    
    # Step 1: Build quasicrystal patch
    print("\n[Step 1] Building E8->H4 quasicrystal patch...")
    Lambda_par, Lambda_perp = build_quasicrystal_patch(n_layers, window_radius)
    n_qc = len(Lambda_par)
    print(f"  Quasicrystal points: {n_qc}")
    results['n_quasicrystal'] = n_qc
    
    if n_qc < 20:
        print("  ERROR: Too few points. Increase n_layers or window_radius.")
        return results
    
    # Step 2: Build Laplacian for quasicrystal
    print("\n[Step 2] Building graph Laplacian on quasicrystal...")
    L_qc = build_delaunay_graph_laplacian(Lambda_par)
    print(f"  Laplacian shape: {L_qc.shape}")
    print(f"  Average degree: {L_qc.diagonal().mean():.2f}")
    
    # Step 3: Compute spectrum
    print("\n[Step 3] Computing quasicrystal spectrum...")
    eigenvalues_qc = compute_laplacian_spectrum(L_qc)
    print(f"  Eigenvalue range: [{eigenvalues_qc.min():.4f}, {eigenvalues_qc.max():.4f}]")
    results['eigenvalues_qc'] = eigenvalues_qc
    
    # Step 4: Compute IDS
    print("\n[Step 4] Computing IDS and fitting Weyl law...")
    E_qc, N_qc = compute_ids(eigenvalues_qc)
    fit_qc = fit_weyl_law(E_qc, N_qc, d=4)
    print(f"  Fitted power: {fit_qc['fitted_power']:.3f} (expected: 2.0)")
    print(f"  sigma_qc = {fit_qc['sigma_fixed_d4']:.6e} +/- {fit_qc['sigma_error']:.2e}")
    results['sigma_qc'] = fit_qc['sigma_fixed_d4']
    results['sigma_qc_error'] = fit_qc['sigma_error']
    
    # Step 5: Reference periodic lattice
    print("\n[Step 5] Computing reference periodic lattice...")
    ref = compute_reference_sigma(n_qc)
    print(f"  sigma_periodic = {ref['sigma']:.6e} +/- {ref['sigma_error']:.2e}")
    results['sigma_periodic'] = ref['sigma']
    
    # Step 6: Compute ratio
    print("\n[Step 6] Computing suppression ratio...")
    sigma_ratio = fit_qc['sigma_fixed_d4'] / ref['sigma'] if ref['sigma'] > 0 else np.nan
    print(f"  sigma_ratio = sigma_qc / sigma_periodic = {sigma_ratio:.6e}")
    print(f"  phi^(-12) = {PHI_M12:.6e}")
    
    if not np.isnan(sigma_ratio) and sigma_ratio > 0:
        match = 100 * min(sigma_ratio, PHI_M12) / max(sigma_ratio, PHI_M12)
        print(f"  Match: {match:.1f}%")
        results['sigma_ratio'] = sigma_ratio
        results['match_percent'] = match
    else:
        results['sigma_ratio'] = np.nan
        results['match_percent'] = 0
    
    # Verdict
    print("\n" + "="*70)
    print("VERDICT")
    print("="*70)
    
    if results.get('match_percent', 0) > 90:
        print(f"\n  [OK] IDS prefactor ratio = {sigma_ratio:.4e}")
        print(f"       phi^(-12) = {PHI_M12:.4e}")
        print(f"       Match: {results['match_percent']:.1f}%")
        print("\n  CONCLUSION: phi^(-12) suppression CONFIRMED by IDS analysis.")
        results['PROVEN'] = True
    elif results.get('match_percent', 0) > 50:
        print(f"\n  [PARTIAL] IDS prefactor ratio = {sigma_ratio:.4e}")
        print(f"            phi^(-12) = {PHI_M12:.4e}")
        print(f"            Match: {results['match_percent']:.1f}%")
        print("\n  CONCLUSION: Results suggestive but not conclusive.")
        print("              Increase patch size for better statistics.")
        results['PROVEN'] = False
    else:
        print(f"\n  [FAIL] IDS prefactor ratio = {sigma_ratio:.4e}")
        print(f"         phi^(-12) = {PHI_M12:.4e}")
        print(f"         Match: {results['match_percent']:.1f}%")
        print("\n  CONCLUSION: phi^(-12) suppression NOT confirmed.")
        results['PROVEN'] = False
    
    print("="*70)
    
    return results


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # Run with moderate patch size (balance speed vs. accuracy)
    results = run_ids_spectral_proof(n_layers=2, window_radius=2.5)
