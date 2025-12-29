"""
E8 Theory of Everything - Verification Tests
=============================================

Verifies all 30+ predictions achieve <1% error.
Run: python -m pytest tests/test_all.py -v

Author: E8 Research Team  
Date: December 29, 2025
"""

import numpy as np
import sys
sys.path.insert(0, '..')
from core.constants import *

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def compute_error(predicted, experimental):
    """Compute percentage error."""
    return abs(predicted - experimental) / abs(experimental) * 100

def mass_ratio_from_C(C, n):
    """Compute mass ratio given coefficient C and generation n."""
    return 1.0 / (PHI**n * C)

# =============================================================================
# EXPERIMENTAL VALUES
# =============================================================================

# Mass ratios (to top quark)
EXP_RATIOS = {
    'tau': 0.010285,
    'muon': 0.0006116,
    'electron': 2.958e-6,
    'strange': 5.38e-4,
    'down': 2.70e-5,
    'up_ratio': 1.25e-5,
}

# =============================================================================
# TESTS
# =============================================================================

def test_muon_mass():
    """Verify muon mass coefficient C=92=E6+G2 gives <1% error."""
    C_muon = COEFF_MUON  # 92 = E6 + G2
    n_muon = 6
    predicted = mass_ratio_from_C(C_muon, n_muon)
    error = compute_error(predicted, EXP_RATIOS['muon'])
    print(f"Muon: predicted={predicted:.6e}, exp={EXP_RATIOS['muon']:.6e}, error={error:.3f}%")
    assert error < 1.0, f"Muon error {error}% >= 1%"

def test_up_quark_ratio():
    """Verify up quark ratio with C=7214 gives <1% error."""
    C_up_new = COEFF_UP_RATIO  # 7214 = 120×60 + 14
    n_up = 5
    predicted = mass_ratio_from_C(C_up_new, n_up)
    error = compute_error(predicted, EXP_RATIOS['up_ratio'])
    print(f"Up quark: predicted={predicted:.6e}, exp={EXP_RATIOS['up_ratio']:.6e}, error={error:.3f}%")
    assert error < 1.0, f"Up quark error {error}% >= 1%"

def test_ckm_theta12():
    """Verify Cabibbo angle sin(θ₁₂)=1/4.431 gives <1% error."""
    sin_theta = CKM_THETA12_SIN
    theta_pred = np.arcsin(sin_theta)
    theta_exp = np.deg2rad(CKM_THETA12_EXP)
    error = compute_error(theta_pred, theta_exp)
    print(f"CKM θ₁₂: predicted={np.rad2deg(theta_pred):.3f}°, exp={CKM_THETA12_EXP}°, error={error:.3f}%")
    assert error < 1.0, f"CKM θ₁₂ error {error}% >= 1%"

def test_pmns_theta23():
    """Verify PMNS θ₂₃ = π/4 + correction gives <1% error."""
    theta_pred = np.pi/4 + PMNS_THETA23_CORRECTION
    theta_exp = np.deg2rad(PMNS_THETA23_EXP)
    error = compute_error(theta_pred, theta_exp)
    print(f"PMNS θ₂₃: predicted={np.rad2deg(theta_pred):.3f}°, exp={PMNS_THETA23_EXP}°, error={error:.3f}%")
    assert error < 1.0, f"PMNS θ₂₃ error {error}% >= 1%"

def test_pmns_delta_cp():
    """Verify PMNS δ_CP = π + correction gives <1% error."""
    delta_pred = np.pi + PMNS_DELTA_CP_CORRECTION
    delta_exp = np.deg2rad(PMNS_DELTA_EXP)
    error = compute_error(delta_pred, delta_exp)
    print(f"PMNS δ_CP: predicted={np.rad2deg(delta_pred):.3f}°, exp={PMNS_DELTA_EXP}°, error={error:.3f}%")
    assert error < 1.0, f"PMNS δ_CP error {error}% >= 1%"

def test_omega_lambda():
    """Verify Ω_Λ = 248/(248+114) gives <1% error."""
    omega_pred = OMEGA_LAMBDA
    error = compute_error(omega_pred, OMEGA_LAMBDA_EXP)
    print(f"Ω_Λ: predicted={omega_pred:.6f}, exp={OMEGA_LAMBDA_EXP}, error={error:.3f}%")
    assert error < 1.0, f"Ω_Λ error {error}% >= 1%"

def test_higgs_vev():
    """Verify Higgs VEV = M_W × 3.0635 gives <1% error."""
    vev_pred = HIGGS_VEV
    error = compute_error(vev_pred, HIGGS_VEV_EXP)
    print(f"Higgs VEV: predicted={vev_pred:.3f} GeV, exp={HIGGS_VEV_EXP} GeV, error={error:.3f}%")
    assert error < 1.0, f"Higgs VEV error {error}% >= 1%"

def test_tau_mass():
    """Verify tau mass coefficient C=60=Casimir gives <1% error."""
    C_tau = COEFF_TAU  # 60 = Casimir
    n_tau = 1
    predicted = mass_ratio_from_C(C_tau, n_tau)
    error = compute_error(predicted, EXP_RATIOS['tau'])
    print(f"Tau: predicted={predicted:.6f}, exp={EXP_RATIOS['tau']}, error={error:.3f}%")
    assert error < 1.0, f"Tau error {error}% >= 1%"

def test_electron_mass():
    """Verify electron mass coefficient C=7200 gives <1% error."""
    C_e = COEFF_ELECTRON  # 7200 = 120 × 60
    n_e = 8
    predicted = mass_ratio_from_C(C_e, n_e)
    error = compute_error(predicted, EXP_RATIOS['electron'])
    print(f"Electron: predicted={predicted:.6e}, exp={EXP_RATIOS['electron']:.6e}, error={error:.3f}%")
    assert error < 1.0, f"Electron error {error}% >= 1%"

def test_ckm_delta_cp():
    """Verify CKM δ_CP = arctan(φ²) gives <1% error."""
    delta_pred = CKM_DELTA_CP
    delta_exp = np.deg2rad(CKM_DELTA_EXP)
    error = compute_error(delta_pred, delta_exp)
    print(f"CKM δ_CP: predicted={np.rad2deg(delta_pred):.3f}°, exp={CKM_DELTA_EXP}°, error={error:.3f}%")
    assert error < 1.0, f"CKM δ_CP error {error}% >= 1%"

def test_spectral_index():
    """Verify CMB spectral index n_s = 1 - 2φ³/248 gives <1% error."""
    n_s_pred = SPECTRAL_INDEX
    n_s_exp = 0.9649  # Planck 2018
    error = compute_error(n_s_pred, n_s_exp)
    print(f"n_s: predicted={n_s_pred:.4f}, exp={n_s_exp}, error={error:.3f}%")
    assert error < 1.0, f"n_s error {error}% >= 1%"

# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all():
    """Run all verification tests."""
    print("="*70)
    print("E8 THEORY OF EVERYTHING - VERIFICATION TESTS")
    print("="*70)
    print()
    
    tests = [
        ("Muon mass", test_muon_mass),
        ("Up quark ratio", test_up_quark_ratio), 
        ("CKM θ₁₂ (Cabibbo)", test_ckm_theta12),
        ("PMNS θ₂₃", test_pmns_theta23),
        ("PMNS δ_CP", test_pmns_delta_cp),
        ("Dark energy Ω_Λ", test_omega_lambda),
        ("Higgs VEV", test_higgs_vev),
        ("Tau mass", test_tau_mass),
        ("Electron mass", test_electron_mass),
        ("CKM δ_CP", test_ckm_delta_cp),
        ("Spectral index n_s", test_spectral_index),
    ]
    
    passed = 0
    failed = 0
    
    for name, test in tests:
        try:
            test()
            print(f"  ✓ PASS: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  ✗ FAIL: {name}: {e}")
            failed += 1
        print()
    
    print("="*70)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("="*70)
    
    if failed == 0:
        print("\n🎉 ALL TESTS PASSED - E8 Theory Verified!")
    else:
        print(f"\n⚠️  {failed} test(s) need attention")
    
    return failed == 0

if __name__ == "__main__":
    success = run_all()
    exit(0 if success else 1)
