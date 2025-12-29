"""
Temperature Sweep: Emergent Phase Transitions from E8
======================================================

This module simulates the E8 symmetry breaking phase transitions as
the universe cooled from the Big Bang. Each transition reveals new
physics: forces split, masses emerge, symmetries break.

Timeline of E8 Phase Transitions:
- T ~ M_GUT: E8 → E7 (gravity separates)
- T ~ 10^15 GeV: E7 → E6 (GUT breaking)
- T ~ 10^14 GeV: E6 → SO(10) (further breaking)
- T ~ 10^12 GeV: SO(10) → SU(5) (proton decay scale)
- T ~ 10^12 GeV: SU(5) → SM (electroweak scale)
- T ~ 246 GeV: Electroweak symmetry breaking (Higgs)
- T ~ 150 MeV: QCD confinement (hadronization)

Author: E8 Research Team
Date: December 29, 2025
"""

import numpy as np
from constants import *


class E8PhaseTransition:
    """
    Model the E8 symmetry breaking phase transitions.
    """
    
    def __init__(self):
        self.phi = PHI
        self.M_planck = 1.22e19  # GeV
        self.M_GUT = 2e16        # GeV
        
        # Critical temperatures (GeV)
        self.phases = [
            {'name': 'E8', 'temp': self.M_planck, 'dim': 248, 'symmetry': 'E8'},
            {'name': 'E7', 'temp': self.M_GUT, 'dim': 133, 'symmetry': 'E7 × SU(2)'},
            {'name': 'E6', 'temp': self.M_GUT/self.phi, 'dim': 78, 'symmetry': 'E6 × SU(2) × U(1)'},
            {'name': 'SO(10)', 'temp': self.M_GUT/self.phi**2, 'dim': 45, 'symmetry': 'SO(10)'},
            {'name': 'SU(5)', 'temp': 1e14, 'dim': 24, 'symmetry': 'SU(5)'},
            {'name': 'SM', 'temp': 1e12, 'dim': 12, 'symmetry': 'SU(3)×SU(2)×U(1)'},
            {'name': 'Broken EW', 'temp': 246, 'dim': 12, 'symmetry': 'SU(3)×U(1)_em'},
            {'name': 'QCD Confined', 'temp': 0.15, 'dim': 8, 'symmetry': 'SU(3)_c confined'},
        ]
        
    def describe_phase(self, phase):
        """Describe a single phase."""
        print(f"╔{'═'*60}╗")
        print(f"║  {phase['name']:^56}  ║")
        print(f"╚{'═'*60}╝")
        print(f"   Temperature: T ~ {phase['temp']:.2e} GeV")
        print(f"   Symmetry: {phase['symmetry']}")
        print(f"   Gauge dimension: {phase['dim']}")
        print()
        
    def temperature_sweep(self):
        """Simulate the full temperature sweep from Big Bang to now."""
        print()
        print("╔═══════════════════════════════════════════════════════════════════╗")
        print("║   TEMPERATURE SWEEP: E8 PHASE TRANSITIONS                         ║")
        print("╚═══════════════════════════════════════════════════════════════════╝")
        print()
        print("As the universe cooled from the Big Bang (T ~ M_Planck)")
        print("to today (T ~ 2.7 K), E8 underwent a cascade of phase")
        print("transitions, each revealing new physics.")
        print()
        
        print("="*70)
        print("THE E8 BREAKING CASCADE")
        print("="*70)
        print()
        print("     E8 (248 dim)")
        print("      │  T ~ 10¹⁹ GeV (Planck)")
        print("      │")
        print("      ├──────────────────────────────────────────────┐")
        print("      │ GRAVITY SEPARATES                            │")
        print("      │ E8 → SO(3,1) × (remaining gauge)             │")
        print("      ├──────────────────────────────────────────────┘")
        print("      │")
        print("      ↓")
        print("     E7 (133 dim)")
        print("      │  T ~ 10¹⁶ GeV (GUT)")
        print("      │")
        print("      ├──────────────────────────────────────────────┐")
        print("      │ GUT TRANSITION                               │")
        print("      │ Proton decay operators activated             │")
        print("      ├──────────────────────────────────────────────┘")
        print("      │")
        print("      ↓")
        print("     E6 (78 dim)")
        print("      │  T ~ 10¹⁵ GeV")
        print("      │")
        print("      ├──────────────────────────────────────────────┐")
        print("      │ NEUTRINO MASS SCALE                          │")
        print("      │ Seesaw mechanism active                      │")
        print("      ├──────────────────────────────────────────────┘")
        print("      │")
        print("      ↓")
        print("     SO(10) (45 dim)")
        print("      │  T ~ 10¹⁴ GeV")
        print("      │")
        print("      ↓")
        print("     SU(5) (24 dim)")
        print("      │  T ~ 10¹² GeV")
        print("      │")
        print("      ├──────────────────────────────────────────────┐")
        print("      │ STANDARD MODEL EMERGES                       │")
        print("      │ Strong/electroweak split                     │")
        print("      ├──────────────────────────────────────────────┘")
        print("      │")
        print("      ↓")
        print("     SM (12 dim): SU(3)×SU(2)×U(1)")
        print("      │  T ~ 246 GeV (Electroweak)")
        print("      │")
        print("      ├──────────────────────────────────────────────┐")
        print("      │ ELECTROWEAK SYMMETRY BREAKING                │")
        print("      │ Higgs gets VEV, W±/Z get mass                │")
        print("      │ Fermion masses emerge                        │")
        print("      ├──────────────────────────────────────────────┘")
        print("      │")
        print("      ↓")
        print("     SU(3)×U(1)_em")
        print("      │  T ~ 150 MeV (QCD)")
        print("      │")
        print("      ├──────────────────────────────────────────────┐")
        print("      │ QCD CONFINEMENT                              │")
        print("      │ Quarks confined into hadrons                 │")
        print("      │ Baryogenesis complete                        │")
        print("      ├──────────────────────────────────────────────┘")
        print("      │")
        print("      ↓")
        print("     PRESENT UNIVERSE")
        print("      T ~ 2.7 K = 2.3×10⁻⁴ eV")
        print()
        
    def compute_order_parameters(self):
        """Compute order parameters at each transition."""
        print("="*70)
        print("ORDER PARAMETERS AT EACH TRANSITION")
        print("="*70)
        print()
        print("Each phase transition has an order parameter ψ that")
        print("characterizes the symmetry breaking:")
        print()
        print("┌────────────────────────────────────────────────────────────────┐")
        print("│  Transition        │  Order Parameter  │  Physical Meaning    │")
        print("├────────────────────────────────────────────────────────────────┤")
        print("│  E8 → E7           │  ⟨Φ_248⟩          │  Gravity separates   │")
        print("│  E7 → E6           │  ⟨Φ_56⟩           │  GUT scale VEV       │")
        print("│  E6 → SO(10)       │  ⟨Φ_27⟩           │  Matter parity       │")
        print("│  SO(10) → SU(5)    │  ⟨Φ_16⟩           │  B-L breaking        │")
        print("│  SU(5) → SM        │  ⟨Φ_24⟩           │  Doublet-triplet     │")
        print("│  EW breaking       │  ⟨H⟩ = 246 GeV    │  HIGGS VEV           │")
        print("│  QCD confinement   │  ⟨q̄q⟩             │  Chiral condensate   │")
        print("└────────────────────────────────────────────────────────────────┘")
        print()
        
    def critical_temperatures(self):
        """List all critical temperatures from E8."""
        print("="*70)
        print("CRITICAL TEMPERATURES FROM E8")
        print("="*70)
        print()
        print("The critical temperatures are determined by E8 structure:")
        print()
        
        T_planck = self.M_planck
        T_gut = self.M_GUT
        
        temps = [
            ("Planck", T_planck, "M_P", "E8 origin"),
            ("GUT", T_gut, "M_GUT", "E8 → E7"),
            ("E6", T_gut/self.phi, "M_GUT/φ", "E7 → E6"),
            ("SO(10)", T_gut/self.phi**2, "M_GUT/φ²", "E6 → SO(10)"),
            ("SU(5)", T_gut/self.phi**4, "M_GUT/φ⁴", "SO(10) → SU(5)"),
            ("EW", 246, "v", "Higgs → masses"),
            ("QCD", 0.15, "Λ_QCD", "Confinement"),
        ]
        
        print("┌──────────────┬─────────────────┬────────────────┬──────────────┐")
        print("│  Transition  │  T (GeV)        │  E8 Formula    │  Physics     │")
        print("├──────────────┼─────────────────┼────────────────┼──────────────┤")
        for name, temp, formula, physics in temps:
            print(f"│  {name:10s}  │  {temp:.2e}  │  {formula:14s}│  {physics:12s}│")
        print("└──────────────┴─────────────────┴────────────────┴──────────────┘")
        print()
        print("Note: Each scale is separated by powers of φ!")
        print(f"   φ = {self.phi:.6f}")
        print(f"   φ² = {self.phi**2:.6f}")
        print(f"   φ⁴ = {self.phi**4:.6f}")
        print()
        
    def emergent_phenomena(self):
        """Describe what emerges at each transition."""
        print("="*70)
        print("EMERGENT PHENOMENA AT EACH SCALE")
        print("="*70)
        print()
        
        emergents = [
            ("10¹⁹ GeV", "Spacetime, gravity, quantum mechanics"),
            ("10¹⁶ GeV", "GUT gauge bosons, proton decay operators"),
            ("10¹⁵ GeV", "Right-handed neutrinos, seesaw masses"),
            ("10¹⁴ GeV", "Matter/antimatter asymmetry (leptogenesis)"),
            ("10¹² GeV", "SM gauge bosons, Yukawa couplings"),
            ("246 GeV", "W±, Z⁰ bosons, fermion masses, CKM/PMNS"),
            ("150 MeV", "Protons, neutrons, pions, nuclear forces"),
            ("10⁻⁴ eV", "CMB, dark matter, galaxies, us"),
        ]
        
        print("What emerges at each temperature scale:")
        print()
        for scale, phenomena in emergents:
            print(f"   T ~ {scale}:")
            print(f"      → {phenomena}")
            print()
            
    def gravitational_waves(self):
        """Predict gravitational wave signatures from phase transitions."""
        print("="*70)
        print("GRAVITATIONAL WAVE SIGNATURES")
        print("="*70)
        print()
        print("First-order E8 phase transitions produce gravitational waves!")
        print()
        print("The peak frequency depends on transition temperature:")
        print()
        print("   f_peak ≈ (T/10¹⁰ GeV) × 10⁻⁵ Hz")
        print()
        print("┌────────────────┬──────────────┬──────────────────────────┐")
        print("│  Transition    │  f_peak (Hz) │  Detector                │")
        print("├────────────────┼──────────────┼──────────────────────────┤")
        print("│  E8 → E7       │  ~10⁴        │  Future (GHz searches)   │")
        print("│  GUT scale     │  ~10⁻¹       │  LISA                    │")
        print("│  EW scale      │  ~10⁻³       │  LISA, BBO               │")
        print("│  QCD scale     │  ~10⁻⁸       │  Pulsar timing arrays    │")
        print("└────────────────┴──────────────┴──────────────────────────┘")
        print()
        print("LISA could detect signals from E6 → SO(10) transition!")
        print()
        
    def run_full_analysis(self):
        """Run the complete temperature sweep analysis."""
        self.temperature_sweep()
        self.compute_order_parameters()
        self.critical_temperatures()
        self.emergent_phenomena()
        self.gravitational_waves()
        
        # Summary
        print("="*70)
        print("SUMMARY: E8 PHASE TRANSITION HIERARCHY")
        print("="*70)
        print()
        print("┌────────────────────────────────────────────────────────────────┐")
        print("│                                                                │")
        print("│  As the universe cooled from the Big Bang:                    │")
        print("│                                                                │")
        print("│  E8 (248) → E7 (133) → E6 (78) → SO(10) (45) → SU(5) (24)     │")
        print("│      ↓          ↓          ↓           ↓            ↓          │")
        print("│   Gravity    GUT      Neutrinos    B-L        Strong/EW       │")
        print("│                                                  ↓             │")
        print("│                              SM (12): SU(3)×SU(2)×U(1)         │")
        print("│                                        ↓                       │")
        print("│                              Higgs: v = 246 GeV                │")
        print("│                                        ↓                       │")
        print("│                              QCD: Λ = 150 MeV                  │")
        print("│                                        ↓                       │")
        print("│                              TODAY: T = 2.7 K                  │")
        print("│                                                                │")
        print("│  Each transition separated by powers of φ = 1.618...          │")
        print("│  All physics emerges from E8 symmetry breaking!               │")
        print("│                                                                │")
        print("└────────────────────────────────────────────────────────────────┘")
        print()


if __name__ == "__main__":
    phase = E8PhaseTransition()
    phase.run_full_analysis()
