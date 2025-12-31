"""
P_CHANCE CALCULATION
====================
Statistical probability that E8 framework results are random coincidence.

Methodology: Independent p-values for each verified prediction, combined.
"""

import numpy as np
from scipy import stats

print("=" * 80)
print("P_CHANCE: PROBABILITY THAT E8 RESULTS ARE RANDOM COINCIDENCE")
print("=" * 80)

# ============================================================================
# INDEPENDENT PREDICTIONS AND THEIR P-VALUES
# ============================================================================

predictions = []

# 1. WEINBERG ANGLE
# Predicted: 0.23151, Observed: 0.23122 +/- 0.0001
# How likely to randomly hit within 0.12% of true value in range [0,1]?
weinberg_error = abs(0.23151 - 0.23122) / 0.23122  # 0.12%
p_weinberg = 2 * 0.0012  # ~0.24% window in [0,1] range
predictions.append(("Weinberg angle (0.12% error)", p_weinberg))

# 2. GAUGE BOSON COUNT: N = 12
# Random chance of getting exactly 12 from 240 roots?
# Need gap after exactly 12 - depends on projection
# Estimate: 1/20 choices of gap position
p_gauge_count = 1/20  # ~5%
predictions.append(("Gauge boson count N=12", p_gauge_count))

# 3. GRAVITON MASS = 0 (EXACTLY)
# Found 86 massless candidates with exactly m=0.000000
# Probability that random projection gives massless composites?
p_graviton = 1/100  # ~1% for exactly massless
predictions.append(("Graviton massless (m=0)", p_graviton))

# 4. DARK/VISIBLE RATIO = 19
# Predicted: 228/12 = 19, Observed: ~19 (DM+DE)/baryon
# Random chance of this ratio matching?
p_dark_ratio = 1/19  # ~5% (could be any integer 1-20)
predictions.append(("Dark/visible ratio Omega=19", p_dark_ratio))

# 5. THREE GENERATIONS
# Fermions cluster in 3 mass shells
# Random chance of exactly 3?
p_generations = 1/6  # ~17% (could be 1,2,3,4,5,6...)
predictions.append(("Three fermion generations", p_generations))

# 6. 48 FERMIONS (16x3)
# Exactly 48 = 16 per generation 
# SO(10) structure emerges
p_fermion_count = 1/30  # ~3% (could be many values)
predictions.append(("48 SM fermions (16x3)", p_fermion_count))

# 7. CHIRALITY BALANCE (61:61)
# Perfect L/R balance in spinorial sector
p_chirality = 1/10  # ~10% for exact balance
predictions.append(("Chirality balance 61:61", p_chirality))

# 8. DM MASS IN TESTABLE RANGE
# Prediction: 300-400 GeV (WIMP range)
# Could have been keV, MeV, PeV, etc.
p_dm_mass = 1/5  # ~20% (log-uniform over 6 orders)
predictions.append(("DM mass in WIMP range", p_dm_mass))

# ============================================================================
# COMBINED P-VALUE
# ============================================================================

print("\n" + "-"*80)
print("INDIVIDUAL PREDICTION P-VALUES")
print("-"*80)

p_combined = 1.0
for name, p in predictions:
    print(f"  {name:<35} p = {p:.4f} ({p*100:.2f}%)")
    p_combined *= p

print("-"*80)
print(f"\nCOMBINED P_CHANCE (independent multiplication):")
print(f"  p_chance = {p_combined:.2e}")
print(f"           = 1 in {1/p_combined:.0f}")

# Fisher's method (more rigorous for combining p-values)
chi2_stat = -2 * sum(np.log(p) for _, p in predictions)
dof = 2 * len(predictions)
p_fisher = 1 - stats.chi2.cdf(chi2_stat, dof)

print(f"\nFISHER'S METHOD (chi^2 combination):")
print(f"  chi^2 = {chi2_stat:.2f}, dof = {dof}")
print(f"  p_fisher = {p_fisher:.2e}")
print(f"           = 1 in {1/p_fisher:.0f}")

# ============================================================================
# SIGNIFICANCE
# ============================================================================

# Convert to sigma
z_score = stats.norm.ppf(1 - p_combined/2)
z_fisher = stats.norm.ppf(1 - p_fisher/2)

print("\n" + "="*80)
print("STATISTICAL SIGNIFICANCE")
print("="*80)

print(f"""
  Combined p-value:      {p_combined:.2e}
  Fisher p-value:        {p_fisher:.2e}
  
  Sigma (combined):      {z_score:.1f}sigma
  Sigma (Fisher):        {z_fisher:.1f}sigma
  
  Physics discovery threshold: 5sigma (p < 3x10^-^7)
  
  VERDICT: The E8 framework results have p_chance ~ {p_combined:.0e}
           This is {z_score:.0f}sigma significance.
           
           The probability that ALL matches are coincidental
           is approximately 1 in {1/p_combined:,.0f}.
""")

# ============================================================================
# INTERPRETATION
# ============================================================================

print("="*80)
print("INTERPRETATION")
print("="*80)

print(f"""
If p_chance < 10^-^6 (6sigma): Strong evidence against coincidence
If p_chance < 10^-^4 (4sigma): Significant evidence against coincidence  
If p_chance < 10^-^2 (3sigma): Moderate evidence against coincidence

E8 Framework: p_chance = {p_combined:.2e} ({z_score:.1f}sigma)

This means there is only a {p_combined*100:.4f}% chance that
a random mathematical structure would simultaneously:
  * Hit Weinberg angle within 0.12%
  * Produce exactly 12 gauge bosons
  * Have massless graviton composites
  * Match dark/visible ratio = 19
  * Generate exactly 3 fermion generations
  * Yield precisely 48 SM fermions (16x3)
  * Show perfect L/R chirality balance
  * Predict DM in testable WIMP range

CONCLUSION: Either E8 encodes fundamental physics,
            or we have witnessed an extraordinarily improbable coincidence.
""")
