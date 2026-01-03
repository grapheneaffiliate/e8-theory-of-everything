"""
GSM SPECTRUM ANALYZER
======================
"Strip-Mining the Riemann Zeros for Geometric Keys"

Purpose: Extract the PRIMARY geometric coefficients from the first 50 Riemann Zeros
         to discover if they map to Lie Algebra dimensions or sporadic group properties.

Hypothesis: The integer coefficients in the geometric formulas (23, 5, 14, ...) 
            encode symmetry groups (E8, Moonshine, etc.)

Author: GSM Research Team
Date: January 2, 2026
"""

import mpmath
from mpmath import sqrt, pi, power, mpf
import json
from datetime import datetime

# Maximum precision
mpmath.mp.dps = 60

# ═══════════════════════════════════════════════════════════════════════════
# KNOWN FIRST 50 RIEMANN ZEROS (Imaginary parts, from Odlyzko tables)
# ═══════════════════════════════════════════════════════════════════════════
RIEMANN_ZEROS = [
    "14.134725141734693790457251983562",
    "21.022039638771554992628479593897",
    "25.010857580145688763213790992563",
    "30.424876125859513210311897530584",
    "32.935061587739189690662368964075",
    "37.586178158825671257217763480705",
    "40.918719012147495187398126914633",
    "43.327073280914999519496122165407",
    "48.005150881167159727942472749428",
    "49.773832477672302181916784678564",
    "52.970321477714460644147296608881",
    "56.446247697063394804367759476706",
    "59.347044002602353079653648674992",
    "60.831778524609809844259901824524",
    "65.112544048081606660875054253184",
    "67.079810529494173714478828896522",
    "69.546401711173979252926857526555",
    "72.067157674481907581883330013489",
    "75.704690699083933168326916762030",
    "77.144840068874805372682664856305",
    "79.337375020249367922763592877116",
    "82.910380854086030183164837494771",
    "84.735492980517050105735311206827",
    "87.425274613125229406531667850919",
    "88.809111207634465423682348079509",
    "92.491899257797855152396663099712",
    "94.651344040519886966597925815199",
    "95.870634228245309758741029219246",
    "98.831194218193692233324420138622",
    "101.31785100573139122878544794027",
    "103.72553804047833941639840810213",
    "105.44662305232609449367083241411",
    "107.16861118427640751512335196308",
    "111.02953554316967452465645030994",
    "111.87465917699263708561207871677",
    "114.32022091545271276589093727619",
    "116.22668032085755438216080431206",
    "118.79078286597621732297913970269",
    "121.37012500242064591894553297837",
    "122.94682929355258820081746033077",
    "124.25681855434576718473200319713",
    "127.51668387959649512427932376691",
    "129.57870419995605098576803390617",
    "131.08768853093265672356163384034",
    "133.49773720299758645013049135890",
    "134.75650975337387133132606415716",
    "138.11604205453344320019155519028",
    "139.73620895212138895045004652061",
    "141.12370740402112376194035381847",
    "143.11184580762063273940512386891",
]


class GSMSpectrumAnalyzer:
    """
    Analyzes Riemann Zeros to extract their dominant geometric structure.
    """
    
    def __init__(self):
        # GSM Constants
        self.phi = (1 + sqrt(5)) / 2
        self.psi = self.phi - 1
        self.pi = pi
        self.Lambda = 16 * sqrt(15)
        
        # Key Lie algebra dimensions and group orders
        self.LIE_DIMENSIONS = {
            1: "U(1)",
            3: "SU(2), SO(3)",
            5: "Pentagonal, Icosahedral A5",
            7: "G2 minimal",
            8: "SU(3)",
            10: "SO(5), Sp(4)",
            14: "G2",
            15: "SU(4)",
            21: "SO(7)",
            23: "Moonshine prime (Happy Family)",
            24: "Leech lattice dim",
            28: "SO(8)",
            35: "SU(6)",
            36: "SO(9)",
            45: "SO(10)",
            52: "F4",
            78: "E6",
            133: "E7",
            248: "E8",
        }
        
        # Fibonacci numbers (for coefficient analysis)
        self.FIBONACCI = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
        
        # Results storage
        self.spectrum = []
    
    def find_best_geometric_base(self, target_zero: mpf, zero_index: int) -> dict:
        """
        For a given Riemann zero, find the best geometric base expression.
        Returns the dominant coefficient structure.
        """
        target = mpf(target_zero)
        
        best_match = {
            "zero_index": zero_index,
            "zero_value": float(target),
            "best_error": float('inf'),
            "best_formula": "",
            "primary_coefficient": None,
            "denominator": None,
            "base_type": "",
            "geometric_structure": "",
        }
        
        # Test various base structures
        # Structure 1: (a/b) * π² * φ^n + (c/d) * π * φ^m - 1/2
        # Structure 2: k * π * φ^n - 1/2
        # Structure 3: Λ * (a/b) / φ^n - 1/2
        
        for num in range(1, 50):  # Numerator
            for den in range(1, 50):  # Denominator
                if num == den:
                    continue
                coeff = mpf(num) / mpf(den)
                
                # Test: coeff * π² * φ² - 1/2
                for phi_pow in range(-4, 6):
                    base = coeff * (self.pi ** 2) * power(self.phi, phi_pow) - mpf('0.5')
                    residual = abs(target - base)
                    
                    if residual < best_match["best_error"] and residual < target * 0.5:
                        best_match["best_error"] = float(residual)
                        best_match["primary_coefficient"] = num
                        best_match["denominator"] = den
                        best_match["base_type"] = f"π²·φ^{phi_pow}"
                        best_match["geometric_structure"] = f"({num}/{den})·π²·φ^{phi_pow} - 1/2"
                
                # Test: coeff * π * φ² - 1/2
                for phi_pow in range(-4, 6):
                    base = coeff * self.pi * power(self.phi, phi_pow) - mpf('0.5')
                    residual = abs(target - base)
                    
                    if residual < best_match["best_error"] and residual < target * 0.5:
                        best_match["best_error"] = float(residual)
                        best_match["primary_coefficient"] = num
                        best_match["denominator"] = den
                        best_match["base_type"] = f"π·φ^{phi_pow}"
                        best_match["geometric_structure"] = f"({num}/{den})·π·φ^{phi_pow} - 1/2"
        
        # Test integer multiples: k * π * φ² - 1/2
        for k in range(1, 50):
            for phi_pow in range(0, 5):
                base = k * self.pi * power(self.phi, phi_pow) - mpf('0.5')
                residual = abs(target - base)
                
                if residual < best_match["best_error"]:
                    best_match["best_error"] = float(residual)
                    best_match["primary_coefficient"] = k
                    best_match["denominator"] = 1
                    best_match["base_type"] = f"π·φ^{phi_pow}"
                    best_match["geometric_structure"] = f"{k}·π·φ^{phi_pow} - 1/2"
        
        # Test Lambda-based: Λ / φ^k - 1/2
        for phi_pow in range(-2, 8):
            base = self.Lambda / power(self.phi, phi_pow) - mpf('0.5')
            residual = abs(target - base)
            
            if residual < best_match["best_error"]:
                best_match["best_error"] = float(residual)
                best_match["primary_coefficient"] = 16  # From Λ = 16√15
                best_match["denominator"] = 1
                best_match["base_type"] = f"Λ/φ^{phi_pow}"
                best_match["geometric_structure"] = f"Λ/φ^{phi_pow} - 1/2"
        
        return best_match
    
    def analyze_all_zeros(self):
        """
        Run spectrum analysis on all 50 zeros.
        """
        print("═" * 80)
        print("         GSM SPECTRUM ANALYZER")
        print("         Extracting Geometric Keys from Riemann Zeros")
        print("═" * 80)
        print(f"\nAnalyzing {len(RIEMANN_ZEROS)} Riemann Zeros...")
        print("-" * 80)
        
        for i, zero_str in enumerate(RIEMANN_ZEROS):
            result = self.find_best_geometric_base(zero_str, i + 1)
            self.spectrum.append(result)
            
            # Print progress
            coeff = result["primary_coefficient"]
            lie_match = self.LIE_DIMENSIONS.get(coeff, "")
            fib_match = "FIB" if coeff in self.FIBONACCI else ""
            
            print(f"Zero #{i+1:2d}: γ = {result['zero_value']:.6f}  →  "
                  f"{result['geometric_structure']:30s}  "
                  f"[Key: {coeff:3d}] {lie_match} {fib_match}")
        
        return self.spectrum
    
    def extract_coefficient_histogram(self):
        """
        Build histogram of primary coefficients.
        """
        histogram = {}
        for result in self.spectrum:
            coeff = result["primary_coefficient"]
            if coeff not in histogram:
                histogram[coeff] = 0
            histogram[coeff] += 1
        return histogram
    
    def generate_report(self):
        """
        Generate comprehensive analysis report.
        """
        print("\n" + "═" * 80)
        print("         SPECTRUM ANALYSIS REPORT")
        print("═" * 80)
        
        # Coefficient histogram
        histogram = self.extract_coefficient_histogram()
        sorted_hist = sorted(histogram.items(), key=lambda x: -x[1])
        
        print("\n┌────────────────────────────────────────────────────────┐")
        print("│         PRIMARY COEFFICIENT FREQUENCY                  │")
        print("├────────────────────────────────────────────────────────┤")
        
        for coeff, count in sorted_hist[:15]:
            lie_info = self.LIE_DIMENSIONS.get(coeff, "Unknown")
            fib = "★FIB" if coeff in self.FIBONACCI else ""
            bar = "█" * count
            print(f"│ {coeff:3d}: {bar:15s} ({count:2d}) → {lie_info:20s} {fib:4s}│")
        
        print("└────────────────────────────────────────────────────────┘")
        
        # Check Lie algebra matches
        print("\n┌────────────────────────────────────────────────────────┐")
        print("│         LIE ALGEBRA / GROUP MATCHES                    │")
        print("├────────────────────────────────────────────────────────┤")
        
        lie_matches = 0
        for coeff in histogram.keys():
            if coeff in self.LIE_DIMENSIONS:
                lie_matches += histogram[coeff]
                print(f"│ {coeff:3d} → {self.LIE_DIMENSIONS[coeff]:40s}     │")
        
        print("└────────────────────────────────────────────────────────┘")
        
        total_zeros = len(self.spectrum)
        print(f"\n  LIE ALGEBRA MATCH RATE: {lie_matches}/{total_zeros} "
              f"= {100*lie_matches/total_zeros:.1f}%")
        
        # Fibonacci coefficient count
        fib_matches = sum(1 for r in self.spectrum 
                        if r["primary_coefficient"] in self.FIBONACCI)
        print(f"  FIBONACCI MATCH RATE:   {fib_matches}/{total_zeros} "
              f"= {100*fib_matches/total_zeros:.1f}%")
        
        # Base type distribution
        base_types = {}
        for r in self.spectrum:
            bt = r["base_type"]
            if bt not in base_types:
                base_types[bt] = 0
            base_types[bt] += 1
        
        print("\n┌────────────────────────────────────────────────────────┐")
        print("│         GEOMETRIC BASE TYPE DISTRIBUTION               │")
        print("├────────────────────────────────────────────────────────┤")
        
        for bt, count in sorted(base_types.items(), key=lambda x: -x[1])[:10]:
            bar = "█" * count
            print(f"│ {bt:15s}: {bar:20s} ({count:2d})                 │")
        
        print("└────────────────────────────────────────────────────────┘")
        
        # Save to JSON
        report = {
            "timestamp": str(datetime.now()),
            "total_zeros_analyzed": total_zeros,
            "coefficient_histogram": histogram,
            "lie_match_rate": lie_matches / total_zeros,
            "fibonacci_match_rate": fib_matches / total_zeros,
            "base_type_distribution": base_types,
            "full_spectrum": self.spectrum,
        }
        
        with open("riemann_spectrum_analysis.json", "w") as f:
            json.dump(report, f, indent=2)
        
        print("\n  Full report saved to: riemann_spectrum_analysis.json")
        
        return report


def main():
    print("""
    ╔════════════════════════════════════════════════════════════════════════╗
    ║                                                                        ║
    ║            G S M   S P E C T R U M   A N A L Y Z E R                  ║
    ║                                                                        ║
    ║   "Strip-Mining the Riemann Zeros for Lie Algebra Dimensions"         ║
    ║                                                                        ║
    ╚════════════════════════════════════════════════════════════════════════╝
    """)
    
    analyzer = GSMSpectrumAnalyzer()
    
    # Run full analysis
    analyzer.analyze_all_zeros()
    
    # Generate report
    analyzer.generate_report()
    
    print("\n" + "═" * 80)
    print("         ANALYSIS COMPLETE")
    print("═" * 80)
    print("\n  KEY QUESTION: Do the coefficients map to Lie Group dimensions?")
    print("  If YES → Riemann Zeros encode the classification of symmetry groups!")
    print("═" * 80)


if __name__ == "__main__":
    main()
