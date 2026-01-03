#!/usr/bin/env python3
"""
GSM Self-Correcting Engine - Two-Stage Discovery System
========================================================

A self-correcting theorem discovery engine that:
1. EXPLORES (Macro-Search): Scans geometric landscape for rough structure (error < 2%)
2. POLISHES (Micro-Refinement): Automatically refines discoveries using golden perturbations

Logic Flow:
  A. Engine finds 199 - 16√15 (Error = 0.05)
  B. Calculates Residual R = target - rough_value
  C. Hunts for R using perturbation terms (φ⁻ᵏ, π⁻ᵏ, fractions)
  D. If solved: combines base + correction
  E. Saves the POLISHED theorem to library

Result: High-precision solutions, not just approximations.

Author: GSI Discovery Framework
Date: 2026-01-01
"""

import sympy as sp
from sympy import sqrt, pi, E, Rational, Integer, fibonacci, lucas, simplify, N
import json
import os
import random
import time
from datetime import datetime
from typing import List, Dict, Tuple, Optional
import hashlib

# ==============================================================================
# CONFIGURATION
# ==============================================================================

LIBRARY_FILE = "gsm_refined_library.json"
SAVE_INTERVAL = 5  # Save after every N polished discoveries
MAX_COMPLEXITY = 100  # Allowed larger for refined theorems
MACRO_TOLERANCE = 0.02  # Explorer triggers refinement if within 2%
MICRO_TOLERANCE = 1e-7  # Refiner considers it a "Lock" if within this error

# Physical constant targets (The "Bounty Board")
TARGETS = {
    'alpha_inverse': 137.035999177,
    'proton_electron': 1836.15267343,
    'muon_electron': 206.7682830,
    'riemann_zero_1': 14.134725142,
    'weinberg_angle': 0.23121,
    'strong_coupling': 0.1180,
}


# ==============================================================================
# REFINED LIBRARY
# ==============================================================================

class RefinedGSMLibrary:
    """
    Enhanced library that stores both rough and refined theorems.
    """
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.knowledge_base: List[Dict] = []
        self.known_hashes: set = set()
        self.load_library()
    
    def load_library(self):
        """Load existing library from disk."""
        if not os.path.exists(self.filepath):
            print(f"[Library] No existing library at {self.filepath}. Starting fresh.")
            return
        
        try:
            with open(self.filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.knowledge_base = data.get('theorems', [])
                self.known_hashes = {entry['hash_id'] for entry in self.knowledge_base if 'hash_id' in entry}
                print(f"[Library] Loaded {len(self.knowledge_base)} theorems.")
        except Exception as e:
            print(f"[Library] Error loading: {e}. Starting fresh.")
    
    def save_library(self):
        """Save library to disk."""
        try:
            data = {
                'metadata': {
                    'last_updated': str(datetime.now()),
                    'total_theorems': len(self.knowledge_base),
                    'version': '2.0_self_correcting'
                },
                'theorems': self.knowledge_base
            }
            with open(self.filepath, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"[Library] Saved {len(self.knowledge_base)} theorems.")
        except Exception as e:
            print(f"[Library] Error saving: {e}")
    
    def add_theorem(self, expr_str: str, value: float, tags: List[str], 
                    error: float = None, base_expr: str = None, correction: str = None):
        """Add a theorem with optional refinement metadata."""
        hash_id = hashlib.md5(expr_str.encode()).hexdigest()[:12]
        
        if hash_id in self.known_hashes:
            return False
        
        entry = {
            'expression': expr_str,
            'value': float(value),
            'error': float(error) if error else 0.0,
            'error_ppm': (error / value * 1e6) if error and value > 0 else 0.0,
            'tags': tags,
            'timestamp': str(datetime.now()),
            'hash_id': hash_id,
            'complexity': len(expr_str)
        }
        
        # Add refinement metadata if present
        if base_expr:
            entry['base_expression'] = base_expr
        if correction:
            entry['correction_term'] = correction
        
        self.knowledge_base.append(entry)
        self.known_hashes.add(hash_id)
        return True
    
    def get_best_for_target(self, target_name: str) -> Optional[Dict]:
        """Get the best (lowest error) theorem for a target."""
        matching = [t for t in self.knowledge_base 
                   if any(target_name in tag for tag in t.get('tags', []))]
        if matching:
            return min(matching, key=lambda x: abs(x.get('error', float('inf'))))
        return None


# ==============================================================================
# SELF-CORRECTING ENGINE
# ==============================================================================

class SelfCorrectingEngine:
    """
    Two-stage discovery engine:
    1. Explorer (Macro): Finds rough geometric structures
    2. Polisher (Micro): Refines using golden perturbations
    """
    
    def __init__(self, library: RefinedGSMLibrary):
        self.lib = library
        
        # ==================================================================
        # 1. MACRO ATOMS (For finding the "Shape" - base structure)
        # ==================================================================
        
        self.phi = Rational(1, 2) + sqrt(5) / 2
        self.psi = Rational(1, 2) - sqrt(5) / 2
        self.Lambda = 16 * sqrt(15)  # Lattice Invariant
        
        self.atoms = [
            # Golden ratio family
            self.phi,
            self.phi**2,
            self.phi**3,
            self.phi**(-1),
            self.phi**(-2),
            sqrt(5),
            
            # Lattice
            self.Lambda,
            sqrt(self.Lambda),
            sqrt(15),
            
            # Transcendentals
            pi,
            pi**2,
            pi**3,
            pi**5,
            E,
            
            # Key integers
            Integer(1),
            Integer(2),
            Integer(3),
            Integer(5),
            Integer(6),
            Integer(7),
            Integer(12),
            Integer(13),
            Integer(23),
            Integer(60),
            Integer(137),
            Integer(360),
            
            # Fibonacci/Lucas
            fibonacci(8),   # 21
            fibonacci(11),  # 89
            fibonacci(13),  # 233
            lucas(7),       # 29
            lucas(11),      # 199
            lucas(18),      # 5778
            
            # Fractions
            Rational(1, 2),
            Rational(1, 3),
            Rational(23, 60),
        ]
        
        # ==================================================================
        # 2. PERTURBATION TOOLKIT (For Polishing) - The Sniper Terms
        # ==================================================================
        
        self.perturbations = []
        
        # Golden powers (key for golden friction corrections)
        for k in range(1, 40):
            term = self.phi**(-k)
            name = f"phi^-{k}"
            self.perturbations.append((term, name))
        
        # Lucas perturbations
        for k in [3, 4, 5, 7, 11]:
            term = lucas(k) * self.phi**(-12)
            name = f"L_{k}*phi^-12"
            self.perturbations.append((term, name))
        
        # Pi powers
        for k in [3, 4, 5, 6]:
            self.perturbations.append((pi**(-k), f"pi^-{k}"))
        
        # Lattice perturbations
        self.perturbations.append((self.Lambda**(-1), "Lambda^-1"))
        self.perturbations.append((sqrt(15)**(-1), "sqrt15^-1"))
        
        # ==================================================================
        # 3. OPERATORS
        # ==================================================================
        
        self.binary_ops = [
            (lambda x, y: x + y, '+'),
            (lambda x, y: x - y, '-'),
            (lambda x, y: x * y, '*'),
            (lambda x, y: x / y, '/'),
        ]
        
        # Statistics
        self.stats = {
            'total_generated': 0,
            'rough_matches': 0,
            'refinements_attempted': 0,
            'refinements_successful': 0,
        }
    
    def evaluate(self, expr: sp.Expr) -> Optional[float]:
        """Safely evaluate expression to float."""
        try:
            val = float(N(expr, 25))  # High precision
            if not (val != val):  # Check for NaN
                return val
            return None
        except:
            return None
    
    def generate_macro(self, depth: int = 3) -> Tuple[sp.Expr, str]:
        """Generate random macro expression for base structure."""
        curr_expr = random.choice(self.atoms)
        curr_str = str(curr_expr)
        
        for _ in range(random.randint(1, depth)):
            op_func, op_sym = random.choice(self.binary_ops)
            other = random.choice(self.atoms)
            
            try:
                # Avoid division by zero
                if op_sym == '/':
                    other_val = self.evaluate(other)
                    if other_val is None or abs(other_val) < 1e-15:
                        continue
                
                curr_expr = simplify(op_func(curr_expr, other))
                curr_str = f"({curr_str} {op_sym} {other})"
            except:
                pass
        
        return curr_expr, curr_str
    
    def refine_discovery(self, base_expr: sp.Expr, base_val: float, 
                        target_val: float, target_name: str) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[str]]:
        """
        The Polisher: Hunts for the residual difference.
        
        Returns: (Refined_Expr_String, New_Value, New_Error, Correction_Description)
        """
        residual = target_val - base_val
        best_correction = None
        best_error = abs(residual)
        correction_str = ""
        correction_desc = ""
        
        self.stats['refinements_attempted'] += 1
        
        # Scan through perturbations
        for p_term, p_name in self.perturbations:
            p_val = self.evaluate(p_term)
            if p_val is None or p_val == 0:
                continue
            
            # Test simple coefficients: (n/d) * perturbation
            for num in range(1, 15):
                for den in range(1, 15):
                    coeff = num / den
                    coeff_rat = Rational(num, den)
                    
                    # Test ADDITION: base + coeff * p_term
                    guess_val = base_val + (coeff * p_val)
                    err = abs(target_val - guess_val)
                    if err < best_error:
                        best_error = err
                        best_correction = coeff_rat * p_term
                        if num == den:
                            correction_str = f" + {p_name}"
                        else:
                            correction_str = f" + ({num}/{den})*{p_name}"
                        correction_desc = f"+{coeff:.4f}*{p_name}"
                    
                    # Test SUBTRACTION: base - coeff * p_term
                    guess_val = base_val - (coeff * p_val)
                    err = abs(target_val - guess_val)
                    if err < best_error:
                        best_error = err
                        best_correction = -coeff_rat * p_term
                        if num == den:
                            correction_str = f" - {p_name}"
                        else:
                            correction_str = f" - ({num}/{den})*{p_name}"
                        correction_desc = f"-{coeff:.4f}*{p_name}"
        
        # Decision Gate: Accept refinement only if significant improvement (10x better)
        original_error = abs(target_val - base_val)
        
        if best_correction is not None and best_error < (original_error * 0.1):
            self.stats['refinements_successful'] += 1
            
            print(f"    >>> POLISHING: Found correction {correction_str}")
            print(f"    >>> Error: {original_error:.6f} -> {best_error:.9f} ({best_error/target_val*1e6:.4f} ppm)")
            
            refined_expr = simplify(base_expr + best_correction)
            refined_str = str(refined_expr)
            
            return refined_str, float(N(refined_expr, 20)), best_error, correction_desc
        
        return None, None, None, None
    
    def run(self, duration_minutes: int = 60):
        """Main discovery loop."""
        start_time = time.time()
        new_polished = 0
        new_rough = 0
        
        print("\n" + "=" * 70)
        print("GSM SELF-CORRECTING ENGINE")
        print("Two-Stage Discovery: Explore + Polish")
        print("=" * 70)
        print(f"Duration: {duration_minutes} minutes")
        print(f"Macro tolerance: {MACRO_TOLERANCE*100:.1f}%")
        print(f"Micro tolerance: {MICRO_TOLERANCE}")
        print(f"Library: {len(self.lib.knowledge_base)} existing theorems")
        print("-" * 70)
        
        print("\nTarget bounty board:")
        for name, val in TARGETS.items():
            print(f"  {name}: {val}")
        
        print("\n" + "-" * 70)
        print("Hunting...\n")
        
        try:
            while (time.time() - start_time) < (duration_minutes * 60):
                # Stage 1: EXPLORE
                expr, expr_str = self.generate_macro()
                self.stats['total_generated'] += 1
                
                val = self.evaluate(expr)
                if val is None or abs(val) > 1e15 or val < 1e-10:
                    continue
                
                # Check against all targets
                for target_name, target_val in TARGETS.items():
                    error = abs(val - target_val)
                    rel_error = error / target_val if target_val != 0 else float('inf')
                    
                    # If within macro tolerance, trigger refinement
                    if rel_error < MACRO_TOLERANCE:
                        self.stats['rough_matches'] += 1
                        
                        print(f"\n[!] ROUGH MATCH for {target_name}")
                        print(f"    Base: {expr_str}")
                        print(f"    Value: {val:.10f} (target: {target_val})")
                        print(f"    Relative error: {rel_error*100:.4f}%")
                        
                        # Stage 2: POLISH
                        refined_str, refined_val, refined_err, correction = self.refine_discovery(
                            expr, val, target_val, target_name
                        )
                        
                        if refined_str:
                            # Save POLISHED theorem
                            tag = f"Target:{target_name}_Refined"
                            added = self.lib.add_theorem(
                                refined_str, refined_val, [tag],
                                error=refined_err,
                                base_expr=expr_str,
                                correction=correction
                            )
                            if added:
                                print(f"  [+] SAVED REFINED: {refined_str}")
                                new_polished += 1
                        else:
                            # Save ROUGH theorem if decent (< 0.5% error)
                            if rel_error < 0.005:
                                tag = f"Target:{target_name}_Rough"
                                added = self.lib.add_theorem(
                                    expr_str, val, [tag], error=error
                                )
                                if added:
                                    print(f"  [~] Saved rough (refinement failed)")
                                    new_rough += 1
                        
                        # Periodic save
                        if (new_polished + new_rough) > 0 and (new_polished + new_rough) % SAVE_INTERVAL == 0:
                            self.lib.save_library()
                
                # Progress report
                if self.stats['total_generated'] % 5000 == 0:
                    elapsed = (time.time() - start_time) / 60
                    print(f"\n[Progress] {elapsed:.1f} min | "
                          f"Generated: {self.stats['total_generated']} | "
                          f"Rough: {self.stats['rough_matches']} | "
                          f"Polished: {new_polished}")
        
        except KeyboardInterrupt:
            print("\n[User Interrupt] Saving...")
        
        # Final save
        self.lib.save_library()
        
        # Summary
        print("\n" + "=" * 70)
        print("ENGINE COMPLETE")
        print("=" * 70)
        print(f"Duration: {(time.time() - start_time) / 60:.1f} minutes")
        print(f"Total generated: {self.stats['total_generated']}")
        print(f"Rough matches found: {self.stats['rough_matches']}")
        print(f"Refinements attempted: {self.stats['refinements_attempted']}")
        print(f"Refinements successful: {self.stats['refinements_successful']}")
        print(f"New polished theorems: {new_polished}")
        print(f"New rough theorems: {new_rough}")
        print(f"Library size: {len(self.lib.knowledge_base)}")
        
        # Show best discoveries
        print("\n--- BEST DISCOVERIES ---")
        for target_name in TARGETS:
            best = self.lib.get_best_for_target(target_name)
            if best:
                ppm = best.get('error_ppm', 0)
                print(f"\n{target_name}:")
                print(f"  Expression: {best['expression']}")
                print(f"  Value: {best['value']:.15f}")
                print(f"  Error: {ppm:.4f} ppm")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 70)
    print("GSM SELF-CORRECTING ENGINE v2.0")
    print("Two-Stage Discovery: Explore + Polish")
    print("=" * 70)
    
    # Initialize library
    library = RefinedGSMLibrary(LIBRARY_FILE)
    
    # Initialize engine
    engine = SelfCorrectingEngine(library)
    
    # Run for specified duration (default 10 minutes for testing)
    # Set to 1440 for full day run
    print("\nStarting self-correcting discovery cycle...")
    print("Press Ctrl+C to stop and save.\n")
    
    engine.run(duration_minutes=10)


if __name__ == "__main__":
    main()
