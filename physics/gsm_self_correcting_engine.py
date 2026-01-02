#!/usr/bin/env python3
"""
GSM Self-Correcting Engine
===========================

Two-Stage Discovery System with automatic refinement.

ARCHITECTURE:
- Stage 1 (Explorer): Macro-search for rough structure (errors < 2%)
- Stage 2 (Polisher): Micro-refinement using perturbation terms

LOGIC FLOW:
  A. Engine finds base ~ target (Error < 2%)
  B. Calculates residual R = target - base
  C. Hunts for R using perturbation terms (n/d) * phi^-k
  D. Combines: refined = base + correction
  E. Saves POLISHED theorem to library

This is a Self-Correcting System. Output is precision-tuned equations,
not rough approximations.

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
SAVE_INTERVAL = 5  # Save after every 5 polished discoveries
MAX_COMPLEXITY = 80  # Allowed to be slightly larger for refined theorems
MACRO_TOLERANCE = 0.02  # Explorer triggers refinement if within 2%
MICRO_TOLERANCE = 1e-7  # Refiner considers it a "Lock" if within this error

# Physical targets (The "Bounty Board")
TARGETS = {
    'alpha_inverse': 137.035999177,
    'proton_electron': 1836.15267343,
    'muon_electron': 206.7682830,
    'riemann_zero_1': 14.134725142,
    'weinberg_angle': 0.23121,
    'tau_electron': 3477.48,
}


# ==============================================================================
# LIBRARY (Persistent Memory)
# ==============================================================================

class RefinedGSMLibrary:
    """Long-term memory with refined theorem storage."""
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.knowledge_base: List[Dict] = []
        self.known_hashes: set = set()
        self.load_library()
    
    def load_library(self):
        """Load existing library from disk."""
        if not os.path.exists(self.filepath):
            print(f"[Library] Creating new library at {self.filepath}")
            self.knowledge_base = []
            return
        
        try:
            with open(self.filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.knowledge_base = data.get('theorems', [])
                self.known_hashes = {
                    entry['hash_id'] for entry in self.knowledge_base 
                    if 'hash_id' in entry
                }
                print(f"[Library] Loaded {len(self.knowledge_base)} theorems")
        except Exception as e:
            print(f"[Library] Error loading: {e}. Starting fresh.")
            self.knowledge_base = []
    
    def save_library(self):
        """Save library to disk."""
        data = {
            'metadata': {
                'last_updated': str(datetime.now()),
                'total_theorems': len(self.knowledge_base),
                'version': '2.0-self-correcting'
            },
            'theorems': self.knowledge_base
        }
        with open(self.filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"[Library] Saved {len(self.knowledge_base)} theorems.")
    
    def add_theorem(self, expr_str: str, value: float, tags: List[str], 
                    error: float = None):
        """Add a theorem with deduplication."""
        hash_id = hashlib.md5(expr_str.encode()).hexdigest()[:12]
        if hash_id in self.known_hashes:
            return False
        
        entry = {
            'expression': expr_str,
            'value': float(value),
            'error': float(error) if error else 0.0,
            'relative_error_ppm': (error / value * 1e6) if error and value else 0.0,
            'tags': tags,
            'timestamp': str(datetime.now()),
            'hash_id': hash_id
        }
        self.knowledge_base.append(entry)
        self.known_hashes.add(hash_id)
        return True
    
    def get_best_by_target(self, target_name: str, n: int = 5) -> List[Dict]:
        """Get the n best (lowest error) theorems for a target."""
        matching = [t for t in self.knowledge_base 
                    if any(target_name in tag for tag in t.get('tags', []))]
        return sorted(matching, key=lambda x: x.get('error', float('inf')))[:n]


# ==============================================================================
# SELF-CORRECTING ENGINE
# ==============================================================================

class SelfCorrectingEngine:
    """
    Two-stage discovery engine with automatic refinement.
    
    Stage 1: Explorer (macro search for ~2% matches)
    Stage 2: Polisher (micro refinement using phi^-k perturbations)
    """
    
    def __init__(self, library: RefinedGSMLibrary):
        self.lib = library
        
        # ==================================================================
        # MACRO ATOMS (For finding the "Shape")
        # ==================================================================
        
        self.phi = Rational(1, 2) + sqrt(5) / 2
        self.Lambda = 16 * sqrt(15)  # 24-cell Lattice Invariant
        
        self.atoms = [
            self.phi,
            self.phi**2,
            self.phi**(-1),
            self.Lambda,
            sqrt(self.Lambda),
            pi,
            pi**2,
            pi**3,
            pi**5,
            E,
            sqrt(5),
            sqrt(15),
            Integer(1),
            Integer(2),
            Integer(3),
            Integer(5),
            Integer(12),
            Integer(137),
            Integer(199),
            lucas(7),       # 29
            lucas(11),      # 199
            fibonacci(12),  # 144
            fibonacci(13),  # 233
        ]
        
        # ==================================================================
        # MICRO ATOMS (Perturbation Toolkit for Polishing)
        # ==================================================================
        
        self.perturbations = []
        # Golden powers down to phi^-35
        for k in range(1, 36):
            self.perturbations.append((self.phi**(-k), f"phi^-{k}"))
        
        # Add pi powers
        self.perturbations.append((pi**(-3), "pi^-3"))
        self.perturbations.append((pi**(-4), "pi^-4"))
        self.perturbations.append((pi**(-5), "pi^-5"))
        
        # Add sqrt(5) terms
        self.perturbations.append((sqrt(5)**(-3), "sqrt(5)^-3"))
        self.perturbations.append((sqrt(5)**(-5), "sqrt(5)^-5"))
        
        # Binary operators
        self.ops = [
            (lambda x, y: x + y, "+"),
            (lambda x, y: x - y, "-"),
            (lambda x, y: x * y, "*"),
            (lambda x, y: x / y, "/"),
        ]
        
        # Statistics
        self.stats = {
            'explored': 0,
            'rough_matches': 0,
            'refined': 0,
            'saved': 0,
        }
    
    def evaluate(self, expr: sp.Expr) -> Optional[float]:
        """Safely evaluate expression to float."""
        try:
            val = float(N(expr, 20))
            if not (val != val) and abs(val) < 1e15:  # Check for NaN and overflow
                return val
            return None
        except:
            return None
    
    def generate_macro(self) -> Tuple[sp.Expr, str]:
        """Generate random macro expression for base structure."""
        curr = random.choice(self.atoms)
        curr_str = str(curr)
        
        for _ in range(random.randint(1, 3)):
            op, sym = random.choice(self.ops)
            other = random.choice(self.atoms)
            try:
                # Avoid division by zero
                if sym == '/':
                    other_val = self.evaluate(other)
                    if other_val is None or abs(other_val) < 1e-10:
                        continue
                
                curr = op(curr, other)
                curr_str = f"({curr_str} {sym} {other})"
            except:
                pass
        
        return curr, curr_str
    
    def refine_discovery(self, base_expr: sp.Expr, base_val: float, 
                         target_val: float, target_name: str) -> Tuple[Optional[str], Optional[float], Optional[float]]:
        """
        The Polisher: Hunts for the residual difference.
        
        Attempts to represent residual as (n/d) * perturbation_term.
        
        Returns:
            (refined_expression_string, new_value, new_error) or (None, None, None)
        """
        residual = target_val - base_val
        best_correction = None
        best_error = abs(residual)
        correction_str = ""
        
        # Scan perturbations with simple rational coefficients
        for p_term, p_name in self.perturbations:
            p_val = float(N(p_term))
            if abs(p_val) < 1e-20:
                continue
            
            # Try simple rational coefficients (n/d) where n,d in [1,12]
            for num in range(1, 13):
                for den in range(1, 13):
                    coeff = num / den
                    
                    # Check Addition: base + (n/d) * p_term
                    guess_val = base_val + (coeff * p_val)
                    err = abs(target_val - guess_val)
                    if err < best_error:
                        best_error = err
                        best_correction = (Rational(num, den) * p_term)
                        correction_str = f" + ({num}/{den})*{p_name}"
                    
                    # Check Subtraction: base - (n/d) * p_term
                    guess_val = base_val - (coeff * p_val)
                    err = abs(target_val - guess_val)
                    if err < best_error:
                        best_error = err
                        best_correction = -(Rational(num, den) * p_term)
                        correction_str = f" - ({num}/{den})*{p_name}"
                    
                    # Also try plain integer coefficients
                    if den == 1:
                        # Try n * p_term
                        for sign in [1, -1]:
                            guess_val = base_val + sign * num * p_val
                            err = abs(target_val - guess_val)
                            if err < best_error:
                                best_error = err
                                if sign == 1:
                                    best_correction = num * p_term
                                    correction_str = f" + {num}*{p_name}"
                                else:
                                    best_correction = -num * p_term
                                    correction_str = f" - {num}*{p_name}"
        
        # Decision Gate: Only accept if error improves by at least 10x
        original_error = abs(target_val - base_val)
        
        if best_correction is not None and best_error < (original_error * 0.1):
            # Significant improvement found!
            print(f"    >>> POLISHING: Found correction{correction_str}")
            print(f"    >>> Error: {original_error:.6f} -> {best_error:.10f}")
            
            refined_expr = base_expr + best_correction
            refined_str = f"{base_expr}{correction_str}"
            
            # Try to simplify
            try:
                simplified = simplify(refined_expr)
                simplified_str = str(simplified)
                if len(simplified_str) < len(refined_str):
                    refined_str = simplified_str
            except:
                pass
            
            return refined_str, float(N(refined_expr)), best_error
        
        return None, None, None
    
    def run(self, minutes: int = 60):
        """
        Main discovery loop.
        
        - Explores random macro expressions
        - When a rough match is found, triggers refinement
        - Saves polished theorems to library
        """
        print("\n" + "=" * 70)
        print("GSM SELF-CORRECTING ENGINE")
        print("=" * 70)
        print(f"Runtime: {minutes} minutes")
        print(f"Macro tolerance: {MACRO_TOLERANCE * 100}%")
        print(f"Library: {len(self.lib.knowledge_base)} existing theorems")
        print("-" * 70)
        
        start_time = time.time()
        
        try:
            while (time.time() - start_time) < (minutes * 60):
                # 1. EXPLORE - Generate random macro expression
                expr, expr_str = self.generate_macro()
                val = self.evaluate(expr)
                if val is None:
                    continue
                
                self.stats['explored'] += 1
                
                # 2. CHECK TARGETS - See if we're close to any bounty
                for name, target in TARGETS.items():
                    if target == 0:
                        continue
                    
                    error = abs(val - target)
                    rel_error = error / abs(target)
                    
                    # If within macro tolerance, trigger refinement
                    if rel_error < MACRO_TOLERANCE:
                        self.stats['rough_matches'] += 1
                        print(f"\n[!] ROUGH MATCH for {name}:")
                        print(f"    Base: {expr_str}")
                        print(f"    Value: {val:.10f} (target: {target:.10f})")
                        print(f"    Error: {rel_error * 100:.4f}%")
                        
                        # 3. TRIGGER REFINEMENT
                        ref_str, ref_val, ref_err = self.refine_discovery(
                            expr, val, target, name
                        )
                        
                        if ref_str:
                            # Save the POLISHED version
                            self.stats['refined'] += 1
                            tag = f"Target:{name}:Refined"
                            added = self.lib.add_theorem(
                                ref_str, ref_val, [tag], error=ref_err
                            )
                            if added:
                                self.stats['saved'] += 1
                                print(f"  [+] SAVED REFINED: {ref_str}")
                                print(f"      Value: {ref_val:.12f}")
                                print(f"      Error: {ref_err:.2e} ({ref_err/target*1e6:.2f} ppm)")
                        else:
                            # Save base if it's still quite good (< 0.1%)
                            if rel_error < 0.001:
                                tag = f"Target:{name}:Rough"
                                added = self.lib.add_theorem(
                                    expr_str, val, [tag], error=error
                                )
                                if added:
                                    self.stats['saved'] += 1
                                    print(f"  [+] Saved rough (refinement failed)")
                        
                        # Periodic save
                        if self.stats['saved'] > 0 and self.stats['saved'] % SAVE_INTERVAL == 0:
                            self.lib.save_library()
                
                # Progress report every 10000 explorations
                if self.stats['explored'] % 10000 == 0:
                    elapsed = (time.time() - start_time) / 60
                    print(f"  [Progress] {elapsed:.1f} min | "
                          f"Explored: {self.stats['explored']} | "
                          f"Rough: {self.stats['rough_matches']} | "
                          f"Refined: {self.stats['refined']}")
        
        except KeyboardInterrupt:
            print("\n[User Interrupt] Saving and exiting...")
        
        # Final save
        self.lib.save_library()
        
        # Summary
        print("\n" + "=" * 70)
        print("SELF-CORRECTING ENGINE COMPLETE")
        print("=" * 70)
        print(f"Duration: {(time.time() - start_time) / 60:.1f} minutes")
        print(f"Explored: {self.stats['explored']}")
        print(f"Rough matches: {self.stats['rough_matches']}")
        print(f"Refined: {self.stats['refined']}")
        print(f"Saved: {self.stats['saved']}")
        print(f"Library size: {len(self.lib.knowledge_base)} theorems")
        
        # Show best results per target
        print("\n--- BEST RESULTS BY TARGET ---")
        for name in TARGETS.keys():
            best = self.lib.get_best_by_target(name, n=1)
            if best:
                t = best[0]
                print(f"\n{name}:")
                print(f"  {t['expression']}")
                print(f"  Value: {t['value']:.12f}")
                print(f"  Error: {t['error']:.2e} ({t.get('relative_error_ppm', 0):.2f} ppm)")


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
    
    # Show existing best if any
    if library.knowledge_base:
        print("\n--- EXISTING BEST DISCOVERIES ---")
        for name in list(TARGETS.keys())[:3]:
            best = library.get_best_by_target(name, n=1)
            if best:
                t = best[0]
                print(f"{name}: {t['expression'][:50]}... (err: {t['error']:.2e})")
    
    # Initialize engine
    engine = SelfCorrectingEngine(library)
    
    # Run discovery
    # Default: 5 minutes for testing. Set to 1440 for 24 hours.
    print("\n" + "-" * 70)
    print("Starting self-correcting discovery (Ctrl+C to stop and save)")
    
    engine.run(minutes=5)


if __name__ == "__main__":
    main()
