#!/usr/bin/env python3
"""
GSM Persistent Discovery System (GSM-PDS)
==========================================

A self-learning geometric theorem discovery engine with long-term memory.

Features:
1. PERSISTENCE: Saves discoveries to gsm_library.json
2. MEMORY: Loads previous knowledge on restart  
3. DEDUPLICATION: Checks mathematical equivalence to avoid duplicates
4. COMPLEXITY SCORING: Prefers elegant, simple formulas
5. RESUME: Stop anytime, restart and continue building knowledge
6. FEEDBACK LOOP: Uses discovered integers as new atoms

This is the start of the GSM Automated Research Lab.

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
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, asdict
import hashlib

# ==============================================================================
# CONFIGURATION
# ==============================================================================

LIBRARY_FILE = "gsm_library.json"
SAVE_INTERVAL = 10  # Save to disk after every N new discoveries
MAX_COMPLEXITY = 50  # Prevent overly massive equations
NUMERIC_TOLERANCE = 1e-9  # For deduplication
TARGET_TOLERANCE = 1e-4  # For physical constant hunting

# Physical constant targets
TARGETS = {
    'alpha_inverse': 137.035999177,    # Fine structure constant inverse
    'proton_electron': 1836.15267343,  # m_p/m_e
    'weinberg_angle': 0.23121,         # sin²θ_W
    'riemann_zero_1': 14.134725142,    # First Riemann zeta zero
    'muon_electron': 206.7682830,      # m_μ/m_e
}


# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass
class Theorem:
    """Represents a discovered mathematical relationship."""
    expression: str
    simplified: str
    value: float
    tags: List[str]
    timestamp: str
    complexity: int
    discovery_method: str
    hash_id: str


class GSMLibrary:
    """
    The Long-Term Memory of the Discovery Engine.
    Manages loading, saving, and checking for duplicates.
    """
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.knowledge_base: List[Dict] = []
        self.known_values: Dict[float, Dict] = {}
        self.known_hashes: set = set()
        self.load_library()
    
    def load_library(self):
        """Load existing library from disk."""
        if not os.path.exists(self.filepath):
            print(f"[Library] No existing library found at {self.filepath}. Creating new.")
            self.knowledge_base = []
            return
        
        try:
            with open(self.filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self.knowledge_base = data.get('theorems', [])
                
                # Build lookup indices
                for entry in self.knowledge_base:
                    self.known_values[entry['value']] = entry
                    if 'hash_id' in entry:
                        self.known_hashes.add(entry['hash_id'])
                
                print(f"[Library] Loaded {len(self.knowledge_base)} theorems from disk.")
                
                # Print summary
                tags_count = {}
                for entry in self.knowledge_base:
                    for tag in entry.get('tags', []):
                        tags_count[tag] = tags_count.get(tag, 0) + 1
                
                print(f"[Library] Categories: {tags_count}")
                
        except Exception as e:
            print(f"[Library] Error loading library: {e}. Starting fresh.")
            self.knowledge_base = []
    
    def save_library(self):
        """Save library to disk."""
        try:
            data = {
                'metadata': {
                    'last_updated': str(datetime.now()),
                    'total_theorems': len(self.knowledge_base),
                    'version': '1.0'
                },
                'theorems': self.knowledge_base
            }
            with open(self.filepath, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"[Library] Saved. Total Theorems: {len(self.knowledge_base)}")
        except Exception as e:
            print(f"[Library] Error saving: {e}")
    
    def is_known_value(self, value: float, tolerance: float = NUMERIC_TOLERANCE) -> bool:
        """Check if a numerical value is already in the library."""
        for known_val in self.known_values.keys():
            if abs(known_val - value) < tolerance:
                return True
        return False
    
    def is_known_expression(self, expr_hash: str) -> bool:
        """Check if an expression hash is already known."""
        return expr_hash in self.known_hashes
    
    def add_theorem(self, expression: str, simplified: str, value: float, 
                    tags: List[str], method: str = "random_search"):
        """Add a new discovery to the database."""
        
        # Create hash for deduplication
        hash_id = hashlib.md5(simplified.encode()).hexdigest()[:12]
        
        if hash_id in self.known_hashes:
            return False  # Already known
        
        entry = {
            'expression': expression,
            'simplified': simplified,
            'value': float(value),
            'tags': tags,
            'timestamp': str(datetime.now()),
            'complexity': len(simplified),
            'discovery_method': method,
            'hash_id': hash_id
        }
        
        self.knowledge_base.append(entry)
        self.known_values[value] = entry
        self.known_hashes.add(hash_id)
        
        return True
    
    def get_best_by_tag(self, tag: str, n: int = 5) -> List[Dict]:
        """Get the n best (lowest complexity) theorems with a given tag."""
        matching = [t for t in self.knowledge_base if tag in t.get('tags', [])]
        return sorted(matching, key=lambda x: x['complexity'])[:n]
    
    def search_near_value(self, target: float, tolerance: float = 0.1) -> List[Dict]:
        """Find theorems with values near a target."""
        return [t for t in self.knowledge_base 
                if abs(t['value'] - target) < tolerance]


# ==============================================================================
# DISCOVERY ENGINE
# ==============================================================================

class RecursiveDiscoveryEngine:
    """
    The main discovery engine using GSM axioms and golden derivative mutation.
    """
    
    def __init__(self, library: GSMLibrary):
        self.lib = library
        
        # ==================================================================
        # GSM AXIOMS (The Geometric Toolkit)
        # ==================================================================
        
        # Golden ratio and conjugate
        self.phi = Rational(1, 2) + sqrt(5) / 2
        self.psi = Rational(1, 2) - sqrt(5) / 2
        
        # Verified Lattice Invariant (24-cell eigenvalue product)
        self.Lambda = 16 * sqrt(15)
        
        # Key constants
        self.sqrt5 = sqrt(5)
        self.sqrt15 = sqrt(15)
        
        # ==================================================================
        # ATOM POOL (Starting expressions for generation)
        # ==================================================================
        
        self.atoms = [
            # Golden ratio family
            self.phi,
            self.phi**2,
            self.phi**(-1),
            self.phi**(-5),
            self.phi**(-12),
            
            # Lattice
            self.Lambda,
            sqrt(self.Lambda),
            
            # Squares roots
            self.sqrt5,
            self.sqrt15,
            sqrt(Integer(2)),
            sqrt(Integer(3)),
            
            # Transcendentals
            pi,
            pi**2,
            pi**3,
            pi**5,
            E,
            
            # Key integers
            Integer(0),
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
            
            # Fibonacci/Lucas
            fibonacci(8),   # 21
            fibonacci(11),  # 89
            fibonacci(13),  # 233
            lucas(7),       # 29
            lucas(11),      # 199
            lucas(18),      # 5778
            
            # Fractions
            Rational(23, 60),
            Rational(1, 137),
        ]
        
        # ==================================================================
        # OPERATORS
        # ==================================================================
        
        self.binary_ops = [
            (lambda x, y: x + y, '+'),
            (lambda x, y: x - y, '-'),
            (lambda x, y: x * y, '*'),
            (lambda x, y: x / y, '/'),
        ]
        
        self.unary_ops = [
            (lambda x: x**2, '²'),
            (lambda x: sqrt(x), '√'),
            (lambda x: 1/x, '⁻¹'),
            (self.golden_derivative, 'D_φ'),
        ]
        
        # Statistics
        self.stats = {
            'total_generated': 0,
            'total_evaluated': 0,
            'integers_found': 0,
            'targets_found': 0,
            'duplicates_skipped': 0,
        }
    
    def golden_derivative(self, expr: sp.Expr) -> sp.Expr:
        """
        The GSM Mutation Operator (Golden Derivative).
        D_φ(f) = φ·f - φ⁻¹·f = f·(φ - φ⁻¹) = f (since φ - φ⁻¹ = 1)
        
        For more interesting transformation, use "Golden Friction":
        F_φ(f) = f · φ⁻¹² (applies 12th inverse golden power)
        """
        return expr * self.phi**(-12)
    
    def evaluate(self, expr: sp.Expr) -> Optional[float]:
        """Safely evaluate expression to float."""
        try:
            val = float(N(expr, 20))
            if not (val != val):  # Check for NaN
                return val
            return None
        except:
            return None
    
    def simplify_symbolic(self, expr: sp.Expr) -> sp.Expr:
        """Try to simplify the expression symbolically."""
        try:
            return simplify(expr)
        except:
            return expr
    
    def generate_candidate(self, depth: int = 3) -> Tuple[sp.Expr, str]:
        """Generate a random geometric expression."""
        
        # Start with a random atom
        curr_expr = random.choice(self.atoms)
        curr_str = str(curr_expr)
        
        for _ in range(random.randint(1, depth)):
            if random.random() < 0.3:  # 30% chance of unary op
                op_func, op_sym = random.choice(self.unary_ops)
                try:
                    # Check for valid sqrt
                    if op_sym == '√':
                        val = self.evaluate(curr_expr)
                        if val is None or val < 0:
                            continue
                    
                    curr_expr = op_func(curr_expr)
                    curr_str = f"{op_sym}({curr_str})"
                except:
                    pass
            else:  # Binary op
                op_func, op_sym = random.choice(self.binary_ops)
                other = random.choice(self.atoms)
                try:
                    # Avoid division by zero
                    if op_sym == '/':
                        other_val = self.evaluate(other)
                        if other_val is None or abs(other_val) < 1e-15:
                            continue
                    
                    curr_expr = op_func(curr_expr, other)
                    curr_str = f"({curr_str} {op_sym} {other})"
                except:
                    pass
        
        return curr_expr, curr_str
    
    def check_significance(self, value: float) -> Tuple[bool, List[str]]:
        """
        Check if a value is significant (integer or near a target).
        Returns (is_significant, list_of_tags).
        """
        tags = []
        
        # Check for integer (quantum number)
        if abs(value - round(value)) < NUMERIC_TOLERANCE and abs(value) > 0.001:
            tags.append('Integer')
            int_val = int(round(value))
            
            # Check if it's a Fibonacci number
            fibs = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987]
            if int_val in fibs:
                tags.append('Fibonacci')
            
            # Check if it's a Lucas number
            lucas_nums = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521]
            if int_val in lucas_nums:
                tags.append('Lucas')
            
            # Check if prime
            if int_val > 1:
                is_prime = all(int_val % i != 0 for i in range(2, int(int_val**0.5) + 1))
                if is_prime:
                    tags.append('Prime')
        
        # Check for physical constants
        for name, target in TARGETS.items():
            if abs(value - target) / target < TARGET_TOLERANCE:
                tags.append(f'Target:{name}')
        
        return len(tags) > 0, tags
    
    def run_cycle(self, duration_minutes: int = 60):
        """
        The Main Loop: Runs for a set time, discovering and saving.
        """
        start_time = time.time()
        new_discoveries = 0
        cycle_start = datetime.now()
        
        print("\n" + "=" * 70)
        print("GSM PERSISTENT DISCOVERY SYSTEM")
        print("=" * 70)
        print(f"Started at: {cycle_start}")
        print(f"Duration: {duration_minutes} minutes")
        print(f"Library contains: {len(self.lib.knowledge_base)} theorems")
        print("-" * 70)
        
        try:
            while (time.time() - start_time) < (duration_minutes * 60):
                # 1. GENERATE
                expr, expr_str = self.generate_candidate()
                self.stats['total_generated'] += 1
                
                # 2. EVALUATE
                val = self.evaluate(expr)
                if val is None or abs(val) > 1e10:
                    continue
                
                self.stats['total_evaluated'] += 1
                
                # 3. CHECK SIGNIFICANCE
                is_significant, tags = self.check_significance(val)
                
                if is_significant:
                    # 4. DEDUPLICATE
                    if self.lib.is_known_value(val):
                        self.stats['duplicates_skipped'] += 1
                        continue
                    
                    # 5. SIMPLIFY
                    simple_expr = self.simplify_symbolic(expr)
                    simple_str = str(simple_expr)
                    
                    # Check complexity
                    if len(simple_str) > MAX_COMPLEXITY:
                        continue
                    
                    # 6. SAVE
                    added = self.lib.add_theorem(
                        expression=expr_str,
                        simplified=simple_str,
                        value=val,
                        tags=tags,
                        method='random_golden_search'
                    )
                    
                    if added:
                        new_discoveries += 1
                        
                        # Print discovery
                        tag_str = ', '.join(tags)
                        print(f"  [!] NEW ({tag_str}): {simple_str} ≈ {val:.10f}")
                        
                        # Update stats
                        if 'Integer' in tags:
                            self.stats['integers_found'] += 1
                        if any('Target:' in t for t in tags):
                            self.stats['targets_found'] += 1
                        
                        # FEEDBACK LOOP: Add small integers to atoms
                        if 'Integer' in tags and abs(val) < 20:
                            int_val = Integer(int(round(val)))
                            if int_val not in self.atoms:
                                self.atoms.append(int_val)
                                print(f"      [+] Added {int_val} to atom pool")
                
                # Periodic save
                if new_discoveries > 0 and new_discoveries % SAVE_INTERVAL == 0:
                    self.lib.save_library()
                
                # Progress report every 1000 evaluations
                if self.stats['total_evaluated'] % 1000 == 0:
                    elapsed = (time.time() - start_time) / 60
                    print(f"  [Progress] {elapsed:.1f} min | "
                          f"Generated: {self.stats['total_generated']} | "
                          f"Discoveries: {new_discoveries}")
        
        except KeyboardInterrupt:
            print("\n[User Interrupt] Saving and exiting...")
        
        # Final save
        self.lib.save_library()
        
        # Summary
        print("\n" + "=" * 70)
        print("CYCLE COMPLETE")
        print("=" * 70)
        print(f"Duration: {(time.time() - start_time) / 60:.1f} minutes")
        print(f"Total generated: {self.stats['total_generated']}")
        print(f"Total evaluated: {self.stats['total_evaluated']}")
        print(f"New discoveries: {new_discoveries}")
        print(f"  - Integers: {self.stats['integers_found']}")
        print(f"  - Targets: {self.stats['targets_found']}")
        print(f"Duplicates skipped: {self.stats['duplicates_skipped']}")
        print(f"Library size: {len(self.lib.knowledge_base)} theorems")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 70)
    print("GSM PERSISTENT DISCOVERY SYSTEM (GSM-PDS)")
    print("Automated Geometric Theorem Discovery with Long-Term Memory")
    print("=" * 70)
    
    # Initialize library
    library = GSMLibrary(LIBRARY_FILE)
    
    # Initialize engine
    engine = RecursiveDiscoveryEngine(library)
    
    # Show current best discoveries by category
    if library.knowledge_base:
        print("\n--- CURRENT BEST DISCOVERIES ---")
        
        for tag in ['Integer', 'Target:alpha_inverse', 'Target:proton_electron']:
            best = library.get_best_by_tag(tag, n=3)
            if best:
                print(f"\n{tag}:")
                for t in best:
                    print(f"  {t['simplified']} ≈ {t['value']:.10f}")
    
    # Run discovery cycle
    # Default: 5 minutes for testing. Set to 1440 for 24 hours.
    print("\n" + "-" * 70)
    print("Starting discovery cycle (Press Ctrl+C to stop and save)")
    
    engine.run_cycle(duration_minutes=5)


if __name__ == "__main__":
    main()
