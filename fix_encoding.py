"""
Fix Unicode encoding issues in E8 modules for Windows compatibility.
Replaces Unicode box characters with ASCII equivalents.
"""

import os
import re

# Files to fix
FILES = [
    "explicit_calculations.py",
    "gauge_boson_assignment.py", 
    "fermion_mapping.py",
    "chirality_triality.py",
    "so10_decomposition.py",
    "dark_matter_candidates.py",
    "cosmology_predictions.py",
    "p_chance_calculation.py",
    "physics/neutrino_sector.py",
    "physics/ckm_matrix.py",
    "deep_simulation.py",
    "physics/e8_unified_engine.py",
    "run_unified_theory.py",
]

# Unicode to ASCII replacements
REPLACEMENTS = {
    '█': '#',
    '─': '-',
    '│': '|',
    '┌': '+',
    '┐': '+',
    '└': '+',
    '┘': '+',
    '├': '+',
    '┤': '+',
    '┬': '+',
    '┴': '+',
    '┼': '+',
    '═': '=',
    '║': '|',
    '╔': '+',
    '╗': '+',
    '╚': '+',
    '╝': '+',
    '▶': '>',
    '✓': '[OK]',
    '✗': '[X]',
    '⚠': '[!]',
    '•': '*',
    '→': '->',
    '←': '<-',
    '↓': 'v',
    '↑': '^',
    '×': 'x',
    '÷': '/',
    '±': '+/-',
    '≈': '~',
    '≠': '!=',
    '≤': '<=',
    '≥': '>=',
    '∑': 'SUM',
    '∫': 'INT',
    '∞': 'inf',
    'π': 'pi',
    'θ': 'theta',
    'λ': 'lambda',
    'μ': 'mu',
    'ν': 'nu',
    'σ': 'sigma',
    'τ': 'tau',
    'φ': 'phi',
    'ψ': 'psi',
    'Ω': 'Omega',
    'α': 'alpha',
    'β': 'beta',
    'γ': 'gamma',
    'δ': 'delta',
    'ε': 'epsilon',
    'ρ': 'rho',
    'χ': 'chi',
    'Ψ': 'Psi',
    'Σ': 'Sigma',
    'Λ': 'Lambda',
    '⁰': '^0',
    '¹': '^1',
    '²': '^2',
    '³': '^3',
    '⁴': '^4',
    '⁵': '^5',
    '⁶': '^6',
    '⁷': '^7',
    '⁸': '^8',
    '⁹': '^9',
    '⁺': '^+',
    '⁻': '^-',
    '₀': '_0',
    '₁': '_1',
    '₂': '_2',
    '₃': '_3',
    '₄': '_4',
    '₅': '_5',
    '₆': '_6',
    '₇': '_7',
    '₈': '_8',
    '₉': '_9',
    '̄': '',  # combining macron
    '–': '-',  # en dash
    '—': '--',  # em dash
    ''': "'",
    ''': "'",
    '"': '"',
    '"': '"',
    '…': '...',
    '∈': 'in',
    '∉': 'not in',
    '⊂': 'subset',
    '⊃': 'superset',
    '∪': 'union',
    '∩': 'intersect',
    'ℏ': 'hbar',
    '½': '1/2',
    '¼': '1/4',
    '¾': '3/4',
    '√': 'sqrt',
    '∂': 'd',
    '∏': 'product',
    '′': "'",
    '″': '"',
    '°': ' deg',
    '∝': ' prop ',
    'η': 'eta',
    'ⁿ': '^n',
    'ₙ': '_n',
    'Δ': 'Delta',
    '🏆': '[ACHIEVEMENT]',
    '🎯': '[TARGET]',
}

def fix_file(filepath):
    """Replace Unicode characters in a file."""
    print(f"Fixing: {filepath}")
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original = content
    
    for unicode_char, ascii_replacement in REPLACEMENTS.items():
        content = content.replace(unicode_char, ascii_replacement)
    
    if content != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"  Fixed {filepath}")
        return True
    else:
        print(f"  No changes needed in {filepath}")
        return False

def main():
    print("=" * 60)
    print("E8 ENCODING FIX - Windows Compatibility")
    print("=" * 60)
    
    fixed = 0
    for filename in FILES:
        if os.path.exists(filename):
            if fix_file(filename):
                fixed += 1
        else:
            print(f"  Skipped: {filename} not found")
    
    print("-" * 60)
    print(f"Fixed {fixed} files")
    print("Run 'python run_unified_theory.py --quick' to test")
    print("=" * 60)

if __name__ == "__main__":
    main()
