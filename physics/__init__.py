"""
E8 Theory of Everything - Rigorous Physics Framework
=====================================================

This package implements a mathematically rigorous E8 grand unified theory.

Unlike the numerology in the 'core' package, this package contains:
- Actual Lie algebra structure with 248 generators, 240 roots
- Explicit Yang-Mills Lagrangian
- Dynamical symmetry breaking mechanism
- Renormalization group running (1-loop and 2-loop)
- Honest treatment of fitted vs derived parameters

WHAT THIS PACKAGE DOES:
1. e8_algebra.py    - Complete E8 Lie algebra with root system
2. yang_mills.py    - The E8 Yang-Mills Lagrangian (real physics foundation)
3. symmetry_breaking.py - How E8 breaks to Standard Model
4. rge.py           - Renormalization group equations (why alpha is not exactly 1/137)
5. fermion_sector.py - Fermion masses and the Yukawa challenge
6. predictions.py   - Testable predictions vs curve fits

HONEST ASSESSMENT:
This framework shows what a real E8 theory requires. It also honestly
acknowledges what remains unsolved:
- Chirality (E8 is vectorlike, SM is chiral)
- Three generations (no E8 explanation)
- Fermion mass hierarchy (Yukawas are fitted, not derived)

Author: Physics Framework Team
Date: December 2025
"""

__all__ = [
    'E8LieAlgebra',
    'E8YangMills', 
    'E8SymmetryBreaking',
    'RenormalizationGroup',
    'FermionSector',
    'E8Predictions'
]


def __getattr__(name):
    """Lazy import to avoid import errors if scipy not available."""
    if name == 'E8LieAlgebra':
        from .e8_algebra import E8LieAlgebra
        return E8LieAlgebra
    elif name == 'E8YangMills':
        from .yang_mills import E8YangMills
        return E8YangMills
    elif name == 'E8SymmetryBreaking':
        from .symmetry_breaking import E8SymmetryBreaking
        return E8SymmetryBreaking
    elif name == 'RenormalizationGroup':
        from .rge import RenormalizationGroup
        return RenormalizationGroup
    elif name == 'FermionSector':
        from .fermion_sector import FermionSector
        return FermionSector
    elif name == 'E8Predictions':
        from .predictions import E8Predictions
        return E8Predictions
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
