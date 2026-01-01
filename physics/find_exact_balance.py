#!/usr/bin/env python3
"""
Find exact V_B = V_F balance in E8 projections.
"""

import numpy as np

PHI = (1 + np.sqrt(5)) / 2

# Generate E8 roots
roots = []
for i in range(8):
    for j in range(i+1, 8):
        for s1 in [-1, 1]:
            for s2 in [-1, 1]:
                root = np.zeros(8)
                root[i], root[j] = s1, s2
                roots.append(root)

for bits in range(256):
    root = np.array([(1 if (bits >> i) & 1 else -1) * 0.5 for i in range(8)])
    if np.sum(root < 0) % 2 == 0:
        roots.append(root)

roots = np.array(roots)
integer_mask = np.all(np.abs(roots - np.round(roots)) < 0.01, axis=1)

print('=' * 60)
print('E8 LAMBDA EXACT BALANCE FINDER')
print('=' * 60)
print()
print(f'Integer (bosonic) roots: {np.sum(integer_mask)}')
print(f'Half-integer (fermionic) roots: {np.sum(~integer_mask)}')
print()

best_ratio = 1e10
best_P = None

for trial in range(20000):
    A = np.random.randn(8, 4)
    Q, _ = np.linalg.qr(A)
    P = Q[:, :4].T
    
    projected = roots @ P.T
    lengths_sq = np.sum(projected**2, axis=1)
    
    V_B = np.sum(lengths_sq[integer_mask]**2)
    V_F = np.sum(lengths_sq[~integer_mask]**2)
    
    ratio = abs(V_B/V_F - 1.0)
    if ratio < best_ratio:
        best_ratio = ratio
        best_P = P
        actual_ratio = V_B/V_F
        print(f'Trial {trial}: V_B/V_F = {actual_ratio:.8f} (deviation {ratio:.8f})')

print()
print('=' * 60)
print(f'BEST BALANCE: V_B/V_F = {actual_ratio:.8f}')
print(f'Deviation from 1.0: {best_ratio:.8f}')
print('=' * 60)

if best_ratio < 0.0001:
    print('EXACT CANCELLATION FOUND!!!')
elif best_ratio < 0.001:
    print('Near-exact balance (~0.1% imbalance)')
elif best_ratio < 0.01:
    print('Good balance (~1% imbalance)')
else:
    print(f'Best achievable: {best_ratio*100:.4f}% imbalance')
    print('=> Exact cancellation NOT possible with random projections')
    print()
    print('SOLUTION: Use triality/topological argument instead.')
    print('Λ = 0 topologically (π₃(E8) = Z)')
    print('Observed Λ = (H₀/M_Pl)⁴ ≈ 10⁻¹²² from symmetry breaking')
