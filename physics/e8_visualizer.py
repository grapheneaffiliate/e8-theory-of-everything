"""
E8 UNIVERSE VISUALIZER
======================
Goal: Visual confirmation of the Standard Model Geometry.
Method: 
1. Load the Universe DNA (4x8 Matrix).
2. Project E8 to 4D -> 3D Shadow.
3. Color code:
   - RED:   Standard Model Bosons (N=12)
   - BLUE:  Generation 27 Fermions (N=16)
   - GREY:  The rest (Dark Sector)
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
from e8_constants import UNIVERSE_MATRIX

# ==========================================
# 1. SETUP
# ==========================================
class E8Crystal:
    """Generate the 240 roots of E8."""
    def __init__(self):
        roots = []
        for i in range(8):
            for j in range(i + 1, 8):
                for s1, s2 in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                    r = np.zeros(8)
                    r[i], r[j] = s1, s2
                    roots.append(r)
        for signs in product([0.5, -0.5], repeat=8):
            if sum(1 for s in signs if s < 0) % 2 == 0:
                roots.append(np.array(signs))
        self.roots = np.array(roots)


E8_ROOTS = E8Crystal().roots


def visualize():
    """
    Generate 3D visualization of the E8 Standard Model structure.
    """
    print("="*70)
    print("GENERATING 3D MODEL OF THE E8 UNIVERSE")
    print("="*70)
    print()
    
    # 1. Project to 4D
    shadows_4d = E8_ROOTS @ UNIVERSE_MATRIX.T
    lengths_sq = np.sum(shadows_4d**2, axis=1)
    lengths = np.sqrt(lengths_sq)
    
    # 2. Identify Groups
    sorted_idx = np.argsort(lengths_sq)
    
    # Bosons: The 12 Lightest
    boson_idx = sorted_idx[:12]
    
    # Fermions: Generation 27 (The "Perfect 16")
    # Find the cluster of 16 around mass ~1.022
    target_mass = 1.022
    fermion_mask = np.abs(lengths - target_mass) < 0.01
    fermion_idx = np.where(fermion_mask)[0]
    
    # If we didn't find exactly 16, expand the search
    if len(fermion_idx) != 16:
        print(f"Found {len(fermion_idx)} fermions at mass ~1.022")
        print("Searching for 16-fold cluster...")
        
        # More robust search: find mass where exactly 16 roots cluster
        for test_mass in np.arange(0.5, 1.5, 0.005):
            test_mask = np.abs(lengths - test_mass) < 0.01
            if np.sum(test_mask) == 16:
                fermion_idx = np.where(test_mask)[0]
                target_mass = test_mass
                print(f"Found 16-fold cluster at mass {target_mass:.3f}")
                break
    
    # Dark Sector: Everything else
    dark_idx = [i for i in range(240) if i not in boson_idx and i not in fermion_idx]
    
    print(f"Bosons (Red):    {len(boson_idx)}")
    print(f"Fermions (Blue): {len(fermion_idx)}")
    print(f"Dark (Grey):     {len(dark_idx)}")
    print()
    
    # 3. Project to 3D for plotting (PCA)
    # Get the best 3D view of the 4D slice
    u, s, vh = np.linalg.svd(shadows_4d.T)
    projection_3d = u[:, :3].T
    roots_3d = shadows_4d @ projection_3d.T
    
    # 4. Plot
    print("Creating 3D visualization...")
    fig = plt.figure(figsize=(14, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot Dark Sector (Faint)
    dark_vecs = roots_3d[dark_idx]
    ax.scatter(dark_vecs[:,0], dark_vecs[:,1], dark_vecs[:,2], 
               c='lightgrey', alpha=0.1, s=10, label=f'Dark Sector ({len(dark_idx)})')
    
    # Plot Fermions (Blue)
    ferm_vecs = roots_3d[fermion_idx]
    ax.scatter(ferm_vecs[:,0], ferm_vecs[:,1], ferm_vecs[:,2], 
               c='blue', alpha=0.8, s=80, label=f'Generation Matter ({len(fermion_idx)})')
    
    # Plot Bosons (Red Stars)
    boson_vecs = roots_3d[boson_idx]
    ax.scatter(boson_vecs[:,0], boson_vecs[:,1], boson_vecs[:,2], 
               c='red', alpha=1.0, s=150, marker='*', label=f'Gauge Bosons ({len(boson_idx)})')
    
    # Draw lines from origin to bosons
    for v in boson_vecs:
        ax.plot([0, v[0]], [0, v[1]], [0, v[2]], 'r-', alpha=0.3, linewidth=0.5)
    
    # Draw lines from origin to fermions
    for v in ferm_vecs:
        ax.plot([0, v[0]], [0, v[1]], [0, v[2]], 'b-', alpha=0.2, linewidth=0.5)
    
    # Mark origin
    ax.scatter([0], [0], [0], c='gold', s=200, marker='o', label='Origin (Vacuum)')
    
    ax.set_xlabel('Principal Axis 1')
    ax.set_ylabel('Principal Axis 2')
    ax.set_zlabel('Principal Axis 3')
    ax.set_title("E8 Standard Model Geometry\nRed★=Forces(12), Blue●=Matter(16), Grey○=Dark(212)")
    ax.legend(loc='upper left')
    
    # Make axes equal
    max_range = np.max(np.abs(roots_3d)) * 0.8
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])
    
    # Save
    filename = "e8_standard_model_geometry.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"\n✓ Image saved to {filename}")
    
    # Also save a rotated view
    ax.view_init(elev=30, azim=45)
    filename2 = "e8_standard_model_geometry_rotated.png"
    plt.savefig(filename2, dpi=300, bbox_inches='tight')
    print(f"✓ Image saved to {filename2}")
    
    print()
    print("="*70)
    print("VISUALIZATION COMPLETE")
    print("="*70)
    print()
    print("Check the PNG files to see:")
    print("  - RED STARS: The 12 gauge bosons (forces)")
    print("  - BLUE DOTS: The 16 fermions (generation of matter)")
    print("  - GREY DOTS: The 212 dark sector particles")
    print()
    
    plt.show()
    
    return {
        'boson_count': len(boson_idx),
        'fermion_count': len(fermion_idx),
        'dark_count': len(dark_idx)
    }


if __name__ == "__main__":
    visualize()
