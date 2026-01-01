#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fix Windows Console Encoding Issues for Unicode Output
Adds UTF-8 encoding configuration to all physics scripts
"""

import os
import sys
import re

ENCODING_FIX = '''import sys

# Fix Windows console encoding for Unicode output
if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except AttributeError:
        import io
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
'''

def fix_file(filepath):
    """Add encoding fix to a Python file if not already present."""
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # Check if already fixed
    if 'sys.stdout.reconfigure' in content or 'Fix Windows console encoding' in content:
        print(f"  [SKIP] {os.path.basename(filepath)} - already fixed")
        return False
    
    # Find the import section
    lines = content.split('\n')
    import_end_idx = 0
    
    for i, line in enumerate(lines):
        if line.startswith('import ') or line.startswith('from '):
            import_end_idx = i + 1
        elif import_end_idx > 0 and line.strip() and not line.startswith('#') and not line.startswith('"""') and not line.startswith("'''"):
            break
    
    # Insert encoding fix after imports
    if import_end_idx > 0:
        lines.insert(import_end_idx, ENCODING_FIX)
        
        # Write back
        new_content = '\n'.join(lines)
        with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
            f.write(new_content)
        
        print(f"  [FIXED] {os.path.basename(filepath)}")
        return True
    else:
        print(f"  [ERROR] {os.path.basename(filepath)} - couldn't find import section")
        return False

def main():
    """Fix all physics scripts with Unicode output."""
    physics_dir = r'C:\Users\atchi\Desktop\e8-theory-of-everything\physics'
    
    files_to_fix = [
        'e8_wave_equation.py',
        'e8_gauge_field.py',
        'e8_dynamical_field_theory.py',
    ]
    
    print("=" * 60)
    print("FIXING WINDOWS CONSOLE ENCODING ISSUES")
    print("=" * 60)
    print()
    
    fixed_count = 0
    for filename in files_to_fix:
        filepath = os.path.join(physics_dir, filename)
        if os.path.exists(filepath):
            if fix_file(filepath):
                fixed_count += 1
        else:
            print(f"  [SKIP] {filename} - file not found")
    
    print()
    print("=" * 60)
    print(f"FIXED {fixed_count}/{len(files_to_fix)} files")
    print("=" * 60)

if __name__ == '__main__':
    main()
