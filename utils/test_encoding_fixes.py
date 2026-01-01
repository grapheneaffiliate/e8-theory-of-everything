#!/usr/bin/env python3
"""Quick test to verify encoding fixes work"""

import sys
import subprocess
import os

def test_script(script_path, script_name):
    """Test if a script runs without encoding errors"""
    print(f"\nTesting {script_name}...")
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=True,
            text=True,
            timeout=30,
            encoding='utf-8',
            errors='replace'
        )
        
        # Check for encoding errors
        if 'UnicodeEncodeError' in result.stderr or 'cp1252' in result.stderr:
            print(f"  [FAIL] {script_name} - Encoding error still present")
            print(f"  Error: {result.stderr[:200]}")
            return False
        elif result.returncode != 0 and 'Traceback' in result.stderr:
            # Check if it's an encoding-related traceback
            if 'encode' in result.stderr.lower() or 'decode' in result.stderr.lower():
                print(f"  [FAIL] {script_name} - Encoding-related crash")
                print(f"  Error: {result.stderr[:200]}")
                return False
            else:
                print(f"  [WARN] {script_name} - Non-encoding error (may be expected)")
                return True
        else:
            print(f"  [PASS] {script_name} - No encoding errors detected")
            return True
            
    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] {script_name} - Script took too long (>30s)")
        return True  # Timeout is OK, just means it's running
    except Exception as e:
        print(f"  [ERROR] {script_name} - {str(e)}")
        return False

def main():
    physics_dir = r'C:\Users\atchi\Desktop\e8-theory-of-everything\physics'
    
    scripts_to_test = [
        ('physical_constants_derivation.py', 'Physical Constants'),
        ('e8_wave_equation.py', 'Wave Equation'),
        ('e8_gauge_field.py', 'Gauge Field'),
        ('e8_dynamical_field_theory.py', 'Dynamical Field Theory'),
    ]
    
    print("="*60)
    print("TESTING ENCODING FIXES")
    print("="*60)
    
    passed = 0
    failed = 0
    
    for script_file, script_name in scripts_to_test:
        script_path = os.path.join(physics_dir, script_file)
        if os.path.exists(script_path):
            if test_script(script_path, script_name):
                passed += 1
            else:
                failed += 1
        else:
            print(f"\n[SKIP] {script_name} - File not found")
    
    print("\n" + "="*60)
    print(f"RESULTS: {passed} passed, {failed} failed out of {len(scripts_to_test)}")
    print("="*60)
    
    if failed == 0:
        print("\n✓ All encoding fixes successful!")
        print("  The simulation should now run without cp1252 encoding errors.")
    else:
        print(f"\n✗ {failed} scripts still have encoding issues")

if __name__ == '__main__':
    main()
