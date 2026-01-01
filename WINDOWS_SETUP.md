# Windows Setup Guide - E8 Theory of Everything

## Encoding Issue Fix

On Windows, the default console encoding (cp1252) cannot display Unicode mathematical symbols used in the simulation output. This has been resolved.

## Solution

### Method 1: Use the Batch File Wrapper (RECOMMENDED)

Simply run the provided batch file instead of calling Python directly:

```batch
run_e8_theory.bat          # Run all modules
run_e8_theory.bat --quick  # Run 3 essential modules
run_e8_theory.bat --full   # Run all 19 modules
run_e8_theory.bat --dynamics  # Run dynamics simulations only
run_e8_theory.bat --proof  # Run topological proof only
```

The batch file automatically sets `PYTHONIOENCODING=utf-8` before running the simulation.

### Method 2: Set Environment Variable Manually

Before running the Python script, set the encoding:

**PowerShell:**
```powershell
$env:PYTHONIOENCODING="utf-8"
python run_unified_theory.py
```

**Command Prompt:**
```cmd
set PYTHONIOENCODING=utf-8
python run_unified_theory.py
```

### Method 3: One-Line Command

```cmd
set PYTHONIOENCODING=utf-8 && python run_unified_theory.py
```

## What Was Fixed

The following files were patched with UTF-8 encoding handling:
- ✓ `run_unified_theory.py` - Main runner script
- ✓ `physics/physical_constants_derivation.py`
- ✓ `physics/e8_wave_equation.py`
- ✓ `physics/e8_gauge_field.py`
- ✓ `physics/e8_dynamical_field_theory.py`

Each file now includes:
```python
if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except AttributeError:
        import io
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
```

However, the environment variable approach is more reliable as it sets encoding before Python starts.

## Verification

To verify the fixes work, run:

```cmd
python test_encoding_fixes.py
```

Expected output:
```
============================================================
TESTING ENCODING FIXES
============================================================

Testing Physical Constants...
  [PASS] Physical Constants - No encoding errors detected

Testing Wave Equation...
  [PASS] Wave Equation - No encoding errors detected

Testing Gauge Field...
  [PASS] Gauge Field - No encoding errors detected

Testing Dynamical Field Theory...
  [PASS] Dynamical Field Theory - No encoding errors detected

============================================================
RESULTS: 4 passed, 0 failed out of 4
============================================================

✓ All encoding fixes successful!
```

## Troubleshooting

If you still see encoding errors:

1. **Check Python version:** Ensure you're using Python 3.7+
   ```cmd
   python --version
   ```

2. **Use the batch file:** Always use `run_e8_theory.bat` instead of calling Python directly

3. **Check console:** Some older Windows configurations may need chcp 65001:
   ```cmd
   chcp 65001
   run_e8_theory.bat
   ```

4. **Redirect output:** If console display fails, redirect to file:
   ```cmd
   run_e8_theory.bat --quick > output.txt 2>&1
   type output.txt
   ```

## Quick Start

**Run the complete simulation:**
```cmd
run_e8_theory.bat
```

**Run a quick test (3 modules):**
```cmd
run_e8_theory.bat --quick
```

## Notes

- Mathematical symbols (∂, θ, φ, α, etc.) will still appear garbled in some console windows, but the program will run successfully
- The actual physics calculations are unaffected - only the display is impacted
- Output files and logs will contain proper UTF-8 encoding

## Files Created

- `run_e8_theory.bat` - Windows launcher with automatic UTF-8 encoding
- `fix_encoding.py` - Script that patched the Python files
- `test_encoding_fixes.py` - Verification script
- `WINDOWS_SETUP.md` - This documentation

## Support

If issues persist, check:
1. Python installation is 64-bit and recent (3.7+)
2. Running from cmd.exe or PowerShell (not IDLE)
3. Have write permissions in the directory
4. No antivirus blocking Python execution

All encoding issues should now be resolved! Enjoy exploring the E8 Theory of Everything.
