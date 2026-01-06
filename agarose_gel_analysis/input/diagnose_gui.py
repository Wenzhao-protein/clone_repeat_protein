#!/usr/bin/env python3
"""Diagnostic script to check GUI status and issues"""

import os
import sys
import subprocess
import glob

print("=" * 80)
print("GUI Diagnostic Tool")
print("=" * 80)

# 1. Check files
print("\n1. Checking files...")
files_to_check = [
    "interactive_gui_new.py",
    "scn_reader.py",
    "2022-12-12_Q6_Q6_Q5_M2_L2AG2_L2BG2_L3AG2_L3BG2.scn"
]

for f in files_to_check:
    exists = "✓" if os.path.exists(f) else "✗"
    size = os.path.getsize(f) if os.path.exists(f) else 0
    print(f"   {exists} {f} ({size} bytes)")

# 2. Check imports
print("\n2. Checking Python imports...")
try:
    import tkinter as tk
    print("   ✓ tkinter")
except Exception as e:
    print(f"   ✗ tkinter: {e}")

try:
    import numpy as np
    print(f"   ✓ numpy {np.__version__}")
except Exception as e:
    print(f"   ✗ numpy: {e}")

try:
    import matplotlib
    print(f"   ✓ matplotlib {matplotlib.__version__}")
except Exception as e:
    print(f"   ✗ matplotlib: {e}")

try:
    import scipy
    print(f"   ✓ scipy {scipy.__version__}")
except Exception as e:
    print(f"   ✗ scipy: {e}")

# 3. Check if GUI can import
print("\n3. Checking interactive_gui_new module...")
try:
    import interactive_gui_new
    print("   ✓ Module imports successfully")
    
    # Check key classes
    if hasattr(interactive_gui_new, 'InteractiveGelAnalyzer'):
        print("   ✓ InteractiveGelAnalyzer class found")
    if hasattr(interactive_gui_new, 'GelAnalysisGUI'):
        print("   ✓ GelAnalysisGUI class found")
    if hasattr(interactive_gui_new, 'main'):
        print("   ✓ main function found")
        
except Exception as e:
    print(f"   ✗ Import error: {e}")
    import traceback
    traceback.print_exc()

# 4. Check DISPLAY
print("\n4. Checking display environment...")
display = os.environ.get('DISPLAY', 'NOT SET')
print(f"   DISPLAY={display}")

if display and display != 'NOT SET':
    try:
        root = tk.Tk()
        root.title("Test")
        root.geometry("200x100")
        root.after(500, root.destroy)
        root.mainloop()
        print("   ✓ Can create Tkinter windows")
    except Exception as e:
        print(f"   ✗ Cannot create Tkinter window: {e}")

# 5. Check running processes
print("\n5. Checking for running GUI processes...")
try:
    result = subprocess.run(
        ["ps", "aux"],
        capture_output=True,
        text=True,
        timeout=5
    )
    gui_processes = [line for line in result.stdout.split('\n') 
                     if 'interactive_gui' in line and 'python' in line]
    if gui_processes:
        print("   Found running GUI processes:")
        for proc in gui_processes:
            print(f"      {proc}")
    else:
        print("   No GUI processes currently running")
except Exception as e:
    print(f"   Could not check processes: {e}")

# 6. Quick syntax check
print("\n6. Checking Python syntax...")
try:
    result = subprocess.run(
        ["python3", "-m", "py_compile", "interactive_gui_new.py"],
        capture_output=True,
        text=True,
        timeout=5
    )
    if result.returncode == 0:
        print("   ✓ No syntax errors")
    else:
        print(f"   ✗ Syntax errors found:\n{result.stderr}")
except Exception as e:
    print(f"   Could not check syntax: {e}")

print("\n" + "=" * 80)
print("Diagnostic complete!")
print("=" * 80)
print("\nTo run the GUI manually, try:")
print("  conda run -n visualization python3 interactive_gui_new.py")
print("\nOr from notebook, execute cell 6 in interactive_gui_launcher.ipynb")
print("=" * 80)
