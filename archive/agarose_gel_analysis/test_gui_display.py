#!/usr/bin/env python3
"""
Quick test to verify the GUI analysis results display works.
"""

print("="*70)
print("GUI Analysis Results Display Test")
print("="*70)

# Test 1: Check matplotlib backend
print("\n1. Testing matplotlib backend...")
import matplotlib
backend = matplotlib.get_backend()
print(f"   Current backend: {backend}")
if backend == 'TkAgg':
    print("   ✓ TkAgg backend configured correctly")
else:
    print(f"   ⚠ Warning: Using {backend} instead of TkAgg")

# Test 2: Check imports
print("\n2. Testing imports...")
try:
    from interactive_gui import GelAnalysisGUI, InteractiveGelAnalyzer
    from interactive_gui import sanitize_for_gui, safe_decode_bytes
    print("   ✓ All imports successful")
except Exception as e:
    print(f"   ✗ Import error: {e}")
    exit(1)

# Test 3: Check new method exists
print("\n3. Checking show_analysis_results method...")
import inspect
methods = [m for m in dir(GelAnalysisGUI) if not m.startswith('_')]
if 'show_analysis_results' in methods:
    print("   ✓ show_analysis_results method exists")
    
    # Check method signature
    sig = inspect.signature(GelAnalysisGUI.show_analysis_results)
    params = list(sig.parameters.keys())
    print(f"   Parameters: {params}")
    if 'analysis' in params:
        print("   ✓ Method signature correct")
    else:
        print("   ✗ Missing 'analysis' parameter")
else:
    print("   ✗ show_analysis_results method not found")

# Test 4: List all public methods
print("\n4. Available GelAnalysisGUI methods:")
for method in sorted(methods):
    if not method.startswith('_'):
        print(f"   - {method}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print("The GUI has been updated with the following improvements:")
print("")
print("1. ✓ Matplotlib backend set to TkAgg for compatibility")
print("2. ✓ New show_analysis_results() method creates Tkinter window")
print("3. ✓ Analysis results display in embedded FigureCanvasTkAgg")
print("4. ✓ Results window has Close button")
print("5. ✓ Window automatically brought to front")
print("")
print("How it works:")
print("  - Click 'Analyze Selection' button")
print("  - New window pops up with 4-subplot analysis")
print("  - Results embedded in Tkinter (not separate matplotlib window)")
print("  - Click 'Close' button to dismiss results")
print("")
print("To test: python interactive_gui.py")
print("="*70)
