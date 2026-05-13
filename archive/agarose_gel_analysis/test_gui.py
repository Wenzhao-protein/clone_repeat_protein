#!/usr/bin/env python3
"""Test script to verify GUI can launch"""

import sys
import traceback

print("=" * 80)
print("Testing GUI Launch")
print("=" * 80)

try:
    print("\n1. Testing imports...")
    import tkinter as tk
    from tkinter import ttk
    import numpy as np
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from scipy.ndimage import rotate as scipy_rotate
    from numpy.polynomial import Polynomial
    print("   ✓ All imports successful")
    
    print("\n2. Testing scn_reader...")
    from scn_reader import read_scn_file_enhanced, get_image_info
    print("   ✓ scn_reader imported")
    
    print("\n3. Testing file read...")
    scn_file = "2022-12-12_Q6_Q6_Q5_M2_L2AG2_L2BG2_L3AG2_L3BG2.scn"
    img_data, dimensions, error = read_scn_file_enhanced(scn_file, prefer_mode='crop')
    if error:
        print(f"   ✗ Error reading file: {error}")
        sys.exit(1)
    print(f"   ✓ File read successful: {img_data.shape}")
    
    print("\n4. Testing interactive_gui_new import...")
    import interactive_gui_new
    print("   ✓ Module imported")
    
    print("\n5. Testing Tkinter window creation...")
    root = tk.Tk()
    root.title("Test Window")
    root.geometry("400x300")
    label = tk.Label(root, text="Test GUI Window", font=("Arial", 16))
    label.pack(pady=50)
    
    # Auto-close after 2 seconds
    def close_window():
        print("   ✓ Window created successfully")
        root.destroy()
        print("\n" + "=" * 80)
        print("ALL TESTS PASSED - GUI should work!")
        print("=" * 80)
    
    root.after(2000, close_window)
    root.mainloop()
    
except Exception as e:
    print(f"\n✗ ERROR: {e}")
    print("\nTraceback:")
    traceback.print_exc()
    sys.exit(1)
