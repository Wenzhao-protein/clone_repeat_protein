#!/usr/bin/env python3
"""Simple test to check if window appears"""

import tkinter as tk
import sys

print("Creating test window...")
print("If you can see this message but no window, there's a display problem.")
print()

try:
    # Create a very simple window
    root = tk.Tk()
    root.title("TEST WINDOW - Can you see this?")
    root.geometry("500x300")
    root.configure(bg='white')
    
    # Add large text
    label = tk.Label(
        root,
        text="✓ SUCCESS!\n\nIf you can see this window,\nTkinter is working correctly!",
        font=("Arial", 18),
        bg='white',
        fg='green',
        pady=50
    )
    label.pack()
    
    # Add button
    btn = tk.Button(
        root,
        text="Close Window",
        command=root.destroy,
        font=("Arial", 14),
        padx=20,
        pady=10
    )
    btn.pack()
    
    # Make sure window appears on top
    root.lift()
    root.attributes('-topmost', True)
    root.after_idle(root.attributes, '-topmost', False)
    root.focus_force()
    
    # Force update
    root.update()
    root.update_idletasks()
    
    print("✓ Window created successfully!")
    print("✓ Window should be visible now")
    print()
    print("Close the window to continue...")
    
    # Run main loop
    root.mainloop()
    
    print("\n✓ Window closed normally")
    print("=" * 60)
    print("Tkinter works! The problem might be:")
    print("  1. Window size too large for your screen")
    print("  2. Matplotlib rendering issue")
    print("  3. Window appearing on different monitor/workspace")
    print("=" * 60)
    
except Exception as e:
    print(f"\n✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
