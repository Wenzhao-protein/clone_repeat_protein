#!/bin/bash
# Complete diagnostic and fix script for GUI display issues

echo "========================================================================"
echo "GUI Display Diagnostic and Fix Script"
echo "========================================================================"
echo ""

# 1. Check DISPLAY
echo "1. Checking DISPLAY environment..."
echo "   DISPLAY=$DISPLAY"
if [ -z "$DISPLAY" ]; then
    echo "   ✗ DISPLAY not set! GUI cannot work without display."
    exit 1
else
    echo "   ✓ DISPLAY is set"
fi
echo ""

# 2. Test X server
echo "2. Testing X server connection..."
if xdpyinfo > /dev/null 2>&1; then
    echo "   ✓ X server is accessible"
else
    echo "   ✗ Cannot connect to X server"
    echo "   Try: export DISPLAY=:0"
    exit 1
fi
echo ""

# 3. Check window manager
echo "3. Checking window manager..."
if pgrep -x "gnome-shell\|kwin\|xfwm4\|openbox\|i3" > /dev/null; then
    WM=$(pgrep -x "gnome-shell\|kwin\|xfwm4\|openbox\|i3" | head -1)
    echo "   ✓ Window manager is running (PID: $WM)"
else
    echo "   ⚠ No common window manager detected"
fi
echo ""

# 4. Check for running GUI processes
echo "4. Checking for running GUI processes..."
GUI_PROCS=$(ps aux | grep "interactive_gui_new\.py" | grep -v grep | wc -l)
if [ $GUI_PROCS -gt 0 ]; then
    echo "   Found $GUI_PROCS GUI process(es) running:"
    ps aux | grep "interactive_gui_new\.py" | grep -v grep | awk '{print "      PID " $2 ": " $11 " " $12 " " $13}'
    echo ""
    read -p "   Kill these processes? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        pkill -f "interactive_gui_new\.py"
        echo "   ✓ Processes killed"
    fi
else
    echo "   No GUI processes currently running"
fi
echo ""

# 5. Run quick Python test
echo "5. Testing Python Tkinter..."
conda run -n visualization python3 << 'PYEOF'
import tkinter as tk
import sys

try:
    root = tk.Tk()
    root.withdraw()  # Don't show window
    root.update()
    print("   ✓ Tkinter works")
    root.destroy()
except Exception as e:
    print(f"   ✗ Tkinter error: {e}")
    sys.exit(1)
PYEOF
echo ""

# 6. Check file
echo "6. Checking interactive_gui_new.py..."
if [ -f "interactive_gui_new.py" ]; then
    SIZE=$(stat -f%z "interactive_gui_new.py" 2>/dev/null || stat -c%s "interactive_gui_new.py" 2>/dev/null)
    echo "   ✓ File exists ($SIZE bytes)"
    
    # Check for syntax errors
    if conda run -n visualization python3 -m py_compile interactive_gui_new.py 2>/dev/null; then
        echo "   ✓ No syntax errors"
    else
        echo "   ✗ Syntax errors found"
        exit 1
    fi
else
    echo "   ✗ File not found"
    exit 1
fi
echo ""

echo "========================================================================"
echo "DIAGNOSTIC COMPLETE"
echo "========================================================================"
echo ""
echo "To run the GUI:"
echo "  cd $(pwd)"
echo "  conda activate visualization"
echo "  python3 interactive_gui_new.py"
echo ""
echo "If window still doesn't appear, try:"
echo "  1. Check if window is minimized or on another workspace"
echo "  2. Alt+Tab to switch between windows"
echo "  3. Check system tray for application icon"
echo "  4. Try smaller window: edit line with geometry(\"1400x1000\")"
echo ""
echo "========================================================================"
