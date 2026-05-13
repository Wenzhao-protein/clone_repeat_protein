#!/usr/bin/env python3
"""Quick launcher for testing GUI"""

import subprocess
import sys

print("Launching interactive_gui_new.py in visualization environment...")
print("=" * 80)

try:
    result = subprocess.run(
        ["conda", "run", "-n", "visualization", "python3", "interactive_gui_new.py"],
        cwd="/home/wenzhao/github_repo/clone_repeat_protein/agarose_gel_analysis/input",
        input="1\n",
        text=True,
        capture_output=False,
        timeout=10
    )
    print("=" * 80)
    print(f"Process completed with exit code: {result.returncode}")
except subprocess.TimeoutExpired:
    print("\n" + "=" * 80)
    print("GUI launched successfully (timeout means it's running)")
    print("=" * 80)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
