#!/usr/bin/env python3
"""
Execute the HURDLER notebook to generate df1 and df2 data files.
"""

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import sys
from pathlib import Path

def run_notebook(notebook_path, timeout=3600):
    """
    Execute a Jupyter notebook
    
    Parameters:
    -----------
    notebook_path : str or Path
        Path to the notebook file
    timeout : int
        Maximum execution time in seconds (default: 3600 = 1 hour)
    """
    notebook_path = Path(notebook_path)
    
    if not notebook_path.exists():
        print(f"Error: Notebook not found at {notebook_path}")
        sys.exit(1)
    
    print(f"Loading notebook: {notebook_path}")
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)
    
    print("Executing notebook...")
    print("This may take 10-30 minutes depending on your system...")
    
    ep = ExecutePreprocessor(timeout=timeout, kernel_name='python3')
    
    try:
        ep.preprocess(nb, {'metadata': {'path': str(notebook_path.parent)}})
        print("\n✓ Notebook execution completed successfully!")
        
        # Save the executed notebook
        output_path = notebook_path.parent / f"{notebook_path.stem}_executed.ipynb"
        with open(output_path, 'w', encoding='utf-8') as f:
            nbformat.write(nb, f)
        print(f"✓ Executed notebook saved to: {output_path}")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Error executing notebook: {e}")
        return False


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Execute HURDLER analysis notebook')
    parser.add_argument('--notebook', type=str, 
                        default='./hurdler_site_combination_analysis.ipynb',
                        help='Path to notebook file')
    parser.add_argument('--timeout', type=int, default=3600,
                        help='Timeout in seconds (default: 3600)')
    
    args = parser.parse_args()
    
    success = run_notebook(args.notebook, timeout=args.timeout)
    sys.exit(0 if success else 1)
