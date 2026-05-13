#!/usr/bin/env python3
"""HURDLER Site Query Tool.

Command-line and library interface for looking up valid HURDLER three-site
combinations from a pre-computed ``df2`` table.

The lookup answers the question: *given two 3-amino-acid windows
(Site I and Site II) within my repeat unit, and a target plasmid, which
restriction-enzyme triples can I use?*

Inputs
------
``./output/hurdler_three_site_combinations_df2.csv`` — produced by
:mod:`hurdler.pipeline` (or its notebook equivalent under
``notebooks/hurdler/``). The path is configurable via ``--df2``.

Examples
--------
Single query::

    PYTHONPATH=src python -m hurdler.query \
        --site-i-aa "NEQ" --site-ii-aa "IQA" --plasmid "pET-28a(+)"

Batch query::

    PYTHONPATH=src python -m hurdler.query \
        --batch data/example_batch_query.csv --output output/results.csv

List the 3-mer AA sequences present in the lookup::

    PYTHONPATH=src python -m hurdler.query --list

See ``docs/workflows/hurdler_site_combinations.md`` for the broader
workflow context.
"""

import pandas as pd
import argparse
import sys
from pathlib import Path


def load_data(df2_path="./output/hurdler_three_site_combinations_df2.csv"):
    """Load the pre-computed ``df2`` table produced by :mod:`hurdler.pipeline`."""
    df2_path = Path(df2_path)

    if not df2_path.exists():
        print(f"Error: {df2_path} not found.")
        print("Please run `python -m hurdler.pipeline` (or the equivalent "
              "notebook under notebooks/hurdler/) to generate it.")
        sys.exit(1)

    return pd.read_csv(df2_path)


def find_hurdler_sites(site_i_3mer_aa, site_ii_3mer_aa, plasmid, df2):
    """
    Find all valid three-site combinations for given 3mer AA sequences and plasmid.
    
    Parameters:
    -----------
    site_i_3mer_aa : str
        3-amino acid sequence for Site I (seamless insert)
    site_ii_3mer_aa : str
        3-amino acid sequence for Site II (silent mutation)
    plasmid : str
        Plasmid name (e.g., 'pET-28a(+)')
    df2 : DataFrame
        The df2 DataFrame containing all combinations
    
    Returns:
    --------
    DataFrame with matching combinations, or None if no matches found
    """
    # Check if plasmid is valid
    plasmid_col = f'{plasmid}_compatible'
    if plasmid_col not in df2.columns:
        available = [col.replace('_compatible', '') for col in df2.columns if col.endswith('_compatible')]
        print(f"Error: Plasmid '{plasmid}' not found.")
        print(f"Available plasmids: {', '.join(available)}")
        return None
    
    # Filter for matching 3mer AA sequences and plasmid compatibility
    results = df2[
        (df2['site_i_3mer_aa'] == site_i_3mer_aa) &
        (df2['site_ii_3mer_aa'] == site_ii_3mer_aa) &
        (df2[plasmid_col] == True)
    ].copy()
    
    return results


def format_results(results, site_i_3mer_aa, site_ii_3mer_aa, plasmid, verbose=True):
    """Format and display query results"""
    if results is None or len(results) == 0:
        print(f"\nNo valid combinations found for:")
        print(f"  Site I 3mer AA:  {site_i_3mer_aa}")
        print(f"  Site II 3mer AA: {site_ii_3mer_aa}")
        print(f"  Plasmid:         {plasmid}")
        return
    
    print(f"\n{'='*80}")
    print(f"Found {len(results)} valid combination(s):")
    print(f"  Site I 3mer AA:  {site_i_3mer_aa}")
    print(f"  Site II 3mer AA: {site_ii_3mer_aa}")
    print(f"  Plasmid:         {plasmid}")
    print(f"{'='*80}")
    
    if verbose:
        for idx, (_, row) in enumerate(results.iterrows(), 1):
            print(f"\nCombination {idx}:")
            print(f"  Site I:   {row['site_i']:<15} (ovhg: {row['ovhg_i']:>3}, frame: {row['site_i_frame']})")
            print(f"            DNA: {row['site_i_dna']}")
            print(f"  Site II:  {row['site_ii']:<15} (ovhg: {row['ovhg_ii']:>3}, frame: {row['site_ii_frame']}, silent mutation)")
            print(f"            DNA: {row['site_ii_dna']}")
            print(f"            Mutated: {row['site_ii_dna_mutated']}")
            print(f"  Site III: {row['site_iii']:<15} (ovhg: {row['ovhg_iii']:>3}, Type IIS)")
    else:
        # Compact format
        print("\n{:<15} {:<15} {:<15} {:<5} {:<5} {:<5}".format(
            "Site I", "Site II", "Site III", "ovhg_I", "ovhg_II", "ovhg_III"
        ))
        print("-" * 80)
        for _, row in results.iterrows():
            print("{:<15} {:<15} {:<15} {:<5} {:<5} {:<5}".format(
                row['site_i'], row['site_ii'], row['site_iii'],
                str(row['ovhg_i']), str(row['ovhg_ii']), str(row['ovhg_iii'])
            ))


def batch_query(input_file, output_file, df2):
    """Process batch queries from CSV file"""
    print(f"Processing batch queries from {input_file}...")
    
    try:
        input_df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    
    # Validate required columns
    required_cols = ['site_i_3mer_aa', 'site_ii_3mer_aa', 'plasmid']
    missing_cols = [col for col in required_cols if col not in input_df.columns]
    if missing_cols:
        print(f"Error: Input file missing required columns: {', '.join(missing_cols)}")
        sys.exit(1)
    
    # Process each query
    all_results = []
    for idx, row in input_df.iterrows():
        site_i_aa = row['site_i_3mer_aa']
        site_ii_aa = row['site_ii_3mer_aa']
        plasmid = row['plasmid']
        
        results = find_hurdler_sites(site_i_aa, site_ii_aa, plasmid, df2)
        
        if results is not None and len(results) > 0:
            results['query_id'] = idx
            results['query_site_i_3mer_aa'] = site_i_aa
            results['query_site_ii_3mer_aa'] = site_ii_aa
            results['query_plasmid'] = plasmid
            all_results.append(results)
            print(f"  Query {idx+1}: Found {len(results)} combination(s)")
        else:
            print(f"  Query {idx+1}: No combinations found")
    
    if all_results:
        final_results = pd.concat(all_results, ignore_index=True)
        final_results.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
        print(f"Total combinations found: {len(final_results)}")
    else:
        print("\nNo valid combinations found for any query.")


def list_available_3mer_aa(df2):
    """List available 3mer AA sequences"""
    print("\n=== Available 3mer AA Sequences ===")
    
    site_i_aa = df2['site_i_3mer_aa'].unique()
    site_ii_aa = df2['site_ii_3mer_aa'].unique()
    
    print(f"\nSite I (Seamless Insert): {len(site_i_aa)} unique sequences")
    print(f"Examples: {', '.join(sorted(site_i_aa)[:20])}")
    
    print(f"\nSite II (Silent Mutation): {len(site_ii_aa)} unique sequences")
    print(f"Examples: {', '.join(sorted(site_ii_aa)[:20])}")
    
    plasmids = [col.replace('_compatible', '') for col in df2.columns if col.endswith('_compatible')]
    print(f"\nAvailable plasmids: {', '.join(plasmids)}")


def main():
    parser = argparse.ArgumentParser(
        description='Query HURDLER three-site combinations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single query
  python hurdler_query.py --site-i-aa "ABC" --site-ii-aa "DEF" --plasmid "pET-28a(+)"
  
  # Batch query from CSV file
  python hurdler_query.py --batch input.csv --output results.csv
  
  # List available 3mer AA sequences
  python hurdler_query.py --list
        """
    )
    
    parser.add_argument('--site-i-aa', type=str, help='3-amino acid sequence for Site I')
    parser.add_argument('--site-ii-aa', type=str, help='3-amino acid sequence for Site II')
    parser.add_argument('--plasmid', type=str, help='Plasmid name')
    parser.add_argument('--batch', type=str, help='Input CSV file for batch queries')
    parser.add_argument('--output', type=str, help='Output CSV file for batch results')
    parser.add_argument('--list', action='store_true', help='List available 3mer AA sequences')
    parser.add_argument('--compact', action='store_true', help='Use compact output format')
    
    args = parser.parse_args()
    
    # Load data
    df2 = load_data()
    
    # List mode
    if args.list:
        list_available_3mer_aa(df2)
        return
    
    # Batch mode
    if args.batch:
        if not args.output:
            print("Error: --output required for batch mode")
            sys.exit(1)
        batch_query(args.batch, args.output, df2)
        return
    
    # Single query mode
    if not all([args.site_i_aa, args.site_ii_aa, args.plasmid]):
        parser.print_help()
        sys.exit(1)
    
    results = find_hurdler_sites(args.site_i_aa, args.site_ii_aa, args.plasmid, df2)
    format_results(results, args.site_i_aa, args.site_ii_aa, args.plasmid, verbose=not args.compact)


if __name__ == '__main__':
    main()
