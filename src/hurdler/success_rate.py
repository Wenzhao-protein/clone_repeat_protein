#!/usr/bin/env python3
"""
HURDLER Success Rate Analysis

Tests HURDLER feasibility for random amino acid sequences of varying lengths.
For each length (4-60), generates 1000 random sequences and checks if valid
three-site combinations exist for each plasmid.

Output: Success rate plot for each plasmid across different sequence lengths.
"""

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from tqdm import tqdm
import sys
import pickle


# Standard amino acids
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

# Plasmids to test
PLASMIDS = [
    'pGEX-4T-1',
    'pMAL-c5X',
    'pET-21a(+)',
    'pET-28a(+)',
    'pET-28a(+)_start_codon',
    'pCold_I',
    'pUC18',
    'pQE-3'
]


def load_data():
    """Load df2 data"""
    df2_path = Path("./output/hurdler_three_site_combinations_df2.csv")
    
    if not df2_path.exists():
        print(f"Error: {df2_path} not found.")
        print("Please run hurdler_site_combination_analysis.ipynb first.")
        sys.exit(1)
    
    print("Loading data...")
    df2 = pd.read_csv(df2_path)
    print(f"Loaded {len(df2):,} combinations")
    
    return df2


def generate_random_sequence(length):
    """Generate a random amino acid sequence"""
    return ''.join(random.choice(AMINO_ACIDS) for _ in range(length))


def extract_3mers(sequence):
    """Extract all unique 3-mer amino acid sequences from a sequence"""
    if len(sequence) < 3:
        return set()
    
    three_mers = set()
    for i in range(len(sequence) - 2):
        three_mers.add(sequence[i:i+3])
    
    return three_mers


def find_3mer_positions(sequence):
    """Return a dict mapping 3mer -> list of start positions in the sequence."""
    pos = {}
    n = len(sequence)
    if n < 3:
        return pos
    for i in range(n - 2):
        triplet = sequence[i:i+3]
        lst = pos.get(triplet)
        if lst is None:
            pos[triplet] = [i]
        else:
            lst.append(i)
    return pos


def check_hurdler_feasibility(three_mers, plasmid, df2):
    """
    Check if HURDLER is feasible for the given 3-mers and plasmid.
    
    Returns True if at least one valid (Site I, Site II) pair exists
    where both 3-mers are in the sequence.
    """
    plasmid_col = f'{plasmid}_compatible'
    
    if plasmid_col not in df2.columns:
        return False
    
    # Filter for this plasmid
    df_plasmid = df2[df2[plasmid_col] == True]
    
    if len(df_plasmid) == 0:
        return False
    
    # Check if any combination has both Site I and Site II 3-mers in our set
    for _, row in df_plasmid.iterrows():
        site_i_aa = row['site_i_3mer_aa']
        site_ii_aa = row['site_ii_3mer_aa']
        
        if site_i_aa in three_mers and site_ii_aa in three_mers:
            return True
    
    return False


def test_sequence_length(length, num_tests, df2, plasmids):
    """
    Test HURDLER feasibility for random sequences of a given length.
    
    Returns a dictionary with success counts for each plasmid.
    """
    results = {plasmid: 0 for plasmid in plasmids}
    
    for _ in range(num_tests):
        # Generate random sequence
        sequence = generate_random_sequence(length)
        
        # Extract 3-mers
        three_mers = extract_3mers(sequence)
        
        # Check feasibility for each plasmid
        for plasmid in plasmids:
            if check_hurdler_feasibility(three_mers, plasmid, df2):
                results[plasmid] += 1
    
    return results


def load_fast_match_package(path='./output/hurdler_fast_match_package.pkl'):
    """Load the optimized fast-match package containing pattern lookup and indices."""
    pkg_path = Path(path)
    if not pkg_path.exists():
        print(f"Error: {pkg_path} not found. Please run hurdler_success_rate_optimized.ipynb to generate it.")
        return None
    with open(pkg_path, 'rb') as f:
        return pickle.load(f)


def build_site_ii_to_site_i_dict(fast_pkg, plasmid):
    """Build mapping: site_ii_3mer_aa -> list[(site_i_3mer_aa, direction)] for a plasmid."""
    mapping = {}
    plook = fast_pkg.get('pattern_lookup', {})
    for pattern, entries in plook.items():
        for entry in entries:
            # Ensure plasmid compatibility
            plasmids_info = entry.get('plasmids', {})
            if not plasmids_info.get(plasmid, False):
                continue
            ii = entry.get('3mer_aa_ii')
            i = entry.get('3mer_aa_i')
            direction = entry.get('search_direction')
            if not ii or not i or not direction:
                continue
            lst = mapping.get(ii)
            item = (i, direction)
            if lst is None:
                mapping[ii] = [item]
            else:
                # deduplicate
                if item not in lst:
                    lst.append(item)
    return mapping


def check_hurdler_success_pattern(sequence, mapping, module_length):
    """
    Check success using pattern-based mapping with direction and index constraints.
    Conditions:
    - Sequence is two repeats of an internal module of length L (provided separately)
    - For a pair (site_ii_3mer -> (site_i_3mer, direction)):
      index difference d must satisfy 5 < d < L
      direction 'right': site_i occurs before site_ii (d = pos_ii - pos_i)
      direction 'left' : site_i occurs after site_ii (d = pos_i - pos_ii)
    Returns True if any valid pair exists; otherwise False.
    """
    pos = find_3mer_positions(sequence)
    # Iterate over all site_ii keys found in sequence
    for site_ii, candidates in mapping.items():
        ii_positions = pos.get(site_ii, [])
        if not ii_positions:
            continue
        for (site_i, direction) in candidates:
            i_positions = pos.get(site_i, [])
            if not i_positions:
                continue
            # Check all combinations of positions
            for pii in ii_positions:
                for pi in i_positions:
                    if direction == 'right':
                        d = pii - pi
                        if d > 5 and d < module_length and pi < pii:
                            return True
                    else:  # 'left'
                        d = pi - pii
                        if d > 5 and d < module_length and pii < pi:
                            return True
    return False


def test_sequence_length_pattern(length, num_tests, plasmids, fast_pkg):
    """
    Test success rate using pattern-based logic on sequences that are two repeats
    of an internal module of the given length.
    """
    results = {plasmid: 0 for plasmid in plasmids}
    # Precompute per-plasmid mapping
    per_plasmid_mapping = {pl: build_site_ii_to_site_i_dict(fast_pkg, pl) for pl in plasmids}
    for _ in range(num_tests):
        module = generate_random_sequence(length)
        sequence = module + module
        for plasmid in plasmids:
            mapping = per_plasmid_mapping.get(plasmid, {})
            if not mapping:
                continue
            if check_hurdler_success_pattern(sequence, mapping, module_length=length):
                results[plasmid] += 1
    return results


def run_analysis(min_length=4, max_length=60, num_tests=1000, use_pattern_based=False):
    """
    Run complete success rate analysis.
    
    Parameters:
    -----------
    min_length : int
        Minimum sequence length to test
    max_length : int
        Maximum sequence length to test
    num_tests : int
        Number of random sequences to test per length
    
    Returns:
    --------
    DataFrame with results
    """
    print("\n" + "="*80)
    print("HURDLER SUCCESS RATE ANALYSIS")
    print("="*80)
    print(f"Sequence length range: {min_length}-{max_length}")
    print(f"Tests per length: {num_tests}")
    print(f"Plasmids: {len(PLASMIDS)}")
    print("="*80 + "\n")
    
    # Load data or fast package depending on mode
    df2 = None
    fast_pkg = None
    if use_pattern_based:
        fast_pkg = load_fast_match_package()
        if fast_pkg is None:
            print("Falling back to df2-based feasibility due to missing fast package.")
            df2 = load_data()
            use_pattern_based = False
    else:
        df2 = load_data()
    
    # Initialize results storage
    results_list = []
    
    # Test each length
    lengths = range(min_length, max_length + 1)
    
    print("\nRunning tests...")
    for length in tqdm(lengths, desc="Testing lengths"):
        # Test this length
        if use_pattern_based:
            length_results = test_sequence_length_pattern(length, num_tests, PLASMIDS, fast_pkg)
        else:
            length_results = test_sequence_length(length, num_tests, df2, PLASMIDS)
        
        # Store results
        for plasmid, success_count in length_results.items():
            success_rate = (success_count / num_tests) * 100
            results_list.append({
                'length': length,
                'plasmid': plasmid,
                'success_count': success_count,
                'total_tests': num_tests,
                'success_rate': success_rate
            })
    
    # Create DataFrame
    results_df = pd.DataFrame(results_list)
    
    return results_df


def plot_results(results_df, output_path='./output/hurdler_success_rate_analysis.png'):
    """
    Create line plot showing success rates across sequence lengths.
    
    Parameters:
    -----------
    results_df : DataFrame
        Results from run_analysis()
    output_path : str
        Path to save the plot
    """
    print("\nGenerating plot...")
    
    # Set style
    sns.set_style('whitegrid')
    plt.figure(figsize=(12, 7))
    
    # Plot each plasmid
    for plasmid in PLASMIDS:
        data = results_df[results_df['plasmid'] == plasmid]
        plt.plot(data['length'], data['success_rate'], 
                marker='o', label=plasmid, linewidth=2, markersize=4)
    
    # Formatting
    plt.xlabel('Sequence Length (amino acids)', fontsize=12)
    plt.ylabel('Success Rate (%)', fontsize=12)
    plt.title('HURDLER Feasibility Success Rate vs. Sequence Length\n(1000 random sequences per length)', 
              fontsize=14, fontweight='bold')
    plt.legend(title='Plasmid', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(results_df['length'].min() - 1, results_df['length'].max() + 1)
    plt.ylim(-5, 105)
    
    # Add horizontal lines at key percentages
    for y in [25, 50, 75, 100]:
        plt.axhline(y=y, color='gray', linestyle='--', alpha=0.3, linewidth=0.5)
    
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    
    # Show
    plt.show()


def print_summary(results_df):
    """Print summary statistics"""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    # Overall statistics by plasmid
    print("\nAverage success rate by plasmid (across all lengths):")
    avg_by_plasmid = results_df.groupby('plasmid')['success_rate'].mean().sort_values(ascending=False)
    for plasmid, rate in avg_by_plasmid.items():
        print(f"  {plasmid:<30}: {rate:>6.2f}%")
    
    # Statistics by length range
    print("\nSuccess rate by sequence length range:")
    
    # Define length ranges
    ranges = [
        (4, 10, "Very short (4-10)"),
        (11, 20, "Short (11-20)"),
        (21, 30, "Medium (21-30)"),
        (31, 40, "Long (31-40)"),
        (41, 50, "Very long (41-50)"),
        (51, 60, "Extra long (51-60)")
    ]
    
    for min_len, max_len, label in ranges:
        data = results_df[(results_df['length'] >= min_len) & (results_df['length'] <= max_len)]
        if len(data) > 0:
            avg_rate = data['success_rate'].mean()
            print(f"  {label:<25}: {avg_rate:>6.2f}%")
    
    # Best and worst cases
    print("\nBest performing combinations:")
    top_5 = results_df.nlargest(5, 'success_rate')[['length', 'plasmid', 'success_rate']]
    for idx, row in top_5.iterrows():
        print(f"  Length {row['length']:>2}, {row['plasmid']:<30}: {row['success_rate']:>6.2f}%")
    
    print("\nWorst performing combinations:")
    bottom_5 = results_df.nsmallest(5, 'success_rate')[['length', 'plasmid', 'success_rate']]
    for idx, row in bottom_5.iterrows():
        print(f"  Length {row['length']:>2}, {row['plasmid']:<30}: {row['success_rate']:>6.2f}%")
    
    # Length at which 50%, 75%, 90% success is achieved
    print("\nSequence length to achieve target success rates:")
    for target in [50, 75, 90]:
        print(f"\n  {target}% success rate:")
        for plasmid in PLASMIDS:
            data = results_df[results_df['plasmid'] == plasmid]
            achieving = data[data['success_rate'] >= target]
            if len(achieving) > 0:
                min_length = achieving['length'].min()
                print(f"    {plasmid:<30}: {min_length:>2} aa")
            else:
                print(f"    {plasmid:<30}: Not achieved")
    
    print("\n" + "="*80)


def export_results(results_df, output_path='./output/hurdler_success_rate_data.csv'):
    """Export results to CSV"""
    results_df.to_csv(output_path, index=False)
    print(f"\nResults exported to: {output_path}")


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze HURDLER success rates')
    parser.add_argument('--min-length', type=int, default=4, help='Minimum sequence length')
    parser.add_argument('--max-length', type=int, default=60, help='Maximum sequence length')
    parser.add_argument('--num-tests', type=int, default=1000, help='Number of tests per length')
    parser.add_argument('--no-plot', action='store_true', help='Skip plotting')
    parser.add_argument('--output', type=str, default='./output/hurdler_success_rate_analysis.png',
                        help='Output plot path')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    parser.add_argument('--pattern-based', action='store_true', help='Use pattern-based direction/index logic with repeated modules')
    
    args = parser.parse_args()
    
    # Set random seed if provided
    if args.seed:
        random.seed(args.seed)
        np.random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")
    
    # Run analysis
    results_df = run_analysis(
        min_length=args.min_length,
        max_length=args.max_length,
        num_tests=args.num_tests,
        use_pattern_based=args.pattern_based
    )
    
    # Export results
    export_results(results_df)
    
    # Print summary
    print_summary(results_df)
    
    # Plot results
    if not args.no_plot:
        plot_results(results_df, output_path=args.output)
    
    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()
