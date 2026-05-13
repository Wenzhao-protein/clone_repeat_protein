#!/usr/bin/env python3
"""
HURDLER Success Rate Test

Tests the success rate of finding valid HURDLER combinations for random amino acid
sequences of varying lengths across different plasmids.

For each sequence length from 4 to 60:
- Generate 1000 random amino acid sequences
- Extract all unique 3mer AA subsequences
- Check if valid HURDLER combinations exist for each plasmid
- Calculate and plot success rates
"""

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from tqdm import tqdm
import sys


# Amino acid alphabet
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def generate_random_aa_sequence(length):
    """Generate a random amino acid sequence of given length"""
    return ''.join(random.choice(AMINO_ACIDS) for _ in range(length))


def extract_unique_3mers(sequence):
    """Extract all unique 3-mer amino acid subsequences from a sequence.
    
    For repeat protein sequences, consider circular boundary conditions:
    - seq[-2:] + seq[0] forms a 3mer
    - seq[-1] + seq[:2] forms a 3mer
    """
    if len(sequence) < 3:
        return set()
    
    three_mers = set()
    
    # Standard 3mers
    for i in range(len(sequence) - 2):
        three_mers.add(sequence[i:i+3])
    
    # Circular boundary 3mers (for repeat proteins)
    if len(sequence) >= 3:
        # seq[-2], seq[-1], seq[0]
        three_mers.add(sequence[-2:] + sequence[0])
        # seq[-1], seq[0], seq[1]
        three_mers.add(sequence[-1] + sequence[:2])
    
    return three_mers


def check_hurdler_feasible(three_mer_set, plasmid, lookup_dict):
    """
    Check if a HURDLER solution exists for the given set of 3mers and plasmid.
    
    Uses lightweight lookup dictionary where:
    - Key: frozenset of 3mer AA pairs
    - Value: dict of compatible plasmids
    
    Returns:
        bool: True if at least one valid combination exists
        int: Number of valid 3mer pairs found
    """
    valid_pairs = 0
    
    # Convert three_mer_set to list for iteration
    three_mer_list = list(three_mer_set)
    n = len(three_mer_list)
    
    # Check all pairs of 3mers in the sequence
    for i in range(n):
        for j in range(i, n):  # Only need to check each pair once due to frozenset
            # Create key
            key = frozenset([three_mer_list[i], three_mer_list[j]])
            
            # Check if this pair exists in lookup and is compatible with plasmid
            if key in lookup_dict:
                plasmid_dict = lookup_dict[key]
                if plasmid in plasmid_dict:
                    valid_pairs += 1
    
    return valid_pairs > 0, valid_pairs


def test_success_rate(min_length=4, max_length=60, n_trials=1000, output_dir='./output'):
    """
    Test HURDLER success rate across different sequence lengths.
    
    Parameters:
    -----------
    min_length : int
        Minimum sequence length to test
    max_length : int
        Maximum sequence length to test
    n_trials : int
        Number of random sequences to test per length
    output_dir : str
        Directory to save results
    """
    print("="*80)
    print("HURDLER SUCCESS RATE TEST")
    print("="*80)
    
    # Load lightweight lookup dictionary (much faster)
    lookup_path = Path(output_dir) / 'hurdler_3mer_lightweight_lookup.pkl'
    if not lookup_path.exists():
        print(f"\nError: {lookup_path} not found.")
        print("Please run create_hurdler_lookup.py first to generate the lookup dictionary.")
        sys.exit(1)
    
    print(f"\nLoading lightweight lookup dictionary from {lookup_path}...")
    import pickle
    with open(lookup_path, 'rb') as f:
        lookup_dict = pickle.load(f)
    print(f"Loaded lookup with {len(lookup_dict):,} unique 3mer AA pairs")
    
    # Get list of all plasmids from the lookup
    all_plasmids = set()
    for plasmid_dict in lookup_dict.values():
        all_plasmids.update(plasmid_dict.keys())
    plasmids = sorted(all_plasmids)
    print(f"\nTesting {len(plasmids)} plasmids: {', '.join(plasmids)}")
    
    # Initialize results storage
    results = []
    
    # Test each length
    print(f"\nTesting sequence lengths {min_length} to {max_length} ({n_trials} trials each)...")
    
    for length in tqdm(range(min_length, max_length + 1), desc="Length"):
        for trial in range(n_trials):
            # Generate random sequence
            sequence = generate_random_aa_sequence(length)
            
            # Extract 3mers
            three_mers = extract_unique_3mers(sequence)
            n_unique_3mers = len(three_mers)
            
            # Test each plasmid
            for plasmid in plasmids:
                feasible, n_pairs = check_hurdler_feasible(three_mers, plasmid, lookup_dict)
                
                results.append({
                    'length': length,
                    'trial': trial,
                    'plasmid': plasmid,
                    'n_unique_3mers': n_unique_3mers,
                    'feasible': feasible,
                    'n_valid_pairs': n_pairs
                })
    
    # Convert to DataFrame
    df_results = pd.DataFrame(results)
    
    # Save raw results
    results_file = Path(output_dir) / 'hurdler_success_rate_raw_results.csv'
    df_results.to_csv(results_file, index=False)
    print(f"\nRaw results saved to: {results_file}")
    
    # Calculate success rates
    df_summary = df_results.groupby(['length', 'plasmid']).agg({
        'feasible': 'mean',  # Success rate
        'n_unique_3mers': 'mean',  # Average number of 3mers
        'n_valid_pairs': 'mean'  # Average number of valid pairs
    }).reset_index()
    
    df_summary.columns = ['length', 'plasmid', 'success_rate', 'avg_n_3mers', 'avg_n_pairs']
    
    # Save summary
    summary_file = Path(output_dir) / 'hurdler_success_rate_summary.csv'
    df_summary.to_csv(summary_file, index=False)
    print(f"Summary saved to: {summary_file}")
    
    return df_summary, df_results


def plot_success_rates(df_summary, output_dir='./output'):
    """
    Plot success rates for each plasmid as a function of sequence length.
    
    Parameters:
    -----------
    df_summary : DataFrame
        Summary DataFrame with columns: length, plasmid, success_rate
    output_dir : str
        Directory to save plots
    """
    print("\nGenerating plots...")
    
    # Set style
    sns.set_style('whitegrid')
    plt.rcParams['figure.figsize'] = (12, 8)
    
    # Create main plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Get unique plasmids
    plasmids = df_summary['plasmid'].unique()
    
    # Color palette
    colors = sns.color_palette('husl', len(plasmids))
    
    # Plot each plasmid
    for idx, plasmid in enumerate(sorted(plasmids)):
        df_plasmid = df_summary[df_summary['plasmid'] == plasmid].sort_values('length')
        
        ax.plot(df_plasmid['length'], 
                df_plasmid['success_rate'] * 100,  # Convert to percentage
                marker='o', 
                linewidth=2, 
                markersize=4,
                label=plasmid,
                color=colors[idx],
                alpha=0.8)
    
    # Formatting
    ax.set_xlabel('Sequence Length (amino acids)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Success Rate (%)', fontsize=14, fontweight='bold')
    ax.set_title('HURDLER Success Rate vs Sequence Length', fontsize=16, fontweight='bold')
    ax.legend(loc='best', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(df_summary['length'].min() - 1, df_summary['length'].max() + 1)
    ax.set_ylim(-5, 105)
    
    # Add reference lines
    ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=90, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = Path(output_dir) / 'hurdler_success_rate_by_length.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_file}")
    
    plt.show()
    plt.close()
    
    # Create additional plot: success rate distribution
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Grouped bar chart for selected lengths
    selected_lengths = [10, 20, 30, 40, 50, 60]
    df_selected = df_summary[df_summary['length'].isin(selected_lengths)]
    
    pivot_data = df_selected.pivot(index='plasmid', columns='length', values='success_rate')
    pivot_data = pivot_data * 100  # Convert to percentage
    pivot_data.plot(kind='bar', ax=ax1, width=0.8)
    
    ax1.set_xlabel('Plasmid', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Success Rate (%)', fontsize=12, fontweight='bold')
    ax1.set_title('Success Rate by Plasmid at Selected Lengths', fontsize=14, fontweight='bold')
    ax1.legend(title='Length', fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    
    # Plot 2: Heatmap
    pivot_data_full = df_summary.pivot(index='plasmid', columns='length', values='success_rate')
    pivot_data_full = pivot_data_full * 100  # Convert to percentage
    
    sns.heatmap(pivot_data_full, 
                annot=False, 
                fmt='.1f', 
                cmap='RdYlGn', 
                ax=ax2,
                cbar_kws={'label': 'Success Rate (%)'},
                vmin=0,
                vmax=100)
    
    ax2.set_xlabel('Sequence Length (amino acids)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Plasmid', fontsize=12, fontweight='bold')
    ax2.set_title('Success Rate Heatmap', fontsize=14, fontweight='bold')
    
    # Show only every 5th length label
    xticks = ax2.get_xticks()
    xticklabels = [int(df_summary['length'].unique()[int(i)]) if i % 5 == 0 and int(i) < len(df_summary['length'].unique()) else '' 
                   for i in xticks]
    ax2.set_xticklabels(xticklabels)
    
    plt.tight_layout()
    
    # Save plot
    plot_file2 = Path(output_dir) / 'hurdler_success_rate_detailed.png'
    plt.savefig(plot_file2, dpi=300, bbox_inches='tight')
    print(f"Detailed plot saved to: {plot_file2}")
    
    plt.show()
    plt.close()


def print_summary_statistics(df_summary, df_results):
    """Print summary statistics"""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    # Overall statistics
    print("\n1. OVERALL SUCCESS RATES")
    for plasmid in sorted(df_summary['plasmid'].unique()):
        df_plasmid = df_summary[df_summary['plasmid'] == plasmid]
        avg_success = df_plasmid['success_rate'].mean() * 100
        min_success = df_plasmid['success_rate'].min() * 100
        max_success = df_plasmid['success_rate'].max() * 100
        
        print(f"\n  {plasmid}:")
        print(f"    Average success rate: {avg_success:.1f}%")
        print(f"    Range: {min_success:.1f}% - {max_success:.1f}%")
    
    # Length-based analysis
    print("\n2. SUCCESS RATE BY LENGTH RANGE")
    length_ranges = [
        (4, 10, "Very Short (4-10)"),
        (11, 20, "Short (11-20)"),
        (21, 30, "Medium (21-30)"),
        (31, 40, "Long (31-40)"),
        (41, 60, "Very Long (41-60)")
    ]
    
    for min_len, max_len, label in length_ranges:
        df_range = df_results[
            (df_results['length'] >= min_len) & 
            (df_results['length'] <= max_len)
        ]
        
        print(f"\n  {label}:")
        for plasmid in sorted(df_summary['plasmid'].unique()):
            df_plasmid = df_range[df_range['plasmid'] == plasmid]
            success_rate = df_plasmid['feasible'].mean() * 100
            print(f"    {plasmid:<25}: {success_rate:>6.1f}%")
    
    # Find critical lengths (where success rate crosses thresholds)
    print("\n3. CRITICAL LENGTHS (50% and 90% success rates)")
    
    thresholds = [0.5, 0.9]
    for threshold in thresholds:
        print(f"\n  {threshold*100:.0f}% Success Rate Achieved At:")
        for plasmid in sorted(df_summary['plasmid'].unique()):
            df_plasmid = df_summary[df_summary['plasmid'] == plasmid].sort_values('length')
            
            # Find first length where success rate >= threshold
            above_threshold = df_plasmid[df_plasmid['success_rate'] >= threshold]
            
            if len(above_threshold) > 0:
                critical_length = above_threshold.iloc[0]['length']
                print(f"    {plasmid:<25}: Length >= {critical_length}")
            else:
                print(f"    {plasmid:<25}: Not achieved in tested range")
    
    # 3mer coverage analysis
    print("\n4. AVERAGE NUMBER OF 3MERS BY LENGTH")
    for length in [10, 20, 30, 40, 50, 60]:
        df_length = df_results[df_results['length'] == length]
        avg_3mers = df_length['n_unique_3mers'].mean()
        print(f"  Length {length:>2}: {avg_3mers:>6.1f} unique 3mers on average")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Test HURDLER success rate across sequence lengths'
    )
    parser.add_argument('--min-length', type=int, default=4,
                        help='Minimum sequence length (default: 4)')
    parser.add_argument('--max-length', type=int, default=60,
                        help='Maximum sequence length (default: 60)')
    parser.add_argument('--n-trials', type=int, default=1000,
                        help='Number of trials per length (default: 1000)')
    parser.add_argument('--output-dir', type=str, default='./output',
                        help='Output directory (default: ./output)')
    parser.add_argument('--skip-test', action='store_true',
                        help='Skip testing and only plot existing results')
    
    args = parser.parse_args()
    
    # Create output directory if needed
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    if args.skip_test:
        # Load existing results
        summary_file = Path(args.output_dir) / 'hurdler_success_rate_summary.csv'
        results_file = Path(args.output_dir) / 'hurdler_success_rate_raw_results.csv'
        
        if not summary_file.exists() or not results_file.exists():
            print("Error: Results files not found. Run without --skip-test first.")
            sys.exit(1)
        
        print("Loading existing results...")
        df_summary = pd.read_csv(summary_file)
        df_results = pd.read_csv(results_file)
    else:
        # Run tests
        df_summary, df_results = test_success_rate(
            min_length=args.min_length,
            max_length=args.max_length,
            n_trials=args.n_trials,
            output_dir=args.output_dir
        )
    
    # Generate plots
    plot_success_rates(df_summary, output_dir=args.output_dir)
    
    # Print statistics
    print_summary_statistics(df_summary, df_results)
    
    print("\n" + "="*80)
    print("Test complete!")
    print("="*80)


if __name__ == '__main__':
    main()
