#!/usr/bin/env python3
"""
Quick test of HURDLER success rate system with a small sample.
"""

import sys
import random
from pathlib import Path

# Test parameters
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
TEST_LENGTHS = [10, 20, 30]
N_TRIALS = 10

def generate_random_aa_sequence(length):
    """Generate a random amino acid sequence"""
    return ''.join(random.choice(AMINO_ACIDS) for _ in range(length))

def extract_unique_3mers_circular(sequence):
    """Extract all unique 3-mer AA including circular boundaries"""
    if len(sequence) < 3:
        return set()
    
    three_mers = set()
    
    # Standard 3mers
    for i in range(len(sequence) - 2):
        three_mers.add(sequence[i:i+3])
    
    # Circular boundary 3mers
    if len(sequence) >= 3:
        three_mers.add(sequence[-2:] + sequence[0])
        three_mers.add(sequence[-1] + sequence[:2])
    
    return three_mers

def quick_test():
    """Run a quick test"""
    print("="*60)
    print("HURDLER QUICK TEST")
    print("="*60)
    
    # Check if lookup exists
    lookup_path = Path('./output/hurdler_3mer_lightweight_lookup.pkl')
    if not lookup_path.exists():
        print(f"\nWaiting for {lookup_path} to be created...")
        print("Please let create_hurdler_lookup.py finish running first.")
        return
    
    print(f"\n✓ Found lookup file: {lookup_path}")
    
    # Load lookup
    import pickle
    print("Loading lookup...")
    with open(lookup_path, 'rb') as f:
        lookup_dict = pickle.load(f)
    
    print(f"✓ Loaded {len(lookup_dict):,} 3mer AA pairs")
    
    # Get plasmids
    all_plasmids = set()
    for plasmid_dict in lookup_dict.values():
        all_plasmids.update(plasmid_dict.keys())
    plasmids = sorted(all_plasmids)
    print(f"✓ Found {len(plasmids)} plasmids: {', '.join(plasmids[:3])}...")
    
    # Run quick test
    print(f"\nRunning quick test with {N_TRIALS} sequences per length...")
    
    for length in TEST_LENGTHS:
        success_count = {p: 0 for p in plasmids}
        
        for _ in range(N_TRIALS):
            seq = generate_random_aa_sequence(length)
            three_mers = extract_unique_3mers_circular(seq)
            
            # Check each plasmid
            for plasmid in plasmids:
                # Check if any 3mer pair is valid for this plasmid
                three_mer_list = list(three_mers)
                found = False
                for i in range(len(three_mer_list)):
                    if found:
                        break
                    for j in range(i, len(three_mer_list)):
                        key = frozenset([three_mer_list[i], three_mer_list[j]])
                        if key in lookup_dict:
                            if plasmid in lookup_dict[key]:
                                success_count[plasmid] += 1
                                found = True
                                break
        
        # Print results
        print(f"\nLength {length}:")
        for plasmid in plasmids:
            rate = success_count[plasmid] / N_TRIALS * 100
            print(f"  {plasmid:<30}: {rate:>5.1f}% ({success_count[plasmid]}/{N_TRIALS})")
    
    print("\n" + "="*60)
    print("Quick test complete!")
    print("System is ready for full success rate testing.")
    print("="*60)

if __name__ == '__main__':
    quick_test()
