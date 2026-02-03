#!/usr/bin/env python3
"""
Demonstration: Old vs New Success Rate Algorithms

This script shows the concrete differences between the two implementations
with actual examples.
"""

import numpy as np

# Example sequence
np.random.seed(42)
aa_alphabet = 'ACDEFGHIKLMNPQRSTVWY'
module = ''.join(np.random.choice(list(aa_alphabet), size=15))
sequence = module + module
L = 15

print("="*80)
print("SUCCESS RATE ALGORITHM COMPARISON")
print("="*80)
print(f"\nModule (L={L}): {module}")
print(f"Sequence (2x repeat, 30 AA): {sequence}")
print()

# Mapping for this test
mapping = {
    'ACD': [('MNP', 'right'), ('QRS', 'left')],
    'EFG': [('TUV', 'right')],
}

print(f"Mapping: {mapping}")
print()

# ============================================================================
# ALGORITHM 1: OLD (INCORRECT)
# ============================================================================
print("="*80)
print("ALGORITHM 1: OLD IMPLEMENTATION (INCORRECT)")
print("="*80)

def old_check(sequence, mapping, L):
    """Old algorithm: only finds first occurrence"""
    for ii_3mer, site_i_list in mapping.items():
        if ii_3mer in sequence:
            ii_pos = sequence.find(ii_3mer)  # ❌ Only first occurrence
            print(f"  Found Site II '{ii_3mer}' at position {ii_pos} (FIRST ONLY)")
            
            for i_3mer, direction in site_i_list:
                for i_pos in range(len(sequence) - 2):
                    if sequence[i_pos:i_pos+3] == i_3mer:
                        d = abs(ii_pos - i_pos)  # ❌ No direction info
                        if 5 < d < L:
                            print(f"    ✓ Site I '{i_3mer}' @ pos {i_pos}, "
                                  f"d=|{ii_pos}-{i_pos}|={d}, direction={direction}")
                            return True
    return False

print("\nProcessing:")
result_old = old_check(sequence, mapping, L)
print(f"\nResult: {'SUCCESS' if result_old else 'FAILURE'}")

# ============================================================================
# ALGORITHM 2: NEW (CORRECT)
# ============================================================================
print("\n" + "="*80)
print("ALGORITHM 2: NEW IMPLEMENTATION (CORRECT)")
print("="*80)

def find_3mer_positions(sequence):
    """Find all positions of each 3mer"""
    pos = {}
    n = len(sequence)
    if n < 3:
        return pos
    for i in range(n - 2):
        triplet = sequence[i:i+3]
        if triplet not in pos:
            pos[triplet] = []
        pos[triplet].append(i)
    return pos

def new_check(sequence, mapping, L):
    """New algorithm: finds all occurrences"""
    pos = find_3mer_positions(sequence)
    
    for site_ii, candidates in mapping.items():
        ii_positions = pos.get(site_ii, [])
        if not ii_positions:
            print(f"  Site II '{site_ii}' not found in sequence")
            continue
        
        print(f"  Found Site II '{site_ii}' at positions {ii_positions} (ALL)")
        
        for (site_i, direction) in candidates:
            i_positions = pos.get(site_i, [])
            if not i_positions:
                print(f"    Site I '{site_i}' not found in sequence")
                continue
            
            print(f"    Site I '{site_i}' at positions {i_positions}, direction='{direction}'")
            
            # Check all combinations
            for pii in ii_positions:
                for pi in i_positions:
                    if direction == 'right':
                        # Site I on left, Site II on right
                        d = pii - pi
                        valid = (d > 5 and d < L and pi < pii)
                        symbol = "✓" if valid else "✗"
                        if valid:
                            print(f"      {symbol} Site I @ {pi}, Site II @ {pii}: "
                                  f"d={pii}-{pi}={d}, 5<{d}<{L}✓, {pi}<{pii}✓ → SUCCESS")
                            return True
                        else:
                            reason = "d out of range" if not (d > 5 and d < L) else "pos order wrong"
                            print(f"      {symbol} Site I @ {pi}, Site II @ {pii}: "
                                  f"d={d} ({reason})")
                    else:  # 'left'
                        # Site I on right, Site II on left
                        d = pi - pii
                        valid = (d > 5 and d < L and pii < pi)
                        symbol = "✓" if valid else "✗"
                        if valid:
                            print(f"      {symbol} Site II @ {pii}, Site I @ {pi}: "
                                  f"d={pi}-{pii}={d}, 5<{d}<{L}✓, {pii}<{pi}✓ → SUCCESS")
                            return True
                        else:
                            reason = "d out of range" if not (d > 5 and d < L) else "pos order wrong"
                            print(f"      {symbol} Site II @ {pii}, Site I @ {pi}: "
                                  f"d={d} ({reason})")
    return False

print("\nProcessing:")
result_new = new_check(sequence, mapping, L)
print(f"\nResult: {'SUCCESS' if result_new else 'FAILURE'}")

# ============================================================================
# COMPARISON
# ============================================================================
print("\n" + "="*80)
print("SUMMARY & COMPARISON")
print("="*80)

print("\nKey Differences:")
print("  1. OLD: Only checks FIRST occurrence of Site II 3mer")
print("     NEW: Checks ALL occurrences")
print()
print("  2. OLD: Uses abs() → loses direction information")
print("     NEW: Uses signed distance → respects direction constraints")
print()
print("  3. OLD: Ignores the 'direction' field in mapping")
print("     NEW: Strictly enforces direction constraints")
print()
print("  4. OLD: Lower success rate (missing valid combinations)")
print("     NEW: Accurate success rate (finds all valid combinations)")

print(f"\nResults for this sequence:")
print(f"  OLD Algorithm: {'SUCCESS ✓' if result_old else 'FAILURE ✗'}")
print(f"  NEW Algorithm: {'SUCCESS ✓' if result_new else 'FAILURE ✗'}")
print()

if result_old != result_new:
    print("⚠️  ALGORITHMS DISAGREE!")
    print("    This demonstrates the bug in the old algorithm.")
else:
    print("✓ Algorithms agree on this example.")
    print("  (Both may miss patterns due to limited test sequence)")

print("\n" + "="*80)
