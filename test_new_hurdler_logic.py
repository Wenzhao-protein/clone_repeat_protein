"""
Test script to verify the new HURDLER success rate calculation logic
"""
import pickle
from collections import defaultdict

# Load lookup dictionary
print("="*80)
print("TESTING NEW HURDLER LOGIC")
print("="*80)

with open('./output/hurdler_lookup_optimized.pkl', 'rb') as f:
    hurdler_lookup = pickle.load(f)

def extract_3mer_pairs_with_positions(sequence):
    """Extract all 3mer AA from doubled sequence with their positions"""
    doubled_seq = sequence + sequence
    
    three_mer_positions = defaultdict(list)
    for i in range(len(doubled_seq) - 2):
        three_mer = doubled_seq[i:i+3]
        three_mer_positions[three_mer].append(i)
    
    return three_mer_positions

def check_hurdler_success_detailed(sequence, hurdler_lookup, plasmid):
    """Check if sequence has valid HURDLER solution with detailed output"""
    module_length = len(sequence)
    three_mer_positions = extract_3mer_pairs_with_positions(sequence)
    
    if plasmid not in hurdler_lookup:
        return False, "Plasmid not in lookup"
    
    plasmid_lookup = hurdler_lookup[plasmid]
    three_mers = list(three_mer_positions.keys())
    
    print(f"\nSequence: {sequence}")
    print(f"Module length: {module_length}")
    print(f"Doubled sequence: {sequence + sequence}")
    print(f"Unique 3mers found: {len(three_mers)}")
    
    # Try all pairs of 3mers
    for i in range(len(three_mers)):
        for j in range(i+1, len(three_mers)):
            site_i_3mer = three_mers[i]
            site_ii_3mer = three_mers[j]
            
            # Check both orders
            for pair in [(site_i_3mer, site_ii_3mer), (site_ii_3mer, site_i_3mer)]:
                if pair not in plasmid_lookup:
                    continue
                
                enzyme_combos = plasmid_lookup[pair]
                
                print(f"\n  Testing pair: {pair}")
                print(f"  Enzyme combinations: {len(enzyme_combos)}")
                
                # Test each enzyme combination
                for enzyme_combo in enzyme_combos:
                    site_i_enzyme, site_ii_enzyme, site_iii_enzyme, search_direction = enzyme_combo
                    
                    site_i_positions = three_mer_positions[pair[0]]
                    site_ii_positions = three_mer_positions[pair[1]]
                    
                    print(f"    Enzyme combo: {site_i_enzyme} + {site_ii_enzyme} + {site_iii_enzyme}")
                    print(f"    Search direction: {search_direction}")
                    print(f"    Site I ({pair[0]}) positions: {site_i_positions}")
                    print(f"    Site II ({pair[1]}) positions: {site_ii_positions}")
                    
                    # Test if any combination of positions satisfies the distance constraint
                    for site_ii_pos in site_ii_positions:
                        if search_direction == 'left':
                            # Search for site_i to the left of site_ii
                            for site_i_pos in site_i_positions:
                                if site_i_pos < site_ii_pos:
                                    distance = site_ii_pos - site_i_pos
                                    valid = 5 < distance < module_length
                                    print(f"      Left search: site_i@{site_i_pos} < site_ii@{site_ii_pos}, distance={distance}, valid={valid}")
                                    if valid:
                                        return True, f"SUCCESS: {pair} with {enzyme_combo}"
                        else:  # search_direction == 'right'
                            # Search for site_i to the right of site_ii
                            for site_i_pos in site_i_positions:
                                if site_i_pos > site_ii_pos:
                                    distance = site_i_pos - site_ii_pos
                                    valid = 5 < distance < module_length
                                    print(f"      Right search: site_i@{site_i_pos} > site_ii@{site_ii_pos}, distance={distance}, valid={valid}")
                                    if valid:
                                        return True, f"SUCCESS: {pair} with {enzyme_combo}"
    
    return False, "No valid enzyme combination found"

# Test case 1: Short sequence that should fail
print("\n" + "="*80)
print("TEST CASE 1: Short sequence (4 AA)")
print("="*80)
test_seq_1 = "ACDE"
plasmid = 'pUC18_compatible'
success, message = check_hurdler_success_detailed(test_seq_1, hurdler_lookup, plasmid)
print(f"\nResult: {success}")
print(f"Message: {message}")

# Test case 2: Medium sequence that should succeed
print("\n" + "="*80)
print("TEST CASE 2: Medium sequence (20 AA)")
print("="*80)
test_seq_2 = "ACDEFGHIKLMNPQRSTVWY"
success, message = check_hurdler_success_detailed(test_seq_2, hurdler_lookup, plasmid)
print(f"\nResult: {success}")
print(f"Message: {message}")

print("\n" + "="*80)
print("TESTING COMPLETE")
print("="*80)
