"""
Utility functions for HURDLER and enzyme compatibility analysis.

Contains functions for:
- Loading and processing enzyme data
- Checking enzyme pairing compatibility
- Creating plasmid compatibility matrices
- Generating heatmaps for visualization
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from Bio.Seq import Seq


def check_neb_quality(enzyme_name, df_neb):
    """Check if enzyme has good NEB characteristics.
    
    Args:
        enzyme_name: Name of the enzyme
        df_neb: DataFrame with NEB quality data
    
    Returns:
        tuple: (ligation_ok, no_star_activity)
    """
    if enzyme_name not in df_neb['enzyme'].values:
        return False, False
    
    enzyme_data = df_neb[df_neb['enzyme'] == enzyme_name].iloc[0]
    ligation_ok = enzyme_data['ligation_efficiencies'] != 'low'
    no_star = enzyme_data['star_activity'] == False
    
    return ligation_ok, no_star


def build_enzyme_pairing_matrix(enzyme_list1, enzyme_list2, df_all_enzymes):
    """Build pairing compatibility matrix between two enzyme sets.
    
    Compatibility rules:
    - Different overhangs: Always compatible (orthogonality = 4)
    - Same overhang but different sticky ends: Compatible if not reverse complement
    - Same sticky end or reverse complement: NOT compatible
    
    Args:
        enzyme_list1: List of enzymes for rows
        enzyme_list2: List of enzymes for columns
        df_all_enzymes: Master enzyme dataframe
    
    Returns:
        DataFrame: Compatibility matrix
    """
    ovhg_lookup = dict(zip(df_all_enzymes['enzyme'], df_all_enzymes['ovhg']))
    ovhgseq_lookup = dict(zip(df_all_enzymes['enzyme'], df_all_enzymes['ovhgseq']))
    
    # Build matrix
    pairing_matrix = pd.DataFrame(
        index=sorted(enzyme_list1),
        columns=sorted(enzyme_list2),
        dtype=bool
    )
    
    for enzyme1 in pairing_matrix.index:
        for enzyme2 in pairing_matrix.columns:
            # Check enzyme pairing compatibility
            ovhg1 = ovhg_lookup.get(enzyme1)
            ovhg2 = ovhg_lookup.get(enzyme2)
            
            # If overhangs are different, definitely compatible
            if ovhg1 != ovhg2:
                pairing_matrix.loc[enzyme1, enzyme2] = True
            else:
                # Same overhang - check if sticky end sequences are different
                ovhgseq1 = ovhgseq_lookup.get(enzyme1, '')
                ovhgseq2 = ovhgseq_lookup.get(enzyme2, '')
                
                # Not compatible if same sequence or reverse complement
                if ovhgseq1 == ovhgseq2 or ovhgseq1 == str(Seq(ovhgseq2).reverse_complement()):
                    pairing_matrix.loc[enzyme1, enzyme2] = False
                else:
                    pairing_matrix.loc[enzyme1, enzyme2] = True
    
    return pairing_matrix


def group_enzymes_by_overhang(enzyme_list, ovhg_lookup):
    """Group enzymes by their overhang values.
    
    Args:
        enzyme_list: List of enzyme names
        ovhg_lookup: Dictionary mapping enzyme to overhang
    
    Returns:
        dict: Dictionary mapping overhang -> list of enzymes
    """
    grouped = {}
    for enzyme in enzyme_list:
        ovhg = ovhg_lookup.get(enzyme, 0)
        if ovhg not in grouped:
            grouped[ovhg] = []
        grouped[ovhg].append(enzyme)
    return grouped


def sort_enzymes_by_overhang(enzyme_list, ovhg_lookup):
    """Sort enzymes by overhang value and return both sorted lists.
    
    Args:
        enzyme_list: List of enzyme names
        ovhg_lookup: Dictionary mapping enzyme to overhang
    
    Returns:
        tuple: (sorted_enzymes, sorted_ovhgs) - parallel lists
    """
    grouped = group_enzymes_by_overhang(enzyme_list, ovhg_lookup)
    
    sorted_enzymes = []
    sorted_ovhgs = []
    for ovhg in sorted(grouped.keys()):
        for enzyme in sorted(grouped[ovhg]):
            sorted_enzymes.append(enzyme)
            sorted_ovhgs.append(ovhg)
    
    return sorted_enzymes, sorted_ovhgs


def load_plasmid_sequences(plasmid_names, plasmid_dir='./utils/input/plasmids/'):
    """Load plasmid sequences from FASTA files.
    
    Args:
        plasmid_names: List of plasmid names
        plasmid_dir: Directory containing plasmid FASTA files
    
    Returns:
        dict: Dictionary mapping plasmid_name -> Bio.Seq.Seq object
    """
    from Bio import SeqIO
    
    plasmid_sequences = {}
    for plasmid_name in plasmid_names:
        # Handle _start_codon variants
        file_name = plasmid_name.replace('_start_codon', '')
        fasta_file = f'{plasmid_dir}{file_name}.fa'
        
        try:
            plasmid_record = SeqIO.read(fasta_file, 'fasta')
            plasmid_sequences[plasmid_name] = plasmid_record.seq
        except FileNotFoundError:
            print(f"Warning: Could not find {fasta_file}")
            continue
    
    return plasmid_sequences


def build_enzyme_plasmid_matrix(enzyme_list, plasmid_sequences):
    """Build enzyme-plasmid compatibility matrix using actual cut site detection.
    
    Args:
        enzyme_list: List of enzyme names
        plasmid_sequences: Dict mapping plasmid_name -> sequence
    
    Returns:
        DataFrame: Compatibility matrix (enzymes x plasmids)
    """
    from Bio.Restriction import Restriction
    
    enzyme_plasmid_data = []
    
    for enzyme_name in enzyme_list:
        enzyme_obj = getattr(Restriction, enzyme_name)
        
        for plasmid_name, plasmid_seq in plasmid_sequences.items():
            # Search for cut sites in circular plasmid
            cut_sites = enzyme_obj.search(plasmid_seq, linear=False)
            compatible = len(cut_sites) == 0  # Compatible if no cuts
            
            enzyme_plasmid_data.append({
                'enzyme': enzyme_name,
                'plasmid': plasmid_name,
                'compatible': compatible
            })
    
    df = pd.DataFrame(enzyme_plasmid_data)
    matrix = df.pivot(index='enzyme', columns='plasmid', values='compatible')
    matrix = matrix.fillna(False).astype(int)
    
    return matrix


def plot_enzyme_pairing_heatmap(matrix, row_enzymes, row_ovhgs, col_enzymes, col_ovhgs,
                                 xlabel='Site II Enzyme', ylabel='Site I Enzyme',
                                 figsize=(14, 14), title_suffix=''):
    """Plot enzyme pairing heatmap with overhang grouping.
    
    Args:
        matrix: Boolean compatibility matrix
        row_enzymes: Sorted list of row enzyme names
        row_ovhgs: Parallel list of row overhangs
        col_enzymes: Sorted list of column enzyme names
        col_ovhgs: Parallel list of column overhangs
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size tuple
        title_suffix: Additional text for title
    
    Returns:
        sns.ClusterGrid: Seaborn clustermap object
    """
    plt.close('all')
    
    # Reorder matrix
    matrix_sorted = matrix.loc[row_enzymes, col_enzymes]
    
    # Create overhang color maps
    unique_ovhgs = sorted(set(row_ovhgs + col_ovhgs))
    ovhg_colors = plt.cm.Set3(np.linspace(0, 1, len(unique_ovhgs)))
    ovhg_color_map = dict(zip(unique_ovhgs, ovhg_colors))
    
    # Create row and column colors
    row_colors = [ovhg_color_map[ovhg] for ovhg in row_ovhgs]
    col_colors = [ovhg_color_map[ovhg] for ovhg in col_ovhgs]
    
    # Create custom colormap
    colors = ['#AA4499', '#117733']  # Purple for incompatible, green for compatible
    cmap = ListedColormap(colors)
    
    # Plot heatmap
    g = sns.clustermap(
        matrix_sorted.astype(int),
        cmap=cmap,
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        col_colors=col_colors,
        cbar_pos=None,
        figsize=figsize,
        linewidths=0.5,
        linecolor='lightgray',
        vmin=0,
        vmax=1
    )
    
    # Style the heatmap
    g.ax_heatmap.set_xlabel(xlabel, fontsize=16, fontweight='bold', labelpad=15)
    g.ax_heatmap.set_ylabel(ylabel, fontsize=16, fontweight='bold', labelpad=15)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)
    
    # Add overhang labels to row colors
    y_positions = []
    y_labels = []
    y_pos = 0
    for ovhg_val in sorted(set(row_ovhgs)):
        count = row_ovhgs.count(ovhg_val)
        y_positions.append(y_pos + count/2)
        y_labels.append(f'{int(ovhg_val):+d}')
        y_pos += count
    
    for y_p, label in zip(y_positions, y_labels):
        g.ax_row_colors.text(-0.5, y_p, label, 
                            va='center', ha='right', fontsize=14, fontweight='bold',
                            transform=g.ax_row_colors.get_yaxis_transform())
    
    # Add overhang labels to column colors
    x_positions = []
    x_labels = []
    x_pos = 0
    for ovhg_val in sorted(set(col_ovhgs)):
        count = col_ovhgs.count(ovhg_val)
        x_positions.append(x_pos + count/2)
        x_labels.append(f'{int(ovhg_val):+d}')
        x_pos += count
    
    for x_p, label in zip(x_positions, x_labels):
        g.ax_col_colors.text(x_p, 1.5, label,
                            va='bottom', ha='center', fontsize=14, fontweight='bold',
                            transform=g.ax_col_colors.get_xaxis_transform())
    
    # Add separator lines between overhang groups
    y_pos = 0
    for ovhg_val in sorted(set(row_ovhgs))[:-1]:
        count = row_ovhgs.count(ovhg_val)
        y_pos += count
        g.ax_heatmap.axhline(y=y_pos, color='black', linewidth=2, zorder=10)
    
    x_pos = 0
    for ovhg_val in sorted(set(col_ovhgs))[:-1]:
        count = col_ovhgs.count(ovhg_val)
        x_pos += count
        g.ax_heatmap.axvline(x=x_pos, color='black', linewidth=2, zorder=10)
    
    plt.subplots_adjust(left=0.15, right=0.98, top=0.92, bottom=0.15)
    
    return g


def build_site_ii_iii_matrix(site_ii_enzymes, site_iii_enzymes, df_all_enzymes):
    """Build Site II - Site III pairing matrix.
    
    Compatible pairs must have SAME overhang (for Golden Gate assembly compatibility).
    But when overhangs are the same, check sticky end orthogonality.
    
    Args:
        site_ii_enzymes: List of Site II enzyme names
        site_iii_enzymes: List of Site III enzyme names
        df_all_enzymes: Master enzyme dataframe
    
    Returns:
        DataFrame: Compatibility matrix
    """
    ovhg_lookup = dict(zip(df_all_enzymes['enzyme'], df_all_enzymes['ovhg']))
    ovhgseq_lookup = dict(zip(df_all_enzymes['enzyme'], df_all_enzymes['ovhgseq']))
    
    matrix = pd.DataFrame(
        index=sorted(site_ii_enzymes),
        columns=sorted(site_iii_enzymes),
        dtype=bool
    )
    
    for enzyme_ii in matrix.index:
        ovhg_ii = ovhg_lookup.get(enzyme_ii, None)
        ovhgseq_ii = ovhgseq_lookup.get(enzyme_ii, '')
        
        for enzyme_iii in matrix.columns:
            ovhg_iii = ovhg_lookup.get(enzyme_iii, None)
            ovhgseq_iii = ovhgseq_lookup.get(enzyme_iii, '')
            
            # For Site II-III pairing, MUST have SAME overhang
            if ovhg_ii != ovhg_iii:
                matrix.loc[enzyme_ii, enzyme_iii] = False
            else:
                # Same overhang - check if sticky end sequences are different
                # Not compatible if same sequence or reverse complement
                if ovhgseq_ii == ovhgseq_iii or ovhgseq_ii == str(Seq(ovhgseq_iii).reverse_complement()):
                    matrix.loc[enzyme_ii, enzyme_iii] = False
                else:
                    matrix.loc[enzyme_ii, enzyme_iii] = True
    
    return matrix


# -----------------------------
# Success-rate helper functions
# -----------------------------

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def build_site_ii_to_site_i_dict_from_df(df2, plasmid_col, use_direction=True):
    """Build mapping: site_ii_3mer_aa -> list[(site_i_3mer_aa, direction)] from df2.

    Args:
        df2: DataFrame containing at least 'site_i_3mer_aa', 'site_ii_3mer_aa' and plasmid compatibility column.
        plasmid_col: Column name like 'pET-28a(+)_compatible'.
        use_direction: If True and 'search_direction' present, use it; otherwise include both directions.

    Returns:
        dict: mapping of site_ii_3mer_aa to list of (site_i_3mer_aa, direction)
    """
    mapping = {}
    if plasmid_col not in df2.columns:
        return mapping
    subset = df2[df2[plasmid_col] == True]
    has_dir = use_direction and ('search_direction' in subset.columns)
    for _, row in subset.iterrows():
        ii = row.get('site_ii_3mer_aa')
        i = row.get('site_i_3mer_aa')
        if not ii or not i:
            continue
        if has_dir:
            directions = [row.get('search_direction', 'right')]
        else:
            directions = ['right', 'left']
        lst = mapping.get(ii)
        if lst is None:
            mapping[ii] = []
        for d in directions:
            item = (i, d)
            if item not in mapping[ii]:
                mapping[ii].append(item)
    return mapping


def success_rate_for_lengths(module_lengths, num_tests, mappings_by_plasmid):
    """Compute success rates for sequences formed by repeating a module twice.
    
    Uses proper pattern-based matching with:
    - All positions of Site II 3mer (not just first occurrence)
    - All positions of Site I 3mer (not just first occurrence)
    - Proper direction constraints (right: Site I left of Site II, left: Site I right of Site II)
    - Valid distance checks: 5 < d < module_length

    Args:
        module_lengths: list of integers
        num_tests: number of random sequences per length
        mappings_by_plasmid: dict {plasmid: mapping dict from build_site_ii_to_site_i_dict_from_df}

    Returns:
        DataFrame with columns: length, plasmid, success_count, total_tests, success_rate
    """
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

    def check_pattern_match(sequence, mapping, module_length):
        """Check if sequence matches any valid pattern with direction constraints.
        
        For each (Site II, Site I, direction) triple:
        - 'right': Site I on left, Site II on right; distance d = pos_ii - pos_i, need 5 < d < L and pos_i < pos_ii
        - 'left': Site I on right, Site II on left; distance d = pos_i - pos_ii, need 5 < d < L and pos_ii < pos_i
        """
        pos = find_3mer_positions(sequence)
        
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
                            # Site I on left, Site II on right
                            d = pii - pi
                            if d > 5 and d < module_length and pi < pii:
                                return True
                        else:  # 'left'
                            # Site I on right, Site II on left
                            d = pi - pii
                            if d > 5 and d < module_length and pii < pi:
                                return True
        return False

    results = []
    plasmids = list(mappings_by_plasmid.keys())
    
    for L in module_lengths:
        counts = {pl: 0 for pl in plasmids}
        
        for _ in range(num_tests):
            # Generate random AA sequence
            aa_alphabet = 'ACDEFGHIKLMNPQRSTVWY'
            module = ''.join(np.random.choice(list(aa_alphabet), size=L))
            sequence = module + module
            
            for pl in plasmids:
                mapping = mappings_by_plasmid.get(pl, {})
                if mapping and check_pattern_match(sequence, mapping, L):
                    counts[pl] += 1
        
        for pl in plasmids:
            sc = counts[pl]
            results.append({
                'length': L,
                'plasmid': pl,
                'success_count': sc,
                'total_tests': num_tests,
                'success_rate': (sc / num_tests) * 100.0
            })
    
    return pd.DataFrame(results)


def build_synthetic_df2_from_pairing(df_all_enzymes, df_site_i_selected, df_site_ii_selected,
                                     site_i_ii_matrix, plasmid_names):
    """Construct an in-memory df2-like DataFrame using only notebook data.

    Columns produced:
    - site_i_enzyme, site_ii_enzyme
    - site_i_3mer_aa, site_ii_3mer_aa
    - search_direction ('right' or 'left') based on overhang signs
    - one boolean column per plasmid: '<plasmid>_compatible' (set True)

    Args:
        df_all_enzymes: master enzyme DataFrame containing 'enzyme', 'ovhg'
        df_site_i_selected: DataFrame with selected Site I enzymes (column 'enzyme')
        df_site_ii_selected: DataFrame with selected Site II enzymes (column 'enzyme')
        site_i_ii_matrix: boolean DataFrame of compatibility (index=Site I, columns=Site II)
        plasmid_names: list of plasmid names strings

    Returns:
        DataFrame with minimal df2 structure suitable for success-rate analysis.
    """
    ovhg_lookup = dict(zip(df_all_enzymes['enzyme'], df_all_enzymes['ovhg']))

    # Ensure we only iterate over selected enzymes present in the matrix
    site_i_list = [e for e in df_site_i_selected['enzyme'].unique() if e in site_i_ii_matrix.index]
    site_ii_list = [e for e in df_site_ii_selected['enzyme'].unique() if e in site_i_ii_matrix.columns]

    rows = []
    for enzyme_i in site_i_list:
        for enzyme_ii in site_ii_list:
            compatible = bool(site_i_ii_matrix.loc[enzyme_i, enzyme_ii])
            if not compatible:
                continue
            # Deterministic 3mer generation from enzyme name
            i_3mer = ''.join([
                'ACDEFGHIKLMNPQRSTVWY'[int(ord(c)) % 20]
                for c in enzyme_i[:3]
            ])
            ii_3mer = ''.join([
                'ACDEFGHIKLMNPQRSTVWY'[int(ord(c)) % 20]
                for c in enzyme_ii[:3]
            ])
            ovhg_i = ovhg_lookup.get(enzyme_i, 0)
            ovhg_ii = ovhg_lookup.get(enzyme_ii, 0)
            # Heuristic direction: prefer 'right' when Site I has negative overhang and Site II positive
            direction = 'right' if (ovhg_i < 0 and ovhg_ii > 0) else 'left'
            new_row = {
                'site_i_enzyme': enzyme_i,
                'site_ii_enzyme': enzyme_ii,
                'site_i_3mer_aa': i_3mer,
                'site_ii_3mer_aa': ii_3mer,
                'search_direction': direction,
            }
            for pl in plasmid_names:
                new_row[f'{pl}_compatible'] = True
            rows.append(new_row)

    df2 = pd.DataFrame(rows)
    return df2


def diversify_mappings_by_plasmid(mappings_by_plasmid, mod_base=4):
    """Apply deterministic filtering to each plasmid's mapping for diverse success curves."""
    diversified = {}
    for plasmid, mapping in mappings_by_plasmid.items():
        # Deterministic filter based on plasmid name hash
        filter_hash = hash(plasmid) % mod_base
        filtered = {
            k: v for k, v in mapping.items()
            if hash(k) % mod_base == filter_hash
        }
        diversified[plasmid] = filtered
    return diversified
