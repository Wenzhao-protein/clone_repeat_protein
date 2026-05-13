#!/usr/bin/env python3
"""
Test script to verify the plot fix
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Create sample data similar to what the notebook generates
np.random.seed(42)
plasmids = ['pGEX-4T-1', 'pMAL-c5X', 'pET-21a(+)', 'pET-28a(+)', 'pET-28a(+)_start_codon', 'pCold_I', 'pUC18', 'pQE-3']
module_lengths = list(range(7, 61))

# Generate synthetic success rates (increasing trend)
data_list = []
for plasmid in plasmids:
    for length in module_lengths:
        # Success rate increases with length (realistic pattern)
        success_rate = 0.1 + (length - 7) * 0.03 + np.random.normal(0, 0.02)
        success_rate = max(0, min(success_rate, 3))  # Clamp to [0, 3]
        data_list.append({'plasmid': plasmid, 'length': length, 'success_rate': success_rate})

results_df = pd.DataFrame(data_list)

# Plot with improved formatting (as per the fix)
sns.set_style('whitegrid')
fig, ax = plt.subplots(figsize=(14, 8))

# Define distinct colors for plasmids
colors = plt.cm.tab10(range(len(results_df['plasmid'].unique())))

for i, plasmid in enumerate(sorted(results_df['plasmid'].unique())):
    data = results_df[results_df['plasmid'] == plasmid].sort_values('length')
    ax.plot(data['length'], data['success_rate'], linewidth=2.5, label=plasmid,
            color=colors[i % len(colors)], alpha=0.85)

ax.set_xlabel('Module Length (amino acids)', fontsize=12, fontweight='bold')
ax.set_ylabel('Success Rate (%)', fontsize=12, fontweight='bold')
ax.set_title('HURDLER Success Rate vs Module Length\n(1000 tests per length, pattern-based repeated modules)',
             fontsize=13, fontweight='bold', pad=20)

# Improve grid appearance
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
ax.set_axisbelow(True)

# Place legend outside, compact
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10, framealpha=0.95)

# Set proper axis limits
ax.set_xlim(5, 62)
ax.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('/tmp/test_plot_improved.png', dpi=100, bbox_inches='tight')
print("✓ Improved plot saved to /tmp/test_plot_improved.png")
print("\nImprovements made:")
print("  1. Removed markers (o) to reduce visual clutter")
print("  2. Increased line width to 2.5 for better visibility")
print("  3. Used distinct colors from tab10 colormap")
print("  4. Improved grid appearance (dashed, lighter, behind plot)")
print("  5. Better legend placement and styling")
print("  6. Larger figure size (14x8) for more space")
print("  7. Bold axis labels and title")
print("  8. Set axis limits explicitly (x: 5-62, y: from 0)")
