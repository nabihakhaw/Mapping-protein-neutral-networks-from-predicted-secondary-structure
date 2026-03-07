"""
Figure 10a — Robustness distributions across strains calculated using site_scanning

Violins: robustness values from first two full site-scanning loops
Points/error bars: mean ± std robustness per strain

Outputs:
- figures/fig10a.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from functions.site_scanning_df import (parse_site_scanning, 
compute_robustness, select_first_n_loops)
from functions.get_labels import (get_labels)

# --- Paths --- 
BASE_DIR = "/path/to/site_scanning_results" 
STRAIN_IDS = "/path/to/strain_ids.txt"
LABELS = "/data/final_df_extra.pkl"
FIG_DIR = "figures"

df = parse_site_scanning(BASE_DIR, STRAIN_IDS)
df = compute_robustness(df)

df_equal = select_first_n_loops(df, n_loops=2, max_pos=565)

# Group by strain
robustness_stats = (
    df_equal.groupby('filename')['robustness']
    .agg(['mean', 'std'])
    .reset_index()
)
robustness_stats.columns = ['strain', 'mean_robustness', 'std_robustness']

# Sort by mean robustness descending
robustness_sorted = robustness_stats.sort_values('mean_robustness', ascending=False).reset_index(drop=True)

# Numeric y positions
y = np.arange(len(robustness_sorted))

# Asymmetric error bars so lower bound doesn't go below 0
lower_err = np.minimum(robustness_sorted['std_robustness'], robustness_sorted['mean_robustness'])
upper_err = robustness_sorted['std_robustness']
xerr = np.array([lower_err, upper_err])
order = robustness_sorted['strain']

# Get labels
labels = get_labels(LABELS, STRAIN_IDS)
strain_labels = [labels.get(s.replace(".txt", ""), s) for s in order]

plt.figure(figsize=(6, 7))

palette = sns.cubehelix_palette(start=0.0, rot=1.5, n_colors=len(order))

# Violin (horizontal density)
sns.violinplot(
    data=df_equal[df_equal['filename'].isin(order)],
    y='filename', x='robustness',
    order=order,
    orient='h',
    inner=None,           # remove internal bars
    cut=0,                # no extrapolation beyond data
    linewidth=0,
    palette=palette       # use generated palette
)

plt.errorbar(
    x=robustness_sorted['mean_robustness'],
    y=y,
    xerr=xerr,
    fmt='o',
    ecolor='grey',
    elinewidth=0.2,
    capsize=0,
    markersize=5,
    color='black',
    linestyle='none'
)


# Replace y-tick labels and remove default axis label
plt.yticks(ticks=y, labels=strain_labels)
plt.ylabel("")  # remove "filename"
plt.xlabel("Robustness", fontsize=13)

# Styling
plt.grid(True, linestyle='--', alpha=0.0)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig9a.png", dpi=300)
