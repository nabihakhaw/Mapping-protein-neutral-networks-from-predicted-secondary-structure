"""
Figure 1A — Rank-ordered neutral component size estimates.

Each point corresponds to one strain.
Y-axis shows log10 of estimated neutral component size.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter, MaxNLocator
import os

from functions.compute_nc_properties import compute_nc_properties

# --- User Configurable Paths ---
BASE_DIR = "path/to/neighbourhood_enumeration_results"
STRAIN_IDS = "path/to/strain_ids.txt"
FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

df_nc = compute_nc_properties(
    base_dir=BASE_DIR,
    strain_id_file=STRAIN_IDS,
    S=20,
    L=566,
)

# Sort by NC size estimate
df_nc = df_nc.sort_values("Log10_Snc_est", ascending=False).reset_index(drop=True)
df_nc["Rank"] = df_nc.index + 1

# --- Plot ---
sns.set_style("white")

plt.figure(figsize=(7, 5))

plt.scatter(df_nc["Rank"], 
            df_nc["Log10_Snc_est"], 
            s=60, 
            c="pink", 
            alpha=0.9, 
            edgecolor="black", 
            linewidth=0.4
            )

plt.xscale("log")

def log_formatter(y, _):
    return r"$10^{%d}$" % int(y)

plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
plt.gca().yaxis.set_major_formatter(FuncFormatter(log_formatter))

plt.xlabel("Rank", fontsize=13)
plt.ylabel("NC size estimate", fontsize=13)

sns.despine()
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig1a", dpi=300)
plt.close()
