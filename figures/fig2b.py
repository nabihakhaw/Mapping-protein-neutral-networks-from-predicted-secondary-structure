"""
Figure 1B — Neutral component size versus robustness.

Each point corresponds to one strain.
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

# --- Plot ---
sns.set_style("white")

plt.figure(figsize=(7, 5))

plt.scatter(
    df_nc["Log10_Snc_est"], 
    df_nc["Robustness"], 
    s=60,
    c="pink", 
    alpha=0.9, 
    edgecolor="black", 
    linewidth=0.4
    )

def log_formatter(x, _):
    return r"$10^{%d}$" % int(x)

plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.gca().xaxis.set_major_formatter(FuncFormatter(log_formatter))

plt.xlabel("NC size estimate", fontsize=13)
plt.ylabel("NC robustness estimate", fontsize=13)

sns.despine()
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig1b.png", dpi=300)
plt.close()
