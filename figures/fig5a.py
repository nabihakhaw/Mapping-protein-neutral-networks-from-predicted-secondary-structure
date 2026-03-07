import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""
Figure 4 - Mean mutational robustness per position averaged across 15 strains. 
Shaded region shows ±1 SD across strains.
"""

from functions.neighbourhood_enu_df import (
    parse_neighbourhood_enum,
    compute_avg_neighbourhood_properties,
)

# --- User Configurable Paths ---
BASE_DIR = "/path/to/neighbourhood_enu_results"
STRAIN_IDS = "/path/to/strain_ids.txt"
FIG_DIR = "figures"

# -------------------------
# Load + compute
# -------------------------
df_long = parse_neighbourhood_enum(
    base_dir=BASE_DIR,
    strain_id_file=STRAIN_IDS,
    S=20,
    L=566,
)

df_avg = compute_avg_neighbourhood_properties(df_long)

mean_robustness = df_avg.groupby("position")["robustness"].agg(["mean", "std"]).reset_index()

plt.figure(figsize=(15, 6))
plt.plot(mean_robustness["position"], mean_robustness["mean"], label="Mean Robustness", lw=1.0, color='cornflowerblue')
plt.fill_between(
    mean_robustness["position"],
    mean_robustness["mean"] - mean_robustness["std"],
    mean_robustness["mean"] + mean_robustness["std"],
    alpha=0.2,
    color='lightsteelblue',
    label="±1 SD"
)

plt.xlabel("Position", fontsize=13)
plt.ylabel("Robustness", fontsize=13)
plt.legend(frameon=False)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.grid(True, linestyle="--", alpha=0.0)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig4.png", dpi=300)
plt.close()
