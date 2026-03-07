import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from functions.neighbourhood_enu_df import (
    parse_neighbourhood_enum,
    compute_avg_neighbourhood_properties,
)
from functions.labels import get_labels 

# --- User Configurable Paths ---
BASE_DIR = "/path/to/neighbourhood_enu_results"
STRAIN_IDS = "/path/to/strain_ids.txt"
LABELS = "/data/final_df_extra.pkl"

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

# -------------------------
# Summary statistics
# -------------------------
robustness_stats = (
    df_avg.groupby("strain")["robustness"]
    .agg(["mean", "std"])
    .reset_index()
)

robustness_stats.columns = ["strain", "mean_robustness", "std_robustness"]

robustness_sorted = robustness_stats.sort_values(
    "mean_robustness", ascending=False
).reset_index(drop=True)

y = np.arange(len(robustness_sorted))
order = robustness_sorted["strain"]

# -------------------------
# Labels
# -------------------------
labels = get_labels(LABELS, STRAIN_IDS)
strain_labels = [labels.get(s, s) for s in order]

# -------------------------
# Error bars (non-negative)
# -------------------------
lower_err = np.minimum(
    robustness_sorted["std_robustness"],
    robustness_sorted["mean_robustness"],
)
upper_err = robustness_sorted["std_robustness"]
xerr = np.array([lower_err, upper_err])

# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(6, 7))
sns.set_style("white")

palette = sns.cubehelix_palette(start=0.0, rot=1.5, n_colors=len(order))

sns.violinplot(
    data=df_avg[df_avg["strain"].isin(order)],
    y="strain",
    x="robustness",
    order=order,
    orient="h",
    inner=None,
    cut=0,
    linewidth=0,
    palette=palette,
)

plt.errorbar(
    x=robustness_sorted["mean_robustness"],
    y=y,
    xerr=xerr,
    fmt="o",
    ecolor="grey",
    elinewidth=0.2,
    capsize=0,
    markersize=5,
    color="black",
    linestyle="none",
)

plt.yticks(y, strain_labels)
plt.xlabel("Robustness", fontsize=13)
plt.ylabel("")

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig2b.png"), dpi=300) 
plt.close()
