import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from functions.neighbourhood_enu_df import (
    parse_neighbourhood_enum,
    compute_avg_neighbourhood_properties,
)

"""
Figure 6 - Robustness vs secondary structure elements
"""

# --- User Configurable Paths ---
dfss = pd.read_pickle("/data/final_df_extra.pkl")
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

merged_list = []

# Load strain IDs
with open(STRAIN_IDS) as f:
    strain_ids = {l.strip().split("\t")[0] for l in f if l.strip()}

for strain in strain_ids:
    SSdf = (
        dfss.loc[dfss["strain_id"] == strain, ["index", "structure"]]
        .rename(columns={"index": "position"})
        .copy()
    )

    SSdf["position"] = SSdf["position"].astype(int)
    SSdf["position"] -= 1

    merged = pd.merge(
        df_avg[df_avg["strain"] == strain],
        SSdf,
        on="position",
        how="inner",
    )

    merged_list.append(merged)

merged_all = pd.concat(merged_list, ignore_index=True)

# -------------------------
# Plot
# -------------------------
order = ["H", "C", "E"]
palette = sns.cubehelix_palette(rot=1.5, light=0.55, gamma=0.6, n_colors=len(order))

plt.figure(figsize=(8, 5))
sns.violinplot(
    data=merged_all,
    x="structure",
    y="robustness",
    order=order,
    palette=palette,
    linewidth=0.0,
)

plt.xlabel("Structural element", fontsize=13)
plt.ylabel("Robustness", fontsize=13)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig6.png", dpi=300)
plt.close()
