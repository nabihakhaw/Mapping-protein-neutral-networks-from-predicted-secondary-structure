import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

"""
Figure 5 - Mean mutational robustness per position per strain across
3 full site scanning loops presented in a diverging heatmap.
"""

# -------------------------
# Paths
# -------------------------
BASE_DIR = "/path/to/site_scanning_results"
STRAIN_IDS = "/path/to/strain_ids.txt"
LABELS = "/data/final_df_extra.pkl"
FIG_DIR = "figures"

# -------------------------
# Load + preprocess data
# -------------------------
df = parse_site_scanning(BASE_DIR, STRAIN_IDS)
df = compute_robustness(df)

# keep all full loops
df_equal = select_first_n_loops(df, n_loops=3, max_pos=565)

# -------------------------
# Mean robustness per strain × position
# -------------------------
robustness_stats = (
    df_equal
    .groupby(["filename", "position"])["robustness"]
    .mean()
    .reset_index()
)

# -------------------------
# Build strain × position matrix
# -------------------------
strains = sorted(robustness_stats["filename"].unique())
positions = np.arange(566)

position_vectors = {}

for strain in strains:
    sub = robustness_stats[robustness_stats["filename"] == strain]
    vec = (
        sub
        .set_index("position")
        .reindex(positions)["robustness"]
        .values
    )
    position_vectors[strain] = vec

matrix = np.vstack([position_vectors[s] for s in strains])

# -------------------------
# Load years & sort strains
# -------------------------
label_map = get_labels(LABELS, STRAIN_IDS)

years = {
    strain: int(label_map[strain.replace(".txt", "")].split("/")[-1])
    for strain in strains
}

order = np.argsort([years[s] for s in strains])

matrix_sorted = matrix[order, :]
strains_sorted = [strains[i] for i in order]
years_sorted = [years[s] for s in strains_sorted]

# -------------------------
# Plot heatmap + year side plot
# -------------------------
fig = plt.figure(figsize=(14, 7))
gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1], wspace=0.02)

# Heatmap
ax0 = fig.add_subplot(gs[0])
im = ax0.imshow(
    matrix_sorted,
    aspect="auto",
    interpolation="nearest",
    cmap="winter"
)
ax0.set_xlabel("Position")
ax0.set_ylabel("Strain")
ax0.set_yticks(range(len(strains_sorted)))
ax0.set_yticklabels([s.replace(".txt", "") for s in strains_sorted])

plt.colorbar(im, ax=ax0, label="Mean robustness")

# Year plot
ax1 = fig.add_subplot(gs[1], sharey=ax0)
ax1.plot(years_sorted, np.arange(len(years_sorted)), "-o", color="gray")
ax1.set_xlabel("Year")
ax1.set_yticks([])
ax1.invert_yaxis()

for spine in ax1.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig5.png", dpi=300)
plt.close()
