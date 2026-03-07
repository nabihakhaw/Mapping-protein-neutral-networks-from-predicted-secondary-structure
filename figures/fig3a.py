import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# --- Paths ---
STRAIN_IDS = "/path/to/strain_ids.txt"
RESULTS_DIR = "/path/to/site_scanning_results"
FIG_DIR = "figures"
STRAINS_TO_PLOT = {
    "ADT78908.txt",
    "ACP41934.txt",
    "ABO38054.txt",
    "ACQ73203.txt"
}

# -------------------------
# Load + process data
# -------------------------
df = parse_site_scanning(
    base_dir=RESULTS_DIR,
    strain_id_file=STRAIN_IDS,
)

# keep only first 2 full loops (as before)
df = select_first_n_loops(df, n_loops=2)

df = df[df["filename"].isin(STRAINS_TO_PLOT)]

# cumulative quantities per strain
df = df.sort_values(["filename", "position"])
df["cum_attempts"] = df.groupby("filename")["attempts"].cumsum()
df["cum_neutrals"] = df.groupby("filename")["is_neutral"].cumsum()

# Prepare color map
colors = sns.cubehelix_palette(start=0.0, rot=1.5, n_colors=len(filenames_to_plot))

# -------------------------
# Labels
# -------------------------
labels = get_labels(LABELS, STRAIN_IDS)

def pretty_label(filename):
    strain = filename.replace(".txt", "")
    return labels.get(strain, strain)

# -------------------------
# Labels (in case they are malformed)
# -------------------------
# end_labels = {
#     "ADT78908.txt": "Netherlands/1953",
#     "ACP41934.txt": "Mexico/2009",
#     "ABO38054.txt": "Marton/1943",
#     "ACQ73203.txt": "Ohio/2007"
# }

# def pretty_label(filename):
#     return end_labels.get(filename, filename.replace(".txt", ""))


# -------------------------
# Plot
# -------------------------
sns.set_style("white")

filenames = df["filename"].unique()
colors = sns.cubehelix_palette(
    start=0.0,
    rot=1.5,
    n_colors=len(filenames),
)

plt.figure(figsize=(8, 5))

for color, fname in zip(colors, filenames):
    sub = df[df["filename"] == fname]

    plt.step(
        sub["cum_attempts"],
        sub["cum_neutrals"],
        where="post",
        linewidth=2.5,
        color=color,
    )

    # label at trajectory end
    plt.text(
        sub["cum_attempts"].iloc[-1] * 1.02,
        sub["cum_neutrals"].iloc[-1],
        pretty_label(fname),
        fontsize=11,
        va="center",
        color="black",
    )

# -------------------------
# Formatting
# -------------------------
plt.xlabel("Cumulative attempts", fontsize=13)
plt.ylabel("Cumulative neutral mutants", fontsize=13)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "fig2a.png"), dpi=300)
plt.close()
