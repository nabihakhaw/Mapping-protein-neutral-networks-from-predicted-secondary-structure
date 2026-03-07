from collections import Counter
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from functions.neighbourhood_enu_df import parse_neighbourhood_enum
from functions.network_construction import build_graphs_from_neighbourhood
from functions.get_labels import get_labels

# =========================
# User-configurable paths
# =========================
BASE_DIR = "/path/to/neighbourhood_enu_results"
SEED_DIR = "/path/to/subsampled"          # contains {strain}_20.txt
STRAIN_IDS = "/path/to/strain_ids.txt"
FIG_DIR = "figures"

# =========================
# Global parameters
# =========================
S = 20
L = 566
ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")

# =========================
# Load neighbourhood enumeration data
# =========================
df_long = parse_neighbourhood_enum(
    base_dir=BASE_DIR,
    strain_id_file=STRAIN_IDS,
    S=S,
    L=L,
)

# =========================
# Build graphs
# =========================
keys = ["AAD17229", "ACP41934", "ABD95350", "ABE11878", "ADT78908"]
strain_dirs = [f"{key}_20" for key in keys]

seq_dfs, graphs = build_graphs_from_neighbourhood(
    df_long=df_long,
    strain_dirs=strain_dirs,
    seed_dir=SEED_DIR,
)

# =========================
# Degree distributions
# =========================
degree_distributions = {}

for strain, G in graphs.items():
    if G.number_of_nodes() == 0:
        continue

    degrees = [deg for _, deg in G.degree()]
    deg_counts = Counter(degrees)

    ks = np.array(sorted(deg_counts.keys()))
    counts = np.array([deg_counts[k] for k in ks])

    degree_distributions[strain] = (ks, counts)

# =========================
# Plot
# =========================
labels = get_labels(STRAIN_IDS, keys)
colors = sns.cubehelix_palette(
    start=0.0, rot=1.5, n_colors=len(degree_distributions)
)

plt.figure(figsize=(10, 6))

for (strain, (ks, counts)), color in zip(degree_distributions.items(), colors):
    plt.plot(
        ks, counts,
        marker="o",
        linestyle="-",
        color=color,
        label=labels.get(strain, strain),
    )

plt.xlabel("Degree", fontsize=14)
plt.ylabel("Count", fontsize=14)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.xlim(0, 30)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig3a.png", dpi=300)
plt.close()
