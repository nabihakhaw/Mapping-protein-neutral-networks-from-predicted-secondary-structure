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
LABELS = "/data/final_df_extra.pkl"
FIG_DIR = "figures"

# =========================
# Global parameters
# =========================
S = 20
L = 566

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
# Plot CCDF of seed degrees
# =========================
colors = sns.cubehelix_palette(start=0.0, rot=1.5, n_colors=len(graphs))
labels = get_labels(LABELS, STRAIN_IDS)

plt.figure(figsize=(10, 6))

for (strain, G), color in zip(graphs.items(), colors):

    seq_df = seq_dfs[strain]
    seed_sequences = seq_df["seed_seq"].unique()

    # filter only seeds that exist in the graph
    seed_nodes_in_graph = [s for s in seed_sequences if s in G]
    if not seed_nodes_in_graph:
        print(f"No seed nodes found in graph for {strain}")
        continue

    # get degree of only seed nodes
    seed_degrees = [G.degree(n) for n in seed_nodes_in_graph]

    # build degree distribution
    deg_counts = Counter(seed_degrees)
    ks = np.array(sorted(deg_counts.keys()))
    counts = np.array([deg_counts[k] for k in ks])

    # compute CCDF
    ccdf = np.cumsum(counts[::-1])[::-1] / np.sum(counts)

    plt.plot(
        ks, ccdf,
        marker='o', linestyle='-',
        color=color,
        label=labels.get(strain, strain)
    )

plt.xlabel("Degree", fontsize=14)
plt.ylabel("CCDF", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(title="Strain", fontsize=10, frameon=False)
plt.grid(False)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig3b.png", dpi=300)
plt.close()
