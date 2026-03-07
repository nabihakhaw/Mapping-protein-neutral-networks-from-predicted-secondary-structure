import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------------
# Paths & parameters
# -------------------------
BASE_ROOT = "/path/to/evolvability_results"
SEED_DIR = "/path/to/subsampled"
dfss = pd.read_pickle("/data/final_df_extra.pkl")
FIG_DIR = "figures"

KEYS = ["ABD95350", "ACP41934"]
S = 20

# -------------------------
# Load wild-type sequences
# -------------------------
wt_seq = {}

for key in KEYS:
    seed_file = os.path.join(SEED_DIR, f"{key}_20.txt")
    with open(seed_file) as f:
        seqs = [line.strip() for line in f if line.strip()]

    if not seqs:
        raise ValueError(f"No seed sequences found for {key}")

    wt_seq[key] = seqs[0]  # assume first = WT
    print(f"{key}: loaded WT sequence")

# -------------------------
# Load wild-type structures
# -------------------------
wt_struct = {}

for key in KEYS:
    ss_df = (
        dfss[dfss["strain_id"] == key]
        .iloc[:, [0, 2]]
        .sort_values("index")
    )

    wt_struct[key] = "".join(ss_df["structure"])
    print(f"{key}: loaded WT structure")

# -------------------------
# Process mutational neighbourhoods
# -------------------------
all_results = []

for key in KEYS:
    dir_path = os.path.join(BASE_ROOT, f"{key}_{S}")
    dfs = []

    for i in range(S):
        file_path = os.path.join(dir_path, f"{i}.txt")

        df = pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=["Position", "Genotype", "Phenotype"],
        )

        df["Neighbourhood"] = i
        df["Position"] = df["Position"].astype(int)

        wt_sequence = wt_seq[key]
        wt_structure = wt_struct[key]

        # amino-acid mutations
        df["Mutations"] = df["Genotype"].apply(
            lambda seq: ";".join(
                f"{wt_sequence[p]}{p+1}{aa}"
                for p, (aa_wt, aa) in enumerate(zip(wt_sequence, seq))
                if aa != aa_wt
            )
        )

        # secondary-structure changes
        df["Structure_changes"] = df["Phenotype"].apply(
            lambda ss: ";".join(
                f"{wt_structure[p]}{p+1}{s}"
                for p, (s_wt, s) in enumerate(zip(wt_structure, ss))
                if s != s_wt
            )
        )

        dfs.append(df)

    all_data = pd.concat(dfs, ignore_index=True)

    phenotype_counts = (
        all_data["Phenotype"]
        .value_counts()
        .reset_index()
        .rename(columns={"index": "Phenotype", "Phenotype": "Count"})
    )

    all_results.append(
        {
            "Key": key,
            "AllData": all_data,
            "PhenotypeCounts": phenotype_counts,
        }
    )

# -------------------------
# Per-seed evolvability metrics
# -------------------------
metrics = []

for res in all_results:
    key = res["Key"]
    data = res["AllData"]
    wt_structure = wt_struct[key]

    for seed_id, df_seed in data.groupby("Neighbourhood"):
        phenos = df_seed["Phenotype"].tolist()

        metrics.append(
            {
                "Key": key,
                "Seed": seed_id,
                "Total": len(phenos),
                "UniqueNonNeutral": len({p for p in phenos if p != wt_structure}),
            }
        )

metrics_df = pd.DataFrame(metrics)


plt.figure(figsize=(8, 6))
sns.set_style('white')

metrics_df["Key_code"] = metrics_df["Key"].astype("category").cat.codes

plt.scatter(
    metrics_df["Total"],
    metrics_df["UniqueNonNeutral"],
    c=metrics_df["Key_code"],
    cmap="PiYG",
    s=100,
    alpha=0.5,
)

plt.xlabel("Total phenotypes accessible per seed", fontsize=14)
plt.ylabel("Unique non-neutral phenotypes per seed", fontsize=14)
plt.grid(False)
sns.despine()
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig7.png", dpi = 300)
plt.close()
