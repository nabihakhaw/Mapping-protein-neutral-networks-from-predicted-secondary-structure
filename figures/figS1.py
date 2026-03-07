import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter

# -------------------------
# Paths & parameters
# -------------------------
BASE_ROOT = "/path/to/neighbourhood_enu_results/ADE75232"
FIG_DIR = "figures"

KEY = "ADE75232"

SAMPLE_SIZES = [10, 20, 50, 100]
N_REPS = 2
L = 566

EXPECTED_POSITIONS = set(range(L))

# -------------------------
# Collect estimates
# -------------------------
results = []

for S in SAMPLE_SIZES:
    for rep in range(1, N_REPS + 1):

        dir_name = f"{KEY}_{S}_{rep}"
        dir_path = os.path.join(BASE_ROOT, dir_name)

        if not os.path.isdir(dir_path):
            print(f"Skipping missing directory: {dir_path}")
            continue

        dfs = []

        for i in range(S):
            file_path = os.path.join(dir_path, f"{i}.txt")

            if not os.path.exists(file_path):
                print(f"  Missing file: {file_path}")
                continue

            df = pd.read_csv(
                file_path,
                sep="\t",
                header=None,
                names=["Position", "Count", "Neighbours"],
            )

            df["Position"] = df["Position"].astype(int)

            # sanity check positions
            missing = EXPECTED_POSITIONS - set(df["Position"])
            if missing:
                print(
                    f"  Missing {len(missing)} positions in {file_path} "
                    f"(e.g. {sorted(missing)[:5]})"
                )

            dfs.append(df)

        if not dfs:
            print(f"No valid data for {dir_name}, skipping.")
            continue

        all_data = pd.concat(dfs, ignore_index=True)

        # mean neutral count per position
        xj = all_data.groupby("Position")["Count"].mean()

        # Nc estimate
        min_terms = np.minimum(1.0 + xj, 20.0)
        log10_snc = np.sum(np.log(min_terms)) / np.log(10.0)

        results.append(
            {
                "SampleSize": S,
                "Replicate": rep,
                "Log10_Snc_est": log10_snc,
            }
        )

# -------------------------
# Results dataframe
# -------------------------
df_res = pd.DataFrame(results)

plt.figure(figsize=(6, 4))

sns.violinplot(
    data=df_res,
    x="SampleSize",
    y="Log10_Snc_est",
    inner=None,
    linewidth=0.01,
    palette="husl",
)

sns.stripplot(
    data=df_res,
    x="SampleSize",
    y="Log10_Snc_est",
    color="k",
    size=2,
    jitter=False,
    alpha=0.4,
)

# median annotations
medians = df_res.groupby("SampleSize")["Log10_Snc_est"].median()
for i, (S, med) in enumerate(medians.items()):
    plt.text(i, med + 0.1, f"{med:.2f}", ha="center", fontsize=9)

# format y-axis as 10^x
plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.gca().yaxis.set_major_formatter(
    FuncFormatter(lambda y, _: rf"$10^{{{int(y)}}}$")
)

plt.xlabel("Sample size (S)")
plt.ylabel(r"$S_{NC,\mathrm{est}}$")
sns.despine()
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig8.png", dpi = 300)
plt.close()
