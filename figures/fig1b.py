import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import FuncFormatter
import seaborn as sns

from functions.site_scanning_df import parse_site_scanning
from functions.compute_nc_properties import compute_nc_properties
from functions.get_labels import get_labels

# -------------------------
# Paths
# -------------------------
RESULTS_DIR = "/path/to/site_scanning_results" 
STRAIN_IDS = "/path/to/strain_ids.txt"
LABELS = "/data/final_df_extra.pkl"
FIG_DIR = "figures"

# -------------------------
# Load + process data
# -------------------------
df = parse_site_scanning(base_dir=RESULTS_DIR, strain_id_file=STRAIN_IDS)
df['strain_id'] = df['filename'].str.replace('.txt', '', regex=False)

# -------------------------
# Select strains
# -------------------------
strain_list = ["AEB21358", "ACQ73203"]

cmap = plt.get_cmap("Set2")
colors = cmap.colors  # This gives a list of 8 RGBA tuples
# Pick color 2, 3
selected_colors = [colors[2], colors[3]] 

# -------------------------
# Estimation functions
# -------------------------
def estimate_k(a, N=ATTEMPTS_PER_SITE):
    """Estimate k from number of attempts"""

    try:
        a = float(a)
    except Exception:
        return 0.0
    if a > 0:
        val = (N + 1) / a - 1.0
    else:
        val = 0.0
    return max(0.0, val)

def compute_estimates_for_strain(df, strain, L=PROTEIN_LENGTH, N=ATTEMPTS_PER_SITE):
    """Compute log10(S_nc_est) as function of number of runs"""
    
    sub_df = df[df['strain_id'] == strain].copy()
    if sub_df.empty:
        raise ValueError(f"No data for strain {strain}")
    
    # Identify runs
    sub_df['run_id'] = (sub_df['position'].diff() < 0).cumsum()
    
    # Drop last run if incomplete
    last_run_id = sub_df['run_id'].max()
    last_run = sub_df[sub_df['run_id'] == last_run_id]
    if L not in last_run['position'].values:
        sub_df = sub_df[sub_df['run_id'] != last_run_id]
    
    # Compute xij
    sub_df['xij'] = sub_df['attempts'].apply(lambda a: estimate_k(a, N))
    
    # Build lookup: (run_id, position) -> xij
    run_pos_xij = sub_df.set_index(['run_id', 'position'])['xij'].to_dict()
    run_ids = sorted(sub_df['run_id'].unique())
    
    if len(run_ids) == 0:
        raise ValueError(f"No complete runs for strain {strain}")
    
    # Compute estimates for increasing number of runs
    results = []
    for S in range(1, len(run_ids) + 1):
        chosen_runs = run_ids[:S]
        mean_xj_per_pos = []

        for pos in range(0, L + 1):
            vals = [run_pos_xij[(r, pos)] for r in chosen_runs if (r, pos) in run_pos_xij]
            if len(vals) > 0:
                mean_xj_per_pos.append(np.mean(vals))

        mean_xj = np.array(mean_xj_per_pos)
        min_terms = np.minimum(1.0 + mean_xj, 20.0)
        log10_snc_samp = np.sum(np.log(min_terms)) / np.log(10.0)

        results.append({
            "Strain": strain,
            "SampleSize": S,
            "log10_Snc_est": log10_snc_samp
        })
    
    return pd.DataFrame(results)

# -------------------------
# Compute site scanning estimates
# -------------------------
df_res = pd.concat([compute_estimates_for_strain(df, s) for s in strain_list], ignore_index=True)

# -------------------------
# Compute absolute values
# -------------------------
df_nc = compute_nc_properties(
    base_dir=BASE_DIR,
    strain_id_file=STRAIN_IDS,
    S=20,
    L=PROTEIN_LENGTH,
)
abs_dict = dict(zip(df_nc["Key"], df_nc["Log10_Snc_est"]))

# -------------------------
# Build labels dict
# -------------------------
label_dict = get_labels(LABELS, STRAIN_IDS)

# -------------------------
# Plotting
# -------------------------
sns.set_style('white')
fig, ax = plt.subplots(figsize=(8, 5))

for i, strain in enumerate(df_res["Strain"].unique()):
    color = selected_colors[i % len(selected_colors)]
    df_strain = df_res[df_res["Strain"] == strain]

    # running average curve
    ax.plot(df_strain["SampleSize"], df_strain["log10_Snc_est"],
            color=color, linewidth=3, linestyle=':', label=f"{strain} running avg")

    # absolute bar
    abs_val = abs_dict.get(strain)
    if abs_val is not None:
        last_x = df_strain["SampleSize"].max()
        ax.hlines(y=abs_val,
                  xmin=last_x - 0.1, xmax=last_x + 0.1,
                  color=color, linewidth=3, linestyle='-',
                  label=f"{strain} absolute (10^{abs_val:.2f})")

    # label at end of curve
    last_x = df_strain["SampleSize"].max()
    last_y = df_strain["log10_Snc_est"].iloc[-1]
    display_label = label_dict.get(strain, strain)
    ax.text(last_x + 0.02, last_y, display_label, color=color, fontsize=10,
            verticalalignment='bottom')

# y-axis formatting: 10^x
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"$10^{{{int(y)}}}$"))
ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.grid(True, linestyle="--", alpha=0.2)

ax.set_xlabel("Number of loops")
ax.set_ylabel(r"$S_{NC, est}$")
sns.despine()
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/fig10.png", dpi = 300)
plt.close()
