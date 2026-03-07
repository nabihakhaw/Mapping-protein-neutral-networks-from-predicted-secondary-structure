import os
import pandas as pd

# -------------------------
# Constants
# -------------------------
PROTEIN_LENGTH = 566
ATTEMPTS_PER_SITE = 19

# -------------------------
# Parsing
# -------------------------
def parse_site_scanning(base_dir, strain_id_file):
    """
    Parse site-scanning output files into a long-form dataframe.

    Parameters
    ----------
    base_dir : str
        Directory containing site-scanning result files (*.txt)
    strain_id_file : str
        Path to file containing strain IDs (first column)

    Returns
    -------
    df : pandas.DataFrame
        Columns:
        ['filename', 'genotype', 'position', 'attempts', 'alphabet', 'is_neutral']
    """

    # Load strain IDs
    with open(strain_id_file) as f:
        strain_ids = {l.strip().split("\t")[0] for l in f if l.strip()}

    records = []

    for filename in os.listdir(base_dir):

        if not filename.endswith(".txt"):
            continue

        strain = filename.replace(".txt", "")
        if strain not in strain_ids:
            continue

        path = os.path.join(base_dir, filename)

        runs, current = [], []
        previous_pos = None

        with open(path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue

                try:
                    genotype, alphabet, pos, attempts = (
                        parts[0],
                        parts[1],
                        int(parts[2]),
                        parts[3],
                    )
                except ValueError:
                    continue

                if genotype.upper() == "NO_MUTANT":
                    attempts = ATTEMPTS_PER_SITE
                    alphabet = "[]"
                    is_neutral = 0
                else:
                    try:
                        attempts = int(attempts)
                    except ValueError:
                        continue
                    is_neutral = 1

                # New run if position resets
                if previous_pos is not None and pos < previous_pos:
                    runs.append(current)
                    current = []

                current.append((pos, attempts, alphabet, genotype, is_neutral))
                previous_pos = pos

        if current:
            runs.append(current)

        for run in runs:
            for pos, attempts, alphabet, genotype, is_neutral in run:
                records.append(
                    {
                        "filename": filename,
                        "genotype": genotype,
                        "position": pos,
                        "attempts": attempts,
                        "alphabet": alphabet,
                        "is_neutral": is_neutral,
                    }
                )

    df = pd.DataFrame.from_records(records)
    return df


# -------------------------
# Derived quantities
# -------------------------
def compute_robustness(df):
    """
    Add robustness column to site-scanning dataframe.
    """
    df = df.copy()
    df["robustness"] = (
        (ATTEMPTS_PER_SITE + 1 - df["attempts"])
        / (ATTEMPTS_PER_SITE * df["attempts"])
    )
    return df


# -------------------------
# Loop selection (helper)
# -------------------------
def select_first_n_loops(df, n_loops=2, max_pos=565):
    """
    Select the first n full site-scanning loops per strain.
    """
    filtered = []

    for strain, sub in df.groupby("filename"):
        sub = sub.sort_values("position").reset_index(drop=True)
        resets = sub.index[sub["position"] == 0].tolist()

        loops = []
        for i in range(len(resets)):
            start = resets[i]
            end = resets[i + 1] if i + 1 < len(resets) else None
            loop_df = sub.iloc[start:end]

            if loop_df["position"].max() >= max_pos:
                loops.append(loop_df)

            if len(loops) >= n_loops:
                break

        if loops:
            filtered.append(pd.concat(loops, ignore_index=True))

    return pd.concat(filtered, ignore_index=True) if filtered else pd.DataFrame()