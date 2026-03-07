import os
import pandas as pd

# -------------------------
# Constants
# -------------------------
ATTEMPTS_PER_SITE = 19


# -------------------------
# Parsing
# -------------------------
def parse_neighbourhood_enum(base_dir, strain_id_file, S=20, L=566):
    """
    Parse neighbourhood enumeration output files into a long-form dataframe.

    Parameters
    ----------
    base_dir : str
        Directory containing neighbourhood enumeration folders (e.g. KEY_S/)
    strain_id_file : str
        Path to file containing strain IDs (first column)
    S : int
        Number of sampled neighbourhoods per strain
    L : int
        Protein length

    Returns
    -------
    df : pandas.DataFrame
        Columns:
        ['strain', 'sequence', 'position', 'neutral_count', 'neutral_neighbours']
    """

    expected_positions = set(range(L))

    # Load strain IDs
    with open(strain_id_file) as f:
        keys = [l.strip().split("\t")[0] for l in f if l.strip()]

    records = []

    for key in keys:
        dir_name = f"{key}_{S}"
        dir_path = os.path.join(base_dir, dir_name)

        if not os.path.isdir(dir_path):
            print(f"Warning: {dir_path} not found, skipping.")
            continue

        for i in range(S):
            file_path = os.path.join(dir_path, f"{i}.txt")

            if not os.path.exists(file_path):
                print(f"Warning: {file_path} not found. Skipping sequence {i}.")
                continue

            df = pd.read_csv(
                file_path,
                sep="\t",
                header=None,
                names=["position", "neutral_count", "neutral_neighbours"],
            )

            df["position"] = df["position"].astype(int)

            # Check completeness
            found_positions = set(df["position"].unique())
            missing = expected_positions - found_positions
            if missing:
                print(
                    f"Missing positions in {file_path}: "
                    f"{sorted(missing)[:10]}"
                    f"{'...' if len(missing) > 10 else ''} "
                    f"({len(missing)} missing total)"
                )

            for _, row in df.iterrows():
                records.append(
                    {
                        "strain": key,
                        "sequence": i,
                        "position": row["position"],
                        "neutral_count": row["neutral_count"],
                        "neutral_neighbours": row["neutral_neighbours"],
                    }
                )

    return pd.DataFrame.from_records(records)


# -------------------------
# Derived quantities
# -------------------------
def compute_avg_neighbourhood_properties(df):
    """
    Compute average neutral counts and robustness per strain and position.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of parse_neighbourhood_enum

    Returns
    -------
    avg_df : pandas.DataFrame
        Columns:
        ['strain', 'position', 'neutral_count', 'robustness']
    """

    avg_df = (
        df.groupby(["strain", "position"], as_index=False)
        .agg({"neutral_count": "mean"})
    )

    avg_df["robustness"] = avg_df["neutral_count"] / ATTEMPTS_PER_SITE

    return avg_df