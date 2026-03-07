import os
import pandas as pd
import numpy as np


def compute_nc_properties(
    base_dir,
    strain_id_file,
    S=20,
    L=566,
):
    """
    Compute neutral component size and robustness estimates from
    sampled mutational neighbourhoods.

    Parameters
    ----------
    base_dir : str
        Base directory containing <strain>_<S>/ folders
    strain_id_file : str
        Path to file containing strain IDs (first column)
    S : int, optional
        Number of sampled mutational neighbourhoods (default: 20)
    L : int, optional
        Protein length (default: 566)

    Returns
    -------
    df_res : pandas.DataFrame
        Columns:
        ['Key', 'Log10_Snc_est', 'Robustness']
    """

    expected_positions = set(range(L))
    results = []

    # Load strain IDs
    with open(strain_id_file) as f:
        keys = [l.strip().split("\t")[0] for l in f if l.strip()]

    for key in keys:
        dir_name = f"{key}_{S}"
        dir_path = os.path.join(base_dir, dir_name)

        dfs = []

        for i in range(S):
            file_path = os.path.join(dir_path, f"{i}.txt")

            if not os.path.exists(file_path):
                print(f"Warning: {file_path} not found. Skipping sequence {i}.")
                continue

            df = pd.read_csv(
                file_path,
                sep="\t",
                header=None,
                names=["Position", "Count", "Neighbours"],
            )

            df["Sequence"] = i
            df["Position"] = df["Position"].astype(int)

            # Check that all positions are present
            found_positions = set(df["Position"].unique())
            missing = expected_positions - found_positions
            if missing:
                print(
                    f"Missing positions in {file_path}: "
                    f"{sorted(missing)[:10]} "
                    f"{'...' if len(missing) > 10 else ''} "
                    f"({len(missing)} missing total)"
                )

            dfs.append(df)

        if not dfs:
            print(f"No data found for {key}, skipping.")
            continue

        all_data = pd.concat(dfs, ignore_index=True)

        # Mean neutral count per position across S sequences
        xj = all_data.groupby("Position")["Count"].mean()

        # Neutral component size estimate
        min_terms = np.minimum(1 + xj, 20)
        log_e = np.sum(np.log(min_terms))
        log10 = log_e / np.log(10)

        # Robustness estimate
        rob = (1 / (19 * L)) * np.sum(xj)

        results.append(
            {
                "Key": key,
                "Log10_Snc_est": log10,
                "Robustness": rob,
            }
        )

    return pd.DataFrame(results)
