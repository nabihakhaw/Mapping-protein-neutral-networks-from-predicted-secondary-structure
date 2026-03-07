import pandas as pd


def get_labels(label_pkl, strain_id_file):
    """
    Build strain labels of the form 'area/year' for a given set of strains.

    Parameters
    ----------
    label_pkl : str
        Path to pickle file containing metadata (must include
        'strain_id', 'area', 'year' columns).
    strain_id_file : str
        Path to text file containing strain IDs (first column).

    Returns
    -------
    labels : dict
        Mapping {strain_id: 'area/year'}
    """

    # Load strain IDs of interest
    with open(strain_id_file) as f:
        strain_ids = {
            l.strip().split("\t")[0]
            for l in f
            if l.strip()
        }

    # Load label dataframe
    label_df = pd.read_pickle(label_pkl)

    # Keep only relevant strains
    label_df = label_df[label_df["strain_id"].isin(strain_ids)]

    # Drop duplicate rows per strain (important!)
    label_df = (
        label_df[["strain_id", "area", "year"]]
        .drop_duplicates(subset="strain_id")
    )

    # Build mapping
    labels = {
        row["strain_id"]: f"{row['area']}/{row['year']}"
        for _, row in label_df.iterrows()
    }

    return labels