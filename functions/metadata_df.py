import os
import re
import glob
import pandas as pd

"""
Parse Porter5 .ss3 files and attach strain metadata (strain ID, year, area)
extracted directly from the hemagglutinin FASTA file.
"""

# -----------------------------
# Paths
# -----------------------------
SS3_GLOB = "/path/to/ss3_files/*.ss3"  # Porter5 predictions
FASTA_FILE = "/data/hemagglutinin.txt"
STRAIN_DIR = "/path/to/strains"

# -----------------------------
# Helper function
# -----------------------------
def extract_year_and_area(description: str):
    """
    Extract (year, area) from a FASTA description line.
    area: first part after strain_id, i.e., between first '/' after strain_id
    year: 4-digit number in the last '/'-separated part
    """
    # Remove leading '>'
    desc = description.lstrip(">")

    parts = desc.split("/")

    # area = second element (if exists)
    area = parts[1].strip() if len(parts) > 1 else None

    # last part = everything after the last '/'
    last_part = parts[-1].strip()

    # Search for 4-digit year only in last part
    match = re.search(r"\b(19\d{2}|20\d{2})\b", last_part)
    year = int(match.group(0)) if match else None

    return year, area

# -----------------------------
# Parse FASTA file to get metadata
# -----------------------------
strain_info = {}

with open(FASTA_FILE) as f:
    for line in f:
        line = line.strip()
        if not line or not line.startswith(">"):
            continue

        # Extract strain ID from the 4th pipe-delimited field
        try:
            strain_id = line.split("|")[3]
        except IndexError:
            continue

        year, area = extract_year_and_area(line)
        strain_info[strain_id] = {"year": year, "area": area}

print(f"Found metadata for {len(strain_info)} strains.")

# -----------------------------
# Read SS3 files and merge metadata
# -----------------------------
dfs = []

for ss3_file in glob.glob(SS3_GLOB):
    # Extract the strain ID from the file name
    strain_id = os.path.basename(ss3_file).split(".")[0]

    # Skip SS3 files not in strain_info
    if strain_id not in strain_info:
        print(f"Warning: No metadata for strain {strain_id}, skipping.")
        continue

    ss3_df = pd.read_csv(
        ss3_file,
        delim_whitespace=True,
        skiprows=1,
        usecols=[0, 1, 2, 3, 4, 5],
        names=["index", "residue", "structure", "Helix", "Sheet", "Coil"],
    )

    meta = strain_info[strain_id]
    ss3_df["strain_id"] = strain_id
    ss3_df["year"] = meta.get("year")
    ss3_df["area"] = meta.get("area")

    dfs.append(ss3_df)

# Concatenate all data into a single DataFrame
final_df_extra = pd.concat(dfs, ignore_index=True)

print(final_df_extra.head())