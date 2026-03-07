"""
Subsample a random set of neutral genotypes from the
site-scanning output and prepare FASTA files plus Porter5 secondary
structure predictions for neighbourhood enumeration.
"""


from __future__ import annotations
import os
import random
import numpy as np
from pathlib import Path
from typing import List


# ------------------- User-configurable paths -------------------

INPUT_DIR = Path("/path/to/site_scanning_results")
STRAIN_IDS = Path("/path/to/strain_ids.txt")
BASE_DIR = Path("/path/to/subsampled")
PORTER5_SCRIPT = Path("/path/to/Porter5/Porter5.py")


NUM_SAMPLES = 20
OVERWRITE = False

# ------------------------- Utilities -------------------------

def load_keys(path: Path) -> List[str]:
    """Load strain identifiers from a tab-delimited file."""
    
    with path.open() as f:
        return [line.strip().split("\t")[0] for line in f if line.strip()]


def subsample_genotypes(input_file: Path, n: int) -> list[str] | None:
    """
    Subsample genotypes from a site-scanning output file.

    Filters:
        - Only sequences of length 566
        - Skip lines containing 'NO_MUTANT'
    """
    if not input_file.exists():
        print(f"[WARN] Input file not found: {input_file}")
        return None

    with input_file.open() as f:
        genotypes = [
            line.strip().split("\t")[0]
            for line in f
            if line.strip() and len(line.strip().split("\t")[0]) == 566 and "NO_MUTANT" not in line
        ]

    if len(genotypes) < n:
        print(f"[WARN] Not enough genotypes ({len(genotypes)} < {n})")
        return None

    return random.sample(genotypes, n)



def run_porter5(fasta_path: Path):
    """Run Porter5 secondary structure prediction on a FASTA file."""
    
    cmd = f"python3 {PORTER5_SCRIPT} -i {fasta_path} --cpu 7 --fast"
    os.system(cmd)


# ------------------------- Main Logic -------------------------
keys = load_keys(STRAIN_IDS)

for key in keys:
    
    input_file = Path(INPUT_DIR) / f"{key}.txt"
    strain_out_dir = Path(BASE_DIR) / key
    strain_out_dir.mkdir(parents=True, exist_ok=True)

    subsampled_file = strain_out_dir / f"{key}.txt"

    # Subsample genotypes
    if subsampled_file.exists() and not OVERWRITE:
        print(f"[{key}] Subsample exists, skipping subsampling.")
    else:
        subsample = subsample_genotypes(input_file, NUM_SAMPLES)
        if subsample is None:
            continue

        with subsampled_file.open("w") as out:
            for g in subsample:
                out.write(g + "\n")

        print(f"[{key}] Subsampled {NUM_SAMPLES} genotypes")


    # Prepare FASTAs + run Porter5
    with subsampled_file.open() as f:
        for i, line in enumerate(f, start=1):
            seq = line.strip()
            if not seq:
                continue

            fasta_path = strain_out_dir / f"seq_{i}.fasta"

            if fasta_path.exists() and not OVERWRITE:
                print(f"[{key}] FASTA exists, skipping seq_{i}")
                continue


            # Write FASTA
            with fasta_path.open("w") as out:
                out.write(f">seq_{i}\n{seq}\n")

            print(f"[{key}] Running Porter5 on seq_{i}")
            run_porter5(fasta_path)
