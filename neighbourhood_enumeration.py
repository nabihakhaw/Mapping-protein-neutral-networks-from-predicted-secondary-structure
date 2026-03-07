"""
Exhaustive neighbourhood enumeration for local exploration.

This script tests all single–amino acid substitutions within a defined
segment ("chunk") of a protein sequence and identifies neutral mutations
based on predicted secondary structure (Porter5).

Intended usage:
    - Parallel execution via SLURM array jobs
    - Each job processes an independent chunk of positions

Input directory is provided via --input-dir and is expected to contain:
    <key>/seq_<index>.fasta
    <key>/seq_<index>.fasta.ss3
"""

from __future__ import annotations
import os
import argparse
import logging
import shutil
import subprocess
import numpy as np
import time
from pathlib import Path
from typing import List, Set, Tuple

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
CHUNK_SIZE = 40 # based on time limit for each job and estimated Porter5 prediction times

# ------------------------- Utilities -------------------------
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    
def load_keys(path: Path) -> List[str]:
    """Load strain identifiers from a tab-delimited file."""
    
    with path.open() as f:
        return [line.strip().split("\t")[0] for line in f if line.strip()]


def load_sequence(fasta_path: Path) -> str:
    """Load the amino acid sequence from a fasta output file. """
    
    with open(fasta_path, 'r') as f:
        return "".join(
            line.strip() for line in f if not line.startswith(">")
        )


def load_secondary_structure(ss3_path: Path) -> str:
    """Load only the secondary structure from a Porter5 SS3 file."""
    
    ss = []
    with open(ss3_path, 'r') as f:
        for line in f.readlines()[1:]: # skip header
            ss.append(line.strip().split("\t")[2])
    return "".join(ss)


def run_porter5(porter5_script: str, fasta_path: str, cpu: int):
    subprocess.run(
        [
            "python3",
            porter5_script,
            "-i",
            fasta_path,
            "--cpu",
            str(cpu),
            "--fast",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )



# ------------------------- Core Logic -------------------------

def get_secondary_structure_for_mutant(
    mutant_seq: str,
    key: str,
    seq_idx: int,
    position: int,
    porter5_script: str,
    base_dir: str,
    cpu: int,
) -> str | None:

    """
    Run Porter5 on a mutant sequence and return its predicted
    secondary structure. Temporary files are cleaned up afterward.
    """

    tmp_dir = f"{base_dir}/{key}/{seq_idx}_{position}"
    os.makedirs(tmp_dir, exist_ok=True)

    fasta_path = f"{tmp_dir}/{key}_{seq_idx}.fasta"
    ss3_path = f"{fasta_path}.ss3"

    with open(fasta_path, "w") as f:
        f.write(f">{key}\n{mutant_seq}")

    run_porter5(porter5_script, fasta_path, cpu)

    if not os.path.exists(ss3_path):
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return None

    ss = load_secondary_structure(ss3_path)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    
    return ss


def enumerate_neutral_mutations(
    seq: str,
    ori_ss: str,
    key: str,
    seq_idx: int,
    chunk_idx: int,
    porter5_script: Path,
    base_dir: Path,
    cpu: int,
):
    """
    Enumerate neutral mutations for a single chunk of positions.
    """
    L = len(seq)
    start = chunk_idx * CHUNK_SIZE
    end = min(start + CHUNK_SIZE, L)

    out_dir = base_dir / key / str(seq_idx)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"chunk_{chunk_idx}.txt"
    
    porter5_script = str(porter5_script)
    base_dir = str(base_dir)

    neutral_mutants = set()

    for i in range(start, end):
        neutral_aas = set()

        # log once per position for speed
        logging.info(f"Position {i}")

        for aa in AMINO_ACIDS:
            if aa == seq[i]:
                continue

            mutant = list(seq)
            mutant[i] = aa
            mutant = "".join(mutant)

            if mutant in neutral_mutants:
                continue

            mutant_ss = get_secondary_structure_for_mutant(
                mutant,
                key,
                seq_idx,
                i,
                porter5_script,
                base_dir,
                cpu,
            )

            if mutant_ss == ori_ss:
                neutral_aas.add(aa)
                neutral_mutants.add(mutant)
                

        with open(out_path, "a") as f:
            f.write(
                f"{i}\t{len(neutral_aas)}\t"
                f"{','.join(sorted(neutral_aas))}\n"
            )



# ------------------------- Main -------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Exhaustive neighbourhood enumeration using Porter5"
    )
    parser.add_argument("key_idx", type=int)
    parser.add_argument("seq_idx", type=int)
    parser.add_argument("chunk_idx", type=int)

    parser.add_argument("--base-dir", required=True)
    parser.add_argument("--strain-ids", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--porter5", required=True)
    parser.add_argument("--cpu", type=int, default=7)

    args = parser.parse_args()

    setup_logging()

    base_dir = Path(args.base_dir)
    strain_ids = Path(args.strain_ids)
    input_dir = Path(args.input_dir)
    porter5_script = Path(args.porter5)

    keys = load_keys(strain_ids)
    key = keys[args.key_idx]

    fasta_path = input_dir / key / f"seq_{args.seq_idx + 1}.fasta"
    ss3_path = fasta_path.with_suffix(".fasta.ss3")

    ori_seq = load_sequence(fasta_path)
    ori_ss = load_secondary_structure(ss3_path)


    enumerate_neutral_mutations(
        ori_seq,
        ori_ss,
        key,
        args.seq_idx,
        args.chunk_idx,
        porter5_script,
        base_dir,
        args.cpu,
    )

    print(
        f"Completed: key={args.key_idx}, "
        f"seq={args.seq_idx}, chunk={args.chunk_idx}"
    )


if __name__ == "__main__":
    main()
