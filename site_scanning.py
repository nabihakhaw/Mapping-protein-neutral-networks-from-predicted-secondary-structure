"""
Site scanning for global exploration via stepwise neutral mutations.

This script mutates all positions in the sequence sequentially with 
single–amino acid substitutions and identifies neutral mutations based 
on predicted secondary structure (Porter5). Each accepted neutral mutant 
becomes the template for subsequent mutations, generating a stepwise 
walk on a neutral network.
"""

from __future__ import annotations
import argparse
import json
import logging
import os
import random
import shutil
import subprocess
import numpy as np
import time
from pathlib import Path
from typing import List, Set, Tuple

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
TARGET_MUTANTS = 5000

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


def load_sequence_and_ss(ss3_path: Path) -> Tuple[str, str]:
    """
    Load the amino acid sequence and predicted secondary structure
    from a Porter5 SS3 output file.
    """
    
    seq, ss = [], []
    with ss3_path.open() as f:
        for line in f.readlines()[1:]:    # Skip header
            fields = line.strip().split("\t")
            seq.append(fields[1])         # Extract sequence
            ss.append(fields[2])          # Extract secondary structure
    return "".join(seq), "".join(ss)


def run_porter5(porter5_script: Path, fasta_path: Path, cpu: int):
    """Run Porter5 secondary structure prediction on a FASTA file."""
    
    cmd = [
        "python3",
        str(porter5_script),
        "-i",
        str(fasta_path),
        "--cpu",
        str(cpu),
        "--fast",
    ]
    subprocess.run(
        cmd, 
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        )


# ------------------------- State Handling -------------------------

def load_failed_positions(path: Path) -> Set[int]:
    """
    Load positions that failed to yield a neutral mutant in the current scan loop.

    These positions are skipped during the ongoing loop to avoid repeatedly
    testing sites where all substitutions have already been exhausted. 
    This prevents the algorithm from becoming stuck in an infinite
    loop testing the same non-viable positions.
    """

    if not path.exists():
        return set()
    with path.open() as f:
        return set(json.load(f).get("failed_positions", []))


def save_failed_positions(path: Path, failed_positions: Set[int]):
    """
    Save positions that failed to yield a neutral mutant during the current scan loop.

    Failed positions are persisted to disk to allow resumption within the same
    loop. At the end of a full pass over the sequence, this file is overwritten
    (cleared), as the mutational context has changed and previously failed
    positions may become viable again.
    """
        
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump({"failed_positions": sorted(failed_positions)}, f)



def get_resume_position(neutral_file: Path, protein_len: int) -> int:
    """
    Determine the next sequence position to mutate based on the last
    recorded entry in the neutral mutant file.

    Returns 0 if the file does not exist or cannot be parsed.
    """
    
    if not neutral_file.exists():
        return 0    # Start at 0 if file doesn't exist

    with neutral_file.open() as f:
        try:
            last = list(f)[-1].strip().split("\t")  # Read only the last line
            return (int(last[-2]) + 1) % protein_len   # Extract position and loop back to 0 if needed. 566 = protein length (HA);
        except (IndexError, ValueError):
            return 0



# ------------------------- Main Logic -------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("key_idx", type=int)
    parser.add_argument("--base-dir", required=True)
    parser.add_argument("--strain-ids", required=True)
    parser.add_argument("--input-dir", required=True) # must have Porter5 predicted structures for original sequences
    parser.add_argument("--porter5", required=True)
    parser.add_argument("--cpu", type=int, default=7)
    args = parser.parse_args()

    setup_logging()
    
    base_dir = Path(args.base_dir)
    keys = load_keys(Path(args.strain_ids))
    key = keys[args.key_idx]

    ss3_path = Path(args.input_dir) / f"{key}.fasta.ss3"
    ori_seq, ori_ss = load_sequence_and_ss(ss3_path)
    protein_len = len(ori_seq)

    # Create files
    neutral_file = base_dir / f"{key}.txt"
    failed_file = base_dir / "failed_positions" / f"{key}.json"
    work_root = base_dir / key
    work_root.mkdir(parents=True, exist_ok=True)

    neutral_mutants: Set[str] = set()
    neutral_seq: List[str] = []

    if neutral_file.exists():
        with neutral_file.open() as f:
            for line in f:
                mutant = line.split("\t")[0]
                neutral_mutants.add(mutant)
                neutral_seq.append(mutant)
    else:
        neutral_mutants.add(ori_seq)
        neutral_seq.append(ori_seq)
        

    failed_positions = load_failed_positions(failed_file)
    start_index = get_resume_position(neutral_file, protein_len)

    logging.info(f"Starting from position {start_index}")
    logging.info(f"Failed positions loaded: {sorted(failed_positions)}")

    while len(neutral_mutants) < TARGET_MUTANTS:

        for i in range(start_index, protein_len):
            if i in failed_positions or len(neutral_mutants) >= TARGET_MUTANTS:
                continue

            seq = list(neutral_seq[-1])
            candidates = list(AMINO_ACIDS - {seq[i]})  # Possible mutations, excluding original residue
            random.shuffle(candidates)  # Shuffle to introduce randomness

            attempts = 0
            tested = set()
            found = False

            for aa in candidates:
                attempts += 1
                tested.add(aa)

                mutant = seq[:] # Copy the current sequence
                mutant[i] = aa
                mutant = "".join(mutant)

                if mutant in neutral_mutants:
                    continue  # Skip already known neutral mutants

                tmp_dir = work_root / f"pos_{i}"
                tmp_dir.mkdir(exist_ok=True)

                fasta = tmp_dir / f"{key}.fasta"
                fasta.write_text(f">{key}\n{mutant}")

                try:
                    run_porter5(Path(args.porter5), fasta, args.cpu)
                except subprocess.CalledProcessError:
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    continue

                # Create SS3 file and read secondary structure as mutant_ss
                ss3 = fasta.with_suffix(".fasta.ss3")
                if not ss3.exists():
                    logging.error(f"SS3 file missing for mutant {mutant} at position {i}")
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    continue

                _, mutant_ss = load_sequence_and_ss(ss3)   #ignore sequence
                shutil.rmtree(tmp_dir, ignore_errors=True)  # Clean up temp mutant directory

                if mutant_ss == ori_ss:
                    neutral_mutants.add(mutant)
                    neutral_seq.append(mutant)
                    found = True

                    with neutral_file.open("a") as f:
                        f.write(
                            f"{mutant}\t{sorted(tested)}\t{i}\t{attempts}\n"
                        )

                    logging.info(
                        f"Neutral mutant accepted at pos {i} after {attempts} tries"
                    )
                    break  # Move to the next position

            if not found:
                failed_positions.add(i)
                save_failed_positions(failed_file, failed_positions)
                with neutral_file.open("a") as f:
                    f.write(f"NO_MUTANT\t[]\t{i}\t{attempts}\n")
                logging.info(f"Position {i} exhausted; no neutral mutant found.")

        start_index = 0  # After finishing one full loop, restart from position 0
        failed_positions.clear()
        save_failed_positions(failed_file, failed_positions)

    logging.info("Neutral walk completed")


if __name__ == "__main__":
    main()
