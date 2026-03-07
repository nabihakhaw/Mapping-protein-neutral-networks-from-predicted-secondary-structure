import os
from Bio import SeqIO

"""
Extract strain IDs from a file, keep only sequences of length 566,
write each sequence to its own directory
"""

# Paths
input_fasta = "/data/hemagglutinin.txt"
output_dir = "/path/to/strains"

TARGET_LENGTH = 566
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

valid_ids = []

# Parse FASTA file
for record in SeqIO.parse(input_fasta, "fasta"):
    # Extract strain ID (assumes pipe-delimited header)
    # Example: ...|...|...|BAA02769|...
    strain_id = record.id.split("|")[3]
    sequence = str(record.seq)

    # Filter by length
    if len(sequence) != TARGET_LENGTH:
        continue

    # Filter non-canonical amino acids
    if not set(sequence).issubset(CANONICAL_AA):
        continue

    valid_ids.append(strain_id)

    # Create output directory
    strain_dir = os.path.join(output_dir, strain_id)
    os.makedirs(strain_dir, exist_ok=True)

    # Write sequence to FASTA
    fasta_path = os.path.join(strain_dir, f"{strain_id}.fasta")
    with open(fasta_path, "w") as f_out:
        f_out.write(f">{strain_id}\n{sequence}\n")

print(
    f"Found {len(valid_ids)} sequences of length {TARGET_LENGTH}. "
)