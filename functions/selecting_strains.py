import os
import re
import pandas as pd

"""
Select representative human influenza strains for downstream analysis.

- Removes non-human strains based on FASTA metadata.
- Identifies top 5 phenotypes per year based on secondary structure.
- Selects earliest strain for each phenotype.
- Optionally includes historically important strains.
- Writes strain IDs to a text file for downstream pipelines.
"""

# -----------------------------
# Paths
# -----------------------------
FASTA_FILE = "/data/hemagglutinin.txt"
OUTPUT_FILE = "/path/to/strain_ids.txt"
METADATA = "/data/final_df_extra.pkl"

final_df_extra = pd.read_pickle(METADATA)

# -----------------------------
# Step 1: Identify non-human accessions
# -----------------------------
nonhuman_accessions = set()

with open(FASTA_FILE) as f:
    for line in f:
        if not line.startswith('>'):
            continue

        # Extract host inside the brackets, e.g., (A/swine/...)
        m = re.search(r'\(A/([^/]+)', line)
        if not m:
            continue
        field = m.group(1).strip()  # e.g., "swine", "American duck", "South Carolina"

        accession = line.split('|')[3]   # Extract accession/strain ID

        # Treat as non-human if the field is not capitalised (if field is area - then human, if fiels is host, then non-human)
        words = field.split()
        if any(not word[0].isupper() for word in words):
            nonhuman_accessions.add(accession)

print(f"Excluded {len(nonhuman_accessions)} non-human accessions.")

# -----------------------------
# Step 2: Compute phenotypes
# -----------------------------
# Concatenate structure strings per strain
phenotypes = final_df_extra.groupby('strain_id')['structure'].apply(''.join).reset_index(name='phenotype')

# Merge phenotype with year (drop duplicates)
phenotype_w_years = pd.merge(
    phenotypes,
    final_df_extra[['strain_id', 'year']],
    on='strain_id'
).drop_duplicates('strain_id').copy()

# Strip any extra spaces
phenotype_w_years['phenotype'] = phenotype_w_years['phenotype'].str.strip()

# -----------------------------
# Step 3: Count phenotypes per year & get top 5
# -----------------------------
phenotype_counts = phenotype_w_years.groupby(['year', 'phenotype']).size().reset_index(name='count')

top_phenotypes_per_year = phenotype_counts.groupby('year').apply(
    lambda x: x.nlargest(5, 'count')
).reset_index(drop=True)

# -----------------------------
# Step 4: Filter phenotype_w_years to include only top phenotypes
# -----------------------------
filtered_phenotypes = phenotype_w_years[
    phenotype_w_years['phenotype'].isin(top_phenotypes_per_year['phenotype'])
]

# -----------------------------
# Step 5: Identify first occurrence (earliest year) for each phenotype
# -----------------------------
first_occurrence_idx = filtered_phenotypes.groupby('phenotype')['year'].idxmin()
rep_genotypes = filtered_phenotypes.loc[first_occurrence_idx, ['phenotype', 'year', 'strain_id']].reset_index(drop=True)

# Merge to keep all top phenotypes
final_rep_genotypes = top_phenotypes_per_year[['phenotype']].merge(
    rep_genotypes,
    on='phenotype',
    how='left'
)

# -----------------------------
# Step 6: Remove non-human strains
# -----------------------------
final_rep_genotypes = final_rep_genotypes[
    ~final_rep_genotypes['strain_id'].isin(nonhuman_accessions)
]

# Optional: manually include historically important strains
manual_includes = ["ABD95350", "ACP41934"]
final_rep_genotypes = pd.concat([
    final_rep_genotypes,
    pd.DataFrame({"strain_id": manual_includes})
]).drop_duplicates()

# -----------------------------
# Step 7: Write strain_ids for downstream analysis
# -----------------------------
final_rep_genotypes['strain_id'].to_csv(OUTPUT_FILE, index=False, header=False)

print(f"Wrote {len(final_rep_genotypes)} human strain IDs to {OUTPUT_FILE}")