# Mapping Protein Neutral Networks from Predicted Secondary Structure

This repository contains two complementary methods to explore a protein genotype–phenotype (GP) map using predicted secondary structure (Porter5) as a coarse-grained phenotype. The goal is to estimate properties of neutral components (NCs) including size, robustness, and local connectivity from computationally feasible samples of amino acid sequence space in a broader aim to understand the organisation of robustness and evolvability in proteins.

## Requirements

- Python 3.8
- Conda environment
- [Porter5](https://github.com/mircare/Porter5) (uses 7 CPUs)
- SLURM scheduler (for large-scale batch execution)

All Python dependencies are listed in `environment.yml`. Create the environment with:
```bash
conda env create -f environment.yml
conda activate porter5
```

### Porter5

This code assumes Porter5 is fully configured and functional, including required database search tools. Porter5 runtime depends on database searches and system setup. In our runs, predictions typically required ~1–2 minutes per sequence.

### Computational Considerations

Site scanning scales linearly with sequence length and neutral expansion depth. Neighbourhood enumeration requires 19L structure predictions per genotype. Therefore, we recommend SLURM array jobs for proteins longer than ~150 residues.

## Methods

### Inputs

This repository expects:

- `--strain-ids <strain_ids_file>`
  List of strain IDs to process.

- `--input-dir <ss3_input_dir>`
  Directory containing predicted structures of selected strain IDs.

### Outputs

- All outputs are written to
  `--base-dir <path_to_output>`


### Method 1: Site scanning (global exploration)

Sequentially mutates each position and retains neutral mutants, defined as those whose predicted secondary structure matches the parent. All neutral mutations are retained, enabling broad exploration of a neutral component.

**Usage example:**
```bash
python3 site_scanning.py \
    <key_index> \
    <array_task_id> \
    --base-dir <path_to_output> \
    --strain-ids <strain_ids_file> \
    --input-dir <ss3_input_dir> \
    --porter5 <path_to_Porter5.py> \
    --cpu 7
```
**Positional arguments**

- `<key_index>` — phenotype/strain NC index  
- `<array_task_id>` — integer identifier for parallel execution

### Pre-processing for Method 2

`choose_random_subset.py` — subsamples neutral genotypes from site-scanning output and prepares FASTA files plus Porter5 predictions for neighbourhood enumeration.

### Method 2: Exhaustive neighbourhood enumeration (local exploration)

For a given genotype, enumerates all 19 substitutions at each of L sites (19L mutants) and predicts secondary structure for each mutant. This provides accurate local versatility estimates for sampled genotypes.

**Usage example:**
```bash
array_idx=$SLURM_ARRAY_TASK_ID
seq_idx=$(( array_idx / <num_chunks> ))
chunk_idx=$(( array_idx % <num_chunks> ))

python3 neighbourhood_enumeration.py \
    <key_index> \
    ${seq_idx} \
    ${chunk_idx} \
    --base-dir <output_dir> \
    --strain-ids <strain_ids_file> \
    --input-dir <ss3_input_dir> \
    --porter5 <path_to_Porter5.py> \
    --cpu 7
```
**Positional arguments**

- `<key_index>` — phenotype/strain NC index  
- `${seq_idx}` — sampled genotype index
- `${chunk_idx}` — slice of sites (for batching long proteins)

### Notes

- Tested on: Linux HPC cluster with SLURM
- Large-scale runs require distributed batch execution.
- Replace placeholder paths (<...>) with your own directory structure.

### Analysis

- `functions/cleaning.py` — sequence filtering (length, canonical amino acids)
- `functions/selecting_strains.py` — removes non-human strains and selects representative phenotypes
- `functions/*_df.py` scripts — intermediate dataframe construction
- `functions/compute_nc_properties.py` — NC size and robustness estimators
- `functions/network_construction.py`, `functions/network_metrics.py` — partial NC graph reconstruction and metrics
- `functions/metadata_df.py`, `functions/get_labels.py` — strain metadata integration


### Reproducibility

Figure generation scripts will be added upon publication. All core algorithms and data-generation pipelines are included. The repository is structured to enable full regeneration of neutral component estimates from raw sequence inputs.

### Citation

Mapping protein neutral networks from predicted secondary structure
Nabiha Khawar, Sebastian E. Ahnert
Preprint: bioRxiv 2026.03.04.709605 https://doi.org/10.64898/2026.03.04.709605
Submitted manuscript: under review at Royal Society Interface. Citation will be updated upon publication.
