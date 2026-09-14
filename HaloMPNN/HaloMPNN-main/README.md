## How to run

Run `submit_nf.sh` to submit a slurm job running the pipeline

```
# Activate conda environment
mamba activate /projects/f_sdk94_1/conda/envs/nextflow

# How to run pipeline with test file
nextflow run main.nf -profile slurm --input_fasta test.fasta

# How to run pipeline resuming previous run
nextflow run main.nf -profile slurm --input_fasta test.fasta -resume

```


## Pipeline overview

```
INPUT.fasta
    │
    ▼
[REMOVE_DUPLICATES]          seqkit rmdup (mmseqs env)
    │
    ▼
[EXTRACT_UNIPROT_IDS]        seqkit seq → uniprot_ids.txt
    │
    ▼
[QUERY_AFDB_API]             afdb_api.py → metadata.csv
    │
    ▼
[FILTER_LENGTH_PLDDT]        filter_length_plddt.py (len>50, pLDDT>80)
    │
[MMSEQS_CLUSTER] ──────────── (runs in parallel on nodups FASTA)
    │                        createdb → cluster → createtsv
    │                        createsubdb → convert2fasta → repr_ids.txt
    ▼
[FILTER_METADATA_REPR]       filter_metadata_repr.py
    │
    ▼
[DOWNLOAD_PDB_STRUCTURES]    afdb_download.py → pdb_files/
    │
    ▼
[PARSE_PDB_TO_PT]            parse_pdb_noX.py  (parallelised per .pdb)
    │
    ▼
[CREATE_CLUSTER_LIST]        ls → awk self-assignment TSV
    │
    ▼
[CREATE_CLUSTER_LISTS_PY]    create_cluster_lists.py (aifold env)
    │
    ▼
[MPNN_RETRAIN]               training.py — GPU job (aifold env)
```

## Requirements

- Nextflow ≥ 23.x
- Conda environments already installed at the paths in `nextflow.config`
- All Python scripts present in `bin/`:
  - `afdb_api.py`
  - `filter_length_plddt.py`
  - `filter_metadata_repr.py`
  - `afdb_download.py`
  - `parse_pdb_noX.py`
  - `create_cluster_lists.py`


## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_fasta` | *(required)* | Input FASTA file |
| `--outdir` | `results` | Output directory |
| `--seqid` | `0.30` | MMseqs2 min-seq-id threshold |
| `--split_memory_limit` | `50G` | MMseqs2 split memory limit |
| `--min_length` | `50` | Minimum sequence length filter |
| `--min_plddt` | `80` | Minimum pLDDT filter |
| `--num_neighbors` | `48` | MPNN num_neighbors |
| `--backbone_noise` | `0.20` | MPNN backbone noise |
| `--num_epochs` | `300` | MPNN training epochs |
| `--batch_size` | `3000` | MPNN batch size |
| `--max_protein_length` | `5000` | MPNN max protein length |
| `--num_examples_per_epoch` | `1000` | MPNN examples per epoch |
| `--mmseqs_conda` | `/projects/.../mmseqs` | Path to mmseqs conda env |
| `--aifold_conda` | `/projects/.../aifold` | Path to aifold conda env |
| `--scripts_dir` | `${projectDir}/bin` | Directory containing Python scripts |

## Output structure

```
results/
├── dedup/
│   ├── input_nodups.fasta
│   ├── dupseqs.fasta
│   ├── numdup.txt
│   └── uniprot_ids.txt
├── metadata/
│   ├── metadata.csv
│   ├── metadata_filtered.csv
│   ├── filter_stats.txt
│   └── repr_metadata.csv
├── clustering/
│   ├── cluster_assignments.tsv
│   ├── repr_seqs.fasta
│   └── repr_ids.txt
├── structures/
│   └── pdb_files/
├── pt_files/
│   └── <id>/
├── cluster/
│   ├── cluster_list.tsv
│   └── data_files/
├── mpnn_training/
│   └── train_output/
└── reports/
    ├── report.html
    ├── timeline.html
    └── trace.txt
```
