# MitoSnake

MitoSnake is a Snakemake workflow for assembling plant mitochondrial genomes,
either from short reads alone (`st`) or from short and long reads together
(`ht`).

Current release: **0.9.0**.

## Requirements

- Linux, or WSL on Windows
- Conda (Miniconda or Anaconda); `install.sh` installs Mamba if it is not found

## Installation

```bash
git clone https://github.com/Liyun-Nie/MitoSnake.git
cd MitoSnake
bash install.sh
```

`install.sh` creates, or updates, the `mitosnake` and `mitosnake_tiara` conda
environments from `envs/`. It then replaces the default Tiara models in
`mitosnake_tiara` with the models in `resources/tiara_models/` (see
[Resources](#resources)). Pass `--no-replace-tiara` to keep the default models.

## Usage

MitoSnake reads samples from a tab-delimited sample list and writes all outputs
relative to the current directory. Activate the environment, change to the
working directory, and start the workflow:

```bash
conda activate mitosnake
cd /path/to/working_dir
python /path/to/MitoSnake/mitosnake.py --module ht \
  --resource-mode high \
  --sample-list sample_list.tsv \
  --jobs 2 \
  --keep-going
```

Run `python /path/to/MitoSnake/mitosnake.py --version` to show the version.
Add `--dry-run` to list the planned Snakemake jobs without running them.

### Sample list

The first line is a header. `sample_id`, `read1` and `read2` are required;
`ht` also requires `long_read` and `longread_type` (`clr`, `hifi` or `ont`).
`genome_size` (Mb) is optional and overrides `--genome-size` for that sample.

```
sample_id	read1	read2	long_read	longread_type	genome_size
Arabidopsis_thaliana	/path/to/R1.fastq.gz	/path/to/R2.fastq.gz	/path/to/long_reads.fastq.gz	hifi	130
```

### Modules

| `--module` | Steps |
|------------|-------|
| `st` | fastp trimming, short-read subsampling, SPAdes assembly, Tiara filtering of contigs |
| `ht` | The `st` steps, then long-read length filtering (≥ 10 kb) and subsampling, mapping to the mitochondrial candidate contigs with minimap2, recruitment filtering, Tiara filtering of recruited reads, and Flye assembly |

### Outputs

Each sample gets a directory named after its `sample_id`:

- `st`: `assembly/contigs.fasta` (SPAdes contigs) and
  `assembly/mitochondrion_contigs.fasta` (contigs kept after Tiara filtering;
  all contigs when `--skip-dl-filter` is set)
- `ht`: the `st` outputs plus `assembly/hybrid_assembly.fasta` and
  `assembly/hybrid_assembly.gfa` (Flye assembly)
- `logs/`: per-step logs and `consolidated_workflow.log`

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `--module` | `ht` | `st` or `ht` |
| `--sample-list` | | Tab-delimited sample list (see above) |
| `--resource-mode` | `default` | Thread, memory and disk presets: `default`, `high` or `ultra` |
| `--cores` | 4 | Passed to Snakemake `--cores` |
| `--jobs` | 1 | Passed to Snakemake `--jobs` (maximum number of jobs run in parallel) |
| `--keep-going` | off | Continue with independent jobs after a job fails |
| `--short-cov` | 2.0 | Target short-read coverage after subsampling |
| `--long-cov` | 10.0 | Target long-read coverage after subsampling |
| `--genome-size` | 400 | Estimated nuclear genome size in Mb, used to compute coverage |
| `--min-match` | 2000 | Minimum matching bases for a long read to be recruited |
| `--min-coverage` | 0.3 | Minimum ratio of matching bases to read length for recruitment |
| `--min-length` | 10000 | Written to `config.yaml`; the `ht` long-read length filter is fixed at 10 kb |
| `--skip-dl-filter` | off | Skip Tiara filtering |
| `--working-dir` | current directory | Where `config.yaml` is written; start MitoSnake from this directory |
| `--conda-path` | detected | Conda installation prefix |
| `--dry-run` | off | Snakemake dry-run |
| `--unlock` | off | Remove a stale Snakemake lock from the working directory |
| `--clean` | off | Remove the `.snakemake` directory |

## Resources

`resources/tiara_models/` contains neural-network and TF-IDF model files that
replace the default Tiara classifiers in the `mitosnake_tiara` environment.
They originate from or are derived from [Tiara](https://github.com/ibe-uw/tiara).
See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for the file list,
provenance and the upstream MIT licence.

We thank the Tiara authors. Please cite: Karlicki M, Antonowicz S,
Karnkowska A. Tiara: deep learning-based classification system for
eukaryotic sequences. Bioinformatics. 2022;38(2):344–350.

## Citation

If you use MitoSnake, please cite the GitHub repository:
https://github.com/Liyun-Nie/MitoSnake (v0.9.0).

## Contact

Liyun Nie (`nieliyun18@163.com`)

## License

MitoSnake source code is released under the MIT License. See `LICENSE`.
The bundled Tiara model files remain subject to the upstream Tiara terms
described in `THIRD_PARTY_NOTICES.md`.
