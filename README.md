# DataSet + Reference ADMIXTURE Pipeline

## Description

This Nextflow DSL2 pipeline is designed to run [ADMIXTURE](https://dalexander.github.io/admixture/) on a user-provided DataSet and a reference dataset.

It performs variant quality control, identifies variants shared between the two datasets, merges their samples, applies LD pruning, and runs ADMIXTURE from K=2 to a user-defined maximum K.

The pipeline is configured for a UNIX environment and, more specifically, the Digital Research Alliance of Canada Narval cluster using Slurm. The configuration can be adapted for other computing environments.

## Pipeline steps

1. Subset the DataSet using a sample keep list.
2. Normalize and filter the DataSet VCFs.
3. Extract the retained DataSet variant positions.
4. Retrieve only those positions from the WGS reference VCFs.
5. Normalize and filter the reference variants.
6. Retain exact CHROM/POS/REF/ALT matches between datasets.
7. Merge the samples chromosome by chromosome.
8. Convert each chromosome directly to PLINK BED/BIM/FAM format.
9. Merge the chromosome-specific PLINK files.
10. Apply genome-wide sample and variant QC.
11. Perform LD pruning.
12. Run ADMIXTURE from K=2 through `max_k`.
13. Generate ADMIXTURE ancestry plots.
14. Generate a cross-validation error plot.

Restricting the WGS reference to positions present in the DataSet avoids processing every variant in the WGS reference.

## Repository structure

```text
ADMIXTURE/
├── assets/
│   └── empty_metadata.tsv
├── bin/
│   ├── plot_admixture.R
│   └── plot_cv.R
├── examples/
│   ├── dataset.keep.example.txt
│   ├── dataset_metadata.example.tsv
│   ├── reference.keep.example.txt
│   └── reference_metadata.example.tsv
├── ADMIXTURE_with_ref.config
├── ADMIXTURE_with_ref.nf
├── README.md
└── Run_ADMIXTURE_with_ref.sh
```

### `assets/`

The `assets/` directory contains supporting files used by the pipeline:

```text
empty_metadata.tsv
```

This file can be supplied as a metadata input when no population labels are available for one of the datasets.

### `bin/`

The `bin/` directory contains the R plotting scripts:

```text
plot_admixture.R
plot_cv.R
```

Nextflow automatically adds executable files from `bin/` to the process `PATH`.


### `examples/`

The `examples/` directory contains templates showing the required formats for:

- DataSet sample keep lists
- Reference sample keep lists
- DataSet population metadata
- Reference population metadata

Copy and modify these files for your analysis or use your pre-existing files as long as they are formatted in the exact same way.

## Requirements

The pipeline requires:

- Nextflow 23.10 or newer
- bcftools
- tabix
- PLINK 1.9
- ADMIXTURE
- R
- R packages `ggplot2` and `patchwork`

The Narval configuration currently uses:

```text
bcftools/1.22
plink/1.9b_6.21-x86_64
r/4.6.1
```

ADMIXTURE must be downloaded separately from its [official repository](https://github.com/NovembreLab/admixture).

## Input requirements

The DataSet and reference VCFs must:

- Use the same genome build.
- Use compatible chromosome names.
- Be split by chromosome.
- Be bgzip-compressed.
- Have a `.tbi` or `.csi` index.

The reference FASTA must have a corresponding `.fai` index.

Sample IDs must be unique between the DataSet and reference datasets.

## Files to prepare

### VCF files

The VCF paths must contain the literal `{chr}` placeholder:

```groovy
dataset_vcfs = '/path/to/dataset/chr{chr}.dataset.vcf.gz'
reference_vcfs = '/path/to/reference/chr{chr}.reference.vcf.gz'
```

For example, `{chr}` will be replaced with `1`, `2`, and so forth.

### Reference FASTA

Provide the reference genome FASTA:

```groovy
reference_fasta = '/path/to/reference/GRCh38.fa'
```

The corresponding index must be present:

```text
GRCh38.fa.fai
```

If needed, create it with:

```bash
samtools faidx GRCh38.fa
```

### Sample keep lists

Prepare one keep list for each dataset. Each file must contain one sample ID per line:

```text
SAMPLE001
SAMPLE002
SAMPLE003
```

Sample IDs must exactly match the corresponding VCF sample IDs.

Example files are provided in:

```text
examples/dataset.keep.example.txt
examples/reference.keep.example.txt
```

Set their paths in the configuration:

```groovy
dataset_keep = '/path/to/dataset.keep.txt'
reference_keep = '/path/to/reference.keep.txt'
```

A sample must not occur in both keep lists.

### Population metadata

The metadata files must be tab-separated and contain the exact columns:

```tsv
sample_id	population
SAMPLE001	Population_1
SAMPLE002	Population_1
```

Additional columns are ignored.

Example files are provided in:

```text
examples/dataset_metadata.example.tsv
examples/reference_metadata.example.tsv
```

Set their paths in the configuration:

```groovy
dataset_metadata = '/path/to/dataset_metadata.tsv'
reference_metadata = '/path/to/reference_metadata.tsv'
```

Population labels determine the facets in the ADMIXTURE plots.

Labels such as:

```text
Cluster_1
Cluster_2
Cluster_10
```

are ordered numerically.

If metadata are unavailable for a dataset, use:

```text
assets/empty_metadata.tsv
```

For example:

```groovy
dataset_metadata = "${projectDir}/assets/empty_metadata.tsv"
```

## Configuration

Before running the pipeline, edit:

```text
ADMIXTURE_with_ref.config
```

At minimum, update:

```groovy
params {
    dataset_vcfs = '/path/to/dataset/chr{chr}.dataset.vcf.gz'
    reference_vcfs = '/path/to/reference/chr{chr}.reference.vcf.gz'
    reference_fasta = '/path/to/reference/GRCh38.fa'

    dataset_keep = '/path/to/dataset.keep.txt'
    reference_keep = '/path/to/reference.keep.txt'

    dataset_metadata = '/path/to/dataset_metadata.tsv'
    reference_metadata = '/path/to/reference_metadata.tsv'

    admixture_executable = '/path/to/admixture'

    outdir = 'results'
    max_k = 10
}
```

### FILTER handling

If the DataSet uses `FILTER=.`, set:

```groovy
dataset_require_pass = false
```

To retain only PASS variants from the reference, set:

```groovy
reference_require_pass = true
```

### Quality-control parameters

```groovy
pre_maf = 0.0
pre_max_missing = 0.02

mind = 0.05
geno = 0.02
maf = 0.01
```

### LD pruning

```groovy
ld_window = 50
ld_step = 5
ld_r2 = 0.1
```

A smaller `ld_r2` value applies more stringent LD pruning.

### ADMIXTURE parameters

```groovy
max_k = 10
admixture_cv = 10
admixture_seed = 2026
```

The pipeline runs ADMIXTURE from K=2 through `max_k`.

The current plotting palette supports a maximum of K=10.

Set:

```groovy
admixture_cv = 0
```

to disable cross-validation.

### Narval allocation

Replace the allocation placeholder in the Narval profile:

```groovy
profiles {
    narval {
        process.executor = 'slurm'
        process.clusterOptions = '--account=your-allocation'
    }
}
```

## Setup on Narval

### Load Nextflow

```bash
module load nextflow
```

### Check the required modules

```bash
module spider bcftools/1.22
module spider plink/1.9b_6.21-x86_64
module spider r/4.6.1
```

The pipeline loads the required modules separately for each process using the configuration file.

### Install ADMIXTURE

Download ADMIXTURE from:

```text
https://github.com/NovembreLab/admixture
```

Make the binary executable:

```bash
chmod +x /path/to/admixture
```

Verify it:

```bash
test -x /path/to/admixture
```

Then set:

```groovy
admixture_executable = '/path/to/admixture'
```

### Check the R packages

```bash
module load StdEnv/2023 r/4.6.1

Rscript -e 'packageVersion("ggplot2")'
Rscript -e 'packageVersion("patchwork")'
```

If necessary, install them in your personal R library:

```bash
Rscript -e 'install.packages(
    c("ggplot2", "patchwork"),
    repos = "https://cloud.r-project.org"
)'
```

### Make the plotting scripts executable

```bash
chmod +x bin/plot_admixture.R bin/plot_cv.R
```

## Running the pipeline

Run the supplied launch script:

```bash
sbatch Run_ADMIXTURE_with_ref.sh
```

Alternatively, run Nextflow directly:

```bash
nextflow run ADMIXTURE_with_ref.nf \
    -c ADMIXTURE_with_ref.config \
    -profile narval \
    -resume
```

Parameters can be overridden from the command line:

```bash
nextflow run ADMIXTURE_with_ref.nf \
    -c ADMIXTURE_with_ref.config \
    -profile narval \
    --max_k 10 \
    --ld_r2 0.1 \
    -resume
```

Using `-resume` allows completed tasks to be reused when their inputs, scripts, and parameters have not changed.

## Output structure

```text
results/
├── 01_qc_vcfs/
│   ├── dataset/
│   └── reference/
├── 02_plink_by_chr/
├── 03_genomewide_plink/
├── 04_plink_qc/
├── 05_admixture_input/
├── 06_admixture/
│   ├── K2/
│   ├── K3/
│   └── ...
├── 07_admixture_plots/
├── 08_cross_validation/
└── pipeline_info/
```

## Main outputs

| Directory | Contents |
|---|---|
| `01_qc_vcfs/dataset` | Normalized and filtered DataSet VCFs |
| `01_qc_vcfs/reference` | Reference VCFs restricted to DataSet positions |
| `02_plink_by_chr` | Chromosome-specific PLINK files and intersection statistics |
| `03_genomewide_plink` | Merged genome-wide PLINK dataset |
| `04_plink_qc` | PLINK dataset after genome-wide QC |
| `05_admixture_input` | LD-pruned dataset used by ADMIXTURE |
| `06_admixture/K*` | ADMIXTURE Q, P, and log files |
| `07_admixture_plots` | ADMIXTURE ancestry plots and sample-order TSVs |
| `08_cross_validation` | Cross-validation error table and plot |
| `pipeline_info` | Nextflow trace, report, timeline, and DAG |

## Plotting outputs

The ADMIXTURE plots:

- Display reference samples on top.
- Display DataSet samples underneath.
- Use population labels as facets.
- Order numbered DataSet clusters numerically.
- Use a fixed 10-color palette.
- Do not display a component legend.

The plotted sample-order TSV contains:

```text
plot_position
sample_id
cohort
population
facet_position
```

The cross-validation directory contains:

```text
admixture_cv_errors.tsv
admixture_cv_plot.png
```
