# DataSet + Reference ADMIXTURE Pipeline

## Description

This Nextflow DSL2 pipeline is designed to run [ADMIXTURE](https://dalexander.github.io/admixture/) analysis on a user-provided DataSet and a reference dataset.

It performs variant QC, identifies variants shared between both datasets, merges their samples, applies LD pruning, and runs ADMIXTURE from K=2 to a user-defined maximum K.

The pipeline is configured for a UNIX environment and, more specifically, the Digital Research Alliance of Canada Narval cluster using Slurm. The configuration can be adapted for other systems.

## Pipeline steps

1. Subset the DataSet using a sample keep list.
2. Normalize and filter the DataSet VCFs.
3. Extract the retained DataSet variant positions.
4. Retrieve only those positions from the WGS reference VCFs.
5. Normalize and filter the reference variants.
6. Retain exact CHROM/POS/REF/ALT matches between datasets.
7. Merge samples separately for each chromosome.
8. Convert each chromosome directly to PLINK BED/BIM/FAM.
9. Merge the chromosome-specific PLINK files.
10. Apply genome-wide sample and variant QC.
11. Perform LD pruning.
12. Run ADMIXTURE from K=2 through `max_k`.
13. Generate ADMIXTURE ancestry and cross-validation plots.

Restricting the WGS reference to DataSet positions avoids processing every WGS variant.

## Directory structure

```text
project/
├── ADMIXTURE_with_ref.nf
├── nextflow.config
├── README.md
└── bin/
    ├── plot_admixture.R
    └── plot_cv.R
```

Nextflow automatically adds executable scripts from `bin/` to the process `PATH`.

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

### Sample keep lists

Prepare one file for each dataset containing one sample ID per line:

```text
SAMPLE001
SAMPLE002
SAMPLE003
```

The sample IDs must exactly match the VCF sample IDs.

### Population metadata

Prepare one tab-separated metadata file for each dataset:

```tsv
sample_id	population
SAMPLE001	Population_1
SAMPLE002	Population_1
```

The column names must be exactly:

```text
sample_id
population
```

Additional columns are ignored.

Population labels are used as facets in the ADMIXTURE plots. Labels such as `Cluster_1`, `Cluster_2`, and `Cluster_10` are ordered numerically.

## Configuration

Copy the example configuration:

```bash
cp nextflow.config.example nextflow.config
```

The following parameters must be adapted:

```groovy
params {
    dataset_vcfs = '/path/to/dataset/chr{chr}.vcf.gz'
    reference_vcfs = '/path/to/reference/chr{chr}.vcf.gz'
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

The VCF patterns must contain the literal `{chr}` placeholder.

### Variant filtering

If the DataSet uses `FILTER=.`, use:

```groovy
dataset_require_pass = false
```

To retain only PASS reference variants, use:

```groovy
reference_require_pass = true
```

The main QC parameters are:

```groovy
pre_maf = 0.0
pre_max_missing = 0.02

mind = 0.05
geno = 0.02
maf = 0.01
```

### LD pruning

The LD-pruning parameters are:

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

The current plotting palette supports up to K=10.

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

Load Nextflow:

```bash
module load nextflow
```

Make the plotting scripts executable:

```bash
chmod +x bin/plot_admixture.R bin/plot_cv.R
```

Confirm that ADMIXTURE is executable:

```bash
test -x /path/to/admixture
```

Confirm that the R packages are available:

```bash
module load StdEnv/2023 r/4.6.1

Rscript -e 'packageVersion("ggplot2")'
Rscript -e 'packageVersion("patchwork")'
```

If either package is unavailable, install it in your personal R library before launching the pipeline.

## Running the pipeline

```bash
nextflow run ADMIXTURE_with_ref.nf \
    -c nextflow.config \
    -profile narval \
    -resume
```

Parameters may also be overridden from the command line:

```bash
nextflow run ADMIXTURE_with_ref.nf \
    -c nextflow.config \
    -profile narval \
    --max_k 10 \
    --ld_r2 0.1 \
    -resume
```

## Outputs

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

### Main outputs

| Directory | Contents |
|---|---|
| `01_qc_vcfs` | Normalized and filtered chromosome VCFs |
| `02_plink_by_chr` | Chromosome-specific PLINK files and intersection statistics |
| `03_genomewide_plink` | Merged genome-wide PLINK dataset |
| `04_plink_qc` | PLINK dataset after genome-wide QC |
| `05_admixture_input` | LD-pruned dataset used by ADMIXTURE |
| `06_admixture/K*` | ADMIXTURE Q, P, and log files |
| `07_admixture_plots` | Ancestry plots and plotted sample-order TSV files |
| `08_cross_validation` | CV-error summary and CV plot |
| `pipeline_info` | Nextflow report, trace, timeline, and DAG |

The ancestry plots display reference samples on top and DataSet samples underneath. Population labels are shown as facets, and every K uses a fixed set of ancestry-component colors.
