# DataSet + Reference ADMIXTURE Pipeline

This Nextflow DSL2 pipeline combines a chromosome-split genotyping dataset with a chromosome-split WGS reference dataset and runs ADMIXTURE.

The pipeline:

1. Filters and normalizes the DataSet VCFs.
2. Restricts the WGS reference to positions present in the DataSet.
3. Retains exact CHROM/POS/REF/ALT matches.
4. Merges samples chromosome by chromosome.
5. Converts each chromosome directly to PLINK format.
6. Merges chromosomes and applies genome-wide QC.
7. Performs LD pruning.
8. Runs ADMIXTURE from K=2 through `max_k`.
9. Produces ADMIXTURE and cross-validation plots.

## Directory structure

```text
project/
├── ADMIXTURE_with_ref.nf
├── ADMIXTURE_with_ref.config
├── README.md
└── bin/
    ├── plot_admixture.R
    └── plot_cv.R
```

Make the plotting scripts executable:

```bash
chmod +x bin/plot_admixture.R bin/plot_cv.R
```

Nextflow automatically adds the `bin/` directory to the process `PATH`.

## Requirements

- Nextflow 23.10 or newer
- bcftools and tabix
- PLINK 1.9
- ADMIXTURE
- R packages:
  - `ggplot2`
  - `patchwork`

Current Narval modules:

```groovy
module = 'StdEnv/2023:bcftools/1.22:plink/1.9b_6.21-x86_64'
```

For plotting:

```groovy
module = 'StdEnv/2023:r/4.6.1'
```

ADMIXTURE is provided using an absolute path in the configuration.

## Input files

Both datasets must:

- Use the same genome build.
- Use compatible chromosome names.
- Be split by chromosome.
- Be bgzip-compressed VCF files.
- Have `.tbi` or `.csi` indexes.

The reference FASTA must have a matching `.fai` index.

Example VCF patterns:

```groovy
dataset_vcfs = '/path/to/dataset/chr{chr}.vcf.gz'
reference_vcfs = '/path/to/reference/chr{chr}.vcf.gz'
```

The literal `{chr}` is replaced with each chromosome number.

## Sample keep lists

Provide one sample ID per line:

```text
SAMPLE001
SAMPLE002
SAMPLE003
```

The IDs must match the corresponding VCF sample IDs exactly.

A sample cannot occur in both keep lists.

```groovy
dataset_keep = '/path/to/dataset.keep.txt'
reference_keep = '/path/to/reference.keep.txt'
```

## Metadata files

Metadata files must be tab-separated and contain:

```tsv
sample_id	population
SAMPLE001	Population_1
SAMPLE002	Population_1
```

Additional columns are ignored.

```groovy
dataset_metadata = '/path/to/dataset_metadata.tsv'
reference_metadata = '/path/to/reference_metadata.tsv'
```

The reference populations are shown in the top row of each ADMIXTURE plot. DataSet populations or clusters are shown underneath.

Labels such as `Cluster_1`, `Cluster_2`, and `Cluster_10` are ordered numerically.

## Main parameters

```groovy
params {
    dataset_vcfs = '/path/to/dataset/chr{chr}.vcf.gz'
    reference_vcfs = '/path/to/reference/chr{chr}.vcf.gz'
    reference_fasta = '/path/to/GRCh38.fa'

    dataset_keep = '/path/to/dataset.keep.txt'
    reference_keep = '/path/to/reference.keep.txt'

    dataset_metadata = '/path/to/dataset_metadata.tsv'
    reference_metadata = '/path/to/reference_metadata.tsv'

    admixture_executable = '/path/to/admixture'

    chromosomes = [
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
        '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'
    ]

    outdir = 'results'

    max_k = 10

    dataset_require_pass = false
    reference_require_pass = true

    pre_maf = 0.0
    pre_max_missing = 0.02

    mind = 0.05
    geno = 0.02
    maf = 0.01

    ld_window = 50
    ld_step = 5
    ld_r2 = 0.1

    admixture_cv = 10
    admixture_seed = 2026
}
```

### FILTER handling

If the DataSet uses `FILTER=.`, set:

```groovy
dataset_require_pass = false
```

If the reference contains PASS and non-PASS variants, use:

```groovy
reference_require_pass = true
```

### ADMIXTURE K values

The pipeline runs from K=2 through `max_k`.

For example:

```groovy
max_k = 10
```

runs K=2, K=3, …, K=10.

The plotting palette currently supports a maximum of K=10.

## Running the pipeline

```bash
module load nextflow

nextflow run ADMIXTURE_with_ref.nf \
    -c ADMIXTURE_with_ref.config \
    -profile narval \
    -resume
```

Parameters can also be overridden:

```bash
nextflow run ADMIXTURE_with_ref.nf \
    -c ADMIXTURE_with_ref.config \
    -profile narval \
    --max_k 8 \
    --ld_r2 0.1 \
    -resume
```

## Outputs

```text
results/
├── 01_qc_vcfs/
├── 02_plink_by_chr/
├── 03_genomewide_plink/
├── 04_plink_qc/
├── 05_admixture_input/
├── 06_admixture/
├── 07_admixture_plots/
├── 08_cross_validation/
└── pipeline_info/
```

Important outputs include:

```text
06_admixture/K*/
    ADMIXTURE Q, P, and log files

07_admixture_plots/
    ADMIXTURE ancestry PNGs
    plotted sample-order TSVs

08_cross_validation/
    admixture_cv_errors.tsv
    admixture_cv_plot.png
```

The ADMIXTURE plots use a fixed 10-color palette and contain no legend. The Reference and DataSet samples are plotted in separate rows so their population facets remain readable.

## Common issues

### `bcftools: command not found`

Make sure the appropriate module is assigned to every process using bcftools:

```groovy
module = 'StdEnv/2023:bcftools/1.22'
```

### ADMIXTURE executable not found

Check that the configured binary exists and is executable:

```bash
test -x /path/to/admixture
```

### Missing reference index

Each reference VCF requires either:

```text
chr1.vcf.gz.tbi
```

or:

```text
chr1.vcf.gz.csi
```

### Existing Nextflow reports

Allow Nextflow to overwrite previous reports:

```groovy
timeline.overwrite = true
report.overwrite = true
trace.overwrite = true
dag.overwrite = true
```