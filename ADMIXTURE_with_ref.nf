/*
 * Author: Justin Pelletier
 * Contact: justin.pelletier2@mcgill.ca
 * Year: 2026
 */

nextflow.enable.dsl = 2


def validateParams() {
    if (!params.cag_vcfs)
        error '--cag_vcfs is required'

    if (!params.reference_vcfs)
        error '--reference_vcfs is required'

    if (!params.reference_fasta)
        error '--reference_fasta is required'

    if (!params.cag_keep)
        error '--cag_keep is required'

    if (!params.reference_keep)
        error '--reference_keep is required'

    if (!params.cag_metadata)
        error '--cag_metadata is required'

    if (!params.reference_metadata)
        error '--reference_metadata is required'

    if (!params.admixture_executable)
        error '--admixture_executable is required'

    if (!file(params.admixture_executable).exists())
        error "ADMIXTURE executable not found: ${params.admixture_executable}"

    if ((params.max_k as int) < 2)
        error '--max_k must be >= 2'

    if ((params.pre_maf as double) < 0 ||
        (params.pre_maf as double) > 0.5)
        error '--pre_maf must be between 0 and 0.5'

    if ((params.pre_max_missing as double) < 0 ||
        (params.pre_max_missing as double) > 1)
        error '--pre_max_missing must be between 0 and 1'

    if (params.run_posthoc) {
        if ((params.admixture_cv as int) < 1) {
            error '--run_posthoc requires --admixture_cv to be greater than zero'
        }

        if ((params.posthoc_top_n as int) < 1) {
            error '--posthoc_top_n must be >= 1'
        }

        if ((params.posthoc_plot_top_n as int) < 1) {
            error '--posthoc_plot_top_n must be >= 1'
        }

        if ((params.posthoc_min_reference_n as int) < 1) {
            error '--posthoc_min_reference_n must be >= 1'
        }
    }
}


/*
 * Normalize and QC the CaG array dataset.
 *
 * A chromosome-position file is also generated so that only potentially
 * shared positions are retrieved from the much larger reference WGS VCF.
 */
process NORMALIZE_CAG {
    tag { "CaG:chr${chr}" }
    publishDir "${params.outdir}/01_qc_vcfs/cag", mode: 'copy'

    input:
    tuple val(chr), path(vcf), path(keep)
    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(chr),
          path("cag.chr${chr}.qc.vcf.gz"),
          path("cag.chr${chr}.qc.vcf.gz.tbi"),
          path("cag.chr${chr}.positions.tsv"),
          emit: qc_vcf

    script:
    def passArg = params.cag_require_pass ? '-f PASS' : ''

    def mafArg = (params.pre_maf as double) > 0
        ? "-i 'MAF>=${params.pre_maf}'"
        : ''

    def missingArg = (params.pre_max_missing as double) < 1
        ? "-e 'F_MISSING>${params.pre_max_missing}'"
        : ''

    """
    set -euo pipefail

    bcftools view \
        --threads ${task.cpus} \
        -S ${keep} \
        ${passArg} \
        -m2 \
        -M2 \
        -v snps \
        ${vcf} \
        -Ou |
        bcftools norm \
            --threads ${task.cpus} \
            -f ${fasta} \
            -Ou |
        bcftools +fill-tags \
            -Ou \
            -- -t MAF,F_MISSING |
        bcftools view \
            --threads ${task.cpus} \
            ${mafArg} \
            ${missingArg} \
            -Oz \
            -o cag.chr${chr}.qc.vcf.gz

    tabix -f -p vcf cag.chr${chr}.qc.vcf.gz

    n_cag=\$(bcftools index -n cag.chr${chr}.qc.vcf.gz)

    if [ "\$n_cag" -eq 0 ]; then
        echo "ERROR: CaG QC retained no variants for chromosome ${chr}" >&2
        exit 1
    fi

    bcftools query \
        -f '%CHROM\\t%POS\\n' \
        cag.chr${chr}.qc.vcf.gz \
        > cag.chr${chr}.positions.tsv

    if [ ! -s cag.chr${chr}.positions.tsv ]; then
        echo "ERROR: CaG position file is empty for chromosome ${chr}" >&2
        exit 1
    fi
    """
}


/*
 * Retrieve only reference variants located at CaG positions before
 * normalization and QC.
 */
process NORMALIZE_REFERENCE {
    tag { "Reference:chr${chr}" }
    publishDir "${params.outdir}/01_qc_vcfs/reference", mode: 'copy'

    input:
    tuple val(chr),
          path(reference_vcf),
          path(reference_index),
          path(reference_keep),
          path(cag_vcf),
          path(cag_tbi),
          path(cag_positions)

    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(chr),
          path(cag_vcf),
          path(cag_tbi),
          path("reference.chr${chr}.qc.vcf.gz"),
          path("reference.chr${chr}.qc.vcf.gz.tbi"),
          emit: paired_qc

    script:
    def passArg = params.reference_require_pass ? '-f PASS' : ''

    def mafArg = (params.pre_maf as double) > 0
        ? "-i 'MAF>=${params.pre_maf}'"
        : ''

    def missingArg = (params.pre_max_missing as double) < 1
        ? "-e 'F_MISSING>${params.pre_max_missing}'"
        : ''

    """
    set -euo pipefail

    bcftools view \
        --threads ${task.cpus} \
        -R ${cag_positions} \
        -S ${reference_keep} \
        ${passArg} \
        -m2 \
        -M2 \
        -v snps \
        ${reference_vcf} \
        -Ou |
        bcftools norm \
            --threads ${task.cpus} \
            -f ${fasta} \
            -Ou |
        bcftools +fill-tags \
            -Ou \
            -- -t MAF,F_MISSING |
        bcftools view \
            --threads ${task.cpus} \
            ${mafArg} \
            ${missingArg} \
            -Oz \
            -o reference.chr${chr}.qc.vcf.gz

    tabix -f -p vcf reference.chr${chr}.qc.vcf.gz

    n_reference=\$(bcftools index -n reference.chr${chr}.qc.vcf.gz)

    if [ "\$n_reference" -eq 0 ]; then
        echo "ERROR: No reference variants retained for chromosome ${chr}" >&2
        echo "Check that both VCFs use matching chromosome names." >&2
        exit 1
    fi
    """
}


/*
 * Identify exact shared alleles, merge the two sample sets and immediately
 * convert the chromosome to PLINK binary format.
 *
 * The merged chromosome VCF is temporary and is not published.
 */
process INTERSECT_MERGE {
    tag { "chr${chr}" }
    publishDir "${params.outdir}/02_plink_by_chr", mode: 'copy'

    input:
    tuple val(chr),
          path(cag),
          path(cag_tbi),
          path(reference),
          path(reference_tbi)

    output:
    tuple val(chr),
          path("chr${chr}.bed"),
          path("chr${chr}.bim"),
          path("chr${chr}.fam"),
          emit: chromosome_bed

    path "chr${chr}.intersection.stats.tsv"
    path "chr${chr}.log"

    script:
    """
    set -euo pipefail

    bcftools isec \
        --threads ${task.cpus} \
        -c none \
        -n=2 \
        -w1 \
        ${cag} \
        ${reference} \
        -Oz \
        -o cag.shared.vcf.gz

    bcftools isec \
        --threads ${task.cpus} \
        -c none \
        -n=2 \
        -w2 \
        ${cag} \
        ${reference} \
        -Oz \
        -o reference.shared.vcf.gz

    tabix -f -p vcf cag.shared.vcf.gz
    tabix -f -p vcf reference.shared.vcf.gz

    n_cag_qc=\$(bcftools index -n ${cag})
    n_reference_qc=\$(bcftools index -n ${reference})
    n_cag_shared=\$(bcftools index -n cag.shared.vcf.gz)
    n_reference_shared=\$(bcftools index -n reference.shared.vcf.gz)

    if [ "\$n_cag_shared" -ne "\$n_reference_shared" ]; then
        echo "ERROR: shared-variant counts differ for chromosome ${chr}" >&2
        echo "CaG: \$n_cag_shared" >&2
        echo "Reference: \$n_reference_shared" >&2
        exit 1
    fi

    if [ "\$n_cag_shared" -eq 0 ]; then
        echo "ERROR: no exact shared variants found for chromosome ${chr}" >&2
        exit 1
    fi

    bcftools merge \
        --threads ${task.cpus} \
        cag.shared.vcf.gz \
        reference.shared.vcf.gz \
        -Ou |
        bcftools annotate \
            --threads ${task.cpus} \
            --set-id '%CHROM:%POS:%REF:%FIRST_ALT' \
            -Oz \
            -o merged.chr${chr}.vcf.gz

    tabix -f -p vcf merged.chr${chr}.vcf.gz

    n_merged=\$(bcftools index -n merged.chr${chr}.vcf.gz)

    if [ "\$n_merged" -ne "\$n_cag_shared" ]; then
        echo "ERROR: merged and shared variant counts differ for chromosome ${chr}" >&2
        echo "Shared: \$n_cag_shared" >&2
        echo "Merged: \$n_merged" >&2
        exit 1
    fi

    plink \
        --vcf merged.chr${chr}.vcf.gz \
        --double-id \
        --allow-extra-chr \
        --make-bed \
        --out chr${chr}

    n_plink=\$(wc -l < chr${chr}.bim)
    n_samples=\$(wc -l < chr${chr}.fam)

    if [ "\$n_plink" -ne "\$n_merged" ]; then
        echo "ERROR: PLINK and VCF variant counts differ for chromosome ${chr}" >&2
        echo "VCF: \$n_merged" >&2
        echo "PLINK: \$n_plink" >&2
        exit 1
    fi

    printf \
        'chromosome\\tcag_qc_variants\\treference_qc_variants\\tshared_exact_variants\\tplink_variants\\tsamples\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        '${chr}' \
        "\$n_cag_qc" \
        "\$n_reference_qc" \
        "\$n_cag_shared" \
        "\$n_plink" \
        "\$n_samples" \
        > chr${chr}.intersection.stats.tsv
    """
}


/*
 * Merge chromosome-specific PLINK binary filesets.
 */
process MERGE_PLINK_CHROMOSOMES {
    tag 'autosomes'
    publishDir "${params.outdir}/03_genomewide_plink", mode: 'copy'

    input:
    tuple val(chromosomes),
          path(beds),
          path(bims),
          path(fams)

    output:
    tuple path('combined_unfiltered.bed'),
          path('combined_unfiltered.bim'),
          path('combined_unfiltered.fam'),
          emit: merged_bed

    path 'chromosome_merge_list.txt'
    path 'combined_unfiltered.log'
    path 'combined_unfiltered.nosex', optional: true
    path 'combined_unfiltered.stats.tsv'

    script:
    def firstChr = chromosomes[0]

    def mergeEntries = (1..<chromosomes.size())
        .collect { i ->
            "'${beds[i].name} ${bims[i].name} ${fams[i].name}'"
        }
        .join(' ')

    """
    set -euo pipefail

    printf '%s\\n' \
        ${mergeEntries} \
        > chromosome_merge_list.txt

    plink \
        --bfile chr${firstChr} \
        --merge-list chromosome_merge_list.txt \
        --allow-extra-chr \
        --make-bed \
        --out combined_unfiltered

    n_expected=\$(awk 'END { print NR }' ${bims})
    n_merged=\$(wc -l < combined_unfiltered.bim)
    n_samples=\$(wc -l < combined_unfiltered.fam)

    if [ "\$n_expected" -ne "\$n_merged" ]; then
        echo "ERROR: chromosome merge changed the total variant count" >&2
        echo "Expected: \$n_expected" >&2
        echo "Merged: \$n_merged" >&2
        exit 1
    fi

    if [ "\$n_samples" -eq 0 ]; then
        echo "ERROR: no samples remain after chromosome merge" >&2
        exit 1
    fi

    printf \
        'dataset\\tsamples\\tvariants\\ncombined_unfiltered\\t%s\\t%s\\n' \
        "\$n_samples" \
        "\$n_merged" \
        > combined_unfiltered.stats.tsv
    """
}


/*
 * Apply sample and variant QC across the complete genome-wide dataset.
 */
process PLINK_QC {
    tag 'combined QC'
    publishDir "${params.outdir}/04_plink_qc", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path('combined_qc.bed'),
          path('combined_qc.bim'),
          path('combined_qc.fam'),
          emit: qc_bed

    path 'combined_qc.log'
    path 'combined_qc.nosex', optional: true
    path 'combined_qc.stats.tsv'

    script:
    """
    set -euo pipefail

    input_samples=\$(wc -l < ${fam})
    input_variants=\$(wc -l < ${bim})

    plink \
        --bfile combined_unfiltered \
        --allow-extra-chr \
        --mind ${params.mind} \
        --geno ${params.geno} \
        --maf ${params.maf} \
        --make-bed \
        --out combined_qc

    output_samples=\$(wc -l < combined_qc.fam)
    output_variants=\$(wc -l < combined_qc.bim)

    printf \
        'stage\\tsamples\\tvariants\\ninput\\t%s\\t%s\\nafter_qc\\t%s\\t%s\\n' \
        "\$input_samples" \
        "\$input_variants" \
        "\$output_samples" \
        "\$output_variants" \
        > combined_qc.stats.tsv
    """
}


/*
 * Perform LD pruning and create the final ADMIXTURE input.
 */
process LD_PRUNE {
    tag 'LD pruning'
    publishDir "${params.outdir}/05_admixture_input", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path('admixture_input.bed'),
          path('admixture_input.bim'),
          path('admixture_input.fam'),
          emit: pruned_bed

    path 'prune.prune.in'
    path 'prune.prune.out'
    path 'prune.log'
    path 'admixture_input.log'
    path 'admixture_input.nosex', optional: true
    path 'admixture_input.stats.tsv'

    script:
    """
    set -euo pipefail

    plink \
        --bfile combined_qc \
        --allow-extra-chr \
        --indep-pairwise \
            ${params.ld_window} \
            ${params.ld_step} \
            ${params.ld_r2} \
        --out prune

    plink \
        --bfile combined_qc \
        --allow-extra-chr \
        --extract prune.prune.in \
        --make-bed \
        --out admixture_input

    n_samples=\$(wc -l < admixture_input.fam)
    n_variants=\$(wc -l < admixture_input.bim)

    printf \
        'dataset\\tsamples\\tvariants\\nadmixture_input\\t%s\\t%s\\n' \
        "\$n_samples" \
        "\$n_variants" \
        > admixture_input.stats.tsv
    """
}


/*
 * Run ADMIXTURE independently for K=1 through max_k.
 */
process RUN_ADMIXTURE {
    tag { "K=${k}" }

    publishDir {
        "${params.outdir}/06_admixture/K${k}"
    }, mode: 'copy'

    input:
    tuple val(k), path(bed), path(bim), path(fam)

    output:
    tuple val(k),
          path("admixture_input.${k}.Q"),
          path("admixture_input.${k}.P"),
          path("admixture.K${k}.log"),
          path(fam),
          emit: results

    script:
    def cvArg = (params.admixture_cv as int) > 0
        ? "--cv=${params.admixture_cv}"
        : ''

    """
    set -euo pipefail

    ${params.admixture_executable} \
        ${cvArg} \
        -j${task.cpus} \
        -s ${params.admixture_seed} \
        admixture_input.bed \
        ${k} \
        |& tee admixture.K${k}.log
    """
}


/*
 * Plot ADMIXTURE ancestry proportions.
 */
process PLOT_ADMIXTURE {
    tag { "K=${k}" }
    publishDir "${params.outdir}/07_admixture_plots", mode: 'copy'

    input:
    tuple val(k), path(q), path(p), path(log), path(fam)
    path cag_metadata, stageAs: 'cag_metadata.tsv'
    path reference_metadata, stageAs: 'reference_metadata.tsv'
    path cag_keep, stageAs: 'cag.keep.txt'
    path reference_keep, stageAs: 'reference.keep.txt'

    output:
    path "admixture.K${k}.png"
    path "admixture.K${k}.plot_order.tsv"

    script:
    """
    set -euo pipefail

    plot_admixture.R \
        ${q} \
        ${fam} \
        ${k} \
        ${cag_metadata} \
        ${reference_metadata} \
        ${cag_keep} \
        ${reference_keep} \
        admixture.K${k}.png \
        admixture.K${k}.plot_order.tsv
    """
}


/*
 * Combine ADMIXTURE cross-validation errors and create the CV plot.
 */
process PLOT_CV {
    tag 'cross-validation'

    publishDir(
        "${params.outdir}/08_cross_validation",
        mode: 'copy',
        overwrite: true
    )

    input:
    path logs

    output:
    path 'admixture_cv_errors.tsv', emit: cv_table
    path 'admixture_cv_plot.png', emit: cv_plot

    script:
    """
    set -euo pipefail

    plot_cv.R \
        admixture_cv_errors.tsv \
        admixture_cv_plot.png \
        ${logs}
    """
}

process RANK_REFERENCE_POPULATIONS_AND_CLUSTERS {
    tag "best CV: K=${k}"

    publishDir(
        "${params.outdir}/09_posthoc",
        mode: 'copy',
        overwrite: true
    )

    input:
    tuple val(k), path(q), path(p), path(admixture_log), path(fam)

    path cag_metadata, stageAs: 'cag_metadata.tsv'
    path reference_metadata, stageAs: 'reference_metadata.tsv'

    output:
    path 'K*.tsv', emit: tables
    path 'K*.png', emit: plots

    script:
    """
    set -euo pipefail

    mkdir -p posthoc_admixture/K${k}

    cp -L ${q} \
        posthoc_admixture/K${k}/admixture_input.${k}.Q

    cp -L ${fam} \
        posthoc_admixture/K${k}/admixture_input.fam

    rank_reference_populations_and_clusters.R \
        --admixture-dir posthoc_admixture \
        --cag-metadata ${cag_metadata} \
        --reference-metadata ${reference_metadata} \
        --output-dir . \
        --k ${k} \
        --top-n ${params.posthoc_top_n} \
        --plot-top-n ${params.posthoc_plot_top_n} \
        --min-reference-n ${params.posthoc_min_reference_n}
    """
}


workflow {
    validateParams()

    /*
     * Reusable input files.
     */
    fasta_ch = Channel.value(
        tuple(
            file(params.reference_fasta, checkIfExists: true),
            file("${params.reference_fasta}.fai", checkIfExists: true)
        )
    )

    chr_values = params.chromosomes as List

    cag_keep_file = file(
        params.cag_keep,
        checkIfExists: true
    )

    reference_keep_file = file(
        params.reference_keep,
        checkIfExists: true
    )

    cag_metadata_ch = Channel.value(
        file(params.cag_metadata, checkIfExists: true)
    )

    reference_metadata_ch = Channel.value(
        file(params.reference_metadata, checkIfExists: true)
    )

    cag_keep_ch = Channel.value(cag_keep_file)
    reference_keep_ch = Channel.value(reference_keep_file)

    /*
     * Construct CaG chromosome inputs.
     */
    cag_ch = Channel
        .fromList(chr_values)
        .map { chr ->
            def cag_path = params.cag_vcfs
                .toString()
                .replace('{chr}', chr.toString())

            tuple(
                chr,
                file(cag_path, checkIfExists: true),
                cag_keep_file
            )
        }

    NORMALIZE_CAG(
        cag_ch,
        fasta_ch
    )

    /*
     * Construct reference chromosome inputs, including TBI or CSI indexes.
     */
    reference_ch = Channel
        .fromList(chr_values)
        .map { chr ->
            def reference_path = params.reference_vcfs
                .toString()
                .replace('{chr}', chr.toString())

            def tbi_path = file("${reference_path}.tbi")
            def csi_path = file("${reference_path}.csi")

            def reference_index = tbi_path.exists()
                ? tbi_path
                : csi_path

            tuple(
                chr,
                file(reference_path, checkIfExists: true),
                file(reference_index, checkIfExists: true),
                reference_keep_file
            )
        }

    /*
     * Join reference inputs with the normalized CaG chromosome files.
     */
    reference_with_cag = reference_ch.join(
        NORMALIZE_CAG.out.qc_vcf,
        by: 0
    )

    NORMALIZE_REFERENCE(
        reference_with_cag,
        fasta_ch
    )

    /*
     * Intersect exact alleles, merge samples and convert immediately to PLINK.
     */
    INTERSECT_MERGE(
        NORMALIZE_REFERENCE.out.paired_qc
    )

    /*
     * Preserve chromosome-specific BED/BIM/FAM groups and order them
     * numerically by chromosome.
     */
    ordered_plink_sets = INTERSECT_MERGE.out.chromosome_bed
        .collect(flat: false)
        .map { rows ->
            def sortedRows = rows.sort { a, b ->
                (a[0] as int) <=> (b[0] as int)
            }

            def chromosomes = sortedRows.collect { row ->
                row[0]
            }

            def beds = sortedRows.collect { row ->
                row[1]
            }

            def bims = sortedRows.collect { row ->
                row[2]
            }

            def fams = sortedRows.collect { row ->
                row[3]
            }

            tuple(chromosomes, beds, bims, fams)
        }

    /*
     * Merge binary chromosome filesets, then apply genome-wide QC.
     */
    MERGE_PLINK_CHROMOSOMES(
        ordered_plink_sets
    )

    PLINK_QC(
        MERGE_PLINK_CHROMOSOMES.out.merged_bed
    )

    LD_PRUNE(
        PLINK_QC.out.qc_bed
    )


    /*
     * Generate one ADMIXTURE input tuple for every K from 2 through max_k.
     */
    admixture_inputs = LD_PRUNE.out.pruned_bed
        .flatMap { bed, bim, fam ->
            (2..(params.max_k as int)).collect { k ->
                tuple(k, bed, bim, fam)
            }
        }

    RUN_ADMIXTURE(
        admixture_inputs
    )

    /*
     * Plot ancestry proportions for every K.
     */
    PLOT_ADMIXTURE(
        RUN_ADMIXTURE.out.results,
        cag_metadata_ch,
        reference_metadata_ch,
        cag_keep_ch,
        reference_keep_ch
    )

    /*
     * Collect ADMIXTURE cross-validation logs from every K.
     */
    cv_logs = RUN_ADMIXTURE.out.results
        .map { k, q, p, admixture_log, fam ->
            admixture_log
        }
        .collect()

    /*
     * Create the CV table and plot.
     */
    if ((params.admixture_cv as int) > 0) {
        PLOT_CV(
            cv_logs
        )

        /*
         * Run the post-hoc analysis using the K with the lowest
         * cross-validation error.
         */
        if (params.run_posthoc) {
            best_k_ch = PLOT_CV.out.cv_table
                .map { cv_file ->
                    def lines = cv_file.text
                        .readLines()
                        .findAll { line ->
                            line.trim()
                        }

                    if (lines.size() < 2) {
                        error(
                            "No CV results were found in ${cv_file}"
                        )
                    }

                    def header = lines[0].split('\t')

                    def k_column = header.findIndexOf { column ->
                        column == 'K'
                    }

                    def cv_column = header.findIndexOf { column ->
                        column == 'CV_error'
                    }

                    if (k_column < 0 || cv_column < 0) {
                        error(
                            "The CV table must contain K and " +
                            "CV_error columns: ${cv_file}"
                        )
                    }

                    def cv_results = lines
                        .drop(1)
                        .collect { line ->
                            def fields = line.split('\t')

                            if (
                                fields.size() <= k_column ||
                                fields.size() <= cv_column
                            ) {
                                error(
                                    "Malformed row in the CV table: ${line}"
                                )
                            }

                            def tested_k = fields[k_column] as int
                            def cv_error = fields[cv_column] as double

                            if (
                                Double.isNaN(cv_error) ||
                                Double.isInfinite(cv_error)
                            ) {
                                error(
                                    "Invalid CV error for K=${tested_k}: " +
                                    "${cv_error}"
                                )
                            }

                            tuple(
                                tested_k,
                                cv_error
                            )
                        }

                    /*
                     * Sort by CV error and then K.
                     * The smaller K is selected in the event of a tie.
                     */
                    def best_result = cv_results.sort { left, right ->
                        def cv_comparison =
                            (left[1] as double) <=>
                            (right[1] as double)

                        if (cv_comparison != 0) {
                            return cv_comparison
                        }

                        return(
                            (left[0] as int) <=>
                            (right[0] as int)
                        )
                    }[0]

                    def best_k = best_result[0] as int
                    def best_cv = best_result[1] as double

                    log.info(
                        "Lowest ADMIXTURE CV error: " +
                        "K=${best_k}, CV_error=${best_cv}"
                    )

                    tuple(best_k)
                }

            /*
             * Join the selected K with its corresponding ADMIXTURE
             * Q, P, log and FAM files.
             */
            best_admixture_result = RUN_ADMIXTURE.out.results
                .join(
                    best_k_ch,
                    by: 0
                )

            RANK_REFERENCE_POPULATIONS_AND_CLUSTERS(
                best_admixture_result,
                cag_metadata_ch,
                reference_metadata_ch
            )
        }
    }
}
