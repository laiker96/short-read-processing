"""Optional condition-consensus and cross-condition ATAC atlas branch."""


rule merge_atac_condition_short_fragments:
    input:
        lambda wc: [
            ATAC_SHORT_BAMS[sample]
            for sample in ATAC_ATLAS_CONDITIONS[wc.condition].samples
        ]
    output:
        bam=(
            f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/bam/{{condition}}.fragments-lt"
            f"{ATAC_FRAGMENT_MAXIMUM}.bam"
        ),
        bai=(
            f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/bam/{{condition}}.fragments-lt"
            f"{ATAC_FRAGMENT_MAXIMUM}.bam.bai"
        )
    wildcard_constraints:
        condition=ATAC_ATLAS_CONDITION_RE
    threads: 6
    resources:
        mem_mb=8000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/logs/bam/merge.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        temporary={output.bam:q}.tmp
        samtools merge -@ {threads} -f "$temporary" {input:q} > {log:q} 2>&1
        samtools quickcheck -v "$temporary" 2>> {log:q}
        samtools index -@ {threads} "$temporary" "$temporary.bai" 2>> {log:q}
        mv "$temporary" {output.bam:q}
        mv "$temporary.bai" {output.bai:q}
        """


rule shift_atac_condition_short_fragments:
    input:
        bam=lambda wc: ATAC_ATLAS_CONDITION_SHORT_BAMS[wc.condition],
        bai=lambda wc: ATAC_ATLAS_CONDITION_SHORT_BAIS[wc.condition]
    output:
        bam=temp(f"{WORK_ROOT}/atac_atlas/{{condition}}.shifted.bam"),
        bai=temp(f"{WORK_ROOT}/atac_atlas/{{condition}}.shifted.bam.bai")
    wildcard_constraints:
        condition=ATAC_ATLAS_CONDITION_RE
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/logs/tracks/alignmentSieve.log"
    shell:
        "mkdir -p $(dirname {output.bam:q}) $(dirname {log:q}) && "
        "rm -f {output.bam:q}.unsorted && "
        "alignmentSieve --ATACshift -b {input.bam:q} -o {output.bam:q}.unsorted "
        "--numberOfProcessors {threads} > {log:q} 2>&1 && "
        "samtools sort -@ {threads} -o {output.bam:q} {output.bam:q}.unsorted "
        "2>> {log:q} && rm -f {output.bam:q}.unsorted && "
        "samtools index -@ {threads} {output.bam:q} {output.bai:q} 2>> {log:q} && "
        "samtools quickcheck -v {output.bam:q} 2>> {log:q}"


rule atac_condition_bigwig:
    input:
        bam=f"{WORK_ROOT}/atac_atlas/{{condition}}.shifted.bam",
        bai=f"{WORK_ROOT}/atac_atlas/{{condition}}.shifted.bam.bai"
    output:
        bigwig=(
            f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/tracks/{{condition}}.fragments-lt"
            f"{ATAC_FRAGMENT_MAXIMUM}.Tn5-shifted.CPM.bw"
        )
    params:
        bin_size=int(ATAC_REFINEMENT["bigwig_bin_size"])
    wildcard_constraints:
        condition=ATAC_ATLAS_CONDITION_RE
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/logs/tracks/bamCoverage.log"
    shell:
        "mkdir -p $(dirname {output.bigwig:q}) $(dirname {log:q}) && "
        "bamCoverage -b {input.bam:q} -o {output.bigwig:q} "
        "--outFileFormat bigwig --normalizeUsing CPM --binSize {params.bin_size} "
        "--numberOfProcessors {threads} > {log:q} 2>&1"


rule call_atac_condition_candidates:
    input:
        bam=lambda wc: ATAC_ATLAS_CONDITION_SHORT_BAMS[wc.condition],
        bai=lambda wc: ATAC_ATLAS_CONDITION_SHORT_BAIS[wc.condition]
    output:
        peaks=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/macs3/{{condition}}_peaks.narrowPeak",
        summits=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/macs3/{{condition}}_summits.bed",
        xls=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/macs3/{{condition}}_peaks.xls",
        treat_bdg=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/macs3/{{condition}}_treat_pileup.bdg",
        control_bdg=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/macs3/{{condition}}_control_lambda.bdg"
    params:
        command=atac_atlas_condition_callpeak_argv
    wildcard_constraints:
        condition=ATAC_ATLAS_CONDITION_RE
    resources:
        mem_mb=8000
    conda:
        "../envs/peaks.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/logs/macs3/callpeak.log"
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"


rule refine_atac_condition_cpm:
    input:
        peaks=lambda wc: ATAC_ATLAS_CONDITION_MACS_PEAKS[wc.condition],
        signal=lambda wc: ATAC_ATLAS_CONDITION_BIGWIGS[wc.condition],
        script=str(REPO_ROOT / "src" / "refine_atac_cpm_peaks.py"),
        implementation=str(
            REPO_ROOT / "src" / "short_read_processing" / "cpm_refinement.py"
        )
    output:
        refined=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/refined/{{condition}}.pooled.CPM-refined.bed",
        excluded=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/refined/{{condition}}.pooled.Excluded.bed",
        stats=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/refined/{{condition}}.pooled.stats.json"
    params:
        minimum_mean_cpm=float(ATAC_REFINEMENT["minimum_mean_cpm"]),
        minimum_mode_prominence=float(ATAC_REFINEMENT["minimum_mode_prominence"]),
        merge_gap_bp=int(ATAC_REFINEMENT["merge_gap_bp"]),
        minimum_length=int(ATAC_REFINEMENT["minimum_length"]),
        maximum_length=int(ATAC_REFINEMENT["maximum_length"])
    wildcard_constraints:
        condition=ATAC_ATLAS_CONDITION_RE
    resources:
        mem_mb=8000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/logs/refined/pooled.log"
    shell:
        "mkdir -p $(dirname {output.refined:q}) $(dirname {log:q}) && "
        "python {input.script:q} --peaks {input.peaks:q} "
        "--signal-bigwig {input.signal:q} --output {output.refined:q} "
        "--excluded {output.excluded:q} --stats {output.stats:q} "
        "--merge-gap-bp {params.merge_gap_bp} --minimum-length {params.minimum_length} "
        "--maximum-length {params.maximum_length} "
        "--minimum-mean-cpm {params.minimum_mean_cpm} "
        "--minimum-mode-prominence {params.minimum_mode_prominence} > {log:q} 2>&1"


rule filter_atac_condition_replicate_support:
    input:
        pooled=lambda wc: ATAC_ATLAS_CONDITION_POOLED_REFINED[wc.condition],
        replicates=lambda wc: [
            ATAC_REFINED_PEAKS[sample]
            for sample in ATAC_ATLAS_CONDITIONS[wc.condition].samples
        ]
    output:
        bed=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/{{condition}}.consensus.bed",
        support=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/{{condition}}.support.tsv",
        stats=f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/{{condition}}.consensus.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "build_atac_atlas.py"),
        replicates=atac_atlas_replicate_arguments,
        minimum_replicates=int(ATAC_ATLAS.get("minimum_replicates", 2)),
        overlap_fraction=float(ATAC_ATLAS.get("replicate_overlap_fraction", 0.5))
    wildcard_constraints:
        condition=ATAC_ATLAS_CONDITION_RE
    resources:
        mem_mb=4000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/conditions/{{condition}}/logs/consensus.log"
    shell:
        "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
        "python {params.script:q} condition --condition-id {wildcards.condition:q} "
        "--pooled-peaks {input.pooled:q} {params.replicates} "
        "--minimum-replicates {params.minimum_replicates} "
        "--overlap-fraction {params.overlap_fraction} --output-bed {output.bed:q} "
        "--support-tsv {output.support:q} --stats-json {output.stats:q} "
        "> {log:q} 2>&1"


rule build_cross_condition_atac_atlas:
    input:
        consensus=list(ATAC_ATLAS_CONSENSUS_PEAKS.values()),
        bigwigs=list(ATAC_ATLAS_CONDITION_BIGWIGS.values())
    output:
        bed=f"{ATAC_ATLAS_ROOT}/atlas.peaks.bed",
        variable_bed=f"{ATAC_ATLAS_ROOT}/atlas.variable.peaks.bed",
        membership=f"{ATAC_ATLAS_ROOT}/atlas.membership.tsv",
        presence=f"{ATAC_ATLAS_ROOT}/atlas.presence.tsv",
        coverage=f"{ATAC_ATLAS_ROOT}/atlas.coverage_fraction.tsv",
        mean_cpm=f"{ATAC_ATLAS_ROOT}/atlas.mean_cpm.tsv",
        maximum_cpm=f"{ATAC_ATLAS_ROOT}/atlas.maximum_cpm.tsv",
        stats=f"{ATAC_ATLAS_ROOT}/atlas.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "build_atac_atlas.py"),
        conditions=atac_atlas_condition_arguments,
        peak_width=int(ATAC_ATLAS.get("peak_width", 250))
    resources:
        mem_mb=12000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/logs/atlas.log"
    shell:
        "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
        "python {params.script:q} atlas {params.conditions} "
        "--grouping-method fixed_window --peak-width {params.peak_width} "
        "--output-bed {output.bed:q} "
        "--variable-bed {output.variable_bed:q} "
        "--membership-tsv {output.membership:q} --presence-tsv {output.presence:q} "
        "--coverage-tsv {output.coverage:q} --mean-cpm-tsv {output.mean_cpm:q} "
        "--maximum-cpm-tsv {output.maximum_cpm:q} --stats-json {output.stats:q} "
        "> {log:q} 2>&1"


rule build_narrow_first_cross_condition_atac_atlas:
    input:
        consensus=list(ATAC_ATLAS_CONSENSUS_PEAKS.values()),
        bigwigs=list(ATAC_ATLAS_CONDITION_BIGWIGS.values())
    output:
        anchors=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.anchors250.bed",
        variable_bed=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.variable.peaks.bed",
        membership=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.membership.tsv",
        presence=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.presence.tsv",
        coverage=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.coverage_fraction.tsv",
        mean_cpm=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.mean_cpm.tsv",
        maximum_cpm=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.maximum_cpm.tsv",
        stats=f"{ATAC_ATLAS_ROOT}/atlas.narrow-first.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "build_atac_atlas.py"),
        conditions=atac_atlas_condition_arguments,
        peak_width=int(ATAC_ATLAS.get("peak_width", 250))
    resources:
        mem_mb=12000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/logs/atlas.narrow-first.log"
    shell:
        "mkdir -p $(dirname {output.anchors:q}) $(dirname {log:q}) && "
        "python {params.script:q} atlas {params.conditions} "
        "--grouping-method fixed_window_narrow_first "
        "--peak-width {params.peak_width} --output-bed {output.anchors:q} "
        "--variable-bed {output.variable_bed:q} "
        "--membership-tsv {output.membership:q} --presence-tsv {output.presence:q} "
        "--coverage-tsv {output.coverage:q} --mean-cpm-tsv {output.mean_cpm:q} "
        "--maximum-cpm-tsv {output.maximum_cpm:q} --stats-json {output.stats:q} "
        "> {log:q} 2>&1"


rule build_atac_dhs_support_fwhm:
    input:
        anchors=f"{ATAC_ATLAS_ROOT}/atlas.peaks.bed",
        consensus=list(ATAC_ATLAS_CONSENSUS_PEAKS.values()),
        chrom_sizes=REFERENCE["chrom_sizes"]
    output:
        bigwig=f"{ATAC_ATLAS_ROOT}/atlas.dhs-support-fraction.bw",
        bed=f"{ATAC_ATLAS_ROOT}/atlas.fwhm-boundaries.bed",
        diagnostics=f"{ATAC_ATLAS_ROOT}/atlas.fwhm-diagnostics.tsv",
        stats=f"{ATAC_ATLAS_ROOT}/atlas.fwhm.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "build_atac_atlas.py"),
        conditions=atac_atlas_condition_dhs_arguments
    resources:
        mem_mb=8000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/logs/atlas.fwhm.log"
    shell:
        "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
        "python {params.script:q} fwhm --anchors-bed {input.anchors:q} "
        "{params.conditions} --chrom-sizes {input.chrom_sizes:q} "
        "--support-bigwig {output.bigwig:q} --output-bed {output.bed:q} "
        "--diagnostics-tsv {output.diagnostics:q} --stats-json {output.stats:q} "
        "> {log:q} 2>&1"


rule build_atac_dhs_center_mode_width:
    input:
        anchors=f"{ATAC_ATLAS_ROOT}/atlas.peaks.bed",
        consensus=list(ATAC_ATLAS_CONSENSUS_PEAKS.values()),
        chrom_sizes=REFERENCE["chrom_sizes"]
    output:
        bed=f"{ATAC_ATLAS_ROOT}/atlas.center-mode-half-prominence-boundaries.bed",
        diagnostics=f"{ATAC_ATLAS_ROOT}/atlas.center-mode-half-prominence-diagnostics.tsv",
        stats=f"{ATAC_ATLAS_ROOT}/atlas.center-mode-half-prominence.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "build_atac_atlas.py"),
        conditions=atac_atlas_condition_dhs_arguments
    resources:
        mem_mb=8000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/logs/atlas.center-mode-half-prominence.log"
    shell:
        "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
        "python {params.script:q} center-mode --anchors-bed {input.anchors:q} "
        "{params.conditions} --chrom-sizes {input.chrom_sizes:q} "
        "--output-bed {output.bed:q} --diagnostics-tsv {output.diagnostics:q} "
        "--stats-json {output.stats:q} > {log:q} 2>&1"


rule build_dhs_driven_atac_atlas:
    input:
        consensus=list(ATAC_ATLAS_CONSENSUS_PEAKS.values()),
        bigwigs=list(ATAC_ATLAS_CONDITION_BIGWIGS.values())
    output:
        anchors=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.anchors250.bed",
        peaks=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.peaks.bed",
        membership=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.membership.tsv",
        presence=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.presence.tsv",
        coverage=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.coverage_fraction.tsv",
        mean_cpm=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.mean_cpm.tsv",
        maximum_cpm=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.maximum_cpm.tsv",
        stats=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "build_atac_atlas.py"),
        conditions=atac_atlas_condition_arguments,
        peak_width=int(ATAC_ATLAS.get("peak_width", 250))
    resources:
        mem_mb=12000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/logs/atlas.dhs-driven.log"
    shell:
        "mkdir -p $(dirname {output.anchors:q}) $(dirname {log:q}) && "
        "python {params.script:q} atlas {params.conditions} "
        "--grouping-method dhs_seed --peak-width {params.peak_width} "
        "--output-bed {output.anchors:q} --variable-bed {output.peaks:q} "
        "--membership-tsv {output.membership:q} --presence-tsv {output.presence:q} "
        "--coverage-tsv {output.coverage:q} --mean-cpm-tsv {output.mean_cpm:q} "
        "--maximum-cpm-tsv {output.maximum_cpm:q} --stats-json {output.stats:q} "
        "> {log:q} 2>&1"


rule shape_dhs_driven_atac_atlas:
    input:
        membership=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.membership.tsv",
        bigwigs=list(ATAC_ATLAS_CONDITION_BIGWIGS.values())
    output:
        bed=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.signal-shaped.peaks.bed",
        bigwig=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.aggregate-shape.bw",
        diagnostics=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.signal-shape.tsv",
        stats=f"{ATAC_ATLAS_ROOT}/atlas.dhs-driven.signal-shape.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "build_atac_atlas.py"),
        bigwigs=atac_atlas_condition_bigwig_arguments,
        minimum_length=int(ATAC_REFINEMENT["minimum_length"]),
        maximum_length=int(ATAC_REFINEMENT["maximum_length"])
    resources:
        mem_mb=12000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_ATLAS_ROOT}/logs/atlas.dhs-driven.signal-shape.log"
    shell:
        "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
        "python {params.script:q} shape --membership-tsv {input.membership:q} "
        "{params.bigwigs} --output-bed {output.bed:q} "
        "--aggregate-bigwig {output.bigwig:q} "
        "--diagnostics-tsv {output.diagnostics:q} --stats-json {output.stats:q} "
        "--window-size 1000 --bin-size 10 --smoothing-bins 3 "
        "--relative-threshold 0.2 --background-mad-multiplier 3 "
        "--minimum-length {params.minimum_length} "
        "--maximum-length {params.maximum_length} > {log:q} 2>&1"
