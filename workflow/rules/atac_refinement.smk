"""Paired-end ATAC short-fragment candidate calling and CPM refinement."""


rule filter_atac_short_fragments:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        bai=lambda wc: FINAL_BAIS[wc.sample]
    output:
        bam=f"{ATAC_REFINEMENT_ROOT}/bam/{{sample}}.fragments-lt{ATAC_FRAGMENT_MAXIMUM}.bam",
        bai=f"{ATAC_REFINEMENT_ROOT}/bam/{{sample}}.fragments-lt{ATAC_FRAGMENT_MAXIMUM}.bam.bai"
    params:
        maximum=ATAC_FRAGMENT_MAXIMUM
    wildcard_constraints:
        sample=ATAC_REFINEMENT_RE
    threads: 4
    resources:
        mem_mb=4000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{ATAC_REFINEMENT_ROOT}/logs/bam/{{sample}}.fragment-filter.log"
    shell:
        "mkdir -p $(dirname {output.bam:q}) $(dirname {log:q}) && "
        "samtools view -@ {threads} -b -f 2 -F 3852 "
        "-e 'tlen != 0 && tlen > -{params.maximum} && tlen < {params.maximum}' "
        "-o {output.bam:q} {input.bam:q} 2> {log:q} && "
        "samtools index -@ {threads} {output.bam:q} {output.bai:q} 2>> {log:q} && "
        "samtools quickcheck -v {output.bam:q} 2>> {log:q}"


rule atac_short_fragment_metrics:
    input:
        source=lambda wc: FINAL_BAMS[wc.sample],
        short=lambda wc: ATAC_SHORT_BAMS[wc.sample]
    output:
        metrics=f"{ATAC_REFINEMENT_ROOT}/qc/{{sample}}.fragment-filter.tsv"
    params:
        maximum=ATAC_FRAGMENT_MAXIMUM
    wildcard_constraints:
        sample=ATAC_REFINEMENT_RE
    threads: 2
    resources:
        mem_mb=2000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{ATAC_REFINEMENT_ROOT}/logs/qc/{{sample}}.fragment-filter.log"
    shell:
        r"""
        mkdir -p $(dirname {output.metrics:q}) $(dirname {log:q})
        source_count=$(samtools view -@ {threads} -c -f 66 -F 3852 {input.source:q} 2> {log:q})
        short_count=$(samtools view -@ {threads} -c -f 66 -F 3852 {input.short:q} 2>> {log:q})
        printf 'metric\tvalue\n' > {output.metrics:q}
        printf 'maximum_template_length_exclusive\t%s\n' {params.maximum} >> {output.metrics:q}
        printf 'source_proper_fragments\t%s\n' "$source_count" >> {output.metrics:q}
        printf 'retained_short_fragments\t%s\n' "$short_count" >> {output.metrics:q}
        awk -v kept="$short_count" -v total="$source_count" \
          'BEGIN {{print "retained_fraction\t" (total ? kept/total : 0)}}' >> {output.metrics:q}
        """


rule shift_atac_short_fragments:
    input:
        bam=lambda wc: ATAC_SHORT_BAMS[wc.sample]
    output:
        bam=temp(f"{WORK_ROOT}/atac_short_fragments/{{sample}}.shifted.bam"),
        bai=temp(f"{WORK_ROOT}/atac_short_fragments/{{sample}}.shifted.bam.bai")
    wildcard_constraints:
        sample=ATAC_REFINEMENT_RE
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_REFINEMENT_ROOT}/logs/tracks/{{sample}}.alignmentSieve.log"
    shell:
        "mkdir -p $(dirname {output.bam:q}) $(dirname {log:q}) && "
        "alignmentSieve --ATACshift -b {input.bam:q} -o {output.bam:q}.unsorted "
        "--numberOfProcessors {threads} > {log:q} 2>&1 && "
        "samtools sort -@ {threads} -o {output.bam:q} {output.bam:q}.unsorted "
        "2>> {log:q} && rm -f {output.bam:q}.unsorted && "
        "samtools index -@ {threads} {output.bam:q} {output.bai:q} 2>> {log:q} && "
        "samtools quickcheck -v {output.bam:q} 2>> {log:q}"


rule atac_short_fragment_bigwig:
    input:
        bam=f"{WORK_ROOT}/atac_short_fragments/{{sample}}.shifted.bam",
        bai=f"{WORK_ROOT}/atac_short_fragments/{{sample}}.shifted.bam.bai"
    output:
        bigwig=(
            f"{ATAC_REFINEMENT_ROOT}/tracks/{{sample}}.fragments-lt{ATAC_FRAGMENT_MAXIMUM}"
            ".Tn5-shifted.CPM.bw"
        )
    params:
        bin_size=int(ATAC_REFINEMENT["bigwig_bin_size"])
    wildcard_constraints:
        sample=ATAC_REFINEMENT_RE
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_REFINEMENT_ROOT}/logs/tracks/{{sample}}.bamCoverage.log"
    shell:
        "mkdir -p $(dirname {output.bigwig:q}) $(dirname {log:q}) && "
        "bamCoverage -b {input.bam:q} -o {output.bigwig:q} "
        "--outFileFormat bigwig --normalizeUsing CPM --binSize {params.bin_size} "
        "--numberOfProcessors {threads} > {log:q} 2>&1"


rule call_lenient_atac_short_fragment_peaks:
    input:
        bam=lambda wc: ATAC_SHORT_BAMS[wc.sample],
        bai=lambda wc: ATAC_SHORT_BAIS[wc.sample]
    output:
        peaks=f"{ATAC_REFINEMENT_ROOT}/macs3/{{sample}}/{{sample}}_peaks.narrowPeak",
        summits=f"{ATAC_REFINEMENT_ROOT}/macs3/{{sample}}/{{sample}}_summits.bed",
        xls=f"{ATAC_REFINEMENT_ROOT}/macs3/{{sample}}/{{sample}}_peaks.xls",
        treat_bdg=f"{ATAC_REFINEMENT_ROOT}/macs3/{{sample}}/{{sample}}_treat_pileup.bdg",
        control_bdg=f"{ATAC_REFINEMENT_ROOT}/macs3/{{sample}}/{{sample}}_control_lambda.bdg"
    params:
        command=atac_refinement_callpeak_argv
    wildcard_constraints:
        sample=ATAC_REFINEMENT_RE
    resources:
        mem_mb=8000
    conda:
        "../envs/peaks.yaml"
    log:
        f"{ATAC_REFINEMENT_ROOT}/logs/macs3/{{sample}}.callpeak.log"
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"


rule refine_atac_short_fragment_cpm:
    input:
        peaks=lambda wc: ATAC_LENIENT_PEAKS[wc.sample],
        signal=lambda wc: ATAC_SHORT_BIGWIGS[wc.sample]
    output:
        refined=f"{ATAC_REFINEMENT_ROOT}/refined/{{sample}}.CPM-refined.bed",
        excluded=f"{ATAC_REFINEMENT_ROOT}/refined/{{sample}}.Excluded.bed",
        stats=f"{ATAC_REFINEMENT_ROOT}/refined/{{sample}}.stats.json"
    params:
        script=str(REPO_ROOT / "src" / "refine_atac_cpm_peaks.py"),
        minimum_mean_cpm=float(ATAC_REFINEMENT["minimum_mean_cpm"]),
        merge_gap_bp=int(ATAC_REFINEMENT["merge_gap_bp"]),
        minimum_length=int(ATAC_REFINEMENT["minimum_length"]),
        maximum_length=int(ATAC_REFINEMENT["maximum_length"])
    wildcard_constraints:
        sample=ATAC_REFINEMENT_RE
    resources:
        mem_mb=8000
    conda:
        "../envs/atac_qc.yaml"
    log:
        f"{ATAC_REFINEMENT_ROOT}/logs/refined/{{sample}}.log"
    shell:
        "mkdir -p $(dirname {output.refined:q}) $(dirname {log:q}) && "
        "python {params.script:q} --peaks {input.peaks:q} "
        "--signal-bigwig {input.signal:q} --output {output.refined:q} "
        "--excluded {output.excluded:q} --stats {output.stats:q} "
        "--merge-gap-bp {params.merge_gap_bp} --minimum-length {params.minimum_length} "
        "--maximum-length {params.maximum_length} "
        "--minimum-mean-cpm {params.minimum_mean_cpm} > {log:q} 2>&1"
