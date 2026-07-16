rule align_and_mark_duplicates:
    input:
        reads=lambda wc: TRIMMED_READS_BY_SAMPLE[wc.sample],
        index=BT2_INDEX
    output:
        bam=f"{WORK_ROOT}/alignment/{{sample}}.marked.bam"
    log:
        f"{RESULT_ROOT}/logs/alignment/{{sample}}.bowtie2.log"
    params:
        reads=bowtie_read_arguments,
        layout=bowtie_layout_arguments,
        preset=bowtie_preset,
        index=lambda wc, input: str(input.index[0]).removesuffix(".1.bt2")
    threads: 8
    resources:
        mem_mb=8000
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        bowtie2 {params.preset} {params.layout} -x {params.index:q} {params.reads} \
          --rg-id {wildcards.sample:q} --rg SM:{wildcards.sample:q} -p {threads} 2> {log:q} \
          | samtools view -u - \
          | samtools sort -n -@ {threads} -o {output.bam:q}.namesort.bam
        samtools fixmate -@ {threads} -m {output.bam:q}.namesort.bam {output.bam:q}.fixmate.bam
        samtools sort -@ {threads} -o {output.bam:q}.coordsort.bam {output.bam:q}.fixmate.bam
        samtools markdup -@ {threads} {output.bam:q}.coordsort.bam {output.bam:q}
        rm -f {output.bam:q}.namesort.bam {output.bam:q}.fixmate.bam {output.bam:q}.coordsort.bam
        """


rule filter_bam:
    input:
        bam=f"{WORK_ROOT}/alignment/{{sample}}.marked.bam",
        blacklist=str(REFERENCE["blacklist_bed"])
    output:
        bam=f"{RESULT_ROOT}/bam/{{sample}}.final.bam",
        bai=f"{RESULT_ROOT}/bam/{{sample}}.final.bam.bai"
    params:
        required=required_flags,
        excluded=excluded_flags,
        mapq=lambda wc: SAMPLES[wc.sample]["parameters"]["alignment"]["mapq_minimum"],
        mitochondrial=str(REFERENCE["mitochondrial_contig"]),
        remove_mito=lambda wc: int(
            SAMPLES[wc.sample]["parameters"]["filtering"]["remove_mitochondrial"]
        )
    threads: 6
    resources:
        mem_mb=6000
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/alignment/{{sample}}.filter.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        samtools view -@ {threads} -h -q {params.mapq} {params.required} -F {params.excluded} {input.bam:q} \
          | awk -v mt={params.mitochondrial:q} -v remove_mt={params.remove_mito} 'BEGIN{{OFS="\t"}} /^@/ {{print; next}} remove_mt == 0 || $3 != mt {{print}}' \
          | samtools view -u - \
          | bedtools intersect -v -abam stdin -b {input.blacklist:q} \
          | samtools sort -@ {threads} -o {output.bam:q} - 2> {log:q}
        samtools index -@ {threads} {output.bam:q} {output.bai:q}
        samtools quickcheck -v {output.bam:q}
        """


rule alignment_stats:
    input:
        bam=f"{RESULT_ROOT}/bam/{{sample}}.final.bam"
    output:
        flagstat=f"{RESULT_ROOT}/qc/alignment/{{sample}}.flagstat.txt",
        stats=f"{RESULT_ROOT}/qc/alignment/{{sample}}.stats.txt",
        idxstats=f"{RESULT_ROOT}/qc/alignment/{{sample}}.idxstats.txt"
    threads: 2
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULT_ROOT}/logs/alignment/{{sample}}.stats.log"
    shell:
        "mkdir -p $(dirname {output.flagstat:q}) $(dirname {log:q}) && "
        "samtools flagstat -@ {threads} {input.bam:q} > {output.flagstat:q} 2> {log:q} && "
        "samtools stats -@ {threads} {input.bam:q} > {output.stats:q} 2>> {log:q} && "
        "samtools idxstats {input.bam:q} > {output.idxstats:q} 2>> {log:q}"
