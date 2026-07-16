rule callpeak_narrow:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        control=peak_control_input
    output:
        peaks=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.narrowPeak",
        summits=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_summits.bed",
        xls=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.xls",
        treat_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_treat_pileup.bdg",
        control_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_control_lambda.bdg"
    params:
        command=callpeak_argv
    wildcard_constraints:
        sample=CALLPEAK_NARROW_RE
    conda:
        "../envs/peaks.yaml"
    resources:
        mem_mb=8000
    log:
        f"{RESULT_ROOT}/logs/peaks/{{sample}}.callpeak.log"
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"


rule callpeak_broad:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        control=peak_control_input
    output:
        peaks=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.broadPeak",
        gapped=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.gappedPeak",
        xls=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.xls",
        treat_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_treat_pileup.bdg",
        control_bdg=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_control_lambda.bdg"
    params:
        command=callpeak_argv
    wildcard_constraints:
        sample=CALLPEAK_BROAD_RE
    conda:
        "../envs/peaks.yaml"
    resources:
        mem_mb=8000
    log:
        f"{RESULT_ROOT}/logs/peaks/{{sample}}.callpeak.log"
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"


rule hmmratac:
    input:
        bam=lambda wc: FINAL_BAMS[wc.sample],
        blacklist=str(REFERENCE["blacklist_bed"])
    output:
        peaks=f"{RESULT_ROOT}/peaks/{{sample}}/{{sample}}_peaks.narrowPeak"
    params:
        lower=lambda wc: SAMPLES[wc.sample]["peak_caller"]["lower"],
        upper=lambda wc: SAMPLES[wc.sample]["peak_caller"]["upper"],
        prescan=lambda wc: SAMPLES[wc.sample]["peak_caller"]["prescan_cutoff"],
        outdir=lambda wc: f"{RESULT_ROOT}/peaks/{wc.sample}"
    wildcard_constraints:
        sample=HMMRATAC_RE
    conda:
        "../envs/peaks.yaml"
    resources:
        mem_mb=12000
    log:
        f"{RESULT_ROOT}/logs/peaks/{{sample}}.hmmratac.log"
    shell:
        "mkdir -p {params.outdir:q} $(dirname {log:q}) && "
        "macs3 hmmratac -i {input.bam:q} -f BAMPE -n {wildcards.sample:q} "
        "--outdir {params.outdir:q} -l {params.lower} -u {params.upper} "
        "-c {params.prescan} -e {input.blacklist:q} > {log:q} 2>&1"
