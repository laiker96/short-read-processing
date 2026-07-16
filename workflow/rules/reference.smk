if AUTO_PREPARE_REFERENCE:
    rule download_reference_fasta:
        output:
            archive=REFERENCE_FASTA_ARCHIVE
        params:
            url=REFERENCE_PREPARATION["fasta"]["url"],
            checksum=aria2_checksum(REFERENCE_PREPARATION["fasta"]),
            directory=lambda wc, output: str(Path(output.archive).parent),
            filename=lambda wc, output: Path(output.archive).name
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/reference/download-fasta.log"
        shell:
            "mkdir -p {params.directory:q} $(dirname {log:q}) && "
            "aria2c --continue=true --allow-overwrite=true --auto-file-renaming=false "
            "--file-allocation=none --max-connection-per-server=8 --split=8 "
            "--min-split-size=5M --check-integrity=true --checksum={params.checksum:q} "
            "--dir={params.directory:q} --out={params.filename:q} {params.url:q} "
            "> {log:q} 2>&1"


    rule prepare_reference_fasta:
        input:
            archive=REFERENCE_FASTA_ARCHIVE
        output:
            fasta=str(REFERENCE["fasta"])
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/reference/prepare-fasta.log"
        shell:
            "mkdir -p $(dirname {output.fasta:q}) $(dirname {log:q}) && "
            "pigz -dc {input.archive:q} > {output.fasta:q}.tmp 2> {log:q} && "
            "mv {output.fasta:q}.tmp {output.fasta:q}"


    rule download_reference_annotation:
        output:
            archive=REFERENCE_ANNOTATION_ARCHIVE
        params:
            url=REFERENCE_PREPARATION["annotation"]["url"],
            checksum=aria2_checksum(REFERENCE_PREPARATION["annotation"]),
            directory=lambda wc, output: str(Path(output.archive).parent),
            filename=lambda wc, output: Path(output.archive).name
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/reference/download-annotation.log"
        shell:
            "mkdir -p {params.directory:q} $(dirname {log:q}) && "
            "aria2c --continue=true --allow-overwrite=true --auto-file-renaming=false "
            "--file-allocation=none --max-connection-per-server=8 --split=8 "
            "--min-split-size=5M --check-integrity=true --checksum={params.checksum:q} "
            "--dir={params.directory:q} --out={params.filename:q} {params.url:q} "
            "> {log:q} 2>&1"


    rule download_reference_blacklist:
        output:
            archive=REFERENCE_BLACKLIST_ARCHIVE
        params:
            url=REFERENCE_PREPARATION["blacklist"]["url"],
            checksum=aria2_checksum(REFERENCE_PREPARATION["blacklist"]),
            directory=lambda wc, output: str(Path(output.archive).parent),
            filename=lambda wc, output: Path(output.archive).name
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/reference/download-blacklist.log"
        shell:
            "mkdir -p {params.directory:q} $(dirname {log:q}) && "
            "aria2c --continue=true --allow-overwrite=true --auto-file-renaming=false "
            "--file-allocation=none --max-connection-per-server=8 --split=8 "
            "--min-split-size=1M --check-integrity=true --checksum={params.checksum:q} "
            "--dir={params.directory:q} --out={params.filename:q} {params.url:q} "
            "> {log:q} 2>&1"


    rule prepare_reference_blacklist:
        input:
            archive=REFERENCE_BLACKLIST_ARCHIVE
        output:
            bed=str(REFERENCE["blacklist_bed"])
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/reference/prepare-blacklist.log"
        shell:
            "mkdir -p $(dirname {output.bed:q}) $(dirname {log:q}) && "
            "pigz -dc {input.archive:q} > {output.bed:q}.tmp 2> {log:q} && "
            "mv {output.bed:q}.tmp {output.bed:q}"


    rule prepare_reference_tss:
        input:
            annotation=REFERENCE_ANNOTATION_ARCHIVE,
            fai=FASTA_INDEX
        output:
            bed=str(REFERENCE["tss_bed"])
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/reference/prepare-tss.log"
        script:
            "../scripts/extract_tss.py"


    rule prepare_reference_autosomes:
        input:
            fai=FASTA_INDEX
        output:
            contigs=str(REFERENCE["autosomes_file"])
        params:
            autosomes=REFERENCE_PREPARATION["autosomes"]
        conda:
            "../envs/reference.yaml"
        log:
            f"{RESULT_ROOT}/logs/reference/prepare-autosomes.log"
        script:
            "../scripts/write_autosomes.py"


rule fasta_index:
    input:
        fasta=str(REFERENCE["fasta"])
    output:
        fai=FASTA_INDEX
    threads: 2
    conda:
        "../envs/reference.yaml"
    log:
        f"{RESULT_ROOT}/logs/reference/samtools-faidx.log"
    shell:
        "mkdir -p $(dirname {log:q}) && samtools faidx {input.fasta:q} > {log:q} 2>&1"


rule chromosome_sizes:
    input:
        fai=FASTA_INDEX
    output:
        sizes=str(REFERENCE["chrom_sizes"])
    conda:
        "../envs/reference.yaml"
    log:
        f"{RESULT_ROOT}/logs/reference/chromosome-sizes.log"
    shell:
        "mkdir -p $(dirname {output.sizes:q}) $(dirname {log:q}) && "
        "cut -f1,2 {input.fai:q} > {output.sizes:q} 2> {log:q}"


rule bowtie2_index:
    input:
        fasta=str(REFERENCE["fasta"])
    output:
        BT2_INDEX
    threads: 8
    resources:
        mem_mb=8000
    params:
        prefix=lambda wc, output: str(output[0]).removesuffix(".1.bt2")
    conda:
        "../envs/reference.yaml"
    log:
        f"{RESULT_ROOT}/logs/reference/bowtie2-build.log"
    shell:
        "mkdir -p $(dirname {params.prefix:q}) $(dirname {log:q}) && "
        "bowtie2-build --threads {threads} {input.fasta:q} {params.prefix:q} > {log:q} 2>&1"
