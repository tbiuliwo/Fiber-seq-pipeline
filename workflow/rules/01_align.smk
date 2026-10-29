if not ALIGNED:
    if not str(REFERENCE).endswith(".mmi"):
        rule index_reference:
            input:
                fasta=REFERENCE
            output:
                mmi=ALIGN_REFERENCE
            threads: config["threads"]["align"]
            log:
                out=f"{LOGS}/01_align/reference/index.out",
                err=f"{LOGS}/01_align/reference/index.err"
            shell:
                r"""
                mkdir -p "$(dirname {output.mmi:q})" "$(dirname {log.out:q})"
                pbmm2 index --preset CCS -j {threads} {input.fasta:q} {output.mmi:q} \
                    > {log.out:q} 2> {log.err:q}
                """

    rule align_hifi:
        input:
            bam=sample_bam,
            reference=ALIGN_REFERENCE
        output:
            bam=f"{RESULTS}/01_align/{{sample}}.aligned.bam"
        threads: config["threads"]["align"]
        log:
            out=f"{LOGS}/01_align/{{sample}}.align.out",
            err=f"{LOGS}/01_align/{{sample}}.align.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log.out:q})"
            pbmm2 align --preset HIFI --sort -j {threads} \
                {input.reference:q} {input.bam:q} {output.bam:q} \
                > {log.out:q} 2> {log.err:q}
            """
else:
    rule sort_aligned:
        input:
            bam=sample_bam,
        output:
            bam=f"{RESULTS}/01_align/{{sample}}.aligned.bam"
        threads: config["threads"]["align"]
        log:
            out=f"{LOGS}/01_align/{{sample}}.sort.out",
            err=f"{LOGS}/01_align/{{sample}}.sort.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log.out:q})"
            samtools sort -@ {threads} -o {output.bam:q} {input.bam:q} \
                > {log.out:q} 2> {log.err:q}
            """


rule filter_aligned_chromosomes:
    input:
        bam=f"{RESULTS}/01_align/{{sample}}.aligned.bam",
        script=str(Path(workflow.basedir).parent / "scripts/01filter_chromosomes.py")
    output:
        bam=f"{RESULTS}/01_align/{{sample}}.chromosomes.bam"
    params:
        chromosomes=ANALYSIS_CHROMOSOMES
    threads: config["threads"]["align"]
    log:
        out=f"{LOGS}/01_align/{{sample}}.chromosomes.out",
        err=f"{LOGS}/01_align/{{sample}}.chromosomes.err"
    shell:
        r"""
        mkdir -p "$(dirname {output.bam:q})" "$(dirname {log.out:q})"
        {PYTHON:q} {input.script:q} --bam {input.bam:q} --output {output.bam:q} \
            --threads {threads} --chromosomes {params.chromosomes:q} \
            > {log.out:q} 2> {log.err:q}
        """


rule index_aligned_chromosomes:
    input:
        bam=f"{RESULTS}/01_align/{{sample}}.chromosomes.bam"
    output:
        bai=f"{RESULTS}/01_align/{{sample}}.chromosomes.bam.bai"
    params:
        samtools=str(Path(config.get("tools", {}).get("bin_dir") or "") / "samtools")
    threads: config["threads"]["align"]
    log:
        out=f"{LOGS}/01_align/{{sample}}.chromosomes_index.out",
        err=f"{LOGS}/01_align/{{sample}}.chromosomes_index.err"
    shell:
        r"""
        mkdir -p "$(dirname {log.out:q})"
        {params.samtools:q} index -@ {threads} {input.bam:q} {output.bai:q} \
            > {log.out:q} 2> {log.err:q}
        """
