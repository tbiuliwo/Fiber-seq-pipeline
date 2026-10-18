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
