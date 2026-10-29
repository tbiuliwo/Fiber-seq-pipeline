# 09. FIBERHMM
if FIBERHMM:

    rule call_footprints:
        input:
            bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/05_fire/{{sample}}.fire.bam.bai",
        output:
            bam=f"{RESULTS}/09_fiberhmm/{{sample}}.footprints.bam",
            qc_json=f"{RESULTS}/09_fiberhmm/qc/{{sample}}.footprints.qc.json",
            qc_tsv=f"{RESULTS}/09_fiberhmm/qc/{{sample}}.footprints.qc.tsv",
            qc_png=f"{RESULTS}/09_fiberhmm/qc/{{sample}}.footprints.qc.png",
            qc_pdf=f"{RESULTS}/09_fiberhmm/qc/{{sample}}.footprints.qc.pdf",
        threads: config["threads"]["fiberhmm"]
        params:
            qc_prefix=f"{RESULTS}/09_fiberhmm/qc/{{sample}}.footprints",
        log:
            out=f"{LOGS}/09_fiberhmm/{{sample}}.out",
            err=f"{LOGS}/09_fiberhmm/{{sample}}.err",
        shell:
            r"""
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log.out:q})"
            fiberhmm-call -i {input.bam:q} -o {output.bam:q} --enzyme hia5 --seq pacbio \
                --cores {threads} --io-threads 1 --qc --qc-output-prefix {params.qc_prefix:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule index_footprints:
        input:
            bam=f"{RESULTS}/09_fiberhmm/{{sample}}.footprints.bam",
        output:
            bai=f"{RESULTS}/09_fiberhmm/{{sample}}.footprints.bam.bai",
        log:
            out=f"{LOGS}/09_fiberhmm/{{sample}}.index.out",
            err=f"{LOGS}/09_fiberhmm/{{sample}}.index.err",
        shell:
            "samtools index {input.bam:q} {output.bai:q} > {log.out:q} 2> {log.err:q}"
