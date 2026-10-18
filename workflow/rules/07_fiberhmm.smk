# 7. FIBERHMM
if FIBERHMM:
    rule call_footprints:
        input:
            bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai"
        output:
            bam=f"{RESULTS}/07_fiberhmm/{{sample}}.footprints.bam"
        threads: config["threads"]["fiberhmm"]
        log:
            out=f"{LOGS}/07_fiberhmm/{{sample}}.out",
            err=f"{LOGS}/07_fiberhmm/{{sample}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log.out:q})"
            fiberhmm-call -i {input.bam:q} -o {output.bam:q} --enzyme hia5 --seq pacbio \
                --cores {threads} --io-threads 1 --no-qc \
                > {log.out:q} 2> {log.err:q}
            """

    rule index_footprints:
        input:
            bam=f"{RESULTS}/07_fiberhmm/{{sample}}.footprints.bam"
        output:
            bai=f"{RESULTS}/07_fiberhmm/{{sample}}.footprints.bam.bai"
        log:
            out=f"{LOGS}/07_fiberhmm/{{sample}}.index.out",
            err=f"{LOGS}/07_fiberhmm/{{sample}}.index.err"
        shell:
            "samtools index {input.bam:q} {output.bai:q} > {log.out:q} 2> {log.err:q}"
