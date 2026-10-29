if MODE != "annotated":
    rule annotate_fibers:
        input:
            bam=downstream_aligned_bam
        output:
            bam=f"{RESULTS}/03_annotate/{{sample}}.fiber.bam"
        params:
            command="predict-m6a" if MODE == "hifi" else "add-nucleosomes"
        threads: config["threads"]["predict_m6a"]
        log:
            out=f"{LOGS}/03_annotate/{{sample}}.out",
            err=f"{LOGS}/03_annotate/{{sample}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log.out:q})"
            ft -t {threads} {params.command} {input.bam:q} {output.bam:q} \
                > {log.out:q} 2> {log.err:q}
            """
