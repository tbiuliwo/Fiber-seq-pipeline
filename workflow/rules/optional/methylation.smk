# Optional DNA methylation analysis (5mC / 5hmC)

rule methylation_pileup:
    input:
        bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
        bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai"

    output:
        mc=f"{RESULTS}/optional/methylation/{{sample}}.5mC.bed.gz",
        mc_idx=f"{RESULTS}/optional/methylation/{{sample}}.5mC.bed.gz.tbi",
        hmc=f"{RESULTS}/optional/methylation/{{sample}}.5hmC.bed.gz",
        hmc_idx=f"{RESULTS}/optional/methylation/{{sample}}.5hmC.bed.gz.tbi",
        summary=f"{RESULTS}/optional/methylation/{{sample}}.summary.json"

    threads:
        config["threads"]["methylation"]

    params:
        prefix=lambda wc: f"{RESULTS}/optional/methylation/{wc.sample}"

    log:
        out=f"{LOGS}/optional/methylation/{{sample}}.out",
        err=f"{LOGS}/optional/methylation/{{sample}}.err"

    shell:
        r"""
        mkdir -p {RESULTS}/optional/methylation
        mkdir -p {LOGS}/optional/methylation

        methbat pileup \
            --threads {threads} \
            --input-bam {input.bam:q} \
            --output-prefix {params.prefix:q} \
            --skip-6mA \
            > {log.out:q} 2> {log.err:q}
        """
