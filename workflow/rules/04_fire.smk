# 04. Single-molecule FIRE calling

rule call_fire:
    input:
        bam=annotated_bam
    output:
        bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam"
    threads:
        config["threads"]["fire"]
    log:
        out=f"{LOGS}/04_fire/{{sample}}.out",
        err=f"{LOGS}/04_fire/{{sample}}.err"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.bam:q})"
        mkdir -p "$(dirname {log.out:q})"

        ft -t {threads} fire \
            {input.bam:q} \
            {output.bam:q} \
            > {log.out:q} \
            2> {log.err:q}
        """

rule index_fire_bam:
    input:
        bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam"
    output:
        bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai"
    threads: config["threads"]["qc"]
    log:
        out=f"{LOGS}/04_fire/{{sample}}.index.out",
        err=f"{LOGS}/04_fire/{{sample}}.index.err"
    shell:
        "samtools index -@ {threads} {input.bam:q} {output.bai:q} > {log.out:q} 2> {log.err:q}"
