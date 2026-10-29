# 05. Single-molecule FIRE calling

rule call_fire:
    input:
        bam=annotated_bam
    output:
        bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam"
    threads:
        config["threads"]["fire"]
    log:
        out=f"{LOGS}/05_fire/{{sample}}.out",
        err=f"{LOGS}/05_fire/{{sample}}.err"
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

rule extract_fire_bed:
    input:
        bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam"
    output:
        bed=f"{RESULTS}/05_fire/{{sample}}.fire.bed"
    threads: config["threads"]["fire"]
    log:
        out=f"{LOGS}/05_fire/{{sample}}.extract.out",
        err=f"{LOGS}/05_fire/{{sample}}.extract.err"
    shell:
        r"""
        set -euo pipefail

        ft -t {threads} fire \
            --extract \
            {input.bam:q} \
            {output.bed:q} \
            > {log.out:q} \
            2> {log.err:q}
        """
rule extract_fire_features:
    input:
        bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam"
    output:
        tsv=f"{RESULTS}/05_fire/{{sample}}.fire.features.tsv"
    threads: config["threads"]["fire"]
    log:
        out=f"{LOGS}/05_fire/{{sample}}.features.out",
        err=f"{LOGS}/05_fire/{{sample}}.features.err"
    shell:
        r"""
        set -euo pipefail

        ft -t {threads} fire \
            --feats-to-text \
            {input.bam:q} \
            {output.tsv:q} \
            > {log.out:q} \
            2> {log.err:q}
        """
        
rule index_fire_bam:
    input:
        bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam"
    output:
        bai=f"{RESULTS}/05_fire/{{sample}}.fire.bam.bai"
    threads: config["threads"]["qc"]
    log:
        out=f"{LOGS}/05_fire/{{sample}}.index.out",
        err=f"{LOGS}/05_fire/{{sample}}.index.err"
    shell:
        "samtools index -@ {threads} {input.bam:q} {output.bai:q} > {log.out:q} 2> {log.err:q}"
