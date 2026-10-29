# 04. Basic Fiber-seq QC

rule fiber_qc:
    input:
        bam=annotated_bam
    output:
        table=f"{RESULTS}/04_qc/{{sample}}.qc.tsv"
    threads:
        config["threads"]["qc"]
    params:
        acf="--acf" if config.get("qc", {}).get("autocorrelation", False) else ""
    log:
        out=f"{LOGS}/04_qc/{{sample}}.out",
        err=f"{LOGS}/04_qc/{{sample}}.err"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.table:q})"
        mkdir -p "$(dirname {log.out:q})"

        ft -t {threads} qc \
            {params.acf} \
            {input.bam:q} \
            {output.table:q} \
            > {log.out:q} \
            2> {log.err:q}
        """

rule alignment_qc:
    input:
        bam=aligned_bam
    output:
        flagstat=f"{RESULTS}/04_qc/{{sample}}.flagstat.txt",
        stats=f"{RESULTS}/04_qc/{{sample}}.alignment.stats.txt"
    threads: config["threads"]["qc"]
    log:
        out=f"{LOGS}/04_qc/{{sample}}.alignment.out",
        err=f"{LOGS}/04_qc/{{sample}}.alignment.err"
    shell:
        r"""
        mkdir -p "$(dirname {output.stats:q})" "$(dirname {log.out:q})"
        samtools flagstat -@ {threads} {input.bam:q} > {output.flagstat:q} 2> {log.err:q}
        samtools stats -@ {threads} {input.bam:q} > {output.stats:q} 2>> {log.err:q}
        echo "Alignment QC complete" > {log.out:q}
        """

rule alignment_pre_filter_qc:
    input:
        bam=f"{RESULTS}/01_align/{{sample}}.aligned.bam"
    output:
        flagstat=f"{RESULTS}/04_qc/{{sample}}.pre_filter.flagstat.txt",
        stats=f"{RESULTS}/04_qc/{{sample}}.pre_filter.stats.txt"
    threads: config["threads"]["qc"]
    log:
        out=f"{LOGS}/04_qc/{{sample}}.pre_filter.out",
        err=f"{LOGS}/04_qc/{{sample}}.pre_filter.err"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.stats:q})" "$(dirname {log.out:q})"

        samtools flagstat -@ {threads} {input.bam:q} \
            > {output.flagstat:q} 2> {log.err:q}

        samtools stats -@ {threads} {input.bam:q} \
            > {output.stats:q} 2>> {log.err:q}

        echo "Pre-filter alignment QC complete" > {log.out:q}
        """
