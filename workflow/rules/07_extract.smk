# 07. Feature extraction

rule extract_features:
    input:
        bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam"
    output:
        m6a=f"{RESULTS}/07_extract/{{sample}}.m6a.bed.gz",
        cpg=f"{RESULTS}/07_extract/{{sample}}.cpg.bed.gz",
        nuc=f"{RESULTS}/07_extract/{{sample}}.nuc.bed.gz",
        msp=f"{RESULTS}/07_extract/{{sample}}.msp.bed.gz"
    threads:
        config["threads"]["extract"]
    params:
        coords=(
            "--reference"
            if COORD == "reference"
            else "--molecular"
        )
    log:
        out=f"{LOGS}/07_extract/{{sample}}.out",
        err=f"{LOGS}/07_extract/{{sample}}.err"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.m6a:q})"
        mkdir -p "$(dirname {log.out:q})"

        ft -t {threads} extract \
            {params.coords} \
            --m6a {output.m6a:q} \
            --cpg {output.cpg:q} \
            --nuc {output.nuc:q} \
            --msp {output.msp:q} \
            {input.bam:q} \
            > {log.out:q} \
            2> {log.err:q}
        """


rule extract_fire:
    input:
        bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam"
    output:
        bed=f"{RESULTS}/07_extract/{{sample}}.fire.bed.gz"
    threads:
        config["threads"]["extract"]
    log:
        out=f"{LOGS}/07_extract/{{sample}}.fire.out",
        err=f"{LOGS}/07_extract/{{sample}}.fire.err"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.bed:q})"
        mkdir -p "$(dirname {log.out:q})"

        ft -t {threads} fire \
            --extract \
            {input.bam:q} \
            {output.bed:q} \
            > {log.out:q} \
            2> {log.err:q}
        """
