if FIBERFOLD:
    rule predict_contacts:
        input:
            script=config["fiberfold"]["script"],
            model=config["fiberfold"]["model"],
            windows=f"{RESULTS}/11_fiberfold_inputs/windows.bed.gz",
            fire=expand(f"{RESULTS}/11_fiberfold_inputs/fire/{{sample}}/{{chrom}}.bw", chrom=config["ctcf"]["chromosomes"], allow_missing=True),
            cpg=expand(f"{RESULTS}/11_fiberfold_inputs/cpg/{{sample}}/{{chrom}}.bw", chrom=config["ctcf"]["chromosomes"], allow_missing=True),
            footprint=expand(f"{RESULTS}/10_ctcf/footprint/{{sample}}/{{chrom}}.bw", chrom=config["ctcf"]["chromosomes"], allow_missing=True),
            direction=expand(f"{RESULTS}/10_ctcf/direction/{{chrom}}.bw", chrom=config["ctcf"]["chromosomes"])
        output:
            predictions=f"{RESULTS}/12_fiberfold/{{sample}}.predictions.pkl"
        threads: config["threads"]["fiberfold"]
        params:
            fire=lambda w: f"{RESULTS}/11_fiberfold_inputs/fire/{w.sample}",
            cpg=lambda w: f"{RESULTS}/11_fiberfold_inputs/cpg/{w.sample}",
            footprint=lambda w: f"{RESULTS}/10_ctcf/footprint/{w.sample}",
            direction=f"{RESULTS}/10_ctcf/direction"
        log:
            out=f"{LOGS}/12_fiberfold/{{sample}}.out",
            err=f"{LOGS}/12_fiberfold/{{sample}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.predictions:q})" "$(dirname {log.out:q})"
            OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
            {PYTHON:q} {input.script:q} --model {input.model:q} --ds {input.windows:q} \
                --mA {params.fire:q} --mC {params.cpg:q} --ctcf_fp {params.footprint:q} \
                --ctcf_dir {params.direction:q} --out {output.predictions:q} \
                > {log.out:q} 2> {log.err:q}
            """
