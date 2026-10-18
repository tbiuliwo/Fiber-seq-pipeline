if FOLD_INPUTS:
    FOLD_ROOT = f"{RESULTS}/09_fiberfold_inputs"
    FOLD_CHROMS = config["ctcf"]["chromosomes"]
    FOLD_SCRIPT = str(Path(workflow.basedir).parent / "scripts/fiberfold_inputs.py")

    rule stage_fold_reference:
        input:
            fasta=config["ctcf"]["reference_fasta"],
            fai=rules.index_ctcf_reference.output.fai
        output:
            fasta=f"{FOLD_ROOT}/reference/genome.fa",
            fai=f"{FOLD_ROOT}/reference/genome.fa.fai"
        params:
            reference=str(Path(config["ctcf"]["reference_fasta"]).resolve())
        log:
            out=f"{LOGS}/09_fiberfold_inputs/reference.out",
            err=f"{LOGS}/09_fiberfold_inputs/reference.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.fasta:q})" "$(dirname {log.out:q})"
            (
                ln -sf {params.reference:q} {output.fasta:q}
                cp {input.fai:q} {output.fai:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_windows:
        input:
            fai=rules.index_ctcf_reference.output.fai,
            script=FOLD_SCRIPT
        output:
            windows=f"{FOLD_ROOT}/windows.bed.gz",
            regions=f"{FOLD_ROOT}/chromosomes.bed"
        params:
            chromosomes=FOLD_CHROMS,
            step=config["fiberfold"].get("window_step", 2097152)
        log:
            out=f"{LOGS}/09_fiberfold_inputs/windows.out",
            err=f"{LOGS}/09_fiberfold_inputs/windows.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.windows:q})" "$(dirname {log.out:q})"
            {PYTHON:q} {input.script:q} windows --fai {input.fai:q} --chromosomes {params.chromosomes:q} \
                --step {params.step} --out {output.windows:q} --regions {output.regions:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_fire_tracks:
        input:
            bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai",
            fai=rules.index_ctcf_reference.output.fai,
            regions=rules.fiberfold_windows.output.regions,
            script=FOLD_SCRIPT
        output:
            tracks=expand(f"{FOLD_ROOT}/fire/{{sample}}/{{chrom}}.bw", chrom=FOLD_CHROMS, allow_missing=True)
        threads: config["threads"].get("fiberfold_inputs", 4)
        params:
            chromosomes=FOLD_CHROMS,
            directory=lambda w: f"{FOLD_ROOT}/fire/{w.sample}"
        log:
            out=f"{LOGS}/09_fiberfold_inputs/{{sample}}.fire.out",
            err=f"{LOGS}/09_fiberfold_inputs/{{sample}}.fire.err"
        shell:
            r"""
            mkdir -p "$(dirname {log.out:q})"
            (
                ft -t {threads} pileup --bed {input.regions:q} --keep-zeros --no-msp --no-nuc {input.bam:q} \
                | {PYTHON:q} {input.script:q} fire --fai {input.fai:q} --chromosomes {params.chromosomes:q} \
                    --out-dir {params.directory:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_cpg_scores:
        input:
            bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai",
            reference=rules.stage_fold_reference.output.fasta,
            fai=rules.stage_fold_reference.output.fai
        output:
            bed=f"{FOLD_ROOT}/cpg_scores/{{sample}}.combined.bed.gz",
            tbi=f"{FOLD_ROOT}/cpg_scores/{{sample}}.combined.bed.gz.tbi",
            bw=f"{FOLD_ROOT}/cpg_scores/{{sample}}.combined.bw"
        threads: config["threads"].get("fiberfold_inputs", 4)
        params:
            prefix=lambda w: f"{FOLD_ROOT}/cpg_scores/{w.sample}"
        log:
            out=f"{LOGS}/09_fiberfold_inputs/{{sample}}.cpg.out",
            err=f"{LOGS}/09_fiberfold_inputs/{{sample}}.cpg.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bw:q})" "$(dirname {log.out:q})"
            aligned_bam_to_cpg_scores --bam {input.bam:q} --ref {input.reference:q} \
                --modsites-mode reference --pileup-mode model --min-coverage 4 \
                --threads {threads} --output-prefix {params.prefix:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_cpg_tracks:
        input:
            bw=rules.fiberfold_cpg_scores.output.bw,
            fai=rules.index_ctcf_reference.output.fai,
            script=FOLD_SCRIPT
        output:
            tracks=expand(f"{FOLD_ROOT}/cpg/{{sample}}/{{chrom}}.bw", chrom=FOLD_CHROMS, allow_missing=True)
        params:
            chromosomes=FOLD_CHROMS,
            directory=lambda w: f"{FOLD_ROOT}/cpg/{w.sample}"
        log:
            out=f"{LOGS}/09_fiberfold_inputs/{{sample}}.cpg_tracks.out",
            err=f"{LOGS}/09_fiberfold_inputs/{{sample}}.cpg_tracks.err"
        shell:
            r"""
            {PYTHON:q} {input.script:q} cpg --bigwig {input.bw:q} --fai {input.fai:q} \
                --chromosomes {params.chromosomes:q} --out-dir {params.directory:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule prepare_fiberfold_inputs:
        input:
            rules.fiberfold_windows.output.windows,
            expand(f"{FOLD_ROOT}/fire/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=FOLD_CHROMS),
            expand(f"{FOLD_ROOT}/cpg/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=FOLD_CHROMS),
            rules.ctcf_direction.output.tracks,
            expand(f"{RESULTS}/08_ctcf/footprint/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=FOLD_CHROMS)
