# MA0139.1 motif preparation and sample-specific CTCF footprints for FiberFold.
if CTCF:
    CTCF_ROOT = f"{RESULTS}/08_ctcf"
    CTCF_REFERENCE = config["ctcf"]["reference_fasta"]
    CTCF_CHROMS = config["ctcf"]["chromosomes"]
    CTCF_SCRIPT = str(Path(workflow.basedir).parent / "scripts/ctcf_tracks.py")
    CTCF_MODULES = str(Path(workflow.basedir).parent / "config/ctcf.modules.yaml")
    CTCF_SOURCE = config["ctcf"].get("motifs") or f"{CTCF_ROOT}/reference/MA0139.1.source.tsv.gz"

    if not config["ctcf"].get("motifs"):
        rule download_ctcf_motifs:
            output:
                source=CTCF_SOURCE
            params:
                url=config["ctcf"]["motif_url"]
            log:
                out=f"{LOGS}/08_ctcf/download.out",
                err=f"{LOGS}/08_ctcf/download.err"
            shell:
                r"""
                mkdir -p "$(dirname {output.source:q})" "$(dirname {log.out:q})"
                curl -fL --retry 3 --output {output.source:q} {params.url:q} \
                    > {log.out:q} 2> {log.err:q}
                """

    rule index_ctcf_reference:
        input:
            fasta=CTCF_REFERENCE
        output:
            fai=f"{CTCF_ROOT}/reference/genome.fa.fai"
        log:
            out=f"{LOGS}/08_ctcf/faidx.out",
            err=f"{LOGS}/08_ctcf/faidx.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.fai:q})" "$(dirname {log.out:q})"
            samtools faidx --fai-idx {output.fai:q} {input.fasta:q} > {log.out:q} 2> {log.err:q}
            """

    rule prepare_ctcf_motifs:
        input:
            source=CTCF_SOURCE,
            fai=rules.index_ctcf_reference.output.fai,
            script=CTCF_SCRIPT
        output:
            bed=f"{CTCF_ROOT}/reference/CTCF.MA0139.1.bed"
        params:
            chromosomes=CTCF_CHROMS,
            min_score=config["ctcf"].get("min_score", 0)
        log:
            out=f"{LOGS}/08_ctcf/motifs.out",
            err=f"{LOGS}/08_ctcf/motifs.err"
        shell:
            r"""
            {PYTHON:q} {input.script:q} prepare --source {input.source:q} --fai {input.fai:q} \
                --chromosomes {params.chromosomes:q} --min-score {params.min_score} \
                --out {output.bed:q} > {log.out:q} 2> {log.err:q}
            """

    rule ctcf_direction:
        input:
            bed=rules.prepare_ctcf_motifs.output.bed,
            fai=rules.index_ctcf_reference.output.fai,
            script=CTCF_SCRIPT
        output:
            tracks=expand(f"{CTCF_ROOT}/direction/{{chrom}}.bw", chrom=CTCF_CHROMS)
        params:
            directory=f"{CTCF_ROOT}/direction",
            chromosomes=CTCF_CHROMS
        log:
            out=f"{LOGS}/08_ctcf/direction.out",
            err=f"{LOGS}/08_ctcf/direction.err"
        shell:
            r"""
            {PYTHON:q} {input.script:q} direction --bed {input.bed:q} --fai {input.fai:q} \
                --chromosomes {params.chromosomes:q} --out-dir {params.directory:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule ctcf_footprint:
        input:
            bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai",
            bed=rules.prepare_ctcf_motifs.output.bed,
            modules=CTCF_MODULES
        output:
            table=f"{CTCF_ROOT}/{{sample}}.footprints.tsv.gz"
        threads: config["threads"].get("ctcf", 4)
        log:
            out=f"{LOGS}/08_ctcf/{{sample}}.footprint.out",
            err=f"{LOGS}/08_ctcf/{{sample}}.footprint.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.table:q})" "$(dirname {log.out:q})"
            ft -t {threads} footprint --bed {input.bed:q} --yaml {input.modules:q} {input.bam:q} \
                2> {log.err:q} | gzip -c > {output.table:q} 2>> {log.err:q}
            echo "CTCF footprint table complete" > {log.out:q}
            """

    rule ctcf_footprint_tracks:
        input:
            table=rules.ctcf_footprint.output.table,
            bed=rules.prepare_ctcf_motifs.output.bed,
            reference=CTCF_REFERENCE,
            fai=rules.index_ctcf_reference.output.fai,
            script=CTCF_SCRIPT
        output:
            tracks=expand(f"{CTCF_ROOT}/footprint/{{sample}}/{{chrom}}.bw", chrom=CTCF_CHROMS, allow_missing=True),
            scores=f"{CTCF_ROOT}/{{sample}}.scores.tsv.gz"
        params:
            directory=lambda w: f"{CTCF_ROOT}/footprint/{w.sample}",
            chromosomes=CTCF_CHROMS
        log:
            out=f"{LOGS}/08_ctcf/{{sample}}.tracks.out",
            err=f"{LOGS}/08_ctcf/{{sample}}.tracks.err"
        shell:
            r"""
            {PYTHON:q} {input.script:q} footprints --table {input.table:q} --bed {input.bed:q} \
                --reference {input.reference:q} --fai {input.fai:q} --chromosomes {params.chromosomes:q} \
                --out-dir {params.directory:q} --scores {output.scores:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule prepare_ctcf:
        input:
            rules.ctcf_direction.output.tracks,
            expand(f"{CTCF_ROOT}/footprint/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=CTCF_CHROMS),
            expand(f"{CTCF_ROOT}/{{sample}}.scores.tsv.gz", sample=SAMPLES)
