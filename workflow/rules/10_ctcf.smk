# Direct official FiberFold commands; conversions have explicit rule outputs.
if CTCF:
    CTCF_ROOT = f"{RESULTS}/10_ctcf"
    CTCF_REFERENCE = config["ctcf"]["reference_fasta"]
    CTCF_CHROMS = config["ctcf"]["chromosomes"]
    CTCF_MISC = Path(config["fiberfold"]["script"]).resolve().parent / "misc"
    CTCF_MODULES = str(Path(workflow.basedir).parent / "config/ctcf.modules.yaml")
    CTCF_SOURCE = config["ctcf"].get("motifs") or f"{CTCF_ROOT}/reference/MA0139.1.source.tsv.gz"

    if not config["ctcf"].get("motifs"):
        rule download_ctcf_motifs:
            output: source=CTCF_SOURCE
            params: url=config["ctcf"]["motif_url"]
            log:
                out=f"{LOGS}/10_ctcf/download.out",
                err=f"{LOGS}/10_ctcf/download.err"
            shell:
                r"""
                mkdir -p "$(dirname {output.source:q})" "$(dirname {log.out:q})"
                curl -fL --retry 3 --output {output.source:q} {params.url:q} \
                    > {log.out:q} 2> {log.err:q}
                """

    rule index_ctcf_reference:
        input: fasta=CTCF_REFERENCE
        output:
            fasta=f"{CTCF_ROOT}/reference/genome.fa",
            fai=f"{CTCF_ROOT}/reference/genome.fa.fai"
        params: reference=str(Path(CTCF_REFERENCE).resolve())
        log:
            out=f"{LOGS}/10_ctcf/faidx.out",
            err=f"{LOGS}/10_ctcf/faidx.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.fai:q})" "$(dirname {log.out:q})"
            (
                ln -sf {params.reference:q} {output.fasta:q}
                samtools faidx {output.fasta:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule filter_ctcf_motifs:
        input:
            source=CTCF_SOURCE,
            fai=rules.index_ctcf_reference.output.fai,
            script=str(CTCF_MISC / "filterOverlapMotifs.py")
        output:
            converted=f"{CTCF_ROOT}/reference/motifs/{{chrom}}.input.bed",
            bed=f"{CTCF_ROOT}/reference/motifs/{{chrom}}.bed",
            chromsizes=f"{CTCF_ROOT}/reference/chromsizes/{{chrom}}.sizes"
        params: min_score=config["ctcf"].get("min_score", 0)
        wildcard_constraints: chrom="|".join(CTCF_CHROMS)
        log:
            out=f"{LOGS}/10_ctcf/motifs.{{chrom}}.out",
            err=f"{LOGS}/10_ctcf/motifs.{{chrom}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bed:q})" "$(dirname {output.chromsizes:q})" "$(dirname {log.out:q})"
            (
                awk -v chrom={wildcards.chrom:q} 'BEGIN {{OFS="\t"}} $1==chrom {{print $1,$2}}' \
                    {input.fai:q} > {output.chromsizes:q}
                test -s {output.chromsizes:q}
                chrom_length=$(cut -f2 {output.chromsizes:q})
                if [[ {input.source:q} == *.gz ]]; then gzip -cd {input.source:q}; else cat {input.source:q}; fi \
                | awk -v chrom={wildcards.chrom:q} -v limit="$chrom_length" -v threshold={params.min_score} \
                    'BEGIN {{OFS="\t"}} /^(#|track|browser)/ {{next}} $1==chrom {{
                        if ((NF!=6 && NF!=7) || $3-$2!=19 || $2<0 || $3>limit || ($NF!="+" && $NF!="-")) {{
                            print "Invalid MA0139.1 BED row: " $0 > "/dev/stderr"; exit 1
                        }}
                        if ($5>=threshold) print $1,$2,$3,$4,$5,$NF
                    }}' \
                | LC_ALL=C sort -s -k2,2n -k3,3n -k5,5nr -k6,6 > {output.converted:q}
                {PYTHON:q} {input.script:q} {output.converted:q} > {output.bed:q}
                wc -l {output.bed:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule prepare_ctcf_motifs:
        input: beds=expand(f"{CTCF_ROOT}/reference/motifs/{{chrom}}.bed", chrom=CTCF_CHROMS)
        output: bed=f"{CTCF_ROOT}/reference/CTCF.MA0139.1.bed"
        log:
            out=f"{LOGS}/10_ctcf/motifs.out",
            err=f"{LOGS}/10_ctcf/motifs.err"
        shell:
            r"""
            mkdir -p "$(dirname {log.out:q})"
            (cat {input.beds:q} > {output.bed:q}; wc -l {output.bed:q}) > {log.out:q} 2> {log.err:q}
            """

    rule ctcf_direction:
        input:
            bed=f"{CTCF_ROOT}/reference/motifs/{{chrom}}.bed",
            chromsizes=f"{CTCF_ROOT}/reference/chromsizes/{{chrom}}.sizes",
            script=str(CTCF_MISC / "encodeMotifDirection.py")
        output: tracks=f"{CTCF_ROOT}/direction/{{chrom}}.bw"
        wildcard_constraints: chrom="|".join(CTCF_CHROMS)
        log:
            out=f"{LOGS}/10_ctcf/direction.{{chrom}}.out",
            err=f"{LOGS}/10_ctcf/direction.{{chrom}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.tracks:q})" "$(dirname {log.out:q})"
            {PYTHON:q} {input.script:q} {input.bed:q} {input.chromsizes:q} {output.tracks:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule ctcf_footprint:
        input:
            bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/05_fire/{{sample}}.fire.bam.bai",
            bed=rules.prepare_ctcf_motifs.output.bed,
            modules=CTCF_MODULES
        output: table=f"{CTCF_ROOT}/{{sample}}.footprints.tsv.gz"
        threads: config["threads"].get("ctcf", 4)
        log:
            out=f"{LOGS}/10_ctcf/{{sample}}.footprint.out",
            err=f"{LOGS}/10_ctcf/{{sample}}.footprint.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.table:q})" "$(dirname {log.out:q})"
            ft -t {threads} footprint --bed {input.bed:q} --yaml {input.modules:q} {input.bam:q} \
                2> {log.err:q} | gzip -n -c > {output.table:q} 2>> {log.err:q}
            echo "CTCF footprint table complete" > {log.out:q}
            """

    rule ctcf_legacy_table:
        input: table=rules.ctcf_footprint.output.table
        output: table=f"{CTCF_ROOT}/{{sample}}.footprints.legacy.tsv.gz"
        log:
            out=f"{LOGS}/10_ctcf/{{sample}}.legacy.out",
            err=f"{LOGS}/10_ctcf/{{sample}}.legacy.err"
        shell:
            r"""
            mkdir -p "$(dirname {log.out:q})"
            (
                gzip -cd {input.table:q} | awk -F '\t' 'BEGIN {{OFS="\t"}}
                    NR==1 {{
                        sub(/^#/, "", $1)
                        split("chrom start end strand n_spanning_fibers n_spanning_msps n_overlapping_nucs module:0-6 module:6-13 module:13-19 footprint_codes fire_quals fiber_names", columns, " ")
                        for (i=1;i<=NF;i++) field[$i]=i
                        for (i=1;i<=13;i++) if (!(columns[i] in field)) {{
                            print "Missing footprint column: " columns[i] > "/dev/stderr"; exit 1
                        }}
                        printf "#"
                    }}
                    NR>1 {{
                        c=$(field["chrom"]); s=$(field["start"]); e=$(field["end"])
                        if (s<last_end[c]) {{print "Unsorted, duplicate or overlapping footprint rows: " c ":" s > "/dev/stderr"; exit 1}}
                        last_end[c]=e
                    }}
                    {{for (i=1;i<=13;i++) printf "%s%s", $(field[columns[i]]), (i==13 ? ORS : OFS)}}' \
                | gzip -n -c > {output.table:q}
                echo "Converted footprint table to official 13-column format"
            ) > {log.out:q} 2> {log.err:q}
            """

    rule ctcf_footprint_tracks:
        input:
            table=rules.ctcf_legacy_table.output.table,
            chromsizes=f"{CTCF_ROOT}/reference/chromsizes/{{chrom}}.sizes",
            reference=rules.index_ctcf_reference.output.fasta,
            fai=rules.index_ctcf_reference.output.fai,
            script=str(CTCF_MISC / "encodeFP.py")
        output: tracks=f"{CTCF_ROOT}/footprint/{{sample}}/{{chrom}}.bw"
        wildcard_constraints: chrom="|".join(CTCF_CHROMS)
        log:
            out=f"{LOGS}/10_ctcf/{{sample}}.tracks.{{chrom}}.out",
            err=f"{LOGS}/10_ctcf/{{sample}}.tracks.{{chrom}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.tracks:q})" "$(dirname {log.out:q})"
            {PYTHON:q} {input.script:q} {input.table:q} {input.chromsizes:q} {input.reference:q} {output.tracks:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule prepare_ctcf:
        input:
            expand(f"{CTCF_ROOT}/direction/{{chrom}}.bw", chrom=CTCF_CHROMS),
            expand(f"{CTCF_ROOT}/footprint/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=CTCF_CHROMS)
