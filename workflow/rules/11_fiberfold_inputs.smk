# Standard tools prepare FiberFold inputs; 
if FOLD_INPUTS:
    FOLD_ROOT = f"{RESULTS}/11_fiberfold_inputs"
    FOLD_CHROMS = config["ctcf"]["chromosomes"]
    FOLD_BIGTOOLS = config.get("tools", {}).get("bigtools", "bigtools")

    rule stage_fold_reference:
        input:
            fasta=config["ctcf"]["reference_fasta"],
            fai=rules.index_ctcf_reference.output.fai
        output:
            fasta=f"{FOLD_ROOT}/reference/genome.fa",
            fai=f"{FOLD_ROOT}/reference/genome.fa.fai"
        params: reference=str(Path(config["ctcf"]["reference_fasta"]).resolve())
        log:
            out=f"{LOGS}/11_fiberfold_inputs/reference.out",
            err=f"{LOGS}/11_fiberfold_inputs/reference.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.fasta:q})" "$(dirname {log.out:q})"
            (ln -sf {params.reference:q} {output.fasta:q}; cp {input.fai:q} {output.fai:q}) \
                > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_chromosome:
        input: fai=rules.index_ctcf_reference.output.fai
        output:
            sizes=f"{FOLD_ROOT}/reference/chromsizes/{{chrom}}.sizes",
            region=f"{FOLD_ROOT}/reference/regions/{{chrom}}.bed"
        wildcard_constraints: chrom="|".join(FOLD_CHROMS)
        log:
            out=f"{LOGS}/11_fiberfold_inputs/reference.{{chrom}}.out",
            err=f"{LOGS}/11_fiberfold_inputs/reference.{{chrom}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.sizes:q})" "$(dirname {output.region:q})" "$(dirname {log.out:q})"
            (
                awk -v chrom={wildcards.chrom:q} 'BEGIN {{OFS="\t"}} $1==chrom {{print $1,$2}}' {input.fai:q} > {output.sizes:q}
                test -s {output.sizes:q}
                awk 'BEGIN {{OFS="\t"}} {{print $1,0,$2}}' {output.sizes:q} > {output.region:q}
                cat {output.sizes:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_windows:
        input: regions=expand(f"{FOLD_ROOT}/reference/regions/{{chrom}}.bed", chrom=FOLD_CHROMS)
        output:
            windows=f"{FOLD_ROOT}/windows.bed.gz",
            regions=f"{FOLD_ROOT}/chromosomes.bed"
        params: step=config["fiberfold"].get("window_step", 2097152)
        log:
            out=f"{LOGS}/11_fiberfold_inputs/windows.out",
            err=f"{LOGS}/11_fiberfold_inputs/windows.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.windows:q})" "$(dirname {log.out:q})"
            (
                cat {input.regions:q} > {output.regions:q}
                bedtools makewindows -b {output.regions:q} -w 2097152 -s {params.step} \
                    | awk '$3-$2==2097152 {{print; n++}} END {{if (!n) {{print "No complete 2097152-bp windows" > "/dev/stderr"; exit 1}}}}' \
                    | gzip -n -c > {output.windows:q}
                echo "Generated complete 2097152-bp prediction windows"
            ) > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_fire_tracks:
        input:
            bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/05_fire/{{sample}}.fire.bam.bai",
            sizes=f"{FOLD_ROOT}/reference/chromsizes/{{chrom}}.sizes",
            region=f"{FOLD_ROOT}/reference/regions/{{chrom}}.bed"
        output:
            bedgraph=temp(f"{FOLD_ROOT}/fire/{{sample}}/{{chrom}}.bedGraph"),
            tracks=f"{FOLD_ROOT}/fire/{{sample}}/{{chrom}}.bw"
        threads: config["threads"].get("fiberfold_inputs", 4)
        params: bigtools=FOLD_BIGTOOLS
        wildcard_constraints: chrom="|".join(FOLD_CHROMS)
        log:
            out=f"{LOGS}/11_fiberfold_inputs/{{sample}}.fire.{{chrom}}.out",
            err=f"{LOGS}/11_fiberfold_inputs/{{sample}}.fire.{{chrom}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.tracks:q})" "$(dirname {log.out:q})"
            (
                {params.bigtools:q} --version
                ft -t {threads} pileup --bed {input.region:q} --keep-zeros --no-msp --no-nuc {input.bam:q} \
                    | awk -F '\t' 'BEGIN {{OFS="\t"}}
                        NR==1 {{
                            sub(/^#/, "", $1); for(i=1;i<=NF;i++) field[$i]=i
                            if (!("coverage" in field) || !("fire_coverage" in field)) {{print "Missing pileup coverage columns" > "/dev/stderr"; exit 1}}
                            next
                        }}
                        {{d=$(field["coverage"]); n=$(field["fire_coverage"])
                          if (d<0 || n<0 || n>d) {{print "Invalid FIRE coverage counts" > "/dev/stderr"; exit 1}}
                          printf "%s\t%d\t%d\t%.17g\n", $1,$2,$3,(d ? n/d : 0)}}' > {output.bedgraph:q}
                if [[ ! -s {output.bedgraph:q} ]]; then
                    awk 'BEGIN {{OFS="\t"}} {{print $1,0,$2,0}}' {input.sizes:q} > {output.bedgraph:q}
                fi
                {params.bigtools:q} bedgraphtobigwig --nthreads {threads} {output.bedgraph:q} {input.sizes:q} {output.tracks:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule fiberfold_cpg_scores:
        input:
            bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/05_fire/{{sample}}.fire.bam.bai",
            reference=rules.stage_fold_reference.output.fasta,
            fai=rules.stage_fold_reference.output.fai
        output:
            bed=f"{FOLD_ROOT}/cpg_scores/{{sample}}.combined.bed.gz",
            tbi=f"{FOLD_ROOT}/cpg_scores/{{sample}}.combined.bed.gz.tbi",
            bw=f"{FOLD_ROOT}/cpg_scores/{{sample}}.combined.bw"
        threads: config["threads"].get("fiberfold_inputs", 4)
        params: prefix=lambda w: f"{FOLD_ROOT}/cpg_scores/{w.sample}"
        log:
            out=f"{LOGS}/11_fiberfold_inputs/{{sample}}.cpg.out",
            err=f"{LOGS}/11_fiberfold_inputs/{{sample}}.cpg.err"
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
            bed=rules.fiberfold_cpg_scores.output.bed,
            tbi=rules.fiberfold_cpg_scores.output.tbi,
            sizes=f"{FOLD_ROOT}/reference/chromsizes/{{chrom}}.sizes"
        output:
            bedgraph=temp(f"{FOLD_ROOT}/cpg/{{sample}}/{{chrom}}.bedGraph"),
            tracks=f"{FOLD_ROOT}/cpg/{{sample}}/{{chrom}}.bw"
        threads: config["threads"].get("fiberfold_inputs", 4)
        params: bigtools=FOLD_BIGTOOLS
        wildcard_constraints: chrom="|".join(FOLD_CHROMS)
        log:
            out=f"{LOGS}/11_fiberfold_inputs/{{sample}}.cpg.{{chrom}}.out",
            err=f"{LOGS}/11_fiberfold_inputs/{{sample}}.cpg.{{chrom}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.tracks:q})" "$(dirname {log.out:q})"
            (
                {params.bigtools:q} --version
                tabix {input.bed:q} {wildcards.chrom:q} | awk -v chrom={wildcards.chrom:q} '$1==chrom {{
                    if ($4<0 || $4>100) {{print "CpG percentage outside 0-100" > "/dev/stderr"; exit 1}}
                    printf "%s\t%d\t%d\t%.17g\n",$1,$2,$3,$4/100
                }}' > {output.bedgraph:q}
                if [[ ! -s {output.bedgraph:q} ]]; then
                    awk 'BEGIN {{OFS="\t"}} {{print $1,0,$2,0}}' {input.sizes:q} > {output.bedgraph:q}
                fi
                {params.bigtools:q} bedgraphtobigwig --nthreads {threads} {output.bedgraph:q} {input.sizes:q} {output.tracks:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule prepare_fiberfold_inputs:
        input:
            rules.fiberfold_windows.output.windows,
            expand(f"{FOLD_ROOT}/fire/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=FOLD_CHROMS),
            expand(f"{FOLD_ROOT}/cpg/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=FOLD_CHROMS),
            expand(f"{RESULTS}/10_ctcf/direction/{{chrom}}.bw", chrom=FOLD_CHROMS),
            expand(f"{RESULTS}/10_ctcf/footprint/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=FOLD_CHROMS)
