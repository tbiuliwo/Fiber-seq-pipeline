# 08. Site-level 5mC pileup and plotting tracks from indexed FIRE BAMs.
if METHYLATION:
    METHYLATION_CONFIG = config.get("methylation", {})
    METHYLATION_ROOT = f"{RESULTS}/08_methylation"
    METHYLATION_TOOLS = {
        name: METHYLATION_CONFIG.get(name) or str(
            Path(config.get("tools", {}).get("bin_dir") or "") / name
        )
        for name in ("methbat", "bgzip", "tabix")
    }

    rule methylation_pileup:
        input:
            bam=f"{RESULTS}/05_fire/{{sample}}.fire.bam",
            bai=f"{RESULTS}/05_fire/{{sample}}.fire.bam.bai"
        output:
            mc=f"{METHYLATION_ROOT}/{{sample}}.5mC.bed.gz",
            mc_idx=f"{METHYLATION_ROOT}/{{sample}}.5mC.bed.gz.tbi",
            summary=f"{METHYLATION_ROOT}/{{sample}}.summary.json"
        threads: config["threads"].get("methylation", 4)
        params:
            methbat=METHYLATION_TOOLS["methbat"],
            prefix=lambda wc: f"{METHYLATION_ROOT}/{wc.sample}",
            min_mapq=METHYLATION_CONFIG.get("min_mapq", 1),
            min_coverage=METHYLATION_CONFIG.get("min_coverage", 4),
            edge_trimming_size=METHYLATION_CONFIG.get("edge_trimming_size", 20)
        log:
            out=f"{LOGS}/08_methylation/{{sample}}.out",
            err=f"{LOGS}/08_methylation/{{sample}}.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.mc:q})" "$(dirname {log.out:q})"
            {params.methbat:q} pileup \
                --threads {threads} \
                --input-bam {input.bam:q} \
                --output-prefix {params.prefix:q} \
                --min-mapq {params.min_mapq} \
                --min-coverage {params.min_coverage} \
                --edge-trimming-size {params.edge_trimming_size} \
                --skip-6mA --skip-5hmC \
                > {log.out:q} 2> {log.err:q}
            """

    rule methylation_bedgraph:
        input:
            bed=rules.methylation_pileup.output.mc
        output:
            bedgraph=f"{METHYLATION_ROOT}/{{sample}}.5mC.bedgraph.gz"
        threads: config["threads"].get("methylation", 4)
        params:
            bgzip=METHYLATION_TOOLS["bgzip"]
        log:
            out=f"{LOGS}/08_methylation/{{sample}}.bedgraph.out",
            err=f"{LOGS}/08_methylation/{{sample}}.bedgraph.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bedgraph:q})" "$(dirname {log.out:q})"
            (
                gzip -cd {input.bed:q} \
                | awk -F '\t' 'BEGIN {{OFS="\t"}}
                    /^##/ {{next}}
                    /^#chrom\t/ {{
                        sub(/^#/, "", $1)
                        for (i=1; i<=NF; i++) field[$i]=i
                        split("chrom start end mod_score type", required, " ")
                        for (i=1; i<=5; i++) if (!(required[i] in field)) {{
                            print "Missing MethBat column: " required[i] > "/dev/stderr"; exit 1
                        }}
                        header=1; next
                    }}
                    /^#/ {{next}}
                    {{
                        if (!header) {{print "Missing MethBat header" > "/dev/stderr"; exit 1}}
                        if ($(field["type"])!="Total") next
                        value=$(field["mod_score"])
                        if (value !~ /^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/ || value<0 || value>100) {{
                            print "Invalid MethBat 5mC percentage: " value > "/dev/stderr"; exit 1
                        }}
                        print $(field["chrom"]),$(field["start"]),$(field["end"]),value
                    }}
                    END {{if (!header) exit 1}}' \
                | {params.bgzip:q} -@ {threads} -c > {output.bedgraph:q}
            ) > {log.out:q} 2> {log.err:q}
            """

    rule index_methylation_bedgraph:
        input:
            bedgraph=rules.methylation_bedgraph.output.bedgraph
        output:
            tbi=f"{METHYLATION_ROOT}/{{sample}}.5mC.bedgraph.gz.tbi"
        params:
            tabix=METHYLATION_TOOLS["tabix"]
        log:
            out=f"{LOGS}/08_methylation/{{sample}}.bedgraph_index.out",
            err=f"{LOGS}/08_methylation/{{sample}}.bedgraph_index.err"
        shell:
            r"""
            mkdir -p "$(dirname {log.out:q})"
            {params.tabix:q} -f -p bed {input.bedgraph:q} > {log.out:q} 2> {log.err:q}
            """

    rule methylation:
        input:
            expand(
                f"{METHYLATION_ROOT}/{{sample}}.{{suffix}}", sample=SAMPLES,
                suffix=["5mC.bed.gz", "5mC.bed.gz.tbi", "summary.json",
                        "5mC.bedgraph.gz", "5mC.bedgraph.gz.tbi"],
            )
