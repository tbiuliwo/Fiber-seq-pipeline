# 02. Optional phasing and haplotagging of the step-01 aligned BAM, before annotation.
if PHASING:
    PHASING_ROOT = f"{RESULTS}/02_phasing"
    PHASING_REFERENCE = PHASING_CONFIG.get("reference_fasta") or REFERENCE
    if str(PHASING_REFERENCE).endswith(".mmi"):
        raise ValueError("phasing.reference_fasta must be a FASTA, not a pbmm2 index")
    PHASING_WHATSHAP = PHASING_CONFIG.get("whatshap", "whatshap")
    PHASING_SAMTOOLS = PHASING_CONFIG.get("samtools") or str(
        Path(config.get("tools", {}).get("bin_dir") or "") / "samtools"
    )
    PHASING_TABIX = PHASING_CONFIG.get("tabix") or str(
        Path(config.get("tools", {}).get("bin_dir") or "") / "tabix"
    )

    rule deepvariant:
        input:
            bam=lambda wc: f"{RESULTS}/01_align/{PHASING_SAMPLE_IDS[wc.phasing_sample]}.chromosomes.bam",
            bai=lambda wc: f"{RESULTS}/01_align/{PHASING_SAMPLE_IDS[wc.phasing_sample]}.chromosomes.bam.bai",
            reference=PHASING_REFERENCE,
            fai=f"{PHASING_REFERENCE}.fai",
            container=PHASING_CONFIG["container"]
        output:
            vcf=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.deepvariant.vcf.gz",
            tbi=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.deepvariant.vcf.gz.tbi"
        wildcard_constraints:
            phasing_sample="|".join(re.escape(alias) for alias in PHASING_SAMPLE_IDS)
        params:
            runner=PHASING_CONFIG.get("singularity", "singularity"),
            binds=" ".join("--bind " + shlex.quote(path) for path in PHASING_CONFIG.get("bind_paths", [])),
            chromosomes=" ".join(ANALYSIS_CHROMOSOMES),
            bam=lambda wc, input: str(Path(input.bam).resolve()),
            reference=str(Path(PHASING_REFERENCE).resolve()),
            intermediate=lambda wc: str(Path(f"{PHASING_ROOT}/{wc.phasing_sample}/intermediate").resolve()),
            internal_logs=lambda wc: str(Path(f"{PHASING_ROOT}/{wc.phasing_sample}/logs").resolve()),
            vcf=lambda wc, output: str(Path(output.vcf).resolve())
        threads: config["threads"].get("deepvariant", 16)
        resources:
            mem_mb=PHASING_CONFIG.get("mem_mb", 64000)
        log:
            out=f"{LOGS}/02_phasing/{{phasing_sample}}.deepvariant.out",
            err=f"{LOGS}/02_phasing/{{phasing_sample}}.deepvariant.err"
        shell:
            r"""
            mkdir -p {params.intermediate:q} {params.internal_logs:q} "$(dirname {log.out:q})"
            run_id="$(date '+%Y%m%d_%H%M%S')_$$"
            for logfile in {log.out:q} {log.err:q}; do
                if [[ -f "$logfile" ]]; then cp -p "$logfile" "$logfile.$run_id.bak"; fi
            done
            {params.runner:q} exec --cleanenv {params.binds} {input.container:q} \
                /opt/deepvariant/bin/run_deepvariant \
                --model_type=PACBIO \
                --ref={params.reference:q} \
                --reads={params.bam:q} \
                --regions={params.chromosomes:q} \
                --output_vcf={params.vcf:q} \
                --intermediate_results_dir={params.intermediate:q} \
                --logging_dir={params.internal_logs:q} \
                --num_shards={threads} \
                > {log.out:q} 2> {log.err:q}
            """

    rule phase_variants:
        input:
            vcf=rules.deepvariant.output.vcf,
            tbi=rules.deepvariant.output.tbi,
            bam=lambda wc: f"{RESULTS}/01_align/{PHASING_SAMPLE_IDS[wc.phasing_sample]}.chromosomes.bam",
            bai=lambda wc: f"{RESULTS}/01_align/{PHASING_SAMPLE_IDS[wc.phasing_sample]}.chromosomes.bam.bai",
            reference=PHASING_REFERENCE,
            fai=f"{PHASING_REFERENCE}.fai"
        output:
            vcf=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.phased.vcf.gz"
        params:
            whatshap=PHASING_WHATSHAP
        # WhatsHap phase is single-threaded; samples can run independently.
        threads: 1
        resources:
            mem_mb=PHASING_CONFIG.get("whatshap_mem_mb", 16000)
        conda: "../envs/phasing.yaml"
        log:
            out=f"{LOGS}/02_phasing/{{phasing_sample}}.phase.out",
            err=f"{LOGS}/02_phasing/{{phasing_sample}}.phase.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.vcf:q})" "$(dirname {log.out:q})"
            {params.whatshap:q} phase --reference {input.reference:q} \
                --output {output.vcf:q} {input.vcf:q} {input.bam:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule index_phased_vcf:
        input:
            vcf=rules.phase_variants.output.vcf
        output:
            tbi=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.phased.vcf.gz.tbi"
        params:
            tabix=PHASING_TABIX
        conda: "../envs/phasing.yaml"
        log:
            out=f"{LOGS}/02_phasing/{{phasing_sample}}.phased_index.out",
            err=f"{LOGS}/02_phasing/{{phasing_sample}}.phased_index.err"
        shell:
            r"""
            mkdir -p "$(dirname {log.out:q})"
            {params.tabix:q} -f -p vcf {input.vcf:q} > {log.out:q} 2> {log.err:q}
            """

    rule phasing_stats:
        input:
            vcf=rules.phase_variants.output.vcf,
            tbi=rules.index_phased_vcf.output.tbi
        output:
            tsv=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.phasing.stats.tsv",
            gtf=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.phasing.blocks.gtf"
        params:
            whatshap=PHASING_WHATSHAP
        threads: 1
        resources:
            mem_mb=PHASING_CONFIG.get("whatshap_mem_mb", 16000)
        conda: "../envs/phasing.yaml"
        log:
            out=f"{LOGS}/02_phasing/{{phasing_sample}}.phasing_stats.out",
            err=f"{LOGS}/02_phasing/{{phasing_sample}}.phasing_stats.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.tsv:q})" "$(dirname {log.out:q})"
            {params.whatshap:q} stats --tsv {output.tsv:q} --gtf {output.gtf:q} \
                {input.vcf:q} > {log.out:q} 2> {log.err:q}
            """

    rule haplotag_aligned_bam:
        input:
            vcf=rules.phase_variants.output.vcf,
            tbi=rules.index_phased_vcf.output.tbi,
            bam=lambda wc: f"{RESULTS}/01_align/{PHASING_SAMPLE_IDS[wc.phasing_sample]}.chromosomes.bam",
            bai=lambda wc: f"{RESULTS}/01_align/{PHASING_SAMPLE_IDS[wc.phasing_sample]}.chromosomes.bam.bai",
            reference=PHASING_REFERENCE,
            fai=f"{PHASING_REFERENCE}.fai"
        output:
            bam=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.haplotagged.bam"
        params:
            whatshap=PHASING_WHATSHAP,
            # One CPU handles allele assignment; the remainder can compress BAM.
            output_threads=lambda wc, threads: max(0, threads - 1)
        threads: config["threads"].get("haplotag", 4)
        resources:
            mem_mb=PHASING_CONFIG.get("whatshap_mem_mb", 16000)
        conda: "../envs/phasing.yaml"
        log:
            out=f"{LOGS}/02_phasing/{{phasing_sample}}.haplotag.out",
            err=f"{LOGS}/02_phasing/{{phasing_sample}}.haplotag.err"
        shell:
            r"""
            mkdir -p "$(dirname {output.bam:q})" "$(dirname {log.out:q})"
            {params.whatshap:q} haplotag --ignore-linked-read \
                --reference {input.reference:q} --output-threads {params.output_threads} \
                --output {output.bam:q} {input.vcf:q} {input.bam:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule index_haplotagged_bam:
        input:
            bam=rules.haplotag_aligned_bam.output.bam
        output:
            bai=f"{PHASING_ROOT}/{{phasing_sample}}/{{phasing_sample}}.haplotagged.bam.bai"
        params:
            samtools=PHASING_SAMTOOLS,
            extra_threads=lambda wc, threads: max(0, threads - 1)
        threads: config["threads"].get("haplotag", 4)
        conda: "../envs/phasing.yaml"
        log:
            out=f"{LOGS}/02_phasing/{{phasing_sample}}.haplotagged_index.out",
            err=f"{LOGS}/02_phasing/{{phasing_sample}}.haplotagged_index.err"
        shell:
            r"""
            mkdir -p "$(dirname {log.out:q})"
            {params.samtools:q} index -@ {params.extra_threads} {input.bam:q} {output.bai:q} \
                > {log.out:q} 2> {log.err:q}
            """

    rule variant_calls:
        input:
            [f"{PHASING_ROOT}/{alias}/{alias}.deepvariant.vcf.gz{suffix}"
             for alias in selected_aliases for suffix in ("", ".tbi")]

    rule phasing:
        input: PHASING_OUTPUTS
