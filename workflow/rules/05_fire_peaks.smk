# 5. Population-level FIRE peaks, with shuffled-fiber FDR calibration.
PEAK_ROOT = f"{RESULTS}/05_fire_peaks"
PEAK_CONFIG = config.get("fire_peaks", {})
PEAK_SCRIPT = str(Path(workflow.basedir).parent / "scripts/fire_peak_background.py")
PEAK_FDR_SCRIPT = str(Path(workflow.basedir).parent / "scripts/fire_peak_fdr.py")


rule fire_peak_background:
    input:
        bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
        bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai",
        script=PEAK_SCRIPT,
        exclude=PEAK_CONFIG.get("exclude", [])
    output:
        shuffled=f"{PEAK_ROOT}/{{sample}}.shuffled.bed.gz",
        excluded=f"{PEAK_ROOT}/{{sample}}.exclude.bed.gz",
        sizes=f"{PEAK_ROOT}/{{sample}}.genome.sizes",
        summary=f"{PEAK_ROOT}/{{sample}}.coverage.json",
        limits=f"{PEAK_ROOT}/{{sample}}.coverage_limits.tsv"
    threads: config["threads"].get("fire_peaks", 4)
    params:
        seed=PEAK_CONFIG.get("seed", 42),
        min_coverage=PEAK_CONFIG.get("min_coverage", 4),
        sd_coverage=PEAK_CONFIG.get("sd_coverage", 5)
    log:
        out=f"{LOGS}/05_fire_peaks/{{sample}}.background.out",
        err=f"{LOGS}/05_fire_peaks/{{sample}}.background.err"
    shell:
        r"""
        mkdir -p "$(dirname {output.shuffled:q})" "$(dirname {log.out:q})"
        {PYTHON:q} {input.script:q} --bam {input.bam:q} --threads {threads} \
            --seed {params.seed} --min-coverage {params.min_coverage} --sd-coverage {params.sd_coverage} \
            --shuffled {output.shuffled:q} --excluded {output.excluded:q} --sizes {output.sizes:q} \
            --summary {output.summary:q} --limits {output.limits:q} --exclude {input.exclude:q} \
            > {log.out:q} 2> {log.err:q}
        """


rule fire_peak_fdr:
    input:
        bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
        bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai",
        shuffled=rules.fire_peak_background.output.shuffled,
        limits=rules.fire_peak_background.output.limits,
        script=PEAK_FDR_SCRIPT
    output:
        fdr=f"{PEAK_ROOT}/{{sample}}.fdr.tsv"
    threads: config["threads"].get("fire_peaks", 4)
    resources:
        mem_mb=PEAK_CONFIG.get("mem_mb", 16000)
    log:
        out=f"{LOGS}/05_fire_peaks/{{sample}}.fdr.out",
        err=f"{LOGS}/05_fire_peaks/{{sample}}.fdr.err"
    shell:
        r"""
        mkdir -p "$(dirname {output.fdr:q})" "$(dirname {log.out:q})"
        (
            ft -t {threads} pileup {input.bam:q} --fire-filter -F 260 \
                --fiber-coverage --no-msp --no-nuc --shuffle {input.shuffled:q} \
            | {PYTHON:q} {input.script:q} --limits {input.limits:q} --out {output.fdr:q}
        ) > {log.out:q} 2> {log.err:q}
        """


rule call_fire_peaks:
    input:
        bam=f"{RESULTS}/04_fire/{{sample}}.fire.bam",
        bai=f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai",
        fdr=rules.fire_peak_fdr.output.fdr,
        limits=rules.fire_peak_background.output.limits
    output:
        peaks=f"{PEAK_ROOT}/{{sample}}.peaks.bed.gz",
        tbi=f"{PEAK_ROOT}/{{sample}}.peaks.bed.gz.tbi"
    threads: config["threads"].get("fire_peaks", 4)
    resources:
        mem_mb=PEAK_CONFIG.get("mem_mb", 16000)
    params:
        max_fdr=PEAK_CONFIG.get("max_fdr", 0.05),
        min_fire_fraction=PEAK_CONFIG.get("min_fire_fraction", 0.1),
        window_size=PEAK_CONFIG.get("window_size", 200)
    log:
        out=f"{LOGS}/05_fire_peaks/{{sample}}.peaks.out",
        err=f"{LOGS}/05_fire_peaks/{{sample}}.peaks.err"
    shell:
        r"""
        mkdir -p "$(dirname {output.peaks:q})" "$(dirname {log.out:q})"
        (
            read -r MIN_COV MAX_COV < {input.limits:q}
            PEAK_TMP=$(mktemp)
            trap 'rm -f "$PEAK_TMP"' EXIT

            
            ft -t {threads} call-peaks {input.bam:q} --fire-filter -F 260 \
                --fdr-table {input.fdr:q} \
                --min-cov "$MIN_COV" --max-cov "$MAX_COV" \
                --max-fdr {params.max_fdr} --min-fire-frac-filter {params.min_fire_fraction} \
                --window-size {params.window_size} \
                --out "$PEAK_TMP"

            awk 'BEGIN {{FS=OFS="\t"}} /^#/ || $NF == "true"' "$PEAK_TMP" \
                | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | bgzip -@ {threads} > {output.peaks:q}
            tabix -f -p bed {output.peaks:q}
        ) > {log.out:q} 2> {log.err:q}
        """


rule fire_peaks:
    input:
        expand(f"{PEAK_ROOT}/{{sample}}.peaks.bed.gz", sample=SAMPLES),
        expand(f"{PEAK_ROOT}/{{sample}}.peaks.bed.gz.tbi", sample=SAMPLES),
        expand(f"{PEAK_ROOT}/{{sample}}.fdr.tsv", sample=SAMPLES)

