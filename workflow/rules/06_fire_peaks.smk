# 06. Run the official FIRE Snakefile; calculations stay in the FIRE checkout.
import json
import re
import shutil
import subprocess
import tomllib

PEAK_CONFIG = config.get("fire_peaks", {})
PEAK_ROOT = f"{RESULTS}/06_fire_peaks"
PEAK_VERSION = PEAK_CONFIG.get("version", "v0.2")
PEAK_FIRE_DIR = Path(
    PEAK_CONFIG.get("fire_dir", "/data/tenghaohan/software/FIRE")
).resolve()
PEAK_PIXI = PEAK_CONFIG.get(
    "pixi", shutil.which("pixi") or str(Path.home() / ".pixi/bin/pixi")
)
PEAK_SAMTOOLS = str(Path(config.get("tools", {}).get("bin_dir") or "") / "samtools")
PEAK_CHROMS = "|".join(re.escape(chrom) for chrom in ANALYSIS_CHROMOSOMES)


# Write the official sample table/config and index a local reference symlink.
rule prepare_fire_peaks:
    input:
        bams=expand(f"{RESULTS}/05_fire/{{sample}}.fire.bam", sample=SAMPLES),
        bais=expand(f"{RESULTS}/05_fire/{{sample}}.fire.bam.bai", sample=SAMPLES),
        reference=PEAK_CONFIG.get("reference_fasta", REFERENCE),
        masks=PEAK_CONFIG.get("exclude", []),
        project=str(PEAK_FIRE_DIR / "pixi.toml"),
    output:
        config=f"{PEAK_ROOT}/config.json",
        manifest=f"{PEAK_ROOT}/manifest.tsv",
        reference=f"{PEAK_ROOT}/reference/genome.fa",
        fai=f"{PEAK_ROOT}/reference/genome.fa.fai",
    params:
        settings=PEAK_CONFIG,
        chromosomes=PEAK_CHROMS,
        ref_name=config.get("reference_name", "hg38"),
        samtools=PEAK_SAMTOOLS,
    log:
        f"{LOGS}/06_fire_peaks/prepare.log",
    run:
        project = tomllib.loads(Path(input.project).read_text())
        version = project.get("workspace", project.get("project", {}))["version"]
        full_version = PEAK_VERSION == f"v{version}"
        if PEAK_VERSION != "v" + ".".join(version.split(".")[:2]) and not full_version:
            raise ValueError(
                f"fire_peaks.version does not match installed FIRE {version}"
            )
        if (
            PEAK_CONFIG.get("seed", 42) != 42
            or PEAK_CONFIG.get("window_size", 200) != 200
        ):
            raise ValueError("Official FIRE uses seed 42 and a 200-bp window.")
        reference = Path(output.reference).absolute()
        reference.unlink(missing_ok=True)
        reference.symlink_to(Path(input.reference).resolve())
        with open(log[0], "w") as handle:
            subprocess.run(
                [params.samtools, "faidx", str(reference)],
                check=True,
                stdout=handle,
                stderr=subprocess.STDOUT,
            )
        with open(output.manifest, "w") as handle:
            handle.write("sample\tbam\n")
            for sample, bam in zip(SAMPLES, input.bams):
                handle.write(f"{sample}\t{Path(bam).resolve()}\n")
        official = {
            "ref": str(reference),
            "ref_name": params.ref_name,
            "manifest": str(Path(output.manifest).absolute()),
            "keep_chromosomes": params.chromosomes,
            "full-version": full_version,
            "excludes": [str(Path(p).resolve()) for p in input.masks],
        }
        for local, remote, default in [
            ("max_fdr", "max_peak_fdr", 0.05),
            ("min_fire_fraction", "min_frac_accessible", 0.10),
            ("min_coverage", "min_coverage", 4),
            ("sd_coverage", "coverage_within_n_sd", 5),
            ("min_fire_fdr", "min_fire_fdr", 0.10),
            ("min_msp", "min_msp", 10),
            ("min_ave_msp_size", "min_ave_msp_size", 10),
        ]:
            official[remote] = PEAK_CONFIG.get(local, default)
        Path(output.config).write_text(json.dumps(official, indent=2) + "\n")


# One outer job; the official Snakemake schedules all FIRE calculations inside it.
rule run_official_fire:
    input:
        config=rules.prepare_fire_peaks.output.config,
        manifest=rules.prepare_fire_peaks.output.manifest,
        reference=rules.prepare_fire_peaks.output.reference,
        fai=rules.prepare_fire_peaks.output.fai,
        project=str(PEAK_FIRE_DIR / "pixi.toml"),
        sources=lambda wc: sorted(
            str(p) for p in (PEAK_FIRE_DIR / "workflow").rglob("*") if p.is_file()
        ),
    output:
        peaks=expand(
            f"{PEAK_ROOT}/results/{{sample}}/{{sample}}-fire-{PEAK_VERSION}-peaks.bed.gz",
            sample=SAMPLES,
        ),
        indexes=expand(
            f"{PEAK_ROOT}/results/{{sample}}/{{sample}}-fire-{PEAK_VERSION}-peaks.bed.gz.tbi",
            sample=SAMPLES,
        ),
        fdr=expand(
            f"{PEAK_ROOT}/results/{{sample}}/additional-outputs-{PEAK_VERSION}/fire-peaks/{{sample}}-{PEAK_VERSION}-fire-score-to-fdr.tbl",
            sample=SAMPLES,
        ),
        complete=f"{PEAK_ROOT}/.complete",
    threads: config["threads"].get("fire_peaks", 4)
    resources:
        mem_mb=PEAK_CONFIG.get("mem_mb", 32000),
    params:
        pixi=PEAK_PIXI,
        directory=str(Path(PEAK_ROOT).absolute()),
        config=str(Path(rules.prepare_fire_peaks.output.config).absolute()),
        conda_prefix=str(
            Path(
                PEAK_CONFIG.get("conda_prefix", f"{PEAK_ROOT}/.snakemake/conda")
            ).absolute()
        ),
    log:
        f"{LOGS}/06_fire_peaks/official.log",
    shell:
        r"""
        {params.pixi:q} run --manifest-path {input.project:q} fire all \
            --directory {params.directory:q} --configfile {params.config:q} \
            --cores {threads} --resources mem_mb={resources.mem_mb} \
            --software-deployment-method conda --conda-prefix {params.conda_prefix:q} \
            > {log:q} 2>&1
        touch {output.complete:q}
        """


rule fire_peaks:
    input:
        rules.run_official_fire.output.complete,
