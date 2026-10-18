"""Fiber-seq workflow. Run from the repository root."""
import shlex
from pathlib import Path


# 1. Load configuration
configfile: "config/config.yaml"

SAMPLES_CONFIG = config["samples"]
SAMPLES = list(SAMPLES_CONFIG)

MODE = config["input_mode"]
ALIGNED = config["input_aligned"]
COORD = config["extract_coordinates"]

RESULTS = config["results_dir"].rstrip("/")
LOGS = config["logs_dir"].rstrip("/")
REFERENCE = config["reference"]
METHYLATION = config.get("methylation", {}).get("enabled", False)

ALIGN_REFERENCE = REFERENCE if str(REFERENCE).endswith(".mmi") else f"{RESULTS}/01_align/reference/reference.hifi.mmi"

## Software
PYTHON = config.get("tools", {}).get("python", "python")
TOOL_BIN = config.get("tools", {}).get("bin_dir")
if TOOL_BIN:
    shell.prefix("export PATH=" + shlex.quote(TOOL_BIN) + ":$PATH; ")
    
## Optional analyses
FIBERHMM = config.get("fiberhmm", {}).get("enabled", False)
FIBERFOLD = config.get("fiberfold", {}).get("enabled", False)
FOLD_INPUTS = config.get("fiberfold", {}).get("prepare_inputs", False) or FIBERFOLD
CTCF = config.get("ctcf", {}).get("enabled", False) or FOLD_INPUTS

# 2. Input BAM paths
def sample_bam(wildcards):
    return SAMPLES_CONFIG[wildcards.sample]

def aligned_bam(wildcards):
    return f"{RESULTS}/01_align/{wildcards.sample}.aligned.bam"

def annotated_bam(wildcards):
    if MODE == "annotated":
        return aligned_bam(wildcards)
    return f"{RESULTS}/02_annotate/{wildcards.sample}.fiber.bam"



# 3. Final outputs
FINAL = []

## Basic quality control
FINAL += expand(
    f"{RESULTS}/03_qc/{{sample}}.{{file}}",
    sample=SAMPLES,
    file=["qc.tsv", "flagstat.txt", "alignment.stats.txt"]
)

## Single-molecule FIRE and feature extraction
FINAL += expand(f"{RESULTS}/04_fire/{{sample}}.fire.bam.bai", sample=SAMPLES)
if config.get("fire_peaks", {}).get("enabled", True):
    FINAL += expand(f"{RESULTS}/05_fire_peaks/{{sample}}.peaks.bed.gz", sample=SAMPLES)
    FINAL += expand(f"{RESULTS}/05_fire_peaks/{{sample}}.peaks.bed.gz.tbi", sample=SAMPLES)
    FINAL += expand(f"{RESULTS}/05_fire_peaks/{{sample}}.fdr.tsv", sample=SAMPLES)
FINAL += expand(
    f"{RESULTS}/06_extract/{{sample}}.{{feature}}.bed.gz",
    sample=SAMPLES,
    feature=["m6a", "cpg", "nuc", "msp", "fire"]
)

## Optional footprint analysis
if FIBERHMM:
    FINAL += expand(f"{RESULTS}/07_fiberhmm/{{sample}}.footprints.bam.bai", sample=SAMPLES)

## CTCF, remaining model inputs, then contact prediction
if CTCF:
    FINAL += expand(f"{RESULTS}/08_ctcf/direction/{{chrom}}.bw", chrom=config["ctcf"]["chromosomes"])
    FINAL += expand(f"{RESULTS}/08_ctcf/footprint/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=config["ctcf"]["chromosomes"])
if FOLD_INPUTS:
    FINAL += [f"{RESULTS}/09_fiberfold_inputs/windows.bed.gz"]
    for feature in ("fire", "cpg"):
        FINAL += expand(f"{RESULTS}/09_fiberfold_inputs/{feature}/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=config["ctcf"]["chromosomes"])
if FIBERFOLD:
    FINAL += expand(f"{RESULTS}/10_fiberfold/{{sample}}.predictions.pkl", sample=SAMPLES)

## Independent methylation analysis, disabled by default
if METHYLATION:
    for suffix in ("5mC.bed.gz", "5mC.bed.gz.tbi", "5hmC.bed.gz", "5hmC.bed.gz.tbi", "summary.json"):
        FINAL += expand(f"{RESULTS}/optional/methylation/{{sample}}.{suffix}", sample=SAMPLES)


# 4. Main rule
rule all:
    input: FINAL

# 5. Load modules in analysis order (execution follows file dependencies).
# Result and log directory names match the module filenames.
include: "rules/01_align.smk"
include: "rules/02_annotate.smk"
include: "rules/03_qc.smk"
include: "rules/04_fire.smk"
include: "rules/05_fire_peaks.smk"
include: "rules/06_extract.smk"
include: "rules/07_fiberhmm.smk"
include: "rules/08_ctcf.smk"
include: "rules/09_fiberfold_inputs.smk"
include: "rules/10_fiberfold.smk"
if METHYLATION:
    include: "rules/optional/methylation.smk"
