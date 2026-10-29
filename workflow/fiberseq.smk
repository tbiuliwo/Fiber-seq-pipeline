"""Fiber-seq workflow. Run from the repository root."""
import re
import shlex
from pathlib import Path


# 1. Load configuration
configfile: "config/config.yaml"

SAMPLES_CONFIG = config["samples"]
SAMPLES = list(SAMPLES_CONFIG)

MODE = config["input_mode"]
ALIGNED = config["input_aligned"]
COORD = config["extract_coordinates"]
ANALYSIS_CHROMOSOMES = config.get("analysis_chromosomes", [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])


RESULTS = config["results_dir"].rstrip("/")
LOGS = config["logs_dir"].rstrip("/")
REFERENCE = config["reference"]
METHYLATION = config.get("methylation", {}).get("enabled", False)

ALIGN_REFERENCE = (REFERENCE
                   if str(REFERENCE).endswith(".mmi")
                   else f"{RESULTS}/01_align/reference/reference.hifi.mmi")

## Software
PYTHON = config.get("tools", {}).get("python", "python")

    
## Optional analyses
FIBERHMM = config.get("fiberhmm", {}).get("enabled", False)
FIBERFOLD = config.get("fiberfold", {}).get("enabled", False)
FOLD_INPUTS = config.get("fiberfold", {}).get("prepare_inputs", False) or FIBERFOLD
CTCF = config.get("ctcf", {}).get("enabled", False) or FOLD_INPUTS
PHASING_CONFIG = config.get("phasing", {})
PHASING = PHASING_CONFIG.get("enabled", False)
PHASING_SAMPLE_IDS = {}
PHASED_SAMPLES = {}
PHASING_OUTPUTS = []
if PHASING:
    aliases = PHASING_CONFIG.get("sample_aliases", {})
    for sample in SAMPLES:
        alias = aliases.get(sample, sample)
        if not isinstance(alias, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", alias):
            raise ValueError(f"Invalid phasing output sample name: {alias!r}")
        if alias in PHASING_SAMPLE_IDS:
            raise ValueError(f"Duplicate phasing output sample name: {alias}")
        PHASING_SAMPLE_IDS[alias] = sample
    selected = PHASING_CONFIG.get("samples")
    if selected is None:
        selected = list(PHASING_SAMPLE_IDS)
    if not isinstance(selected, list) or not selected:
        raise ValueError("phasing.samples must be a non-empty list, or null for all samples")
    selected_aliases = []
    for sample in selected:
        alias = aliases.get(sample, sample)
        if alias not in PHASING_SAMPLE_IDS:
            raise ValueError(f"Unknown phasing sample: {sample}")
        if alias not in selected_aliases:
            selected_aliases.append(alias)
    PHASED_SAMPLES = {PHASING_SAMPLE_IDS[alias]: alias for alias in selected_aliases}
    PHASING_OUTPUTS = [
        f"{RESULTS}/02_phasing/{alias}/{alias}.{suffix}"
        for alias in selected_aliases
        for suffix in (
            "deepvariant.vcf.gz", "deepvariant.vcf.gz.tbi",
            "phased.vcf.gz", "phased.vcf.gz.tbi",
            "haplotagged.bam", "haplotagged.bam.bai",
            "phasing.stats.tsv", "phasing.blocks.gtf",
        )
    ]
if CTCF and not set(config["ctcf"]["chromosomes"]).issubset(ANALYSIS_CHROMOSOMES):
    raise ValueError("ctcf.chromosomes must be a subset of analysis_chromosomes")

# 2. Input BAM paths
def sample_bam(wildcards):
    return SAMPLES_CONFIG[wildcards.sample]

def aligned_bam(wildcards):
    return f"{RESULTS}/01_align/{wildcards.sample}.chromosomes.bam"

def downstream_aligned_bam(wildcards):
    """Use haplotagged alignments only for samples selected for optional phasing."""
    alias = PHASED_SAMPLES.get(wildcards.sample)
    if alias is not None:
        return f"{RESULTS}/02_phasing/{alias}/{alias}.haplotagged.bam"
    return aligned_bam(wildcards)

def annotated_bam(wildcards):
    if MODE == "annotated":
        return downstream_aligned_bam(wildcards)
    return f"{RESULTS}/03_annotate/{wildcards.sample}.fiber.bam"



# 3. Final outputs
FINAL = list(PHASING_OUTPUTS)

## Basic quality control
FINAL += expand(
    f"{RESULTS}/04_qc/{{sample}}.{{file}}",
    sample=SAMPLES,
    file=["qc.tsv", "flagstat.txt", "alignment.stats.txt"]
)

## Single-molecule FIRE and feature extraction
FINAL += expand(
    f"{RESULTS}/05_fire/{{sample}}.{{suffix}}",
    sample=SAMPLES,
    suffix=[
        "fire.bam.bai",
        "fire.bed",
        "fire.features.tsv"
    ]
)
if config.get("fire_peaks", {}).get("enabled", True):
    FINAL += [f"{RESULTS}/06_fire_peaks/.complete"]
FINAL += expand(
    f"{RESULTS}/07_extract/{{sample}}.{{feature}}.bed.gz",
    sample=SAMPLES,
    feature=["m6a", "cpg", "nuc", "msp", "fire"]
)

## Optional footprint analysis
if FIBERHMM:
    FINAL += expand(f"{RESULTS}/09_fiberhmm/{{sample}}.footprints.bam.bai", sample=SAMPLES)
    FINAL += expand(
        f"{RESULTS}/09_fiberhmm/qc/{{sample}}.footprints.qc.{{extension}}",
        sample=SAMPLES,
        extension=["json", "tsv", "png", "pdf"],
    )

## CTCF, remaining model inputs, then contact prediction
if CTCF:
    FINAL += expand(f"{RESULTS}/10_ctcf/direction/{{chrom}}.bw", chrom=config["ctcf"]["chromosomes"])
    FINAL += expand(f"{RESULTS}/10_ctcf/footprint/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=config["ctcf"]["chromosomes"])
if FOLD_INPUTS:
    FINAL += [f"{RESULTS}/11_fiberfold_inputs/windows.bed.gz"]
    for feature in ("fire", "cpg"):
        FINAL += expand(f"{RESULTS}/11_fiberfold_inputs/{feature}/{{sample}}/{{chrom}}.bw", sample=SAMPLES, chrom=config["ctcf"]["chromosomes"])
if FIBERFOLD:
    FINAL += expand(f"{RESULTS}/12_fiberfold/{{sample}}.predictions.pkl", sample=SAMPLES)

## Step 08: independent site-level 5mC analysis and plotting tracks
if METHYLATION:
    for suffix in ("5mC.bed.gz", "5mC.bed.gz.tbi", "summary.json", "5mC.bedgraph.gz", "5mC.bedgraph.gz.tbi"):
        FINAL += expand(f"{RESULTS}/08_methylation/{{sample}}.{suffix}", sample=SAMPLES)

# 4. Main rule
rule all:
    input: FINAL

# 5. Load modules in analysis order (execution follows file dependencies).
# Result and log directory names match the module filenames.
include: "rules/01_align.smk"
include: "rules/02_phasing.smk"
include: "rules/03_annotate.smk"
include: "rules/04_qc.smk"
include: "rules/05_fire.smk"
include: "rules/06_fire_peaks.smk"
include: "rules/07_extract.smk"
include: "rules/08_methylation.smk"
include: "rules/09_fiberhmm.smk"
include: "rules/10_ctcf.smk"
include: "rules/11_fiberfold_inputs.smk"
include: "rules/12_fiberfold.smk"
