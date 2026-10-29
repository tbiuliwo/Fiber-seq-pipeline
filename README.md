# Fiber-seq Snakemake pipeline

A modular workflow for PacBio Fiber-seq BAMs. It covers alignment, annotation,
QC, FIRE peaks and optional phasing, methylation and FiberFold analyses.

## Workflow

Rules are in `workflow/rules/`. Numbers organize modules; execution follows
input/output dependencies.

| Step | Module | Function |
|---|---|---|
| 01 | align | Alignment, chromosome filtering and indexing |
| 02 | phasing | Optional variant calling, phasing and BAM haplotagging |
| 03 | annotate | m6A, nucleosome and MSP annotation |
| 04 | qc | Fiber-seq QC, flagstat and alignment statistics |
| 05 | fire | Single-molecule FIRE calls |
| 06 | fire_peaks | Run the official FIRE peak-calling workflow |
| 07 | extract | Export m6A, CpG, nucleosome, MSP and FIRE BED files |
| 08 | methylation | Optional 5mC pileup and indexed bedGraph tracks |
| 09 | fiberhmm | Optional footprint analysis |
| 10 | ctcf | CTCF footprint and motif-direction tracks |
| 11 | fiberfold_inputs | Prepare feature tracks and prediction windows |
| 12 | fiberfold | Optional contact prediction |

## Requirements

- Core: Python 3.11+, Snakemake, pbmm2, samtools, fibertools-rs (`ft`) and pysam.
- FIRE peaks: a separate FIRE checkout, Pixi and its Conda environments.
- Phasing: a DeepVariant PACBIO container, Singularity/Apptainer and WhatsHap.
  `workflow/envs/phasing.yaml` provides the WhatsHap environment.
- Methylation: MethBat, bgzip and tabix. FiberHMM: `fiberhmm-call` with QC support.
- CTCF/FiberFold: the external FiberFold scripts, their Python dependencies,
  bedtools, bigtools and pb-CpG-tools; prediction also needs the trained model.

For CTCF, apply the bundled compatibility patch to a compatible FiberFold checkout
before running its scripts:

```bash
git -C /path/to/fiberFold apply /path/to/Fiber-seq-pipeline/workflow/patches/fiberfold-ctcf.patch
```

## Configuration

Run from the repository root. On a fresh checkout, create a local configuration:

```bash
cp config/example.yaml config/config.yaml
```

Edit sample BAM paths, reference FASTA, software paths, threads and optional
analysis switches in `config/config.yaml`. This local file is ignored by Git.
Reference data, input BAMs, external software and model weights are supplied separately.

- `input_mode: modified`: existing m6A calls; add nucleosome/MSP annotations.
- `input_mode: hifi`: predict m6A from HiFi kinetics, then infer annotations.
- `input_mode: annotated`: reuse existing annotations with aligned input.
- `input_aligned`: choose alignment or sorting of already aligned BAMs.
- `analysis_chromosomes`: chromosomes retained in the BAM and its reference header;
  the helper is `scripts/01filter_chromosomes.py`.

Enable optional analyses through their configuration sections. Phasing is disabled
by default; when enabled, selected samples use haplotagged BAMs before annotation.
Phasing requires a reference FASTA with an existing `.fai` index.
Methylation exports Total 5mC and skips 6mA/5hmC.
Set both `fiberfold.enabled` and `fiberfold.prepare_inputs` to false when not using
the model branch; input preparation also enables CTCF and requires FiberFold scripts.
Configure `fire_peaks.fire_dir` and `fire_peaks.pixi` for the official FIRE checkout.

## Run

Activate an environment containing the required tools, then preview or execute:

```bash
# Preview only.
snakemake -s workflow/fiberseq.smk -n --cores 8

# Run configured targets.
snakemake -s workflow/fiberseq.smk --cores 16

# Select a branch and its required upstream inputs.
snakemake -s workflow/fiberseq.smk fire_peaks -n --cores 8
```

Other aggregate targets are `phasing`, `variant_calls`, `methylation`,
`prepare_ctcf` and `prepare_fiberfold_inputs`; enable their configuration switches
before selecting them. Add `--sdm conda` to use the bundled WhatsHap environment.

Results and logs are written to `results/` and `logs/`, with numbered module
subdirectories. The official FIRE layout is retained under `results/06_fire_peaks/`.
Scientific parameters and input changes may require recomputing downstream results.
