# Fiber-seq Snakemake pipeline

A modular workflow for PacBio Fiber-seq BAMs, including instrument-called
m6A/5mC input. Run commands from this repository's root.

## Analysis route

```text
shared reference indexing + pbmm2 HIFI alignment + sorting,
     or samtools sorting for aligned input
  -> m6A prediction / nucleosome + MSP annotation / reuse existing annotations
  -> basic ft qc (optional autocorrelation)
  -> single-molecule FIRE -> indexed BAM
                         -> shuffled-fiber FDR calibration -> population FIRE peaks
                         -> m6A, CpG 5mC, nucleosome, MSP and FIRE BEDs
                         -> optional FiberHMM footprint BAM + index

alignment -> samtools flagstat + stats
indexed FIRE BAM -> automatic FIRE/CpG tracks + prediction windows
all four feature tracks + model -> optional FiberFold predictions
indexed FIRE BAM + CTCF motifs -> footprint fractions + shared direction bigWigs
indexed FIRE BAMs -> optional prepare_fire_peaks -> official FIRE config + manifest
```

Rule filenames are numbered in analysis order. Snakemake schedules actual
execution from input/output dependencies; basic QC can run alongside FIRE.
Result and log directory names match the corresponding rule filenames.

## Rule organization

| Rule file | Function |
|---|---|
| `01_align.smk` | Reference indexing, alignment and sorting |
| `02_annotate.smk` | m6A / nucleosome / MSP annotation |
| `03_qc.smk` | Basic QC on aligned and annotated BAMs |
| `04_fire.smk` | Single-molecule FIRE and BAM indexing |
| `05_fire_peaks.smk` | Shuffled-background FDR and population FIRE peaks; optional official-input export |
| `06_extract.smk` | Export feature BEDs |
| `07_fiberhmm.smk` | Optional general footprint analysis |
| `08_ctcf.smk` | CTCF-specific feature tracks |
| `09_fiberfold_inputs.smk` | FIRE/CpG tracks and model windows |
| `10_fiberfold.smk` | Contact matrix prediction |
| `optional/methylation.smk` | Independent 5mC/5hmC analysis |

Peak calling is enabled by default and can also be selected with the `fire_peaks`
target. Official-pipeline input export remains an explicit `prepare_fire_peaks` target. QC reports
statistics and does not act as a pass/fail gate for downstream analyses.
The optional methylation module defaults to off and is not a FiberFold
dependency; enabling it requires `methbat` and `threads.methylation`.

## Software

Activate an analysis environment with `pbmm2`, `samtools`, `ft`, Python and
`pysam`. FIRE peak calling additionally requires `bedtools`, `bgzip`, `tabix`
and an `ft` version supporting `call-peaks --fdr-table` (the local version is
0.13.1). Background construction and FDR accumulation use streaming Python
scripts and do not require the full official FIRE workflow or its R environment.
FiberFold input preparation also requires `pyBigWig` and
`pb-CpG-tools` (3.0.0 is installed in the local fiberseq environment). Optional branches require `fiberhmm-call` or the FiberFold Python
stack and model. Snakemake may be installed separately.

On this machine:

```bash
conda activate fiberseq
snakemake -s workflow/fiberseq.smk --list-rules
snakemake -s workflow/fiberseq.smk -n -p --cores 8
snakemake -s workflow/fiberseq.smk -p --cores 16
```

`tools.bin_dir` in the config pins the analysis executables to the existing
`fiberseq` environment, even when the Snakemake Pixi launcher changes PATH.
Change this path on other machines, or set it to null and activate your own
environment. `tools.python` can select a different Python executable.

## Configuration

Edit `config/config.yaml`. The local defaults are the configured WT, R354L and 4A downlinked BAMs
and the supplied GRCh38.p14 FASTA. These files are external inputs and are
not distributed with this repository. Use your own paths to reproduce.

| input_mode | Required input | Annotation step |
|---|---|---|
| `modified` | m6A MM/ML calls already present | `ft add-nucleosomes` |
| `hifi` | HiFi kinetics fi/fp/ri/rp | `ft predict-m6a`, including nucleosome/MSP inference |
| `annotated` | m6A and Ma/MA or legacy nucleosome/MSP tags | Reuse annotations |

Set `input_aligned: false` to run `pbmm2 align --preset HIFI --sort` using
`reference` (FASTA or a HIFI-compatible pbmm2 .mmi). Already aligned input
is sorted into a new BAM in the results directory. Use already aligned input for `annotated` mode to avoid remapping interval annotations.

All samples in one run share an input mode. Use separate configs for mixed
input types. FASTQ and subreads are not substitutes for HiFi kinetics BAMs.
The default `modified` mode uses existing m6A calls rather than re-predicting
them. This version runs directly from the supplied configuration without
custom parameter validation or BAM prechecks.

`extract_coordinates` controls m6A/CpG/nucleosome/MSP BED coordinates;
FIRE extraction uses reference coordinates. CpG BED is per-molecule 5mC,
not a depth-aggregated methylation bigWig. 5hmC is retained in BAM tags but
is not separately exported by this workflow.

## Outputs

- `01_align/`: sorted aligned BAM; `reference/` contains the reusable pbmm2 index.
- `02_annotate/`: annotated BAM, for hifi/modified modes.
- `03_qc/`: ft QC table (with autocorrelation by default), flagstat and stats.
- `04_fire/`: single-molecule FIRE BAM and BAI.
- `05_fire_peaks/`: per-sample `.peaks.bed.gz` + `.tbi`, `.fdr.tsv`, shuffled
  fibers, exclusion masks and coverage summaries; optional official config/manifest.
- `06_extract/`: m6a/cpg/nuc/msp/fire `.bed.gz` files.
- `07_fiberhmm/`: optional footprint BAM and BAI.
- `08_ctcf/`: CTCF motif BED, footprint tables/fractions and chromosome bigWigs.
- `09_fiberfold_inputs/`: FIRE/CpG bigWigs, raw CpG scores and prediction windows.
- `10_fiberfold/`: optional prediction pickle.
- `optional/methylation/`: independent optional 5mC/5hmC analysis.

All paths are relative to `results_dir`. Commands have stdout/stderr logs
under `logs_dir`. QC gathers statistics; it does not filter reads or create
a separate “QC-passed BAM”. Full runs produce large BAM and BED files.

## Optional branches

### FiberHMM

Set `fiberhmm.enabled: true`. The independent branch uses bundled Hia5
PacBio models and `fiberhmm-call`. It does not replace core ft annotations
or create the CTCF-specific bigWig needed by FiberFold. Automatic secondary
QC is disabled for this branch; core ft QC remains available.

### FiberFold

All four model features and prediction windows are generated from the BAM
pipeline automatically. No `fiberfold.inputs` path mapping is needed.

```yaml
fiberfold:
  enabled: false          # true to execute model prediction too
  prepare_inputs: true    # prepare all model inputs without prediction
  window_step: 2097152
  script: /path/to/fiberFold/predict.py
  model: /path/to/model.ckpt
```

To prepare only the inputs and their upstream dependencies:

```bash
snakemake -s workflow/fiberseq.smk prepare_fiberfold_inputs -n --cores 8
snakemake -s workflow/fiberseq.smk prepare_fiberfold_inputs --cores 8
```

To run prediction, set `fiberfold.enabled: true` and run the main workflow.
Enabling prediction automatically enables input preparation and CTCF rules.
FIRE/CpG tracks are generated per sample; motif direction and windows are shared.
The selected chromosomes come from `ctcf.chromosomes` and must exist in the
FASTA used for alignment. See `docs/fiberfold_inputs.md` for calculation details.

### CTCF features

`ctcf.enabled: true` adds CTCF outputs to the default targets. Enabling
FiberFold also activates these rules. To run only the CTCF branch and its
upstream BAM dependencies:

```bash
snakemake -s workflow/fiberseq.smk prepare_ctcf -n --cores 8
snakemake -s workflow/fiberseq.smk prepare_ctcf --cores 8
```

The default source is the hg38 JASPAR 2022 MA0139.1 track used by the
[FiberFold authors](https://github.com/altemoselab/fiberFold). You can set
`ctcf.motifs` to a local BED6 or JASPAR seven-column TSV(.gz). Coordinates
must match `ctcf.reference_fasta`; the bundled module definition requires
19-bp motifs. Default chromosomes are chr1–chr22 and chrX. Update the list
if predicting other chromosomes. No de novo motif scanning is performed.

Rules sort motifs and retain the higher-scoring instance when overlapping
19-bp motifs compete, separately on each chromosome. They then run
`ft footprint` on each indexed FIRE BAM, with modules [0,6), [6,13), [13,19).
The encoder reads the current table by column name rather than old fixed
column positions. It follows the author's encodeFP.py reference-order A/T
module-eligibility convention on both strands, for compatibility with the
released model. A molecule qualifies if the motif is MSP-spanning (bit 0)
and every A/T-containing module has its footprint bit set. The score is
qualifying molecules divided by all spanning molecules, not just MSPs.
Do not divide this fraction by coverage a second time.

Outputs:

```text
08_ctcf/reference/CTCF.MA0139.1.bed
08_ctcf/<sample>.footprints.tsv.gz
08_ctcf/<sample>.scores.tsv.gz
08_ctcf/footprint/<sample>/chr1.bw ...
08_ctcf/direction/chr1.bw ...
```

Direction is +1 on forward motifs, -1 on reverse motifs and 0 elsewhere;
it is shared across samples. Footprint tracks range from 0 to 1 and are
sample-specific. Background and uncovered sites have value 0; the score
table retains spanning counts to distinguish zero coverage from no footprint.
These motif-based scores are not experimental confirmation of CTCF binding.

`scripts/ctcf_tracks.py` writes zero background without chromosome-sized
arrays. It does not arbitrarily discard the last chromosomes. FASTA indexes
are written inside results, leaving the source FASTA directory untouched.

### Whole-genome FIRE peaks

The default workflow includes population-level peak calling. To run this branch
and its upstream dependencies explicitly:

```bash
snakemake -s workflow/fiberseq.smk fire_peaks -n --cores 8
snakemake -s workflow/fiberseq.smk fire_peaks --cores 8
```

`fire_peaks.enabled: false` removes peak calling from the default targets;
the explicit target remains available. Outputs are:

```text
05_fire_peaks/<sample>.peaks.bed.gz
05_fire_peaks/<sample>.peaks.bed.gz.tbi
05_fire_peaks/<sample>.fdr.tsv
05_fire_peaks/<sample>.coverage.json
```

The background script extracts eligible fiber spans with `--fire-filter -F 260`
and calculates a bp-weighted median coverage. Fibers overlapping low/high-depth
regions by at least 20% are excluded from the shuffle. The retained fibers are
shuffled on their original chromosomes using seed 42, avoiding uncovered regions
and the configured assembly-gap, blacklist and segmental-duplication BEDs.

`ft pileup --fiber-coverage --no-msp --no-nuc --shuffle` then calculates real
and shuffled FIRE scores with the same read filters. A streaming script builds
score histograms weighted by genomic interval length, using 0.01-wide score bins;
at each threshold, FDR is
cumulative shuffled bp divided by cumulative real bp, capped at 1. Thresholds
are retained without the lossy table simplification seen in the local native
caller. Finally, `ft call-peaks --fdr-table` calls and merges peaks using this
table. Defaults are maximum FDR 0.05, minimum FIRE fraction 0.1 and a 200-bp
local-maxima window. This fraction is an additional filter, not a replacement
for FDR calibration. Peaks flagged outside the coverage limits are removed from
the final BED. Logs are written under `logs_dir/05_fire_peaks/`.

Whole-chromosome pileups and peak calling use substantial memory. Each job
declares `fire_peaks.mem_mb` (128000 MB by default); use
`--resources mem_mb=256000` to cap the total memory scheduled for these jobs
at 256000 MB, adjusting the value to your machine.

The local `ft` 0.13.1 `call-peaks` implementation ignores its `--shuffled`
argument and generates an internal background instead. Supplying our externally
calibrated `--fdr-table` avoids that behavior; see `docs/validation.md` for the
reproducing fixture. The standalone single-molecule FIRE BED from extraction
remains a different output from these population peaks.

This integrated branch does not reproduce every output or implementation detail
of the full official FIRE workflow: wide peaks, track hubs and its additional
reports are not generated here. To export inputs for that workflow:

```bash
snakemake -s workflow/fiberseq.smk prepare_fire_peaks --cores 8
```

This optional target creates an absolute-path `manifest.tsv` with `sample` and
`bam` columns and a FIRE `config.yaml` inside `05_fire_peaks/`.
Set `fire_peaks.reference_fasta` to a FASTA matching the alignment reference;
an .mmi is insufficient for the official export. Run the official pipeline in
a separate working directory after installing its per-rule dependencies:

```bash
mkdir -p results/fire_peaks_run
cd results/fire_peaks_run
pixi run --manifest-path /data/tenghaohan/software/FIRE/pixi.toml fire \
  --configfile /absolute/path/to/results/05_fire_peaks/config.yaml \
  --software-deployment-method conda --cores 16
```

Whole-genome peak calling needs adequate effective coverage. Synthetic fixtures
can validate file dependencies, null-background handling and the numerical
calling path; they do not establish biological accuracy, full-sample peak
quality or equivalence to the official workflow. Exporting configuration does
not mean the full official FIRE pipeline is installed or executed.

## Validation

See `docs/validation.md` for execution evidence and remaining limitations.
