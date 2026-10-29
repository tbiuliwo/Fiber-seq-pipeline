# Fiber-seq Snakemake pipeline

A modular workflow for PacBio Fiber-seq BAMs, including instrument-called
m6A/5mC input. Run commands from this repository's root.

## Analysis route

```text
shared reference indexing + pbmm2 HIFI alignment + sorting,
     or samtools sorting for aligned input
  -> keep chr1–chr22 + chrX + chrY (reads and BAM reference header)
  -> optional step-02 DeepVariant -> WhatsHap phase -> aligned-BAM haplotagging
  -> m6A prediction / nucleosome + MSP annotation / reuse existing annotations
  -> basic ft qc (optional autocorrelation)
  -> single-molecule FIRE -> indexed BAM
                         -> shuffled-fiber FDR calibration -> population FIRE peaks
                         -> m6A, CpG 5mC, nucleosome, MSP and FIRE BEDs
                         -> optional FiberHMM footprint BAM + index

alignment -> samtools flagstat + stats
phased VCF + indexed aligned BAM -> haplotagged aligned BAM + phasing statistics
indexed FIRE BAM -> MethBat 5mC pileup + plotting tracks (step 08)
indexed FIRE BAM -> automatic FIRE/CpG tracks + prediction windows
all four feature tracks + model -> optional FiberFold predictions
indexed FIRE BAM + CTCF motifs -> footprint fractions + shared direction bigWigs
indexed FIRE BAM -> official FIRE v0.2.0 peak-analysis rules
```

Rule filenames are numbered in analysis order. Snakemake schedules actual
execution from input/output dependencies; basic QC can run alongside FIRE.
Result and log directory names match the corresponding rule filenames.

## Rule organization

| Rule file | Function |
|---|---|
| `01_align.smk` | Reference indexing, alignment and sorting |
| `02_phasing.smk` | Optional variant calling, phasing and aligned-BAM haplotagging |
| `03_annotate.smk` | m6A / nucleosome / MSP annotation |
| `04_qc.smk` | Basic QC on aligned and annotated BAMs |
| `05_fire.smk` | Single-molecule FIRE and BAM indexing |
| `06_fire_peaks.smk` | Direct invocation of the official FIRE Snakefile |
| `07_extract.smk` | Export feature BEDs |
| `08_methylation.smk` | MethBat 5mC pileup, summary and indexed bedGraph tracks |
| `09_fiberhmm.smk` | Optional general footprint analysis |
| `10_ctcf.smk` | CTCF-specific feature tracks |
| `11_fiberfold_inputs.smk` | FIRE/CpG tracks and model windows |
| `12_fiberfold.smk` | Contact matrix prediction |

Peak calling is enabled by default and can also be selected with the `fire_peaks`
target. QC reports statistics and does not act as a pass/fail gate for downstream analyses.
Step 08 is independent of FiberFold. It is enabled in the local config and
disabled in the portable example; it requires MethBat, bgzip and tabix.

## Software

Activate an analysis environment with `pbmm2`, `samtools`, `ft`, Python and
`pysam`. Step 06 launches the official FIRE workflow through Pixi. Configure
`fire_peaks.fire_dir`, `fire_peaks.pixi` and, optionally, `fire_peaks.conda_prefix`.
The official workflow manages its own Conda environments, including R for its
reports. The outer workflow does not need `--sdm conda` for step 06. Preparation
uses `samtools` from `tools.bin_dir`, or the active PATH when that setting is null.
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

`analysis_chromosomes` is the exact post-alignment chromosome whitelist, defaulting
to chr1–chr22, chrX and chrY. The complete `.aligned.bam` is retained; downstream
analyses and QC use `.chromosomes.bam`, which excludes unmapped records and
all other references from both reads and the BAM header. Missing selected
references cause an error. Changing this list requires recomputing downstream
results; the existing alignment can be reused. Keep `ctcf.chromosomes` within
this whitelist.

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
- `02_phasing/<alias>/`: optional variant calls, phased VCF, aligned BAM with
  haplotype labels and index, phasing statistics and retained intermediates.
- `03_annotate/`: annotated BAM, for hifi/modified modes.
- `04_qc/`: ft QC table (with autocorrelation by default), flagstat and stats.
- `05_fire/`: single-molecule FIRE BAM and BAI.
- `06_fire_peaks/`: official FIRE work directory; its `results/` contains per-sample
  CRAM, peaks, pileup, FDR, TrackHub and figures.
- `07_extract/`: m6a/cpg/nuc/msp/fire `.bed.gz` files.
- `08_methylation/`: 5mC site BED and index, summary JSON, and 0–100% bedGraph
  plotting track and index. 5hmC outputs are not required or generated.
- `09_fiberhmm/`: optional footprint BAM and BAI; `qc/` contains JSON/TSV/PNG/PDF reports.
- `10_ctcf/`: CTCF motif BED, footprint tables/fractions and chromosome bigWigs.
- `11_fiberfold_inputs/`: FIRE/CpG bigWigs, raw CpG scores and prediction windows.
- `12_fiberfold/`: optional prediction pickle.

All paths are relative to `results_dir`. Commands have stdout/stderr logs
under `logs_dir`. QC gathers statistics; it does not filter reads or create
a separate “QC-passed BAM”. Full runs produce large BAM and BED files.

## Optional branches

### Step 02: haplotype phasing

Set `phasing.enabled: true` to include DeepVariant PACBIO calling, WhatsHap
diploid phasing, phasing statistics and aligned-BAM haplotagging in the main target.
Both the local config and portable example leave it disabled. All three stages
use step-01 chromosome-filtered BAMs. Selected samples then use the new
haplotagged BAM as input to step-03 annotation, or directly to downstream
analysis in annotated mode. Disabled or unselected samples bypass phasing.

```bash
# Preview only step 02 and its upstream dependencies.
snakemake -s workflow/fiberseq.smk phasing -n -p --cores 16

# Variant calling only, preserving the earlier target scope.
snakemake -s workflow/fiberseq.smk variant_calls -n --cores 16

# One output alias; phasing.enabled must first be true.
snakemake -s workflow/fiberseq.smk \
  results/02_phasing/4A/4A.haplotagged.bam.bai -n --cores 16
```

Sample IDs come from the main `samples` mapping. Configure the container,
container bind paths, memory and optional output aliases under `phasing`, and
threads under `threads.deepvariant` and `threads.haplotag`. WhatsHap phase uses
one CPU per sample. The FASTA must have an existing `.fai`;
set `phasing.reference_fasta` explicitly if the alignment reference is `.mmi`.
Existing DeepVariant outputs were moved from the old step-13 directory without
recalculation. Enabling phasing in a later run changes downstream BAM inputs;
existing unphased results do not become phased results by renaming directories.
WhatsHap is not currently installed in the local analysis environment. Install
it there, or pass `--sdm conda` to use `workflow/envs/phasing.yaml` for the new
rules. This environment has not been downloaded or created during integration.
See [step-02 configuration, commands and current scope](docs/phasing.md).

### Step 08: 5mC methylation

`08_methylation.smk` uses indexed step-05 FIRE BAMs and runs MethBat with
`--skip-6mA --skip-5hmC`. It requires only the 5mC site BED, tabix index and
summary JSON, then exports and indexes bedGraph tracks for plotting. The
bedGraph uses `Total` rows and the named `mod_score` column without rescaling
its 0–100 values. `plotting/05_figureA/figureA.ini` points to these tracks.

```bash
snakemake -s workflow/fiberseq.smk methylation -n -p --cores 4
```

The completed results and logs have been relocated to `08_methylation/` with
their original contents and modification times. Historical paths recorded
inside summary JSON and logs remain unchanged as execution provenance.
See [step-08 settings and migration notes](docs/methylation.md).

The workflow now loads all 12 rule modules; optional branches are controlled
by configuration. Numbering, results/log paths and script references were
updated without executing any analysis or dry-run. See
[the renumbering record](docs/workflow_renumbering.md).

### FiberHMM

Set `fiberhmm.enabled: true`. The independent branch uses bundled Hia5
PacBio models and `fiberhmm-call`. It does not replace core ft annotations
or create the CTCF-specific bigWig needed by FiberFold. Automatic FiberHMM QC
is enabled with an explicit output prefix. JSON, TSV, PNG and PDF reports are
tracked under `09_fiberhmm/qc/<sample>.footprints.qc.*` and included in the
default targets when FiberHMM is enabled. Install plotting support in the same
environment as `fiberhmm-call`:

```bash
conda install -n fiberseq -c conda-forge matplotlib
```

The four reports are required outputs: if automatic QC fails or plotting support
is missing, Snakemake reports missing outputs rather than accepting only the BAM.

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

`11_fiberfold_inputs.smk` directly calls bedtools, ft, pb-CpG-tools, tabix and
bigtools; it no longer calls `scripts/fiberfold_inputs.py` (kept for reference).
Set `tools.bigtools` to the converter executable; this machine reuses bigtools
0.5.8 from the installed FIRE environment. Simple AWK conversions retain FIRE
coverage fractions and convert CpG percentages to 0–1. Independent chromosome
jobs write the final bigWigs; normal missing intervals are left sparse and the
model reader replaces them with zero. CpG text-to-float conversion may differ
from the previous bigWig-based conversion by float32 rounding.

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
19-bp motifs. Default chromosomes are chr1–chr22, chrX and chrY. Update the list
if predicting other chromosomes. No de novo motif scanning is performed.

`10_ctcf.smk` directly runs the official scripts from the `misc/` directory
beside `fiberfold.script`: `filterOverlapMotifs.py`,
`encodeMotifDirection.py`, and `encodeFP.py`. The old `ctcf_tracks.py` is kept
for reference but is no longer called. There is no Python launch wrapper or
runtime source patching.

Two explicit conversions are needed: source motifs become validated, sorted
BED6 per chromosome; the current `ft footprint` table becomes the official
13-column format by selecting named columns. Duplicate, overlapping or
unsorted footprint rows fail instead of being silently pooled. Official
motif filtering retains its greedy and touching-interval behavior and runs
separately on each chromosome to avoid cross-chromosome comparisons.

Both installed official encoders have a minimal static `import sys` fix.
`encodeFP.py` also follows the supplied chromosome list instead of skipping
the last three entries or hardcoding chrY/chrM/chrEBV exclusions. Its scoring
functions are unchanged. For entirely empty chromosomes it writes one zero
interval, since header-only bigWigs return empty vectors to the model reader.
The original files, exact patch, and upstream commit are stored in
`backups/ctcf_direct_official_20261026_220102/`.

Outputs:

```text
10_ctcf/reference/CTCF.MA0139.1.bed
10_ctcf/<sample>.footprints.tsv.gz
10_ctcf/<sample>.footprints.legacy.tsv.gz
10_ctcf/footprint/<sample>/chr1.bw ...
10_ctcf/direction/chr1.bw ...
```

Direction is +1 on forward motifs, -1 on reverse motifs and 0 elsewhere;
it is shared across samples. Footprint tracks range from 0 to 1 and are
sample-specific. Footprint background and uncovered sites remain missing in
the official sparse output; the official model reader replaces these with 0.
The fraction already includes depth normalization; do not divide it again.
These motif-based scores are not experimental confirmation of CTCF binding.

Each chromosome rule writes the official bigWig directly to its final path;
there is no subsequent split or rewrite. FASTA links, indexes and converted
inputs are explicit Snakemake outputs under `10_ctcf/`. Source reference
directories are untouched. The official direction encoder uses
chromosome-sized arrays. Logs retain `.out`/`.err` pairs, with chromosome
names added for independent jobs.

The custom `<sample>.scores.tsv.gz` is no longer a required output: the
official scripts do not produce it and FiberFold prediction does not use it.
Existing score files are retained as historical results and are not refreshed.
See `docs/ctcf.md` for commands, compatibility changes and validation.

### Whole-genome FIRE peaks

Step 06 directly invokes `/data/tenghaohan/software/FIRE/workflow/Snakefile`
through the official Pixi `fire` task. No coverage, background, FDR or merging
rules are copied into the parent workflow. `06_fire_peaks.smk` only prepares
inputs and launches the author's `all` target for the configured samples.

```bash
snakemake -s workflow/fiberseq.smk fire_peaks -n --cores 8
snakemake -s workflow/fiberseq.smk fire_peaks --cores 8 \
  --resources mem_mb=16000
```

The active analysis environment is still required for upstream steps; preparation
uses its `samtools`. On this machine `tools.bin_dir` locates that executable
explicitly. `fire_peaks.pixi` locates Pixi, and `fire_peaks.fire_dir` locates the
FIRE checkout. `fire_peaks.conda_prefix` optionally reuses the installed official
Conda environments. The local configuration reuses
`/data/tenghaohan/software/FIRE/fire-test-data/.snakemake/conda`.

The two calculation-facing rules are:

1. `prepare_fire_peaks`: wait for the indexed step-05 BAMs, create an absolute-path
   sample table and official configuration, and index a local reference symlink.
2. `run_official_fire`: run `pixi run --manifest-path ... fire all` with that
   configuration and an isolated working directory. Official Snakemake schedules
   all internal jobs. The outer job's threads and memory budget are forwarded.

The `fire_peaks` target selects this branch. The parent `all` target includes
its `.complete` marker when `fire_peaks.enabled` is true. The marker is written
only after official `all` succeeds. Peak BEDs, indexes and FDR tables are also
explicit outputs of the wrapper. One official invocation handles all configured
samples, avoiding competing Snakemake locks in the shared working directory.

The official workflow applies its own single-molecule FIRE model and filters
to step-05 BAMs, producing the author's filtered CRAM files. Step 05 remains
unchanged. The complete official output set includes standard/wide/1%-FDR peaks,
pileup, QC, TrackHub files and R figures. Peak rows retain the official
`pass_coverage` flag. The author's existing local `awk` coordinate-rounding fix
is used because this is the same FIRE checkout.

Official work directory and representative outputs:

```text
06_fire_peaks/
  config.json
  manifest.tsv
  reference/genome.fa + .fai
  .complete
  .snakemake/
  temp/
  results/<sample>/
    <sample>-fire-v0.2-filtered.cram + .crai
    <sample>-fire-v0.2-peaks.bed.gz + .tbi
    <sample>-fire-v0.2-pileup.bed.gz + .tbi
    additional-outputs-v0.2/
      fire-peaks/<sample>-v0.2-fire-score-to-fdr.tbl
      fire-peaks/<sample>-v0.2-fire-elements.bed.gz + .tbi
      coverage/...
      analysis/...
    trackHub-v0.2/...
```

The inner `results/` directory is the official layout, now preserved as-is. The
old copied-rule outputs outside this directory are not moved or deleted.
`fire_peaks.version` must match the installed FIRE version: `v0.2` uses its normal
major/minor label, while `v0.2.0` enables the author's `full-version` setting.

Configuration translation retains `max_fdr` -> `max_peak_fdr`,
`min_fire_fraction` -> `min_frac_accessible`, `sd_coverage` ->
`coverage_within_n_sd`, and optional element/fiber filters. `analysis_chromosomes`
becomes the official chromosome regular expression. hg38 exclusion masks are
added by FIRE itself, along with any `fire_peaks.exclude` files. The checkout
fixes shuffle seed at 42 and the local-maxima window at 200 bp; unsupported
values are rejected during preparation instead of silently ignored.

The outer dry-run previews preparation and the launch job, not FIRE's internal
DAG. After preparing inputs, inspect the inner plan with the same command shown
in `docs/fire_peaks.md` plus `-n`. Inner progress and errors are captured in
`logs_dir/06_fire_peaks/official.log`; reference preparation logs are in
`prepare.log`. The author manages intermediate files and restart state inside
its own work directory.

## Validation

See `docs/validation.md` for execution evidence and remaining limitations.

## Independent plotting

Plotting scripts and configurations live in [`plotting/`](plotting/README.md),
with numbered topic directories: `01_alignment`, `02_qc`, `03_fire_compare`,
`04_nucleosome`, `05_figureA` and `06_tss`. Plots are run directly, outside
Snakemake. Alignment comparisons have three self-contained numbered scripts
for read length, depth and median summaries; edit the settings at the top of
each script and run it directly against existing samtools statistics.

The former `results/11_visualization/` and `results/12_fire_compare/` directories
have moved intact to `analysis/11_visualization/` and `analysis/12_fire_compare/`.
The latter retains its comparison tables, BEDs, rankings, chunk state and plots.
Other numbered analysis results remain under `results/`. Generated files under
`analysis/` are ignored by Git.

```bash
# Inspect a self-contained plotting entry point.
/data/tenghaohan/miniconda3/envs/fiberseq/bin/python \
  plotting/04_nucleosome/plot_nucleosome.py --help

# Figure A uses the configuration beside its plotting command.
bash plotting/05_figureA/plot.sh
```

The legacy single-sample QC and FIRE comparison scripts still depend
on a previously removed `_common.py`; see the plotting README before using them.
The standalone numbered alignment and QC scripts do not need this shared code.
The remaining plotting scripts and configurations were relocated without
changing their scientific parameters. See the numbered directory guide for
inputs, outputs and current dependencies.
