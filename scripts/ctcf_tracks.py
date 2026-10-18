"""Prepare MA0139.1 motifs and FiberFold-compatible CTCF bigWigs."""
import argparse
from collections import defaultdict
import csv
import gzip
from pathlib import Path

import pyBigWig
import pysam

MODULES = ((0, 6), (6, 13), (13, 19))


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def sizes(path):
    with open(path) as handle:
        return {f[0]: int(f[1]) for line in handle if (f := line.split())}


def bed_rows(path):
    with open(path) as handle:
        return [(f[0], int(f[1]), int(f[2]), f[3], float(f[4]), f[5])
                for line in handle if (f := line.split())]


def prepare(args):
    lengths = sizes(args.fai)
    motifs = defaultdict(list)
    with text_open(args.source) as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.split()
            # JASPAR track: chrom/start/end/name/score/pval/strand; BED6 also supported.
            if f[0] not in args.chromosomes or f[0] not in lengths:
                continue
            chrom, start, end, name, score, strand = f[0], int(f[1]), int(f[2]), f[3], float(f[4]), f[-1]
            if end - start != 19 or start < 0 or end > lengths[chrom] or strand not in ("+", "-"):
                raise ValueError(f"Expected 19-bp MA0139.1 BED coordinates: {line.strip()}")
            if score >= args.min_score:
                motifs[chrom].append((chrom, start, end, name, score, strand))
    retained = []
    for chrom in args.chromosomes:
        selected = []
        for row in sorted(motifs[chrom], key=lambda r: (r[1], r[2], -r[4], r[5])):
            if selected and row[1] < selected[-1][2]:
                if row[4] > selected[-1][4]:
                    selected[-1] = row
            else:
                selected.append(row)
        retained.extend(selected)
    with open(args.out, "w") as out:
        for chrom, start, end, name, score, strand in retained:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score:g}\t{strand}\n")
    print(f"Retained {len(retained)} non-overlapping motifs")


def write_tracks(directory, lengths, chromosomes, rows):
    """Sparse features plus explicit zero background, without chromosome-sized arrays."""
    Path(directory).mkdir(parents=True, exist_ok=True)
    for chrom in chromosomes:
        length = lengths[chrom]
        starts, ends, values = [], [], []
        cursor = 0
        for start, end, value in sorted(rows.get(chrom, [])):
            if start > cursor:
                starts.append(cursor); ends.append(start); values.append(0.0)
            starts.append(start); ends.append(end); values.append(float(value))
            cursor = end
        if cursor < length:
            starts.append(cursor); ends.append(length); values.append(0.0)
        with pyBigWig.open(str(Path(directory) / f"{chrom}.bw"), "w") as out:
            out.addHeader([(chrom, length)])
            out.addEntries([chrom] * len(starts), starts, ends=ends, values=values)


def direction(args):
    rows = defaultdict(list)
    for chrom, start, end, name, score, strand in bed_rows(args.bed):
        rows[chrom].append((start, end, 1 if strand == "+" else -1))
    write_tracks(args.out_dir, sizes(args.fai), args.chromosomes, rows)


def footprints(args):
    counts = defaultdict(lambda: [0, 0])
    with pysam.FastaFile(args.reference, filepath_index=args.fai) as reference, text_open(args.table) as handle:
        header = handle.readline().lstrip("#").rstrip().split("\t")
        for row in csv.DictReader(handle, fieldnames=header, delimiter="\t"):
            chrom, start, end = row["chrom"], int(row["start"]), int(row["end"])
            if chrom not in args.chromosomes:
                continue
            depth = int(row["n_spanning_fibers"])
            if not depth:
                continue
            seq = reference.fetch(chrom, start, end).upper()
            # Match the author's encodeFP.py: A/T eligibility in reference order,
            # on both strands. Keep this convention for the released model.
            required = 1
            for i, (left, right) in enumerate(MODULES, 1):
                if any(base in seq[left:right] for base in "AT"):
                    required |= 1 << i
            codes = [int(code) for code in row["footprint_codes"].split(",") if code not in ("", ".")]
            positive = sum(code >= 0 and (code & required) == required for code in codes)
            if positive > depth:
                raise ValueError("Footprint count exceeds spanning fibers")
            counts[(chrom, start, end, row["strand"])][0] += positive
            counts[(chrom, start, end, row["strand"])][1] += depth
    rows = defaultdict(list)
    with gzip.open(args.scores, "wt") as scores:
        scores.write("chrom\tstart\tend\tstrand\tfootprinted_fibers\tspanning_fibers\tfraction\n")
        for chrom, start, end, name, score, strand in bed_rows(args.bed):
            positive, depth = counts[(chrom, start, end, strand)]
            value = positive / depth if depth else 0.0
            rows[chrom].append((start, end, value))
            scores.write(f"{chrom}\t{start}\t{end}\t{strand}\t{positive}\t{depth}\t{value}\n")
    write_tracks(args.out_dir, sizes(args.fai), args.chromosomes, rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, function in (("prepare", prepare), ("direction", direction), ("footprints", footprints)):
        cmd = sub.add_parser(name)
        cmd.set_defaults(function=function)
        cmd.add_argument("--fai", required=True)
        cmd.add_argument("--chromosomes", nargs="+", required=True)
        if name == "prepare":
            cmd.add_argument("--source", required=True)
            cmd.add_argument("--out", required=True)
            cmd.add_argument("--min-score", type=float, default=0)
        else:
            cmd.add_argument("--bed", required=True)
            cmd.add_argument("--out-dir", required=True)
        if name == "footprints":
            cmd.add_argument("--table", required=True)
            cmd.add_argument("--reference", required=True)
            cmd.add_argument("--scores", required=True)
    args = parser.parse_args()
    args.function(args)
