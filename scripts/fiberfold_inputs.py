"""Generate FiberFold FIRE/CpG bigWigs and full-length prediction windows."""
import argparse
import csv
import gzip
import math
from pathlib import Path
import sys

import pyBigWig

WINDOW_SIZE = 2097152


def read_sizes(path, chromosomes):
    with open(path) as handle:
        lengths = {f[0]: int(f[1]) for line in handle if (f := line.split())}
    return {chrom: lengths[chrom] for chrom in chromosomes}


class TrackWriter:
    """Write sorted intervals, merge adjacent equal values and fill gaps with zero."""
    def __init__(self, path, chrom, length):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.handle = pyBigWig.open(str(path), "w")
        self.handle.addHeader([(chrom, length)])
        self.chrom, self.length, self.cursor = chrom, length, 0
        self.pending = None
        self.starts, self.ends, self.values = [], [], []

    def flush(self):
        if self.starts:
            self.handle.addEntries([self.chrom] * len(self.starts), self.starts,
                                   ends=self.ends, values=self.values)
            self.starts, self.ends, self.values = [], [], []

    def append(self, start, end, value):
        if self.pending and self.pending[1] == start and self.pending[2] == value:
            self.pending = (self.pending[0], end, value)
        else:
            if self.pending:
                left, right, score = self.pending
                self.starts.append(left); self.ends.append(right); self.values.append(float(score))
                if len(self.starts) >= 10000:
                    self.flush()
            self.pending = (start, end, value)

    def add(self, start, end, value):
        if not (self.cursor <= start < end <= self.length) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"Invalid/out-of-order interval or fraction: {self.chrom}:{start}-{end}={value}")
        if start > self.cursor:
            self.append(self.cursor, start, 0.0)
        self.append(start, end, value)
        self.cursor = end

    def close(self):
        if self.cursor < self.length:
            self.append(self.cursor, self.length, 0.0)
        if self.pending:
            left, right, value = self.pending
            self.starts.append(left); self.ends.append(right); self.values.append(float(value))
        self.flush()
        self.handle.close()


def windows(args):
    lengths = read_sizes(args.fai, args.chromosomes)
    count = 0
    with gzip.open(args.out, "wt") as out, open(args.regions, "w") as regions:
        for chrom, length in lengths.items():
            regions.write(f"{chrom}\t0\t{length}\n")
            for start in range(0, length - WINDOW_SIZE + 1, args.step):
                out.write(f"{chrom}\t{start}\t{start + WINDOW_SIZE}\n")
                count += 1
    if not count:
        raise ValueError("No selected chromosome is long enough for a 2097152-bp prediction window")
    print(f"Wrote {count} prediction windows")


def fire(args):
    lengths = read_sizes(args.fai, args.chromosomes)
    writers = {chrom: TrackWriter(Path(args.out_dir) / f"{chrom}.bw", chrom, length)
               for chrom, length in lengths.items()}
    handle = sys.stdin if args.table == "-" else open(args.table)
    try:
        header = handle.readline().lstrip("#").rstrip().split("\t")
        for row in csv.DictReader(handle, fieldnames=header, delimiter="\t"):
            chrom = row["chrom"]
            if chrom not in writers:
                continue
            depth, positives = int(row["coverage"]), int(row["fire_coverage"])
            fraction = positives / depth if depth else 0.0
            writers[chrom].add(int(row["start"]), int(row["end"]), fraction)
    finally:
        if handle is not sys.stdin:
            handle.close()
        for writer in writers.values():
            writer.close()


def cpg(args):
    lengths = read_sizes(args.fai, args.chromosomes)
    with pyBigWig.open(args.bigwig) as source:
        available = source.chroms()
        for chrom, length in lengths.items():
            out = TrackWriter(Path(args.out_dir) / f"{chrom}.bw", chrom, length)
            try:
                if chrom not in available:
                    continue
                # pb-CpG-tools bigWig scores are percentages (0..100).
                for left in range(0, min(length, available[chrom]), 1000000):
                    right = min(left + 1000000, length, available[chrom])
                    for start, end, percent in source.intervals(chrom, left, right) or ():
                        out.add(max(start, left), min(end, right), percent / 100.0)
            finally:
                out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, function in (("windows", windows), ("fire", fire), ("cpg", cpg)):
        cmd = commands.add_parser(name)
        cmd.set_defaults(function=function)
        cmd.add_argument("--fai", required=True)
        cmd.add_argument("--chromosomes", nargs="+", required=True)
        if name == "windows":
            cmd.add_argument("--step", type=int, default=WINDOW_SIZE)
            cmd.add_argument("--out", required=True)
            cmd.add_argument("--regions", required=True)
        else:
            cmd.add_argument("--out-dir", required=True)
        if name == "fire":
            cmd.add_argument("--table", default="-")
        if name == "cpg":
            cmd.add_argument("--bigwig", required=True)
    args = parser.parse_args()
    args.function(args)
