"""Build a chromosome-preserving shuffled-fiber background for ft call-peaks."""
import argparse
from collections import Counter
import gzip
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

import pysam


def run(command, output=None):
    if output is None:
        subprocess.run(command, check=True)
    else:
        with open(output, "w") as handle:
            subprocess.run(command, stdout=handle, check=True)


def coverage_limits(bedgraph, minimum=4, sd=5):
    """Use the bp-weighted median of nonzero eligible-fiber coverage."""
    histogram = Counter()
    with open(bedgraph) as handle:
        for line in handle:
            _, start, end, depth = line.split()
            depth = int(depth)
            if depth:
                histogram[depth] += int(end) - int(start)
    total = sum(histogram.values())
    if not total:
        raise ValueError("No eligible fibers remain for FIRE peak calling.")
    cumulative = 0
    for median, bases in sorted(histogram.items()):
        cumulative += bases
        if cumulative * 2 >= total:
            break
    low = round(max(minimum, median - sd * math.sqrt(median)))
    high = round(median + sd * math.sqrt(median))
    return median, low, high


def open_bed(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def count_lines(path):
    with open(path) as handle:
        return sum(1 for _ in handle)


def background(args):
    root = Path(args.shuffled).parent
    root.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(args.bam) as bam:
        sizes = dict(zip(bam.references, bam.lengths))
    with open(args.sizes, "w") as handle:
        for chrom in sorted(sizes):
            handle.write(f"{chrom}\t{sizes[chrom]}\n")

    with tempfile.TemporaryDirectory(prefix=".peak-background-", dir=root) as directory:
        tmp = Path(directory)
        fibers = tmp / "fibers.bed"
        # Keep precisely the molecules used by call-peaks. Column 5 records the
        # original alignment start; ft uses it with chromosome/read name as key.
        command = ["ft", "-t", str(args.threads), "extract", "--all", "-", "--simplify",
                   "--fire-filter", "-F", "260", args.bam]
        with open(tmp / "unsorted.bed", "w") as out:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
            try:
                for line in process.stdout:
                    if line.startswith("#"):
                        continue
                    chrom, start, end, name = line.split("\t", 4)[:4]
                    out.write(f"{chrom}\t{start}\t{end}\t{name}\t{start}\n")
            finally:
                process.stdout.close()
            if process.wait():
                raise RuntimeError("ft extract failed while preparing FIRE peak background.")
        sort_command = ["sort", "-S", "1G", "-T", directory,
                        "-k1,1", "-k2,2n", "-k3,3n", "-k4,4"]
        os.environ["LC_ALL"] = "C"
        run(sort_command + [str(tmp / "unsorted.bed")], fibers)
        count = count_lines(fibers)
        coverage = tmp / "coverage.bedgraph"
        run(["bedtools", "genomecov", "-bga", "-i", str(fibers), "-g", args.sizes], coverage)
        median, low, high = coverage_limits(coverage, args.min_coverage, args.sd_coverage)
        bad = tmp / "unreliable.bed"
        with open(coverage) as handle, open(bad, "w") as out:
            for line in handle:
                chrom, start, end, depth = line.split()
                if int(depth) <= low or int(depth) >= high:
                    out.write(f"{chrom}\t{start}\t{end}\n")
        filtered = tmp / "filtered.bed"
        run(["bedtools", "intersect", "-sorted", "-v", "-f", "0.2",
             "-a", str(fibers), "-b", str(bad)], filtered)
        retained = count_lines(filtered)
        if not retained:
            raise ValueError(f"No fibers remain for shuffling (median coverage {median}, "
                             f"limits {low}–{high}). Inspect sample coverage.")

        # Do not move fibers into unobserved sequence, assembly gaps or configured
        # blacklists. Read masks as streams rather than genome-sized arrays.
        run(["bedtools", "genomecov", "-bga", "-i", str(filtered), "-g", args.sizes], coverage)
        masks = tmp / "masks.bed"
        with open(masks, "w") as out:
            with open(coverage) as handle:
                for line in handle:
                    chrom, start, end, depth = line.split()
                    if int(depth) == 0:
                        out.write(f"{chrom}\t{start}\t{end}\n")
            for path in args.exclude:
                with open_bed(path) as handle:
                    for line in handle:
                        if not line.strip() or line.startswith(("#", "track", "browser")):
                            continue
                        chrom, start, end = line.split()[:3]
                        if chrom in sizes:
                            start, end = max(0, int(start)), min(sizes[chrom], int(end))
                            if end > start:
                                out.write(f"{chrom}\t{start}\t{end}\n")
        sorted_masks = tmp / "sorted-masks.bed"
        run(sort_command + [str(masks)], sorted_masks)
        excluded = tmp / "exclude.bed"
        run(["bedtools", "merge", "-i", str(sorted_masks)], excluded)
        shuffled = tmp / "shuffled.bed"
        run(["bedtools", "shuffle", "-chrom", "-seed", str(args.seed),
             "-maxTries", "10000", "-excl", str(excluded), "-g", args.sizes,
             "-i", str(filtered)], shuffled)
        shuffled_count = count_lines(shuffled)
        if shuffled_count != retained:
            raise ValueError(f"Only {shuffled_count}/{retained} fibers could be shuffled. "
                             "Inspect coverage and exclusion masks; an incomplete null is invalid.")
        sorted_shuffle = tmp / "sorted-shuffled.bed"
        run(sort_command + [str(shuffled)], sorted_shuffle)
        pysam.tabix_compress(str(sorted_shuffle), args.shuffled, force=True)
        pysam.tabix_compress(str(excluded), args.excluded, force=True)

    summary = {"median_coverage": median, "min_coverage": low, "max_coverage": high,
               "eligible_fibers": count, "shuffled_fibers": retained,
               "seed": args.seed, "exclude": args.exclude,
               "coverage_source": "FIRE-filtered aligned fiber spans"}
    Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n")
    Path(args.limits).write_text(f"{low}\t{high}\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--shuffled", required=True)
    parser.add_argument("--excluded", required=True)
    parser.add_argument("--sizes", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--limits", required=True)
    parser.add_argument("--exclude", nargs="*", default=[])
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-coverage", type=int, default=4)
    parser.add_argument("--sd-coverage", type=float, default=5)
    background(parser.parse_args())
