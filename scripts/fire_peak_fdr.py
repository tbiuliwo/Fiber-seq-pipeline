"""Estimate empirical FIRE-score FDR from a streamed real/shuffled ft pileup."""
import argparse
from collections import Counter
import math
from pathlib import Path
import sys


def score_histograms(lines, low, high):
    real, shuffled = Counter(), Counter()
    columns = None
    covered_bp = 0
    for line in lines:
        if line.startswith("#chrom"):
            columns = {name: index for index, name in enumerate(line.strip().split("\t"))}
            continue
        if line.startswith("#") or not line.strip():
            continue
        if columns is None:
            raise ValueError("Missing ft pileup header.")
        row = line.rstrip().split("\t")
        if not (low <= int(row[columns["coverage"]]) <= high
                and low <= int(row[columns["coverage_shuffled"]]) <= high):
            continue
        bp = int(row[columns["end"]]) - int(row[columns["start"]])
        covered_bp += bp
        for name, histogram in (("score", real), ("score_shuffled", shuffled)):
            score = float(row[columns[name]])
            if score >= 0 and math.isfinite(score):
                # 0.01 score bins keep memory bounded. Thresholds are the lower
                # bin edges, so all bp counted at t have a binned score >= t.
                histogram[math.floor(score * 100)] += bp
    if not covered_bp or not real or not shuffled:
        raise ValueError("No usable real/shuffled FIRE scores within coverage limits. "
                         "Inspect coverage and FIRE annotations.")
    return real, shuffled


def fdr_rows(real, shuffled):
    real_bp = shuffled_bp = 0
    rows = [(-1.0, 1.0, 1, 1)]
    for score in sorted(real.keys() | shuffled.keys(), reverse=True):
        real_bp += real[score]
        shuffled_bp += shuffled[score]
        if real_bp:
            rows.append((score / 100, min(1.0, shuffled_bp / real_bp), shuffled_bp, real_bp))
    return sorted(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limits", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    low, high = map(int, Path(args.limits).read_text().split())
    real, shuffled = score_histograms(sys.stdin, low, high)
    rows = fdr_rows(real, shuffled)
    with open(args.out, "w") as out:
        out.write("threshold\tFDR\tshuffled_bp\treal_bp\n")
        for threshold, fdr, null_bp, real_bp in rows:
            out.write(f"{threshold:.2f}\t{fdr:.12g}\t{null_bp}\t{real_bp}\n")
    print(f"FDR table: {len(rows)} thresholds; real bp={sum(real.values())}; "
          f"shuffled bp={sum(shuffled.values())}")


if __name__ == "__main__":
    main()
