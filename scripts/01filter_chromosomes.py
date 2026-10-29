"""Keep selected mapped chromosomes and rebuild BAM reference IDs/header."""
import argparse
import json

import pysam


def filter_chromosomes(bam, output, chromosomes, threads=1):
    with pysam.AlignmentFile(bam, "rb", threads=threads) as source:
        missing = set(chromosomes) - set(source.references)
        if missing:
            raise ValueError(f"Chromosomes absent from BAM header: {sorted(missing)}")
        selected = set(chromosomes)
        header = source.header.to_dict()
        # Preserve the input coordinate order, even when the config order differs.
        header["SQ"] = [sq for sq in header["SQ"] if sq["SN"] in selected]
        new_ids = {sq["SN"]: i for i, sq in enumerate(header["SQ"])}
        mapping = {i: new_ids[name] for i, name in enumerate(source.references)
                   if name in selected}
        total = kept = 0
        with pysam.AlignmentFile(output, "wb", header=header, threads=threads) as dest:
            for read in source.fetch(until_eof=True):
                total += 1
                if read.is_unmapped or read.reference_id not in mapping:
                    continue
                read.reference_id = mapping[read.reference_id]

                dest.write(read)
                kept += 1
    return {"total_records": total, "kept_records": kept,
            "discarded_records": total - kept, "chromosomes": chromosomes}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    print(json.dumps(filter_chromosomes(args.bam, args.output, args.chromosomes,
                                       args.threads), indent=2))
