import gzip
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

import pyBigWig
import pysam

spec = importlib.util.spec_from_file_location("ctcf", Path(__file__).parents[1] / "scripts/ctcf_tracks.py")
ctcf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctcf)


class TracksTest(unittest.TestCase):
    def test_motif_selection_and_tracks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 100 + "\n>chr3\n" + "A" * 100 + "\n")
            pysam.faidx(str(reference))
            fai = str(reference) + ".fai"
            source = root / "jaspar.tsv.gz"
            with gzip.open(source, "wt") as out:
                out.write("chr1\t10\t29\tCTCF\t900\t400\t-\n")
                out.write("chr1\t12\t31\tCTCF\t800\t300\t+\n")
                out.write("chr1\t29\t48\tCTCF\t700\t300\t+\n")
                out.write("chr2\t10\t29\tCTCF\t900\t400\t+\n")
            bed = root / "motifs.bed"
            common = dict(fai=fai, chromosomes=["chr1", "chr2", "chr3"])
            ctcf.prepare(SimpleNamespace(**common, source=source, out=bed, min_score=0))
            self.assertEqual(len(ctcf.bed_rows(bed)), 3)
            ctcf.direction(SimpleNamespace(**common, bed=bed, out_dir=root / "direction"))
            with pyBigWig.open(str(root / "direction/chr1.bw")) as bw:
                self.assertEqual(bw.values("chr1", 9, 11), [0.0, -1.0])
                self.assertEqual(bw.values("chr1", 29, 30), [1.0])
            table = root / "footprints.tsv"
            table.write_text(
                "#chrom\tstart\tend\tstrand\tn_spanning_fibers\tfootprint_codes\n"
                "chr1\t10\t29\t-\t2\t15,1\n"
                "chr1\t29\t48\t+\t0\t.\n"
                "chr2\t10\t29\t+\t1\t1\n")
            ctcf.footprints(SimpleNamespace(**common, reference=str(reference), table=table,
                                           bed=bed, scores=root / "scores.tsv.gz", out_dir=root / "footprint"))
            with pyBigWig.open(str(root / "footprint/chr1.bw")) as bw:
                self.assertEqual(bw.values("chr1", 10, 11), [0.5])
                self.assertEqual(bw.values("chr1", 29, 30), [0.0])
            with pyBigWig.open(str(root / "footprint/chr2.bw")) as bw:
                self.assertEqual(bw.values("chr2", 10, 11), [1.0])
            with pyBigWig.open(str(root / "footprint/chr3.bw")) as bw:
                self.assertEqual(bw.values("chr3", 0, 100), [0.0] * 100)


if __name__ == "__main__":
    unittest.main()
