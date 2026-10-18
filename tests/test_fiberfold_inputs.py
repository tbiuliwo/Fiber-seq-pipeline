import gzip
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

import pyBigWig

spec = importlib.util.spec_from_file_location("fold_inputs", Path(__file__).parents[1] / "scripts/fiberfold_inputs.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class InputsTest(unittest.TestCase):
    def test_fractions_and_percentage_conversion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fai = root / "ref.fai"
            fai.write_text("chr1\t20\t0\t0\t0\nchr2\t20\t0\t0\t0\n")
            table = root / "pileup.tsv"
            table.write_text("#chrom\tstart\tend\tcoverage\tfire_coverage\tscore\n"
                             "chr1\t2\t5\t4\t2\t99\nchr1\t5\t7\t4\t2\t99\n"
                             "chr1\t10\t12\t0\t0\t-1\n")
            module.fire(SimpleNamespace(fai=fai, chromosomes=["chr1", "chr2"], table=table, out_dir=root / "fire"))
            with pyBigWig.open(str(root / "fire/chr1.bw")) as bw:
                self.assertEqual(bw.values("chr1", 0, 8), [0, 0, .5, .5, .5, .5, .5, 0])
                self.assertEqual(bw.values("chr1", 10, 12), [0, 0])
            source = root / "combined.bw"
            with pyBigWig.open(str(source), "w") as bw:
                bw.addHeader([("chr1", 20)])
                bw.addEntries(["chr1", "chr1", "chr1"], [2, 4, 6], ends=[3, 5, 7], values=[0., 50., 100.])
            module.cpg(SimpleNamespace(fai=fai, chromosomes=["chr1", "chr2"], bigwig=str(source), out_dir=root / "cpg"))
            with pyBigWig.open(str(root / "cpg/chr1.bw")) as bw:
                self.assertEqual(bw.values("chr1", 4, 7), [.5, 0, 1])
            with pyBigWig.open(str(root / "cpg/chr2.bw")) as bw:
                self.assertEqual(bw.values("chr2", 0, 20), [0] * 20)

    def test_windows_keep_fixed_length_and_skip_tail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fai = root / "ref.fai"
            n = module.WINDOW_SIZE
            fai.write_text(f"chr1\t{2*n+100}\t0\t0\t0\nchr2\t10\t0\t0\t0\n")
            module.windows(SimpleNamespace(fai=fai, chromosomes=["chr1", "chr2"], step=n,
                                           out=root / "windows.gz", regions=root / "regions.bed"))
            with gzip.open(root / "windows.gz", "rt") as handle:
                rows = [line.split() for line in handle]
            self.assertEqual(len(rows), 2)
            self.assertEqual([int(row[2]) - int(row[1]) for row in rows], [n, n])
            self.assertEqual((root / "regions.bed").read_text().splitlines()[-1], "chr2\t0\t10")


if __name__ == "__main__":
    unittest.main()
