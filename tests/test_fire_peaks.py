import importlib.util
from io import StringIO
from pathlib import Path
import tempfile
import unittest


def load_script(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parents[1] / "scripts" / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


background = load_script("fire_peak_background")
fdr = load_script("fire_peak_fdr")

HEADER = "#chrom\tstart\tend\tcoverage\tscore\tcoverage_shuffled\tscore_shuffled\n"


class FirePeakFdrTest(unittest.TestCase):
    def test_intervals_contribute_base_pairs_rather_than_rows(self):
        pileup = StringIO(HEADER +
                         "chr1\t0\t100\t8\t90\t8\t20\n"
                         "chr1\t100\t101\t8\t20\t8\t90\n")
        real, null = fdr.score_histograms(pileup, 4, 20)
        self.assertEqual(real, {9000: 100, 2000: 1})
        self.assertEqual(null, {2000: 100, 9000: 1})
        rows = {threshold: (rate, null_bp, real_bp)
                for threshold, rate, null_bp, real_bp in fdr.fdr_rows(real, null)}
        self.assertEqual(rows[90], (0.01, 1, 100))
        self.assertEqual(rows[20], (1.0, 101, 101))

    def test_real_and_null_use_the_same_coverage_domain(self):
        pileup = StringIO(HEADER +
                         "chr1\t0\t10\t4\t40\t12\t30\n"
                         "chr1\t10\t20\t8\t80\t3\t70\n"
                         "chr1\t20\t30\t13\t90\t8\t60\n"
                         "chr1\t30\t40\t8\t85\t13\t65\n"
                         "chr1\t40\t50\t3\t95\t8\t75\n")
        real, null = fdr.score_histograms(pileup, 4, 12)
        # Both boundary depths are included; either track outside the bounds
        # removes the same genomic interval from both distributions.
        self.assertEqual(real, {4000: 10})
        self.assertEqual(null, {3000: 10})

    def test_equal_scores_are_combined_before_computing_fdr(self):
        pileup = StringIO(HEADER +
                         "chr1\t0\t3\t8\t50\t8\t50\n"
                         "chr1\t3\t10\t8\t50\t8\t50\n"
                         "chr1\t10\t30\t8\t80\t8\t-1\n")
        rows = fdr.fdr_rows(*fdr.score_histograms(pileup, 4, 20))
        self.assertEqual([row[0] for row in rows], [-1.0, 50.0, 80.0])
        self.assertEqual(rows[1], (50.0, 1 / 3, 10, 30))
        self.assertEqual(rows[2], (80.0, 0.0, 0, 20))

    def test_empty_or_unusable_null_fails_instead_of_creating_zero_fdr(self):
        examples = (
            HEADER,
            HEADER + "chr1\t0\t20\t8\t80\t8\t-1\n",
            HEADER + "chr1\t0\t20\t8\t80\t8\tnan\n",
            HEADER + "chr1\t0\t20\t8\t80\t0\t70\n",
        )
        for pileup in examples:
            with self.subTest(pileup=pileup):
                with self.assertRaisesRegex(ValueError, "No usable real/shuffled"):
                    fdr.score_histograms(StringIO(pileup), 4, 20)


class FirePeakCoverageTest(unittest.TestCase):
    def test_median_weights_covered_bases_and_ignores_uncovered_sequence(self):
        with tempfile.TemporaryDirectory() as directory:
            bedgraph = Path(directory) / "coverage.bedgraph"
            bedgraph.write_text("chr1\t0\t10000\t0\n"
                                "chr1\t10000\t10001\t1\n"
                                "chr1\t10001\t10002\t2\n"
                                "chr1\t10002\t10102\t20\n")
            self.assertEqual(background.coverage_limits(bedgraph), (20, 4, 42))
            self.assertEqual(background.coverage_limits(bedgraph, minimum=10, sd=2),
                             (20, 11, 29))

    def test_no_covered_bases_cannot_define_background(self):
        with tempfile.TemporaryDirectory() as directory:
            bedgraph = Path(directory) / "coverage.bedgraph"
            bedgraph.write_text("chr1\t0\t1000\t0\n")
            with self.assertRaisesRegex(ValueError, "No eligible fibers"):
                background.coverage_limits(bedgraph)


if __name__ == "__main__":
    unittest.main()
