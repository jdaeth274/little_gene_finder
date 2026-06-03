import importlib.util
import sys
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "cps_gene_pipeline.py"
spec = importlib.util.spec_from_file_location("cps_gene_pipeline", MODULE_PATH)
cps = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = cps
spec.loader.exec_module(cps)


def make_hit(**overrides):
    values = dict(
        isolate="ERR001",
        gene="wzg",
        source="blastn",
        contig="contig1",
        qstart=1,
        qend=100,
        sstart=10,
        send=109,
        pident=99.0,
        length=100,
        evalue=1e-50,
        bitscore=200.0,
        qlen=100,
        slen=1000,
    )
    values.update(overrides)
    return cps.Hit(**values)


class CpsGenePipelineTests(unittest.TestCase):
    def test_complete_single_hit_classification(self):
        call = cps.classify_hits(
            [make_hit()],
            complete_qcov=0.90,
            partial_qcov=0.20,
            fragmented_gap=30,
            max_query_overlap=10,
        )
        self.assertEqual(call.status, "complete")
        self.assertEqual(call.qcov, 1.0)
        self.assertEqual(call.n_hsps, 1)

    def test_fragmented_high_coverage_split_across_contigs(self):
        call = cps.classify_hits(
            [
                make_hit(contig="contig1", qstart=1, qend=50, sstart=10, send=59, length=50, bitscore=100.0),
                make_hit(contig="contig2", qstart=51, qend=100, sstart=5, send=54, length=50, bitscore=100.0),
            ],
            complete_qcov=0.90,
            partial_qcov=0.20,
            fragmented_gap=30,
            max_query_overlap=10,
        )
        self.assertEqual(call.status, "fragmented")
        self.assertEqual(call.qcov, 1.0)
        self.assertEqual(call.n_contigs, 2)

    def test_reverse_strand_extraction_is_oriented_to_query(self):
        hit = make_hit(sstart=8, send=3, qstart=1, qend=6)
        call = cps.classify_hits([hit], 0.9, 0.2, 30, 10)
        seq = cps.extract_call_sequence(call, {"contig1": "AACCCGGGTT"}, gap_n=5)
        self.assertEqual(seq, cps.revcomp("CCCGGG"))


if __name__ == "__main__":
    unittest.main()
