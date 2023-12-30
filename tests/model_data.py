import unittest
from . import test_config

bp_vector_schema: list[str] = [
    "A",
    "C",
    "G",
    "T",
    "PROMOTOR_MOTIF",
    "ORF",
    "exon",
    "mRNA",
    "miRNA",
    "rRNA",
    "CDS",
    "POLY_ADENYL",
    "gene",
]


class TestAddFunction(unittest.TestCase):
    def test_add_positive_numbers(self):
        pass
