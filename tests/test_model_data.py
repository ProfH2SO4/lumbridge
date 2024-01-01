import unittest, os
from enum import Enum

from .test_config import INPUT_FASTA, OUTPUT_FOLDER


class WriteRule(Enum):
    START = 1
    MIDDLE = 2
    END = 3


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


def parse_line(line):
    parts = line.strip().split("\t")
    start = int(parts[0])
    end = int(parts[1])
    # The rest of the parts can be processed as needed, depending on their expected types
    return start, end, parts[2:]


def get_vector_for_base(base: str, vector_schema: list[str]) -> int:
    """
    Get the position of the base in the vector_schema.

    :param base: A string representing a base pair (e.g., 'A', 'C', 'G', 'T').
    :param vector_schema: A list of base pair representations.
    :return: The position (index) of the base in the vector_schema.
    """
    try:
        return vector_schema.index(base)
    except ValueError:
        # Handle the case where the base is not in vector_schema
        return -1  # Or any other default value or error handling as needed


def process_model_data_line(line):
    # Removing leading and trailing brackets and splitting
    return [eval(vector) for vector in line.strip().split("\t")]


class TestModelData(unittest.TestCase):
    def test_brackets(self):
        model_data_folder: str = f"{OUTPUT_FOLDER}/model_data"
        model_data_file = os.path.join(model_data_folder, "arabidopsis_test.txt")

        with open(model_data_file, "r") as file:
            line_number = 0
            for line in file:
                line_number += 1
                if line.startswith("#"):
                    continue  # Skip header lines

                # Count opening and closing brackets
                open_brackets = line.count("[")
                close_brackets = line.count("]")

                assert (
                    open_brackets == close_brackets
                ), f"Unmatched brackets in line {line_number}: {line}"

    def test_fasta_data(self):
        model_data_folder: str = f"{OUTPUT_FOLDER}/model_data"

        # Load the FASTA file
        fasta_file = os.path.join(
            INPUT_FASTA, "arabidopsis_test.fasta"
        )  # Replace with your actual file name
        with open(fasta_file, "r") as file:
            fasta_data = file.read().splitlines()
        fasta_data = [
            line for line in fasta_data if not line.startswith(">")
        ]  # Remove header lines
        fasta_sequence = "".join(fasta_data)  # Concatenate to a single string

        # Load the model data file
        model_data_file = os.path.join(
            model_data_folder, "arabidopsis_test.txt"
        )  # Replace with your actual file name
        with open(model_data_file, "r") as file:
            model_data = file.read().splitlines()
        model_data = [
            line for line in model_data if not line.startswith("#")
        ]  # Remove header lines
        model_data_sequence = []
        for line in model_data:
            model_data_sequence.extend(process_model_data_line(line))

        # Iterate over each row in the FASTA data and check the corresponding vectors in model data
        for i, base in enumerate(fasta_sequence):
            expected_pos = get_vector_for_base(base, bp_vector_schema)
            vector_ = model_data_sequence[i]
            # Assert that the expected position is set (i.e., has a value of 1)
            assert (
                vector_[expected_pos] == 1
            ), f"Mismatch at position {i} for base {base}, expected position {expected_pos}"

    def test_orf(self):
        model_data_folder: str = f"{OUTPUT_FOLDER}/model_data"
        orf_output_folder: str = f"{OUTPUT_FOLDER}/orf_in_gff3"

        # Load the ORF file
        orf_output_file = os.path.join(orf_output_folder, "arabidopsis_test.txt")
        with open(orf_output_file, "r") as file:
            orf_data = file.readlines()[1:]

        # Parse the ORF file to extract ORF start and end positions where 'In_Gene' is not 0
        orf_positions = []
        for line in orf_data:
            parts = line.strip().split("\t")
            if parts[-1] != "0":
                orf_start, orf_end = int(parts[0]), int(parts[1])
                orf_positions.append((orf_start, orf_end))

        # Load the model data file
        model_data_file = os.path.join(model_data_folder, "arabidopsis_test.txt")
        with open(model_data_file, "r") as file:
            model_data = file.read().splitlines()
        model_data = [line for line in model_data if not line.startswith("#")]
        model_data_sequence = []
        for line in model_data:
            model_data_sequence.extend(process_model_data_line(line))

        # Iterate over each sequence in the ORF data and check the corresponding vectors in model data
        for start, end in orf_positions:
            for i in range(start, end + 1):
                vector_ = model_data_sequence[i - 1]
                position_to_look: int = bp_vector_schema.index("ORF")

                actual_rules = set(vector_[position_to_look])

                # Check for start and end rules
                if i in [start_pos for start_pos, _ in orf_positions]:
                    assert (
                        1 in actual_rules
                    ), f"ORF start rule not marked at position {i}"
                if i in [end_pos for _, end_pos in orf_positions]:
                    assert 3 in actual_rules, f"ORF end rule not marked at position {i}"

                # Check for middle rule, considering overlapping ORFs
                is_middle_position = any(
                    start_pos < i < end_pos for start_pos, end_pos in orf_positions
                )
                if 2 in actual_rules:
                    assert (
                        is_middle_position
                    ), f"ORF middle rule incorrectly marked at position {i}"
                elif is_middle_position:
                    assert (
                        2 not in actual_rules
                    ), f"Missing ORF middle rule at position {i}"

    def test_gff3(self):
        model_data_folder: str = f"{OUTPUT_FOLDER}/model_data"
        gff3_output_folder: str = f"{OUTPUT_FOLDER}/gff3_positive_strand"

        features_to_check: list[str] = [
            "ORF",
            "exon",
            "mRNA",
            "miRNA",
            "rRNA",
            "CDS",
            "gene",
        ]

        # Load the GFF3 file
        gff3_file: str = os.path.join(gff3_output_folder, "arabidopsis_test.gff3")
        feature_positions = {feature: [] for feature in features_to_check}

        with open(gff3_file, "r") as file:
            for line in file:
                if line.startswith("#"):
                    continue  # Skip header lines
                parts = line.strip().split("\t")
                feature_type, start, end = parts[2], int(parts[3]), int(parts[4])
                if feature_type in features_to_check:
                    feature_positions[feature_type].append((start, end))

        # Load the model data file
        model_data_file = os.path.join(model_data_folder, "arabidopsis_test.txt")
        with open(model_data_file, "r") as file:
            model_data = file.read().splitlines()
        model_data = [line for line in model_data if not line.startswith("#")]
        model_data_sequence = []
        for line in model_data:
            model_data_sequence.extend(process_model_data_line(line))

        # Iterate over each feature type and their positions and test against model data
        for feature in features_to_check:
            position_to_look: int = bp_vector_schema.index(feature)
            for start, end in feature_positions[feature]:
                for i in range(start, end + 1):
                    vector_ = model_data_sequence[
                        i - 1
                    ]  # Adjust for zero-based indexing
                    actual_rules = vector_[position_to_look]

                    # Check if the actual rules list contains the expected rule
                    expected_rule = (
                        WriteRule.START.value
                        if i == start
                        else WriteRule.END.value
                        if i == end
                        else WriteRule.MIDDLE.value
                    )
                    assert (
                        expected_rule in actual_rules
                    ), f"Feature {feature} at position {i} expected rule {expected_rule}, found {actual_rules}, at vector_pos {position_to_look}, vector: {vector_}"

    def test_promotor_motif(self):
        model_data_folder: str = f"{OUTPUT_FOLDER}/model_data"
        promotor_output_folder: str = f"{OUTPUT_FOLDER}/homer2_annotation"

        # Load the promoter data
        promotor_file: str = os.path.join(
            promotor_output_folder, "arabidopsis_test.txt"
        )
        promoter_positions = {}
        with open(promotor_file, "r") as pf:
            next(pf)  # Skip header
            for line in pf:
                start, end, others = parse_line(line)
                promoter_positions.setdefault(start, []).append(WriteRule.START.value)
                for pos in range(start + 1, end):
                    promoter_positions.setdefault(pos, []).append(
                        WriteRule.MIDDLE.value
                    )
                promoter_positions.setdefault(end, []).append(WriteRule.END.value)

        # Load the model data file
        model_data_file = os.path.join(model_data_folder, "arabidopsis_test.txt")
        with open(model_data_file, "r") as file:
            model_data = [line for line in file if not line.startswith("#")]

        # Process the model data
        model_data_sequence = []
        for line in model_data:
            vectors = [eval(vector) for vector in line.strip().split("\t")]
            model_data_sequence.extend(vectors)

        position_to_check = bp_vector_schema.index("PROMOTOR_MOTIF")

        for index, vector in enumerate(
            model_data_sequence, start=1
        ):  # Assuming 1-based indexing
            actual_types = vector[position_to_check]
            expected_types = promoter_positions.get(index, 0)

            self.assertEqual(
                actual_types,
                expected_types,
                f"Mismatch at position {index}: vector: {vector} Expected {expected_types}, found {actual_types}",
            )

    def test_poly_adenyl(self):
        model_data_folder: str = f"{OUTPUT_FOLDER}/model_data"
        poly_adenylation_folder: str = f"{OUTPUT_FOLDER}/poly_adenylation"

        # Load the polyadenylation data
        poly_adenyl_file = os.path.join(poly_adenylation_folder, "arabidopsis_test.txt")
        poly_adenyl_positions = {}
        with open(poly_adenyl_file, "r") as pf:
            next(pf)  # Skip header
            for line in pf:
                start, end, _ = parse_line(line)
                poly_adenyl_positions.setdefault(start, []).append(
                    WriteRule.START.value
                )
                for pos in range(start + 1, end):
                    poly_adenyl_positions.setdefault(pos, []).append(
                        WriteRule.MIDDLE.value
                    )
                poly_adenyl_positions.setdefault(end, []).append(WriteRule.END.value)

        # Load the model data file
        model_data_file = os.path.join(model_data_folder, "arabidopsis_test.txt")
        with open(model_data_file, "r") as file:
            model_data = [line for line in file if not line.startswith("#")]

        # Process the model data
        model_data_sequence = []
        for line in model_data:
            vectors = [eval(vector) for vector in line.strip().split("\t")]
            model_data_sequence.extend(vectors)

        position_to_check = bp_vector_schema.index("POLY_ADENYL")

        for index, vector in enumerate(
            model_data_sequence, start=1
        ):  # Assuming 1-based indexing
            actual_signal = vector[position_to_check]
            expected_signal = poly_adenyl_positions.get(index, 0)

            self.assertEqual(
                actual_signal,
                expected_signal,
                f"Mismatch at position {index}: Expected {expected_signal}, found {actual_signal}",
            )
