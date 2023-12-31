import unittest, os

from .test_config import INPUT_FASTA, INPUT_GFF3, OUTPUT_FOLDER

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


class TestModelData(unittest.TestCase):
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
        model_data_sequence = [
            eval(vector)
            for row in model_data
            for vector in row.strip("[]").split("],[")
        ]

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
        model_data_sequence = [
            eval(vector)
            for row in model_data
            for vector in row.strip("[]").split("],[")
        ]

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
