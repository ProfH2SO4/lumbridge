import ast
import unittest, os
from enum import Enum


class WriteRule(Enum):
    START = 1
    MIDDLE = 2
    END = 3


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


def extract_bp_vector_schema(file_path: str) -> list[str]:
    """
    Extracts the bp_vector_schema from a given file.
    """
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#bp_vector_schema="):
                # Extracting the list part of the line and evaluating it
                schema_str = line.split("=", 1)[1].strip()
                # Using eval to convert the string representation of the list to an actual list
                bp_vector_schema = eval(schema_str)
                return bp_vector_schema


def load_feature_data(path_promotor_file: str) -> dict:
    promoter_positions = {}
    with open(path_promotor_file, "r") as pf:
        next(pf)  # Skip header
        for line in pf:
            start, end, others = parse_line(line)
            promoter_positions.setdefault(start, []).append(WriteRule.START.value)
            for pos in range(start + 1, end):
                promoter_positions.setdefault(pos, []).append(WriteRule.MIDDLE.value)
            promoter_positions.setdefault(end, []).append(WriteRule.END.value)
    return promoter_positions


def load_orf_data(orf_file_path: str) -> str:
    with open(orf_file_path, "r") as file:
        orf_data = file.readlines()[1:]

        # Parse the ORF file to extract ORF start and end positions where 'In_Gene' is not 0
        orf_positions = []
        for line in orf_data:
            parts = line.strip().split("\t")
            if parts[-1] != "0":
                orf_start, orf_end = int(parts[0]), int(parts[1])
                orf_positions.append((orf_start, orf_end))

        return orf_positions


def load_fasta_data(fasta_file_path: str) -> str:
    """Load and concatenate sequence data from a FASTA file, excluding headers."""
    with open(fasta_file_path, "r") as file:
        fasta_data = [line.strip() for line in file if not line.startswith(">")]
    return "".join(fasta_data)


def process_model_data_line(line: str) -> list:
    """Process a single line of model data, splitting by tabs and processing each part."""
    return [process_vector(part) for part in line.split("\t")]


def process_vector(vector: str) -> list[int]:
    processed_vector = [int(element) for element in ast.literal_eval(vector)]
    return processed_vector


def load_model_data(model_data_file_path: str) -> list:
    """Load model data from a file, excluding header lines, and split each line by tabs."""
    vectors = []
    with open(model_data_file_path, "r") as file:
        # Exclude header lines and strip whitespace
        for line in file:
            if not line.startswith("#"):
                [vectors.append(i) for i in process_model_data_line(line)]
    # Process each line by splitting by tabs and processing each part
    return vectors


def compare_sequences(fasta_sequence: str, model_data_sequence: list, bp_vector_schema: list):
    """Assert that the FASTA sequence and model data sequence match according to the bp_vector_schema."""
    for i, base in enumerate(fasta_sequence):
        expected_pos = get_vector_for_base(base, bp_vector_schema)
        vector_ = model_data_sequence[i]
        assert vector_[expected_pos] == 1, f"Mismatch at position {i} for base {base}, expected position {expected_pos}"


def compare_orf_sequences(
    orf_positions: list[tuple[int, int]], model_data_sequence: list[list[int]], bp_vector_schema: list[str]
):
    # Iterate over each sequence in the ORF data and check the corresponding vectors in model data
    position_to_look: int = bp_vector_schema.index("ORF")
    for start, end in orf_positions:
        for i in range(start, end + 1):
            vector_ = model_data_sequence[i]  # Access the vector at the current position

            # Check the start of the ORF
            if i == start:
                assert (
                    vector_[position_to_look] == WriteRule.START
                ), f"Expected WriteRule.START at position {i}, got {vector_[position_to_look]}"
            # Check the end of the ORF
            elif i == end:
                assert (
                    vector_[position_to_look] == WriteRule.END
                ), f"Expected WriteRule.END at position {i}, got {vector_[position_to_look]}"
            # Check the middle of the ORF
            else:
                assert (
                    vector_[position_to_look] == WriteRule.MIDDLE
                ), f"Expected WriteRule.MIDDLE at position {i}, got {vector_[position_to_look]}"


def compare_feature_sequences(promotor_positions, model_data_sequence, bp_vector_schema, position_to_check):

    for index, vector in enumerate(model_data_sequence, start=1):  # Assuming 1-based indexing
        actual_types = vector[position_to_check]
        expected_types = promotor_positions.get(index, 0)

        assert (
            actual_types,
            expected_types,
            f"Mismatch at position {index}: vector: {vector} Expected {expected_types}, found {actual_types}",
        )


def test_model_data_factory(input_dir_path_fasta, input_dir_path_gff3, input_dir_path_orf, output_dir_path):
    def _init(self, methodName="runTest"):
        TestModelData.__init__(
            self, methodName, input_dir_path_fasta, input_dir_path_gff3, input_dir_path_orf, output_dir_path
        )

    return type(f"TestModelData", (TestModelData,), {"__init__": _init})


class TestModelData(unittest.TestCase):
    def __init__(
        self, methodName: str, input_dir_path_fasta, input_dir_path_gff3, input_dir_path_orf, output_dir_path: str
    ):
        super(TestModelData, self).__init__(methodName)
        self.input_dir_path_fasta = input_dir_path_fasta
        self.input_dir_path_gff3 = input_dir_path_gff3
        self.input_dir_path_orf = input_dir_path_orf
        self.output_dir_path = output_dir_path

    def test_brackets(self):
        model_data_folder: str = f"{self.output_dir_path}/model_data"
        for filename in os.listdir(model_data_folder):
            filepath = os.path.join(model_data_folder, filename)

            with open(filepath, "r") as file:
                line_number = 0
                for line in file:
                    line_number += 1
                    if line.startswith("#"):
                        continue  # Skip header lines

                    # Count opening and closing brackets
                    open_brackets = line.count("[")
                    close_brackets = line.count("]")

                    assert open_brackets == close_brackets, f"Unmatched brackets in line {line_number}: {line}"

    def test_fasta_data(self):
        model_data_folder = f"{self.output_dir_path}/model_data"

        for filename in os.listdir(model_data_folder):
            model_data_file_path = os.path.join(model_data_folder, filename)
            file_base = filename[:-4]  # Assuming the extension is always 4 characters long, e.g., .txt

            fasta_file_path = os.path.join(self.input_dir_path_fasta, f"{file_base}.fasta")
            fasta_sequence = load_fasta_data(fasta_file_path)

            bp_vector_schema = extract_bp_vector_schema(model_data_file_path)
            model_data_sequence = load_model_data(model_data_file_path)

            compare_sequences(fasta_sequence, model_data_sequence, bp_vector_schema)

    def test_orf(self):
        model_data_folder = f"{self.output_dir_path}/model_data"
        for filename in os.listdir(model_data_folder):
            model_data_file_path = os.path.join(model_data_folder, filename)
            file_base = filename[:-4]  # Assuming the extension is always 4 characters long, e.g., .txt

            orf_file_path = os.path.join(f"{self.output_dir_path}/orf_in_gff3", f"{file_base}.txt")
            orf_positions: list[tuple] = load_orf_data(orf_file_path)

            bp_vector_schema = extract_bp_vector_schema(model_data_file_path)
            model_data_sequence = load_model_data(model_data_file_path)

            compare_orf_sequences(orf_positions, model_data_sequence, bp_vector_schema)

    def test_promotor_motif(self):
        model_data_folder = f"{self.output_dir_path}/model_data"
        promotor_output_folder: str = f"{self.output_dir_path}/homer2_annotation"
        for filename in os.listdir(model_data_folder):
            model_data_file_path = os.path.join(model_data_folder, filename)
            file_base = filename[:-4]  # Assuming the extension is always 4 characters long, e.g., .txt

            promotor_file_path = os.path.join(promotor_output_folder, f"{file_base}.txt")
            promotor_positions: list[tuple] = load_feature_data(promotor_file_path)

            bp_vector_schema = extract_bp_vector_schema(model_data_file_path)
            model_data_sequence = load_model_data(model_data_file_path)
            position_to_check = bp_vector_schema.index("PROMOTOR_MOTIF")

            compare_feature_sequences(promotor_positions, model_data_sequence, bp_vector_schema, position_to_check)

    def test_poly_adenyl(self):
        model_data_folder = f"{self.output_dir_path}/model_data"
        poly_adenylation_folder: str = f"{self.output_dir_path}/poly_adenylation"

        for filename in os.listdir(model_data_folder):
            model_data_file_path = os.path.join(model_data_folder, filename)
            file_base = filename[:-4]  # Assuming the extension is always 4 characters long, e.g., .txt

            poly_adenyl_file_path = os.path.join(poly_adenylation_folder, f"{file_base}.txt")
            poly_adenyl_positions: dict[int, list[int]] = load_feature_data(poly_adenyl_file_path)

            bp_vector_schema = extract_bp_vector_schema(model_data_file_path)
            model_data_sequence = load_model_data(model_data_file_path)
            position_to_check = bp_vector_schema.index("POLY_ADENYL")

            compare_feature_sequences(poly_adenyl_positions, model_data_sequence, bp_vector_schema, position_to_check)

    # def test_gff3(self):
    #     model_data_folder: str = f"{OUTPUT_FOLDER}/model_data"
    #     gff3_output_folder: str = f"{OUTPUT_FOLDER}/gff3_positive_strand"
    #
    #     features_to_check: list[str] = [
    #         "ORF",
    #         "exon",
    #         "mRNA",
    #         "miRNA",
    #         "rRNA",
    #         "CDS",
    #         "gene",
    #     ]
    #
    #     # Load the GFF3 file
    #     gff3_file: str = os.path.join(gff3_output_folder, "arabidopsis_test.gff3")
    #     feature_positions = {feature: [] for feature in features_to_check}
    #
    #     with open(gff3_file, "r") as file:
    #         for line in file:
    #             if line.startswith("#"):
    #                 continue  # Skip header lines
    #             parts = line.strip().split("\t")
    #             feature_type, start, end = parts[2], int(parts[3]), int(parts[4])
    #             if feature_type in features_to_check:
    #                 feature_positions[feature_type].append((start, end))
    #
    #     # Load the model data file
    #     model_data_file = os.path.join(model_data_folder, "arabidopsis_test.txt")
    #     with open(model_data_file, "r") as file:
    #         model_data = file.read().splitlines()
    #     model_data = [line for line in model_data if not line.startswith("#")]
    #     model_data_sequence = []
    #     for line in model_data:
    #         model_data_sequence.extend(process_model_data_line(line))
    #
    #     # Iterate over each feature type and their positions and test against model data
    #     for feature in features_to_check:
    #         position_to_look: int = bp_vector_schema.index(feature)
    #         for start, end in feature_positions[feature]:
    #             for i in range(start, end + 1):
    #                 vector_ = model_data_sequence[
    #                     i - 1
    #                 ]  # Adjust for zero-based indexing
    #                 actual_rules = vector_[position_to_look]
    #
    #                 # Check if the actual rules list contains the expected rule
    #                 expected_rule = (
    #                     WriteRule.START.value
    #                     if i == start
    #                     else WriteRule.END.value
    #                     if i == end
    #                     else WriteRule.MIDDLE.value
    #                 )
    #                 assert (
    #                     expected_rule in actual_rules
    #                 ), f"Feature {feature} at position {i} expected rule {expected_rule}, found {actual_rules}, at vector_pos {position_to_look}, vector: {vector_}"
    #
