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


def load_gff3_data(gff3_file: str, bp_vector_schema: list[str], max_overlap: int) -> dict:
    feature_pos: dict[str, list] = {}

    def is_overlap_allowed(
        existing_intervals: list[tuple[int, int]], new_interval: tuple[int, int], max_overlap: int
    ) -> bool:
        """Check if the new interval overlaps with existing intervals more than the allowed max_overlap."""
        start_new, end_new = new_interval
        for start_existing, end_existing in existing_intervals:
            # Calculate overlap
            overlap = min(end_existing, end_new) - max(start_existing, start_new) + 1
            # If overlap is more than max_overlap, return False
            if overlap > max_overlap:
                return False
        return True

    with open(gff3_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip header lines
            parts = line.strip().split("\t")
            feature_type, start, end = parts[2], int(parts[3]), int(parts[4])
            if feature_type in bp_vector_schema:
                new_interval = (start, end)
                # Only append the new interval if it's allowed
                if is_overlap_allowed(feature_pos.setdefault(feature_type, []), new_interval, max_overlap):
                    feature_pos[feature_type].append(new_interval)

    return feature_pos


def load_orf_data(orf_file_path: str) -> dict[int, list[int]]:
    orf_positions = {}
    with open(orf_file_path, "r") as file:
        next(file)  # Skip header
        # Parse the ORF file to extract ORF start and end positions where 'In_Gene' is not 0
        for line in file:
            parts = line.strip().split("\t")
            if parts[-1] != "0":
                start, end = int(parts[0]), int(parts[1])
                orf_positions.setdefault(start, []).append(WriteRule.START.value)
                for pos in range(start + 1, end):
                    orf_positions.setdefault(pos, []).append(WriteRule.MIDDLE.value)
                orf_positions.setdefault(end, []).append(WriteRule.END.value)
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


def compare_feature_sequences(positions: dict[int, list[int]], model_data_sequence, position_to_check):
    for index, vector in enumerate(model_data_sequence, start=0):  # Assuming 1-based indexing
        actual_types = vector[position_to_check]
        expected_types = positions.get(index, 0)
        if isinstance(expected_types, list):  # when found in positions
            expected_types = expected_types[0]
        assert (
            actual_types == expected_types
        ), f"Mismatch at position {index}: vector: {vector} Expected {expected_types}, found {actual_types}"


def compare_gff3_sequences(features_pos: dict[str, list], model_data_sequence, bp_vector_schema: list[str]):
    # Iterate over each feature type and their positions and test against model data
    for feature, pos in features_pos.items():
        position_to_look: int = bp_vector_schema.index(feature)
        for start, end in pos:
            for i in range(start, end + 1):
                vector_ = model_data_sequence[i]
                actual_rules = vector_[position_to_look]

                # Check if the actual rules list contains the expected rule
                expected_rule = (
                    WriteRule.START.value if i == start else WriteRule.END.value if i == end else WriteRule.MIDDLE.value
                )
                if isinstance(expected_rule, list):
                    expected_rule = expected_rule[0]
                assert (
                    expected_rule == actual_rules
                ), f"Feature {feature} at position {i} expected rule {expected_rule}, found {actual_rules}, at vector_pos {position_to_look}, vector: {vector_}"


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
            orf_positions: dict[int, list[int]] = load_orf_data(orf_file_path)

            bp_vector_schema = extract_bp_vector_schema(model_data_file_path)
            model_data_sequence = load_model_data(model_data_file_path)
            position_to_check = bp_vector_schema.index("ORF")
            compare_feature_sequences(orf_positions, model_data_sequence, position_to_check)

    def test_promotor_motif(self):
        model_data_folder = f"{self.output_dir_path}/model_data"
        promotor_output_folder: str = f"{self.output_dir_path}/homer2_annotation"
        for filename in os.listdir(model_data_folder):
            model_data_file_path = os.path.join(model_data_folder, filename)
            file_base = filename[:-4]  # Assuming the extension is always 4 characters long, e.g., .txt

            promotor_file_path = os.path.join(promotor_output_folder, f"{file_base}.txt")
            promotor_positions: dict = load_feature_data(promotor_file_path)

            bp_vector_schema = extract_bp_vector_schema(model_data_file_path)
            model_data_sequence = load_model_data(model_data_file_path)
            position_to_check = bp_vector_schema.index("PROMOTOR_MOTIF")

            compare_feature_sequences(promotor_positions, model_data_sequence, position_to_check)

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

            compare_feature_sequences(poly_adenyl_positions, model_data_sequence, position_to_check)

    def test_gff3(self):
        model_data_folder = f"{self.output_dir_path}/model_data"
        gff3_output_folder: str = f"{self.output_dir_path}/gff3_positive_strand"

        for filename in os.listdir(model_data_folder):
            model_data_file_path = os.path.join(model_data_folder, filename)
            file_base = filename[:-4]  # Assuming the extension is always 4 characters long, e.g., .txt

            bp_vector_schema = extract_bp_vector_schema(model_data_file_path)

            gff3_output_file_path = os.path.join(gff3_output_folder, f"{file_base}.gff3")
            feature_positions: dict[str, list[tuple]] = load_gff3_data(gff3_output_file_path, bp_vector_schema, 0)

            model_data_sequence = load_model_data(model_data_file_path)
            compare_gff3_sequences(feature_positions, model_data_sequence, bp_vector_schema)
