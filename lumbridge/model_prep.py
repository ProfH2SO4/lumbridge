import os
from enum import IntEnum
from datetime import datetime

from .common import create_folder

__all__ = ["transform_data_to_vectors"]


class WriteRule(IntEnum):
    START = 1
    MIDDLE = 2
    END = 3


def calculate_position(index_feature: int, max_feature_overlap: int) -> tuple[int, int]:
    feature_len = 1  # Length of a single feature representation

    start_position = 4 + (index_feature - 4) * feature_len * (max_feature_overlap + 1)
    end_position = start_position + feature_len - 1
    return (start_position, end_position)


def parse_line(line):
    parts = line.strip().split("\t")
    start = int(parts[0])
    end = int(parts[1])
    # The rest of the parts can be processed as needed, depending on their expected types
    return start, end, parts[2:]


def create_file_header(
    path_to_file: str,
    bp_vector_schema: list[str],
    max_feature_overlap: int,
    version: list[int],
) -> None:
    header_: str = (
        "#HEADER#\n"
        f"#DATE={datetime.utcnow().isoformat()}\n"
        f"#pre_processing_version={version}\n"
        f"#bp_vector_schema={bp_vector_schema}\n"
        "#description of nucleotide:A=[1, 0, 0, 0], C=[0, 1, 0, 0], G=[0, 0, 1, 0], T=[0, 0, 0, 1]\n"
        "#description of feature:0=no_present, 1=start, 2=continuation/ongoing, 3=end\n"
        f"#max_feature_overlap={max_feature_overlap}\n"
        "####END####\n"
    )
    with open(path_to_file, "w") as f:
        f.write(header_)


def write_fasta_file(
    fasta_file: str,
    model_file: str,
    bp_vector_schema: list[str],
    max_feature_overlap: int,
):
    num_of_features: int = len(bp_vector_schema) - 4
    feature_len: int = 1  # Length of a single feature representation
    total_feature_len: int = feature_len * num_of_features * (max_feature_overlap + 1)

    # Mapping nucleotides to vectors
    nucleotide_to_vector = {
        "A": [1, 0, 0, 0] + [0] * total_feature_len,
        "C": [0, 1, 0, 0] + [0] * total_feature_len,
        "G": [0, 0, 1, 0] + [0] * total_feature_len,
        "T": [0, 0, 0, 1] + [0] * total_feature_len,
    }

    # Read the FASTA file and write to the model file
    with open(fasta_file, "r") as fasta, open(model_file, "a") as model:
        vector_count = 0  # Counter for the number of vectors written to the current line
        for line in fasta:
            if line.startswith(">"):
                continue  # Skip header lines

            for char in line.strip():
                vector = nucleotide_to_vector.get(char.upper(), None)
                if vector:
                    vector_count += 1

                    # Write the vector to the model file
                    if vector_count == 2:
                        # At the fifth vector, write it followed by a newline and reset the counter
                        model.write(f"{vector}\n")
                        vector_count = 0
                    else:
                        # Otherwise, write the vector followed by a tab
                        model.write(f"{vector}\t")


def parse_promoter_file(promoter_file: str) -> dict:
    """
    Parses the promoter file to identify the positions of promoters.

    :param promoter_file: The path to the promoter file.
    :return: A dictionary mapping promoter positions to their respective WriteRule.
    """
    promoter_positions = {}
    with open(promoter_file, "r") as file:
        next(file)  # Skip header
        for line in file:
            start, end, _ = parse_line(line)
            promoter_positions.setdefault(start, []).append(WriteRule.START)
            for pos in range(start + 1, end):
                promoter_positions.setdefault(pos, []).append(WriteRule.MIDDLE)
            promoter_positions.setdefault(end, []).append(WriteRule.END)
    return promoter_positions


def parse_gff3_file(gff3_file_path: str, feature_positions: dict[str, dict]) -> None:
    """
    Parses the gff3 file to identify the positions in feature_positions.

    :param gff3_file_path: The path to the gff3 file.
    :param feature_positions: features to extract
    :return: None.
    """
    with open(gff3_file_path, "r") as gff3:
        for line in gff3:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            feature_type = parts[2]
            if feature_type in feature_positions:
                start, end = int(parts[3]), int(parts[4])
                # Append feature status to lists at each position
                feature_positions[feature_type].setdefault(start, []).append(WriteRule.START)
                for pos in range(start + 1, end):
                    feature_positions[feature_type].setdefault(pos, []).append(WriteRule.MIDDLE)
                feature_positions[feature_type].setdefault(end, []).append(WriteRule.END)


def update_model_with_features(
    model_file: str, promoter_positions: dict, position_to_write: tuple[int, int], max_overlap: int
):
    """
    Updates the model file with promoter positions.

    :param model_file: The path to the model file.
    :param promoter_positions: A dictionary of promoter positions and their types.
    :param position_to_write: The start and end positions to write the promoters.
    :param max_overlap: The maximum allowed overlap for features.
    :return:
    """
    temp_file_path = model_file + ".tmp"
    with open(model_file, "r") as model, open(temp_file_path, "w") as temp_file:
        model_vector_index = -1
        for line in model:
            if line.startswith("#"):
                temp_file.write(line)
                continue

            updated_line, model_vector_index = update_line(
                line, promoter_positions, position_to_write, max_overlap, model_vector_index
            )
            temp_file.write(updated_line)

    # Replace the old model file with the updated one
    os.replace(temp_file_path, model_file)


def update_model_with_features_gff3(
    model_file: str, feature_positions: dict, full_bp_vector_schema: list[str], max_overlap: int
) -> None:
    temp_file_path = model_file + ".tmp"
    index = -1
    with open(model_file, "r") as mf, open(temp_file_path, "w") as temp_file:
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)
                continue

            vectors = [eval(vector) for vector in line.strip().split("\t")]
            updated_line = []
            for _, vector in enumerate(vectors):
                index += 1
                for feature, positions in feature_positions.items():
                    if index in positions:
                        feature_index = full_bp_vector_schema.index(feature)
                        for offset in range(max_overlap + 1):
                            pos = calculate_position(feature_index, max_overlap)[0] + offset
                            if vector[pos] == 0:  # Write only the first matching rule
                                vector[pos] = positions[index][offset].value
                                break  # Break after finding a space to write

                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")
    os.rename(temp_file_path, model_file)


def update_line(
    line: str, promoter_positions: dict, position_to_write: tuple[int, int], max_overlap: int, vector_index: int
) -> str:
    """
    Updates a line of the model file based on promoter positions.

    :param vector_index: Index of the vector in the model file
    :param line: A line from the model file.
    :param promoter_positions: A dictionary of promoter positions and their types.
    :param position_to_write: The start and end positions to write the promoters.
    :param max_overlap: The maximum allowed overlap for features.
    :return: The updated line for the model file.
    """
    vectors = [eval(vector) for vector in line.strip().split("\t")]
    updated_vectors = []
    for vector in vectors:
        vector_index += 1
        if vector_index in promoter_positions:
            promoter_type = promoter_positions[vector_index]
            for offset in range(max_overlap + 1):
                pos = position_to_write[0] + offset
                if vector[pos] == 0:  # Write only the first matching rule
                    vector[pos] = promoter_type[offset].value
                    break  # Break after finding a space to write
        updated_vectors.append(str(vector))
    return "\t".join(updated_vectors) + "\n", vector_index


def write_promotor_motifs_file(
    promoter_file: str, model_file: str, position_in_schema: int, max_feature_overlap: int
) -> None:
    """
    Updates the model file with promoter motifs based on the promoter file.

    :param promoter_file: The path to the promoter file.
    :param model_file: The path to the model file.
    :param position_in_schema: The position in the schema for writing promoters.
    :param max_feature_overlap: The maximum allowed overlap for features.
    :return:
    """
    position_to_write = calculate_position(position_in_schema, max_feature_overlap)
    promoter_positions = parse_promoter_file(promoter_file)
    update_model_with_features(model_file, promoter_positions, position_to_write, max_feature_overlap)


def write_gff3_file(
    gff3_file: str,
    model_file: str,
    core_bp_vector_schema: list[str],
    full_bp_vector_schema: list[str],
    max_feature_overlap: int,
) -> None:
    """


    :param gff3_file:
    :param model_file:
    :param core_bp_vector_schema:
    :param full_bp_vector_schema:
    :param max_feature_overlap:
    :return:
    """
    feature_positions: dict[str, dict] = {
        feature: {} for feature in full_bp_vector_schema if feature not in core_bp_vector_schema
    }
    parse_gff3_file(gff3_file, feature_positions)
    update_model_with_features_gff3(model_file, feature_positions, full_bp_vector_schema, max_feature_overlap)


def parse_orf_file(path_orf_file: str) -> dict:
    orf_positions = {}
    with open(path_orf_file, "r") as of:
        next(of)  # Skip header
        for line in of:
            parts = line.strip().split("\t")
            if parts[4] != "0":  # Check if 'In_Gene' is not '0'
                start, end = int(parts[0]), int(parts[1])
                # Add ORF rules as a list for each position
                for pos in range(start, end + 1):
                    orf_rule = WriteRule.START if pos == start else WriteRule.END if pos == end else WriteRule.MIDDLE
                    orf_positions.setdefault(pos, []).append(orf_rule)
    return orf_positions


def parse_poly_adenyl_file(path_poly_a_file: str) -> dict[int, list[int]]:
    poly_adenyl_positions = {}
    with open(path_poly_a_file, "r") as of:
        next(of)  # Skip header
        for line in of:
            parts = line.strip().split("\t")
            start, end = int(parts[0]), int(parts[1])
            poly_adenyl_positions.setdefault(start, []).append(WriteRule.START)
            for pos in range(start + 1, end):
                poly_adenyl_positions.setdefault(pos, []).append(WriteRule.MIDDLE)
            poly_adenyl_positions.setdefault(end, []).append(WriteRule.END)
    return poly_adenyl_positions


def write_orf_file(orf_file: str, model_file: str, position_in_schema: int, max_feature_overlap: int) -> None:
    """
    Updates the model file with opening reading frame based on the orf file.

    :param orf_file: The path to the orf file.
    :param model_file: The path to the model file.
    :param position_in_schema: The position in the schema for writing promoters.
    :param max_feature_overlap: The maximum allowed overlap for features.
    :return:
    """
    position_to_write = calculate_position(position_in_schema, max_feature_overlap)
    promoter_positions = parse_orf_file(orf_file)
    update_model_with_features(model_file, promoter_positions, position_to_write, max_feature_overlap)


def write_poly_a_file(poly_a_file: str, model_file: str, position_in_schema: int, max_feature_overlap: int) -> None:
    """
    Updates the model file with poly_a  based on the poly_a  file.

    :param poly_a_file: The path to the poly_a  file.
    :param model_file: The path to the model file.
    :param position_in_schema: The position in the schema for writing promoters.
    :param max_feature_overlap: The maximum allowed overlap for features.
    :return:
    """
    position_to_write = calculate_position(position_in_schema, max_feature_overlap)
    poly_a_positions = parse_poly_adenyl_file(poly_a_file)
    update_model_with_features(model_file, poly_a_positions, position_to_write, max_feature_overlap)


def transform_data_to_vectors(
    fasta_folder: str,
    gff3_folder: str,
    orf_folder: str,
    promotor_motifs_folder: str,
    poly_a_folder: str,
    output_folder: str,
    max_feature_overlap: int,
    gff3_features: list[str],
    version: list[int],
) -> None:
    """
    Transforms genomic data from various input files into vector representations, placing features at their respective
    positions within the vectors. This process facilitates the analysis and interpretation of genomic elements,
    such as open reading frames (ORFs), promoter motifs, and polyadenylation signals,
    by converting them into a numerical format suitable for computational models.

    :param fasta_folder: Directory containing FASTA files with nucleotide sequences.
    :param gff3_folder: Directory containing GFF3 files with annotations for genomic features.
    :param orf_folder: Directory containing files with identified open reading frames (ORFs).
    :param promotor_motifs_folder: Directory containing files with promoter motif sequences.
    :param poly_a_folder: Directory containing files with polyadenylation signal sequences.
    :param output_folder: Directory where the output vector files will be saved.
    :param max_feature_overlap: Maximum allowed overlap between features in the vector representation.
    :param gff3_features: List of feature types from GFF3 files to include in the vector representation.
    :param version: Version of the vector representation format, represented as a list of integers.
    :return: None. The function does not return a value but writes the vector representations to the specified output directory.
    """
    core_bp_vector_schema: list[str] = [
        "A",
        "C",
        "G",
        "T",
        "PROMOTOR_MOTIF",
        "ORF",
        "POLY_ADENYL",
    ]
    full_bp_vector_schema: list[str] = core_bp_vector_schema + gff3_features
    create_folder(output_folder)
    for filename in os.listdir(fasta_folder):
        ending = ".fasta"
        if not filename.endswith(ending):
            continue
        base_part, extension = os.path.splitext(filename)
        output_file: str = f"{output_folder}/{base_part}.txt"
        create_file_header(output_file, full_bp_vector_schema, max_feature_overlap, version)
        write_fasta_file(
            f"{fasta_folder}/{filename}",
            output_file,
            full_bp_vector_schema,
            max_feature_overlap,
        )
        write_promotor_motifs_file(
            f"{promotor_motifs_folder}/{base_part}.txt",
            output_file,
            full_bp_vector_schema.index("PROMOTOR_MOTIF"),
            max_feature_overlap,
        )
        write_orf_file(
            f"{orf_folder}/{base_part}.txt",
            output_file,
            full_bp_vector_schema.index("ORF"),
            max_feature_overlap,
        )
        write_poly_a_file(
            f"{poly_a_folder}/{base_part}.txt",
            output_file,
            full_bp_vector_schema.index("POLY_ADENYL"),
            max_feature_overlap,
        )
        write_gff3_file(
            f"{gff3_folder}/{base_part}.gff3",
            output_file,
            core_bp_vector_schema,
            full_bp_vector_schema,
            max_feature_overlap,
        )
