import os, tempfile
from enum import Enum
from datetime import datetime

from .common import create_folder

__all__ = ["transform_data_to_vectors"]


class WriteRule(Enum):
    START = [1, 0, 0]
    MIDDLE = [0, 1, 0]
    END = [0, 0, 1]


def calculate_position(index_feature: int, max_feature_overlap: int) -> tuple[int, int]:
    feature_len = 3  # Length of a single feature representation

    start_position = 4 + (index_feature - 4) * feature_len * max_feature_overlap
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
        "#description of feature:[0, 0, 0]=no_present, [1, 0, 0]=start, [0, 1, 0]=continuation/ongoing, [0, 0, 1]=end\n"
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
    feature_len: int = 3  # Length of a single feature representation
    total_feature_len: int = feature_len * num_of_features * max_feature_overlap

    # Mapping nucleotides to vectors
    nucleotide_to_vector = {
        "A": [1, 0, 0, 0] + [0] * total_feature_len,
        "C": [0, 1, 0, 0] + [0] * total_feature_len,
        "G": [0, 0, 1, 0] + [0] * total_feature_len,
        "T": [0, 0, 0, 1] + [0] * total_feature_len,
    }

    # Read the FASTA file and write to the model file
    with open(fasta_file, "r") as fasta, open(model_file, "a") as model:
        vector_count = (
            0  # Counter for the number of vectors written to the current line
        )
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


def write_promotor_motifs_file(
    promotor_file: str,
    model_file: str,
    position_in_schema: int,
    max_feature_overlap: int,
) -> None:

    position_to_write: tuple[int, int] = calculate_position(
        position_in_schema, max_feature_overlap
    )

    # Read promoter file
    promoter_positions = {}
    with open(promotor_file, "r") as pf:
        next(pf)  # Skip header
        for line in pf:
            start, end, other_parts = parse_line(
                line
            )  # Assuming parse_line is defined elsewhere
            promoter_positions.setdefault(start, []).append(WriteRule.START)
            for pos in range(start + 1, end):
                promoter_positions.setdefault(pos, []).append(WriteRule.MIDDLE)
            promoter_positions.setdefault(end, []).append(WriteRule.END)

    # Update the model file
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
                if index in promoter_positions:
                    promoter_types = promoter_positions[index]
                    for offset in range(0, max_feature_overlap * 3, 3):
                        pos = position_to_write[0] + offset
                        if pos + 2 < len(vector) and vector[pos : pos + 3] == [0, 0, 0]:
                            rule_to_write = promoter_types[
                                0
                            ]  # Write only the first matching rule
                            vector[pos : pos + 3] = rule_to_write.value
                            break  # Break after finding a space to write

                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")

    # Replace the old model file with the updated one
    os.rename(temp_file_path, model_file)


def write_gff3_file(
    gff3_file: str,
    model_file: str,
    bp_vector_schema: list[str],
    max_feature_overlap: int,
) -> None:
    # Initialize dictionaries to store positions for each feature
    feature_positions = {
        feature: {}
        for feature in bp_vector_schema
        if feature not in ["A", "C", "G", "T", "PROMOTOR_MOTIF", "ORF", "POLY_ADENYL"]
    }

    # Read GFF3 file and get positions for each feature
    with open(gff3_file, "r") as gff3:
        for line in gff3:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            feature_type = parts[2]
            if feature_type in feature_positions:
                start, end = int(parts[3]), int(parts[4])
                # Append feature status to lists at each position
                for pos in range(start, end + 1):
                    feature_positions[feature_type].setdefault(pos, []).append(
                        WriteRule.MIDDLE.value
                        if pos != start and pos != end
                        else WriteRule.START.value
                        if pos == start
                        else WriteRule.END.value
                    )

    # Open the model file for reading and a temporary file for writing
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
                        feature_index = bp_vector_schema.index(feature)
                        for offset in range(0, max_feature_overlap * 3, 3):
                            pos = (
                                calculate_position(feature_index, max_feature_overlap)[
                                    0
                                ]
                                + offset
                            )
                            if pos + 2 < len(vector) and vector[pos : pos + 3] == [
                                0,
                                0,
                                0,
                            ]:
                                vector[pos : pos + 3] = positions[index][
                                    0
                                ]  # Write only the first matching rule
                                break  # Break after finding a space to write

                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")
    os.rename(temp_file_path, model_file)


def write_orf_file(
    orf_file: str, model_file: str, position_in_schema: int, max_feature_overlap: int
) -> None:
    # Calculate positions to write ORF information
    position_to_write: tuple[int, int] = calculate_position(
        position_in_schema, max_feature_overlap
    )

    # Step 1: Read ORF file and get start and end positions for ORFs within genes
    orf_positions = {}
    with open(orf_file, "r") as of:
        next(of)  # Skip header
        for line in of:
            parts = line.strip().split("\t")
            if parts[4] != "0":  # Check if 'In_Gene' is not '0'
                start, end = int(parts[0]), int(parts[1])
                # Add ORF rules as a list for each position
                for pos in range(start, end + 1):
                    orf_rule = (
                        WriteRule.START
                        if pos == start
                        else WriteRule.END
                        if pos == end
                        else WriteRule.MIDDLE
                    )
                    orf_positions.setdefault(pos, []).append(orf_rule.value)

    # Step 2: Update the model file based on ORF positions
    temp_file_path = model_file + ".tmp"
    index = -1
    with open(model_file, "r") as mf, open(temp_file_path, "w") as temp_file:
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue

            vectors = [eval(vector) for vector in line.strip().split("\t")]
            updated_line = []
            for _, vector in enumerate(vectors):
                index += 1
                if index in orf_positions:
                    orf_types = orf_positions[index]
                    for offset in range(0, max_feature_overlap * 3, 3):
                        pos = position_to_write[0] + offset
                        if pos + 2 < len(vector) and vector[pos : pos + 3] == [0, 0, 0]:
                            rule_to_write = orf_types[
                                0
                            ]  # Write only the first matching rule
                            vector[pos : pos + 3] = rule_to_write
                            break  # Break after finding a space to write

                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")

    # Replace the old model file with the updated temporary file
    os.rename(temp_file_path, model_file)


def write_poly_adenyl_file(
    poly_adenyl_file: str,
    model_file: str,
    position_in_schema: int,
    max_feature_overlap: int,
) -> None:
    # Calculate positions to write polyadenylation information
    position_to_write: tuple[int, int] = calculate_position(
        position_in_schema, max_feature_overlap
    )

    # Step 1: Read polyadenylation file and get start and end positions
    poly_adenyl_positions = {}
    with open(poly_adenyl_file, "r") as of:
        next(of)  # Skip header
        for line in of:
            parts = line.strip().split("\t")
            start, end = int(parts[0]), int(parts[1])
            poly_adenyl_positions.setdefault(start, []).append(WriteRule.START.value)
            for pos in range(start + 1, end):
                poly_adenyl_positions.setdefault(pos, []).append(WriteRule.MIDDLE.value)
            poly_adenyl_positions.setdefault(end, []).append(WriteRule.END.value)

    # Step 2: Update the model file based on polyadenylation positions
    temp_file_path = model_file + ".tmp"
    index = -1
    with open(model_file, "r") as mf, open(temp_file_path, "w") as temp_file:
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue

            vectors = [eval(vector) for vector in line.strip().split("\t")]
            updated_line = []
            for _, vector in enumerate(vectors):
                index += 1
                if index in poly_adenyl_positions:
                    poly_types = poly_adenyl_positions[index]
                    for offset in range(0, max_feature_overlap * 3, 3):
                        pos = position_to_write[0] + offset
                        if pos + 2 < len(vector) and vector[pos : pos + 3] == [0, 0, 0]:
                            vector[pos : pos + 3] = poly_types[
                                0
                            ]  # Write only the first matching rule
                            break  # Break after finding a space to write
                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")
    os.rename(temp_file_path, model_file)


def transform_data_to_vectors(
    fasta_folder: str,
    gff3_folder: str,
    orf_folder: str,
    promotor_motifs_folder: str,
    poly_adenyl_folder: str,
    output_folder: str,
    max_feature_overlap: int,
    version: list[int],
) -> None:
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
    create_folder(output_folder)
    for filename in os.listdir(fasta_folder):
        ending = ".fasta"
        if not filename.endswith(ending):
            continue
        base_part, extension = os.path.splitext(filename)
        output_file: str = f"{output_folder}/{base_part}.txt"
        create_file_header(output_file, bp_vector_schema, max_feature_overlap, version)
        write_fasta_file(
            f"{fasta_folder}/{filename}",
            output_file,
            bp_vector_schema,
            max_feature_overlap,
        )
        write_promotor_motifs_file(
            f"{promotor_motifs_folder}/{base_part}.txt",
            output_file,
            bp_vector_schema.index("PROMOTOR_MOTIF"),
            max_feature_overlap,
        )
        write_gff3_file(
            f"{gff3_folder}/{base_part}.gff3",
            output_file,
            bp_vector_schema,
            max_feature_overlap,
        )
        write_orf_file(
            f"{orf_folder}/{base_part}.txt",
            output_file,
            bp_vector_schema.index("ORF"),
            max_feature_overlap,
        )
        write_poly_adenyl_file(
            f"{poly_adenyl_folder}/{base_part}.txt",
            output_file,
            bp_vector_schema.index("POLY_ADENYL"),
            max_feature_overlap,
        )
