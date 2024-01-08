import os, tempfile
from enum import Enum
from datetime import datetime

from .common import create_folder

__all__ = ["transform_data_to_vectors"]


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


def create_file_header(path_to_file: str, bp_vector_schema: list[str]) -> None:
    header_: str = (
        "#HEADER#\n"
        f"#DATE: {datetime.utcnow().date()}\n"
        f"#bp_vector_schema: {bp_vector_schema}\n"
        "#description of feature: 0 = no present, 1 = start, 2 = continuation/ongoing, 3 = end\n"
        "####END####\n"
    )
    with open(path_to_file, "w") as f:
        f.write(header_)


def write_fasta_file(fasta_file: str, model_file: str, vector_base: list[int]):
    vec_len: int = len(vector_base)
    # Mapping nucleotides to vectors
    nucleotide_to_vector = {
        "A": [1] + [0] * (vec_len - 1),
        "C": [0] + [1] + [0] * (vec_len - 2),
        "G": [0] * 2 + [1] + [0] * (vec_len - 3),
        "T": [0] * 3 + [1] + [0] * (vec_len - 4),
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
                    if vector_count == 5:
                        # At the fifth vector, write it followed by a newline and reset the counter
                        model.write(f"{vector}\n")
                        vector_count = 0
                    else:
                        # Otherwise, write the vector followed by a tab
                        model.write(f"{vector}\t")


def write_promotor_motifs_file(
    promotor_file: str, model_file: str, position_to_write: int
) -> None:
    # Step 1: Read promoter file and get start and end positions
    promoter_positions = {}
    with open(promotor_file, "r") as pf:
        next(pf)  # Skip header
        for line in pf:
            start, end, other_parts = parse_line(line)
            promoter_positions.setdefault(start, []).append(WriteRule.START)
            for pos in range(start + 1, end):
                promoter_positions.setdefault(pos, []).append(WriteRule.MIDDLE)
            promoter_positions.setdefault(end, []).append(WriteRule.END)

    # Step 2: Read and update the model file line by line
    temp_file_path = model_file + ".tmp"
    with open(model_file, "r") as mf, open(temp_file_path, "w") as temp_file:
        vector_count = 0
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue

            vectors = [eval(vector) for vector in line.strip().split("\t")]
            updated_line = []
            for index, vector in enumerate(vectors):
                vector_count += 1
                if vector_count in promoter_positions:
                    # Combine all promoter types for this position
                    promoter_types = promoter_positions[vector_count]
                    vector[position_to_write] = [rule.value for rule in promoter_types]
                updated_line.append(str(vector))

            temp_file.write("\t".join(updated_line) + "\n")

    # Step 3: Replace the old model file with the updated temporary file
    os.rename(temp_file_path, model_file)


def write_gff3_file(
    gff3_file: str, model_file: str, bp_vector_schema: list[str]
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
                    )
                feature_positions[feature_type].setdefault(start, []).append(
                    WriteRule.START.value
                )
                feature_positions[feature_type].setdefault(end, []).append(
                    WriteRule.END.value
                )

    # Open the model file for reading and a temporary file for writing
    temp_file_path = model_file + ".tmp"
    with open(model_file, "r") as mf, open(temp_file_path, "w") as temp_file:
        vector_count = 0
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)
                continue

            vectors = [eval(vector) for vector in line.strip().split("\t")]
            updated_line = []
            for index, vector in enumerate(vectors):
                vector_count += 1
                for feature in feature_positions:
                    if vector_count in feature_positions[feature]:
                        feature_index = bp_vector_schema.index(feature)
                        # Combine all feature statuses for this position
                        feature_statuses = feature_positions[feature][vector_count]
                        vector[feature_index] = feature_statuses
                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")

    # Replace the old model file with the updated temporary file
    os.rename(temp_file_path, model_file)


def write_orf_file(orf_file: str, model_file: str, position_to_write: int) -> None:
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
    with open(model_file, "r") as mf, open(temp_file_path, "w") as temp_file:
        vector_count = 0  # Counter for the position of each vector
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue
            vectors = [eval(vector) for vector in line.strip().split("\t")]
            updated_line = []
            for index, vector in enumerate(vectors):
                vector_count += 1
                if vector_count in orf_positions:
                    # Assign a list of ORF rules to the position
                    vector[position_to_write] = orf_positions[vector_count]
                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")

    # Replace the old model file with the updated temporary file
    os.rename(temp_file_path, model_file)


def write_poly_adenyl_file(
    poly_adenyl_file: str, model_file: str, position_to_write: int
) -> None:
    # Step 1: Read polyadenylation file and get start and end positions
    poly_adenyl = {}
    with open(poly_adenyl_file, "r") as of:
        next(of)  # Skip header
        for line in of:
            parts = line.strip().split("\t")
            start, end = int(parts[0]), int(parts[1])
            poly_adenyl.setdefault(start, []).append(WriteRule.START.value)
            for pos in range(start + 1, end):
                poly_adenyl.setdefault(pos, []).append(WriteRule.MIDDLE.value)
            poly_adenyl.setdefault(end, []).append(WriteRule.END.value)

    # Step 2: Update the model file based on polyadenylation positions
    temp_file_path = model_file + ".tmp"
    with open(model_file, "r") as mf, open(temp_file_path, "w") as temp_file:
        vector_count = 0  # Counter for the position of each vector
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue
            vectors = [eval(vector) for vector in line.strip().split("\t")]
            updated_line = []
            for index, vector in enumerate(vectors):
                vector_count += 1
                if vector_count in poly_adenyl:
                    # Combine all poly-adenylation types for this position
                    poly_types = poly_adenyl[vector_count]
                    vector[position_to_write] = [rule for rule in poly_types]
                else:
                    vector[position_to_write] = 0
                updated_line.append(str(vector))
            temp_file.write("\t".join(updated_line) + "\n")

    # Replace the old model file with the updated temporary file
    os.rename(temp_file_path, model_file)


def transform_data_to_vectors(
    fasta_folder: str,
    gff3_folder: str,
    orf_folder: str,
    promotor_motifs_folder: str,
    poly_adenyl_folder: str,
    output_folder: str,
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
    bp_vector: list[int] = [0 for i in bp_vector_schema]
    create_folder(output_folder)
    for filename in os.listdir(fasta_folder):
        ending = ".fasta"
        if not filename.endswith(ending):
            continue
        base_part, extension = os.path.splitext(filename)
        output_file: str = f"{output_folder}/{base_part}.txt"
        create_file_header(output_file, bp_vector_schema)
        write_fasta_file(f"{fasta_folder}/{filename}", output_file, bp_vector)
        write_promotor_motifs_file(
            f"{promotor_motifs_folder}/{base_part}.txt",
            output_file,
            bp_vector_schema.index("PROMOTOR_MOTIF"),
        )
        write_gff3_file(
            f"{gff3_folder}/{base_part}.gff3", output_file, bp_vector_schema
        )
        write_orf_file(
            f"{orf_folder}/{base_part}.txt", output_file, bp_vector_schema.index("ORF")
        )
        write_poly_adenyl_file(
            f"{poly_adenyl_folder}/{base_part}.txt",
            output_file,
            bp_vector_schema.index("POLY_ADENYL"),
        )
