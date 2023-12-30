import os, tempfile
from datetime import datetime

from .common import create_folder

__all__ = ["transform_data_to_vectors"]


def create_file_header(path_to_file: str, bp_vector_schema: list[str]) -> None:
    header_: str = (
        "#HEADER#\n"
        f"#DATE: {datetime.utcnow().date()}\n"
        f"#bp_vector_schema: {bp_vector_schema}\n"
        "#0 = no present, 1 = present\n"
        "#last element gene has own description: 0 = no present, 1 = gene_start, 2 = present, 3 = gene_end\n"
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
                    # Write the vector to the model file
                    model.write(f"{vector},")
                    vector_count += 1

                    # Check if 5 vectors have been written, then start a new line
                    if vector_count == 5:
                        model.write("\n")  # Start a new line
                        vector_count = 0  # Reset the counter


def write_promotor_motifs_file(
    promotor_file: str, model_file: str, position_to_write: int
) -> None:
    # Step 1: Read promoter file and get start and end positions
    promoter_positions = set()
    with open(promotor_file, "r") as pf:
        next(pf)  # Skip header
        for line in pf:
            start, end, _, _ = line.strip().split("\t")
            promoter_positions.update(range(int(start), int(end) + 1))

    # Step 2: Update the model file based on ORF positions
    with open(model_file, "r") as mf, tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as temp_file:
        vector_count = 0  # Counter for the position of each vector
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue
            vectors = eval(line.strip())
            updated_line = []
            for index, vector in enumerate(vectors):
                vector_count += 1
                if vector_count in promoter_positions:
                    vector[position_to_write] = 1
                updated_line.append(str(vector))
            temp_file.write(",".join(updated_line) + "\n")

    # Step 3: Replace the old model file with the updated temporary file
    os.replace(temp_file.name, model_file)


def write_gff3_file(gff3_file: str, model_file: str) -> None:
    bp_vector_schema: list[str] = [
        "A",
        "C",
        "G",
        "T",
        "PROMOTOR_MOTIF",
        "exon",
        "mRNA",
        "miRNA",
        "rRNA",
        "CDS",
        "POLY_ADENYL",
        "gene",
    ]

    # Initialize dictionaries to store positions for each feature
    feature_positions = {
        feature: {}
        for feature in bp_vector_schema
        if feature not in ["A", "C", "G", "T", "PROMOTOR_MOTIF", "POLY_ADENYL"]
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
                for pos in range(start, end + 1):
                    feature_positions[feature_type][pos] = 2  # Mark as within feature
                feature_positions[feature_type][start] = 1  # Mark start of feature
                feature_positions[feature_type][end] = 3  # Mark end of feature

        # Open the model file for reading and a temporary file for writing
        with open(model_file, "r") as mf, tempfile.NamedTemporaryFile(
            mode="w", delete=False
        ) as temp_file:
            vector_count = 0
            for line in mf:
                if line.startswith("#"):
                    temp_file.write(line)  # Write header lines as is
                    continue

                vectors = eval(line.strip())
                updated_line = []
                for index, vector in enumerate(vectors):
                    vector_count += 1
                    for feature in feature_positions:
                        if vector_count in feature_positions[feature]:
                            feature_index = bp_vector_schema.index(feature)
                            feature_status = feature_positions[feature][vector_count]
                            vector[feature_index] = feature_status
                    updated_line.append(str(vector))
                temp_file.write(",".join(updated_line) + "\n")

        # Replace the old model file with the updated temporary file
        os.replace(temp_file.name, model_file)


def write_orf_file(orf_file: str, model_file: str, position_to_write: int) -> None:
    # Step 1: Read ORF file and get start and end positions for ORFs within genes
    orf_positions = set()
    with open(orf_file, "r") as of:
        next(of)  # Skip header
        for line in of:
            parts = line.strip().split("\t")
            if parts[4] != "0":  # Check if 'In_Gene' is not '0'
                start, end = int(parts[0]), int(parts[1])
                orf_positions.update(range(start, end + 1))

    # Step 2: Update the model file based on ORF positions
    with open(model_file, "r") as mf, tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as temp_file:
        vector_count = 0  # Counter for the position of each vector
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue
            vectors = eval(line.strip())
            updated_line = []
            for index, vector in enumerate(vectors):
                vector_count += 1
                if vector_count in orf_positions:
                    vector[position_to_write] = 1
                updated_line.append(str(vector))
            temp_file.write(",".join(updated_line) + "\n")

    # Step 3: Replace the old model file with the updated temporary file
    os.replace(temp_file.name, model_file)


def write_poly_adenyl_file(
    poly_adenyl_file: str, model_file: str, position_to_write: int
) -> None:
    poly_adenyl = set()
    with open(poly_adenyl_file, "r") as of:
        next(of)  # Skip header
        for line in of:
            parts = line.strip().split("\t")
            start, end = int(parts[0]), int(parts[1])
            poly_adenyl.update(range(start, end + 1))

    with open(model_file, "r") as mf, tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as temp_file:
        vector_count = 0  # Counter for the position of each vector
        for line in mf:
            if line.startswith("#"):
                temp_file.write(line)  # Write header lines as is
                continue
            vectors = eval(line.strip())
            updated_line = []
            for index, vector in enumerate(vectors):
                vector_count += 1
                if vector_count in poly_adenyl:
                    vector[position_to_write] = 1
                updated_line.append(str(vector))
            temp_file.write(",".join(updated_line) + "\n")

    # Step 3: Replace the old model file with the updated temporary file
    os.replace(temp_file.name, model_file)


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
            f"{promotor_motifs_folder}/{base_part}.txt", output_file, 4
        )
        write_gff3_file(f"{gff3_folder}/{base_part}.gff3", output_file)
        write_orf_file(f"{orf_folder}/{base_part}.txt", output_file, 5)
        write_poly_adenyl_file(f"{poly_adenyl_folder}/{base_part}.txt", output_file, 11)
