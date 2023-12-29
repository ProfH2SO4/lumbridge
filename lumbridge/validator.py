import os

from .exc import WrongInputPath, WrongFile, MissingFile

__all__ = ["check_input_structure", "check_strands_in_orf_file"]


def find_missing_file(
    file_name_fasta_no_ending: list[str],
    file_name_gff3_no_ending: list[str],
    file_name_orf_no_ending: list[str],
) -> dict[str, set[str]]:
    # Convert lists to sets for faster lookup
    set_fasta = set(file_name_fasta_no_ending)
    set_gff3 = set(file_name_gff3_no_ending)
    set_orf = set(file_name_orf_no_ending)

    # Find missing elements
    missing_in_fasta = (set_gff3.union(set_orf)) - set_fasta
    missing_in_gff3 = (set_fasta.union(set_orf)) - set_gff3
    missing_in_orf = (set_fasta.union(set_gff3)) - set_orf

    # Return the result as a dictionary
    return {
        "missing_in_fasta": missing_in_fasta,
        "missing_in_gff3": missing_in_gff3,
        "missing_in_orf": missing_in_orf,
    }


def check_input_structure(fasta_folder: str, gff3_folder: str, orf_folder: str) -> None:
    """
    Check content of fasta_folder, gff3_folder, orf_folder.
    All base names (arabidopsis_thaliana_ch_1) in folders must be same.
    e.g
    In fasta_folder arabidopsis_thaliana_ch_1.fasta
    In gff3_folder arabidopsis_thaliana_ch_1.gff3
    In orf_folder arabidopsis_thaliana_ch_1.cds

    :fasta_folder: A path to fasta_folder
    :gff3_folder: A path to gff3_folder
    :orf_folder: A path to orf_folder
    :return: None If ok else raise Error
    """
    if not (os.path.exists(fasta_folder) and os.path.isdir(fasta_folder)):
        raise WrongInputPath(fasta_folder)
    if not (os.path.exists(gff3_folder) and os.path.isdir(gff3_folder)):
        raise WrongInputPath(gff3_folder)
    if not (os.path.exists(orf_folder) and os.path.isdir(orf_folder)):
        raise WrongInputPath(orf_folder)

    file_name_fasta_no_ending: list[str] = []
    file_name_gff3_no_ending: list[str] = []
    file_name_orf_no_ending: list[str] = []

    # Enter fasta Folder
    for filename in os.listdir(fasta_folder):
        ending = ".fasta"
        fai_ending = ".fai"
        if not os.path.isfile(os.path.join(fasta_folder, filename)):
            continue
        if not (filename.endswith(ending) or filename.endswith(fai_ending)):
            raise WrongFile(
                path=os.path.join(fasta_folder, filename), expected_ending=ending
            )
        if filename.endswith(ending):
            file_name_fasta_no_ending.append(filename[: -len(ending)])

    for filename in os.listdir(gff3_folder):
        ending = ".gff3"
        if not os.path.isfile(os.path.join(gff3_folder, filename)):
            continue
        if not filename.endswith(ending):
            raise WrongFile(
                path=os.path.join(gff3_folder, filename), expected_ending=ending
            )
        file_name_gff3_no_ending.append(filename[: -len(ending)])

    for filename in os.listdir(orf_folder):
        ending = ".cds"
        if not os.path.isfile(os.path.join(orf_folder, filename)):
            continue
        if not filename.endswith(ending):
            raise WrongFile(
                path=os.path.join(orf_folder, filename), expected_ending=ending
            )
        file_name_orf_no_ending.append(filename[: -len(ending)])

    if not (
        len(file_name_fasta_no_ending)
        == len(file_name_gff3_no_ending)
        == len(file_name_orf_no_ending)
    ):
        ret_missing: dict[str, set[str]] = find_missing_file(
            file_name_fasta_no_ending, file_name_gff3_no_ending, file_name_orf_no_ending
        )

        raise MissingFile(missing_files=ret_missing)


def check_strands_in_orf_file(file_path) -> bool:
    """Check if ORF file contains ORFs from both strands."""
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith(">"):
                if ":c" in line:
                    return True  # Found an ORF from the reverse strand
    return False  # No ORF from the reverse strand found
