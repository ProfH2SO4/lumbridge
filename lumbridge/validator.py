import os

from .exc import WrongInputPath, WrongFile, MissingFile

__all__ = ["check_input_structure"]


def find_missing_file(file_name_fasta_no_ending: list[str],
                      file_name_gff3_no_ending: list[str],
                      file_name_orf_no_ending: list[str]) -> dict[str, set[str]]:
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


def check_input_structure(fasta_path: str, gff3_path: str, orf_path: str) -> None:
    """
    Check content of fasta_folder, gff3_folder, orf_folder
    :return:
    """
    if not (os.path.exists(fasta_path) and os.path.isdir(fasta_path)):
        raise WrongInputPath(fasta_path)
    if not (os.path.exists(gff3_path) and os.path.isdir(gff3_path)):
        raise WrongInputPath(gff3_path)
    if not (os.path.exists(orf_path) and os.path.isdir(orf_path)):
        raise WrongInputPath(orf_path)

    file_name_fasta_no_ending: list[str] = []
    file_name_gff3_no_ending: list[str] = []
    file_name_orf_no_ending: list[str] = []

    #Enter fasta Folder
    for filename in os.listdir(fasta_path):
        ending = ".fasta"
        if not os.path.isfile(os.path.join(fasta_path, filename)):
            continue
        if not filename.endswith(ending):
            raise WrongFile(path=os.path.join(fasta_path, filename), expected_ending=ending)
        file_name_fasta_no_ending.append(filename[:-len(ending)])

    for filename in os.listdir(gff3_path):
        ending = ".gff3"
        if not os.path.isfile(os.path.join(gff3_path, filename)):
            continue
        if not filename.endswith(ending):
            raise WrongFile(path=os.path.join(gff3_path, filename), expected_ending=ending)
        file_name_gff3_no_ending.append(filename[:-len(ending)])

    for filename in os.listdir(orf_path):
        ending = ".cds"
        if not os.path.isfile(os.path.join(orf_path, filename)):
            continue
        if not filename.endswith(ending):
            raise WrongFile(path=os.path.join(orf_path, filename), expected_ending=ending)
        file_name_orf_no_ending.append(filename[:-len(ending)])

    if not (len(file_name_fasta_no_ending) == len(file_name_gff3_no_ending) == len(file_name_orf_no_ending)):
        ret_missing: dict[str, set[str]] = find_missing_file(file_name_fasta_no_ending,
                                                     file_name_gff3_no_ending,
                                                     file_name_orf_no_ending)

        raise MissingFile(missing_files=ret_missing)

def check_orf_presents(fasta_path: str, orf_path: str):
    pass



