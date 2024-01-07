from types import ModuleType
from os.path import isfile
import os

from .common import (
    create_output_orf_forward_strand,
    find_gff3_and_orf_intervals,
    find_poly_adi_sequences,
    extract_positive_elements_gff3,
)
from .validator import check_input_structure, check_strands_in_orf_file
from .homer2 import annotate_homer2_motifs, make_homer2_output
from .model_prep import transform_data_to_vectors

import config
import log

__version__ = "0.0.1"
__last_update__ = "2023-10-25T20:00:00"


def compile_config(app_config: ModuleType, path: str):
    try:
        with open(path, "rb") as rnf:
            exec(compile(rnf.read(), "config.py", "exec"), app_config.__dict__)
    except OSError as e:
        print(f"File at {path} could not be loaded because of error: {e}")
        raise e from e


def load_config(test_config_path: None | str) -> ModuleType:
    """
    Load local config.py.
    If exists config.py in /etc/lumbridge/ then overrides parameters in local config.py.
    @return: configuration file
    """
    app_config: ModuleType = config
    path: str = "/etc/lumbridge/config.py"

    if not isfile(path) and not test_config_path:
        return app_config
    if isfile(path):
        compile_config(app_config, path)
    if test_config_path:
        compile_config(app_config, test_config_path)
    return app_config


def parse_namespace(config_: ModuleType) -> dict[str, any]:
    """
    Parse configuration file file to dict.
    @param config_: configuration file
    @return: parsed configuration file
    """
    parsed: dict[str, any] = {}
    for key, value in config_.__dict__.items():
        if not key.startswith("__"):
            parsed[key] = value
    return parsed


def create_file_if_not_exists(path_to_file: str) -> None:
    # Check if the directory exists, and create it if it doesn't
    directory = os.path.dirname(path_to_file)
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Check if the file exists, and create it if it doesn't
    if not os.path.exists(path_to_file):
        with open(path_to_file, "w") as file:
            pass  # Create an empty file


def run(test_config: str | None = None) -> None:
    # Load Config
    config_: ModuleType = load_config(test_config)
    parsed_config: dict[str, any] = parse_namespace(config_)

    print("============ Setting Up Logger ============")
    if parsed_config["LOG_CONFIG"]["handlers"].get("file", None):
        file_path: str = parsed_config["LOG_CONFIG"]["handlers"]["file"].get("filename")
        create_file_if_not_exists(file_path)
    log.set_up_logger(parsed_config["LOG_CONFIG"])

    log.info("------ Check if all files provided -------")
    check_input_structure(
        parsed_config["INPUT_FASTA"],
        parsed_config["INPUT_GFF3"],
        parsed_config["INPUT_ORF"],
    )

    log.info("------ Extract gff3 positive elements  -------")
    gff3_folder_pos_strand: str = (
        f"{parsed_config['OUTPUT_FOLDER']}/gff3_positive_strand"
    )
    extract_positive_elements_gff3(parsed_config["INPUT_GFF3"], gff3_folder_pos_strand)

    log.info("------ Create file with ORF only on forward strand  -------")
    orf_output_folder: str = (
        f"{parsed_config['OUTPUT_FOLDER']}/orf_folder_positive_strand"
    )
    create_output_orf_forward_strand(parsed_config["INPUT_ORF"], orf_output_folder)

    log.info("------ Find ORf overlaps in GFF3 file  -------")
    orf_in_gff3_folder: str = f"{parsed_config['OUTPUT_FOLDER']}/orf_in_gff3"
    find_gff3_and_orf_intervals(
        orf_output_folder, gff3_folder_pos_strand, orf_in_gff3_folder
    )

    log.info("------ Find Polyadenylation sequences in Fasta file  -------")
    find_poly_adi_sequences(
        parsed_config["INPUT_FASTA"], parsed_config["OUTPUT_FOLDER"]
    )

    log.info("------ Make homer2 output -------")
    make_homer2_output(
        parsed_config["HOMER2_BIN_PATH"],
        parsed_config["INPUT_FASTA"],
        gff3_folder_pos_strand,
        parsed_config["HOMER2_OUTPUT_FOLDER"],
        parsed_config["HOMER2_UPSTREAM_GEN_SEQ_LENGTH"],
        parsed_config["HOMER2_DOWNSTREAM_GEN_SEQ_LENGTH"],
        cpu_cores=parsed_config["HOMER2_CPU_CORES"],
    )

    log.info("------ Annotate homer2 output to fasta -------")
    annotate_homer2_motifs(
        parsed_config["INPUT_FASTA"],
        gff3_folder_pos_strand,
        parsed_config["HOMER2_OUTPUT_FOLDER"],
        parsed_config["HOMER2_UPSTREAM_GEN_SEQ_LENGTH"],
        parsed_config["HOMER2_DOWNSTREAM_GEN_SEQ_LENGTH"],
        parsed_config["HOMER2_P_THRESHOLD"],
        output_folder=parsed_config["OUTPUT_FOLDER"],
    )

    log.info("------ Prepare data for model -------")
    transform_data_to_vectors(
        fasta_folder=parsed_config["INPUT_FASTA"],
        gff3_folder=gff3_folder_pos_strand,
        orf_folder=orf_in_gff3_folder,
        promotor_motifs_folder=f"{parsed_config['OUTPUT_FOLDER']}/homer2_annotation",
        poly_adenyl_folder=f"{parsed_config['OUTPUT_FOLDER']}/poly_adenylation",
        output_folder=f"{parsed_config['OUTPUT_FOLDER']}/model_data",
    )

    log.info("------ Done  -------")
    print("------ Done  -------")
