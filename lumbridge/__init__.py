from types import ModuleType
from os.path import isfile

from .common import create_output_orf_forward_strand, find_gff3_and_orf_intervals, find_poly_adi_sequences
from .validator import check_input_structure, check_strands_in_orf_file
from .homer2 import annotate_homer2_motifs, make_homer2_output
import config

__version__ = "0.0.1"
__started_at__ = "2023-10-25T20:00:00"


def load_config() -> ModuleType:
    """
    Load local config.py.
    If exists config.py in /etc/lumbridge/ then overrides parameters in local config.py.
    """
    app_config: ModuleType = config
    path: str = "/etc/lumbridge/config.py"

    if not isfile(path):
        return app_config
    try:
        with open(path, "rb") as rnf:
            exec(compile(rnf.read(), "config.py", "exec"), app_config.__dict__)
    except OSError as e:
        print(f"File at {path} could not be loaded because of error: {e}")
        raise e from e
    return app_config


def parse_namespace(config: ModuleType) -> dict[str, any]:
    """
    Parse namespace to list.
    """
    parsed: dict[str, any] = {}
    for key, value in config.__dict__.items():
        if not key.startswith('__'):
            parsed[key] = value
    return parsed


def run():
    # Load Config
    config_: ModuleType = load_config()
    parsed_config: dict[str, any] = parse_namespace(config_)
    print("------ Check if all files provided -------")
    check_input_structure(parsed_config["INPUT_FASTA"], parsed_config["INPUT_GFF3"], parsed_config["INPUT_ORF"])

    print("------ Create file with ORF only on forward strand  -------")
    create_output_orf_forward_strand(parsed_config["INPUT_ORF"], parsed_config["OUTPUT_FOLDER"])
    print("------ Find ORf overlaps in GFF3 file  -------")
    find_gff3_and_orf_intervals(f"{parsed_config['OUTPUT_FOLDER']}/orf_folder_forward_strand",
                                parsed_config["INPUT_GFF3"],
                                f"{parsed_config['OUTPUT_FOLDER']}")
    print("------ Find Polyadenylation sequences in Fasta file  -------")
    find_poly_adi_sequences(parsed_config["INPUT_FASTA"], parsed_config['OUTPUT_FOLDER'])
    print("------ Make homer2 output -------")
    make_homer2_output(parsed_config["INPUT_FASTA"],
                       parsed_config["INPUT_GFF3"],
                       parsed_config["HOMER2_OUTPUT_FOLDER"],
                       parsed_config["HOMER2_SEQUENCE_LEN_BEFORE_GENE"],
                       )
    # print("------ Annotate homer2 output to fasta -------")
    # annotate_homer2_motifs(parsed_config["INPUT_FASTA"],
    #                        parsed_config["INPUT_GFF3"],
    #                        parsed_config["HOMER2_OUTPUT_FOLDER"],
    #                        parsed_config["HOMER2_SEQUENCE_LEN_BEFORE_GENE"],
    #                        parsed_config["HOMER2_P_THRESHOLD"],
    #                        output_folder=parsed_config['OUTPUT_FOLDER'])
    print("------ Done  -------")
