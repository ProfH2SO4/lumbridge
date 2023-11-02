from types import ModuleType
from os.path import isfile

from .validator import check_input_structure
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
    check_input_structure(parsed_config["INPUT_FASTA"], parsed_config["INPUT_GFF3"], parsed_config["INPUT_ORF"])
    print("------ All files provided -------")








