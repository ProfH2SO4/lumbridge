import unittest
import os

from tests import test_config
from tests.test_model_data import test_model_data_factory


def check_paths_exist(paths: list[str]) -> None:
    """
    Checks if the given paths exist. Prints a message for each path indicating whether it exists.
    """
    for path in paths:
        if not os.path.exists(path) and os.path.isdir(path):
            raise f"Tests cannot be started cuz {path} is not dir"


if __name__ == "__main__":
    suite = unittest.TestSuite()

    input_dir_path_fasta = test_config.input_dir_path_fasta
    input_dir_path_gff3 = test_config.input_dir_path_gff3
    input_dir_path_orf = test_config.input_dir_path_orf
    output_dir_path = test_config.output_dir_path

    check_paths_exist([input_dir_path_fasta, input_dir_path_gff3, input_dir_path_orf, output_dir_path])

    test_class = test_model_data_factory(input_dir_path_fasta, input_dir_path_gff3, input_dir_path_orf, output_dir_path)
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
