import unittest


from tests.test_model_data import test_model_data_factory


if __name__ == "__main__":
    suite = unittest.TestSuite()

    input_dir_path_fasta = "/home/matej/git/lumbridge/test_data/fasta_folder"
    input_dir_path_gff3 = "/home/matej/git/lumbridge/test_data/gff3_folder"
    input_dir_path_orf = "/home/matej/git/lumbridge/test_data/orf_folder"
    output_dir_path = "/home/matej/git/lumbridge/lumbridge_output"

    # Create an instance of the test case class with specific parameters
    test_class = test_model_data_factory(input_dir_path_fasta, input_dir_path_gff3, input_dir_path_orf, output_dir_path)

    # Add tests from the customized test case class
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
