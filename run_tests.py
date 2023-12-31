import unittest


from tests import model_data
from lumbridge import run

if __name__ == "__main__":
    # run("./tests/test_config.py")
    # Create a test loader
    suite = unittest.TestSuite()

    # Discover and load all tests in the 'tests' directory
    suite.addTests(unittest.TestLoader().loadTestsFromModule(model_data))

    # Create a test runner that will display the results to the console
    runner = unittest.TextTestRunner(verbosity=2)

    # Run the tests
    runner.run(suite)
