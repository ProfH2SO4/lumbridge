import unittest

from lumbridge import run

if __name__ == "__main__":
    # run("./tests/test_config.py")
    # Create a test loader
    loader = unittest.TestLoader()

    # Discover and load all tests in the 'tests' directory
    suite = loader.discover(start_dir="./tests", pattern="test_*.py")

    # Create a test runner that will display the results to the console
    runner = unittest.TextTestRunner()

    # Run the tests
    runner.run(suite)
