# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- function `transform_data_to_vectors`
    - the logic were change from  os.replace(temp_file.name, model_file) to
os.rename(temp_file_path, model_file) => speed
