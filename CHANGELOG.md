# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- function `transform_data_to_vectors`
    - the logic were change from  os.replace(temp_file.name, model_file) to
os.rename(temp_file_path, model_file) => speed

- In features is better to have [0 or 1 or 2 or 3] compared to
one hot encoded [0, 0, 0]. Pros: The file is more readable.
The features are in Model embedded.

- In case for `gene` prediction there cannot be all features from gff3
features such as `mRNA` does not make sense cuz where `gene` there is `mRNA`
also feature `exon` => The user specifies in `config.py` which features want
from `.gff3` file. e.g `GFF3_FEATURES = ["miRNA", "rRNA", "gene"]`
