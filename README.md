
# Lumbridge

## Introduction
Lumbridge is a project focused on Plant DNA Annotation and Preparation for Neural Networks.
It aims to bridge the gap between annotations of different bioinformatics tools.

## Table of Contents
- [Project Goal](#project-goal)
- [Requirements and Dependencies](#requirements-and-dependencies)
- [Execution Methods](#execution-methods)
- [Installation and Setup](#installation-and-setup)
- [Usage Instructions](#usage-instructions)
- [Project Workflow](#project-workflow)
- [Structure of Model File](#structure-of-model-file)
- [License Information](#license-information)
- [Contact Information](#contact-information)
- [Acknowledgments](#acknowledgments)

## Project Goal

#### Objective: Plant DNA Annotation and Preparation for Modeling

The primary objective of this project is to meticulously annotate and prepare plant DNA sequences
for subsequent computational modeling.  This involves a comprehensive analysis and labeling of
various genetic elements within the DNA sequences of plants,  followed by a systematic preparation of
these annotated sequences to be utilized effectively in predictive modeling and simulation studies.

Key aspects of this process include:
- Utilizing gff3 and .cds files as primary inputs for foundational genomic information.
- Employing advanced bioinformatics tools to identify and annotate additional key genomic features that are
not explicitly detailed in the gff3 and .cds files.
- This includes the identification of promoter motifs, polyadenylation sequences, and other regulatory regions
that play a crucial role in gene expression and regulation.
- Integrating comprehensive data from various sources and formats to build a complete and accurate
representation of the plant DNA for modeling purposes.
- Ensuring the integrity and accuracy of the annotated data, setting a strong foundation for predictive modeling
and simulation studies.

This meticulous preparation is crucial for the success of modeling efforts, as it lays the
foundational groundwork for accurate and reliable biological simulations and analyses.



## Requirements and Dependencies

Before beginning the setup, ensure your system meets the following requirements:
- **Operating System:** Linux
- **Python Version:** Python 3.10 or higher
- **Additional Tools:** Homer2, Bedtools

### Execution Methods

The project can be executed using two primary methods:

1. **Docker**: Utilizes a Docker container to encapsulate the environment and dependencies, ensuring a consistent and isolated execution across different systems.

2. **Linux System**: Direct execution on a Linux-based system, where dependencies and environment configurations are managed locally.

## Installation and Setup

### Docker

1. **Build image:**
   ```bash
   docker build -t lumbridge-app .
   ```

1. **Run container:**
   ```bash
   docker run -p 4000:80 lumbridge-app
   ```


### Setting Up Python Virtual Environment
A Python virtual environment is recommended for managing the project's dependencies.
Follow these steps to create and activate your virtual environment:

1. **Create a Virtual Environment:**
   ```bash
   python3.10 -m venv venv
   ```

2. **Activate the Virtual Environment:**
   ```bash
   source venv/bin/activate
   ```

3. **Install Required Python Packages:**
   ```bash
   pip install -Ur requirements.txt
   ```


### Homer2 installation
Homer2 is a key tool required for this project.
Follow the detailed installation guide to set up Homer2 on your system:

1. **Download Installation Script:**
   ```bash
   wget http://homer.ucsd.edu/homer/configureHomer.pl
   ```

2. **Run Installation Script:**
   ```bash
   perl configureHomer.pl -install
   ```

3. **Add to PATH in ~/.bash_profile:**
   ```text
   # Add the following line to your ~/.bash_profile
   PATH=$PATH:/home/matej/homer2/bin
   # Verify installation
   which findMotifs.pl # Should output path in this case `/home/matej/homer2/bin/findMotifs.pl`
   ```

4. **Update Configuration File:**
   ```text
   # Update HOMER2_BIN_PATH in the configuration file
   HOMER2_BIN_PATH="/home/matej/homer2/bin"
   ```

Refer to the [Homer2 Installation Guide](http://homer.ucsd.edu/homer/introduction/install.html) for more details.

### Install Bedtools
1. **Run Installation Script:**
   ```bash
   sudo apt-get install bedtools
   ```

## Usage Instructions

#### Configuration Options:
To customize the application's behavior, you have two options for configuration:

1. **Create a Custom Configuration File:**
   You can create your own `config.py` file in the `/etc/lumbridge/config.py` directory.
This custom configuration file allows you to specify settings tailored to your needs.

2. **Modify the Local Configuration:**
   Alternatively, you can directly modify the local `config.py` file in the project's root directory.
This is a straightforward way to change settings if you don't require a separate configuration file.

#### Default Data:
By default, the application uses test data from *Arabidopsis Thaliana*. This dataset serves as a standard reference
for initial runs and testing purposes. You can replace it with your specific data in the configuration settings.


## Project Workflow

The workflow is predicated on the following assumptions and requirements for the input files:

1. **Fasta Files:**
   - The fasta files are assumed to be from the forward (positive) strand.
   - They contain the nucleotide sequences that will be used for further analysis and annotation.

2. **GFF3 and ORF Files:**
   - These files may contain annotations from both the forward and reverse strands of DNA.
   - However, for the purpose of this project, only annotations from the forward strand are considered.
   - The GFF3 (General Feature Format version 3) files provide detailed annotations about genomic features.
   - ORF (Open Reading Frame) files detail regions of the genome that potentially code for proteins.

3. **File Naming Conventions:**
   - It's crucial that the files which are related (fasta, gff3, and orf) have the same base name
in their filenames to be recognized as belonging together.
   - For example, files for one chromosome of *Arabidopsis thaliana* should be named consistently
like `arabidopsis_chr1.fasta`, `arabidopsis_chr1.gff3`, and `arabidopsis_chr1.cds`.

This structured approach ensures that the data is correctly aligned and processed in subsequent stages of the project, allowing for accurate analysis and modeling based on the genomic information provided in these files.

```plaintext
fasta_folder    gff3_folder    orf_folder
      |             |             |
      v             v             v
    ---------------------------------
                 Lumbridge
    ---------------------------------
      |                             |
      v                             v
  model_data                     others
```


### Structure of Model File

The data file is structured into two distinct sections: the header and the data content.

#### Header
The header section encapsulates critical metadata about the file, including the date of creation and the schema of the base pair (bp) vector. Each element in the bp vector represents a specific genomic feature, encoded as follows:

```plaintext
#HEADER#
#DATE: 2024-01-07
#bp_vector_schema: ['A', 'C', 'G', 'T', 'PROMOTOR_MOTIF', 'ORF', 'exon', 'mRNA', 'miRNA', 'rRNA', 'CDS', 'POLY_ADENYL', 'gene']
#description of feature: 0 = not present, 1 = start, 2 = continuation/ongoing, 3 = end
####END####
```

#### Data Content
In the data section, each base pair (bp) from the FASTA file is represented by a vector.
While the first four positions of the vector represent the nucleotides (A, C, G, T),
subsequent positions may contain lists of integers, indicating the occurrence of multiple genomic features
at the same position. For example:

- [1, 0, 0, 0, 0, 0, [2], [2, 1], 0, 0, 0, 0, [1, 1]]

In this representation, [2] signifies the ongoing presence of a feature, whereas [2, 1]
indicates overlapping features at the same genomic location.



## License Information
MIT

## Contact Information
forgac.matej@gmail.com


## Acknowledgments
In this project were used additional programs: [Homer2 ](http://homer.ucsd.edu/homer/) and
[Bedtools2](https://github.com/arq5x/bedtools2). Thanks all people for contributions to these projects.

### Citation
Please cite the following article if you use BEDTools in your research:

Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.
#### Also, if you use pybedtools, please cite the following.

Dale RK, Pedersen BS, and Quinlan AR. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics (2011). doi:10.1093/bioinformatics/btr539
