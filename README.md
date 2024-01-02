# lumbridge


## Project Goal

### Objective: Plant DNA Annotation and Preparation for Modeling

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



## Project Setup and Requirements

This section of document outlines the setup and installation requirements for the project.
Ensure that you have the necessary environment and dependencies installed on your system.


### System Requirements
Before beginning the setup, ensure your system meets the following requirements:
- **Operating System:** Linux
- **Python Version:** Python 3.10 or higher
- **Additional Tools:** Homer2

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

- [Installation Guide](http://homer.ucsd.edu/homer/introduction/install.html)


### Running the Application

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
