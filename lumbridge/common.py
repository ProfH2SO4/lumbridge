import os


def create_folder(new_folder_path: str) -> None:
    if not os.path.isdir(new_folder_path):
        try:
            os.mkdir(new_folder_path)
            print(f"Directory {new_folder_path} created successfully.")
        except OSError as error:
            print(f"Creation of the directory {new_folder_path} failed due to: {error}")
            raise error

# Function to parse ORF data
def parse_orf_file(orf_file):
    orfs = []
    with open(orf_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                parts = line.split()
                location_info = parts[0]
                start, end = [int(x) for x in location_info.split(':')[1].split('-')]
                orfs.append((start, end))
    return orfs


# Function to parse GFF3 data
def parse_gff3_file(gff3_file):
    genes = []
    with open(gff3_file, 'r') as file:
        for line in file:
            if not line.startswith('#') and '\tgene\t' in line:
                parts = line.split('\t')
                start = int(parts[3])
                end = int(parts[4])
                genes.append((start, end))
    return genes


def get_overlap_type(orf_start: int, orf_end: int, gene_start: int, gene_end: int) -> int:
    """
    None = 0
    Full = 1
    Partial = 2
    :return:
    """
    if orf_start >= gene_start and orf_end <= gene_end:
        return 1
    if (orf_start < gene_end and orf_end > gene_start) or \
       (orf_end > gene_start and orf_start < gene_end):
        return 2
    return 0


def is_orf_in_gene(orf, genes) -> int | tuple[int, int, int]:
    for gene in genes:
        overlap_type: int = get_overlap_type(orf[0], orf[1], gene[0], gene[1])
        if overlap_type != 0:
            return (gene[0], gene[1], overlap_type)
    return 0


# Main function to write intervals to a file
def write_orf_gff3_intervals(orf_file: str, gff3_file: str, output_file: str) -> None:
    orfs: list[tuple[int, int]] = parse_orf_file(orf_file)
    genes: list[tuple[int, int]] = parse_gff3_file(gff3_file)

    with open(output_file, 'w') as file:
        file.write("ORF_Start\tORF_End\tGFF3_Start\tGFF3_End\tIn_Gene\n")
        for orf in orfs:
            in_gene: int | tuple = is_orf_in_gene(orf, genes)
            if isinstance(in_gene, int):
                file.write(f"{orf[0]}\t{orf[1]}\t-\t-\t{in_gene}\n")
            else:
                file.write(f"{orf[0]}\t{orf[1]}\t{in_gene[0]}\t{in_gene[1]}\t{in_gene[2]}\n")


def find_gff3_and_orf_intervals(path_orf_folder: str, path_gff3_folder: str, path_output_dir: str) -> None:
    """ Find orf in gff3's genes. """
    for filename in os.listdir(path_gff3_folder):
        ending_gff3: str = ".gff3"
        ending_cds: str = ".cds"
        filename_body: str = filename[:-len(ending_gff3)]
        filename_orf: str = f"{path_orf_folder}/{filename_body}{ending_cds}"
        file_name_output: str = f"{filename_body}_orf_in_gff3.txt"
        write_orf_gff3_intervals(filename_orf,
                                 f"{path_gff3_folder}/{filename}",
                                 f"{path_output_dir}/{file_name_output}")


def create_one_strand_orf_file(orf_path: str, output_dir: str) -> None:

    # Determine the new file name
    dir_name, file_name = os.path.split(orf_path)
    name_part, extension = os.path.splitext(file_name)
    new_file_name = f"{name_part}{extension}"
    new_file_path = os.path.join(output_dir, new_file_name)

    # Open the original file and the new file for writing
    with open(orf_path, 'r') as original_file, open(new_file_path, 'w') as new_file:
        # Initialize a variable to keep track of whether to write lines
        write_lines = False

        for line in original_file:
            # Check if the line is a header
            if line.startswith('>'):
                # Determine if it's a forward strand ORF
                write_lines = ':c' not in line
            # Write lines that belong to the forward strand ORF
            if write_lines:
                new_file.write(line)


def create_output_orf_forward_strand(input_orf_folder_path: str, output_folder_path: str) -> None:
    orf_output_folder: str = "orf_folder_forward_strand"
    new_folder_path: str = f"{output_folder_path}/{orf_output_folder}"
    if not os.path.isdir(new_folder_path):
        try:
            os.mkdir(new_folder_path)
            print(f"Directory {new_folder_path} created successfully.")
        except OSError as error:
            print(f"Creation of the directory {new_folder_path} failed due to: {error}")
            raise error

    for filename in os.listdir(input_orf_folder_path):
        ending = ".cds"
        if filename.endswith(ending):
            create_one_strand_orf_file(f"{input_orf_folder_path}/{filename}", new_folder_path)


def write_poly_adi_seq(fasta_file_path: str, poly_adil_path: str) -> None:
    pol_adi_seq = ["ATTAAA", "TATAAA", "AGTAAA", "AAGAAA", "ACTAAA", "ATATAA", "ATTACA", "ATTGAA", "ATTATT"]

    with open(fasta_file_path, 'r') as file, open(poly_adil_path, 'w') as out_file:
        data = file.read()
        out_file.write(f"start\tend\tseq\n")
        for seq in pol_adi_seq:
            pos = -1
            while True:
                pos = data.find(seq, pos + 1)
                if pos == -1:
                    break
                out_file.write(f"{pos}\t{pos+len(seq)}\t{seq}\t\n")


def find_poly_adi_sequences(fasta_folder_path: str, output_folder_path: str) -> None:
    """ Find Polyadenylation sequences """
    poly_adi_output_folder_path: str = f"{output_folder_path}/termination"
    create_folder(poly_adi_output_folder_path)
    for filename in os.listdir(fasta_folder_path):
        ending_fasta: str = ".fasta"
        filename_body: str = filename[:-len(ending_fasta)]
        poly_adil_file_path = f"{poly_adi_output_folder_path}/{filename_body}.txt"
        write_poly_adi_seq(f"{fasta_folder_path}/{filename}", poly_adil_file_path)
