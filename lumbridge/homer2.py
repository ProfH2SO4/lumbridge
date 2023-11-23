import os
from Bio import SeqIO

from .common import parse_gff3_file, create_folder, extract_promoters_to_bed, get_fasta_from_bed, create_genome_file


__all__ = ["annotate_homer2_motifs", "make_homer2_output"]


def is_value_lower(value_to_check: float, upper_bound: float) -> bool:
    return value_to_check <= upper_bound


def get_motif_files_with_tolerance(motif_path_folder: str,
                                   max_tolerated_ambiguity: int = 0,
                                   p_threshold: float = 0.05) -> list[tuple[str, str]]:
    file_ending: str = ".motif"
    # Define ambiguous bases
    dna_bases: set = {'A', 'T', 'C', 'G'}
    motif_file_paths: list[tuple[str, str]] = []
    for file in os.listdir(motif_path_folder):
        if file.endswith(file_ending):
            file_path = os.path.join(motif_path_folder, file)
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        sequence_part: str = line.split('\t', 1)[0]
                        sequence_part = sequence_part[1:]
                        p_part: str = line.rsplit(',', 1)[-1].strip()
                        if not is_value_lower(float(p_part[2:]), p_threshold):
                            continue
                        # Count ambiguous bases in the motif
                        count = sum(base not in dna_bases for base in sequence_part)
                        if count <= max_tolerated_ambiguity:
                            motif_file_paths.append((sequence_part, p_part[2:]))
    return motif_file_paths


def get_position_in_fasta(fasta_file_path: str, start_intervals: list[int], end_intervals: list[int],
                          sequences: list[tuple[str, str]]) -> list[tuple[int, int, str, str]]:
    ret: list[tuple[int, int, str, str]] = []
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        for start_interval, end_interval in zip(start_intervals, end_intervals):
            subsequence = str(record.seq[start_interval - 1:end_interval])
            for seq in sequences:
                position = subsequence.find(seq[0])
                if position != -1:
                    ret.append((start_interval + (position -1), start_interval + (position -1) + len(seq[0]) - 1, seq[0], seq[1]))
    return ret


def write_to_output_folder(output_folder: str,
                           file_body: str,
                           lines: list[tuple[int, int, str, str]],
                           ):
    homer2_output_path: str = f"{output_folder}/homer2_annotation"
    create_folder(homer2_output_path)
    with open(f"{homer2_output_path}/{file_body}.txt", 'w') as file:
        file.write("start\tend\tp_value\tseq\n")
        for interval in lines:
            file.write(f"{interval[0]}\t{interval[1]}\t{interval[3]}\t{interval[2]}\n")


def annotate_homer2_motifs(fasta_folder_path: str,
                           gff3_folder_path: str,
                           homer2_folder_path: str,
                           sequence_len_before_gene: int,
                           p_threshold: float,
                           output_folder: str
                           ):
    homer2_folder_post_fix: str = "_homer2"
    # /motifResults/knownResults
    for file in os.listdir(gff3_folder_path):
        gene_boundary: list[tuple[int, int]] = parse_gff3_file(f"{gff3_folder_path}/{file}")  # all gene intervals
        name_part, _ = os.path.splitext(file)
        # open homer2
        homer_folder_path: str = f"{homer2_folder_path}/{name_part}{homer2_folder_post_fix}"
        homer_known_motif_path: str = f"{homer_folder_path}/motifResults/knownResults"
        ret_homer_motifs: list[tuple[str, str]] = get_motif_files_with_tolerance(homer_known_motif_path,
                                                                                      p_threshold=p_threshold)
        # open fasta
        fasta_file_path: str = f"{fasta_folder_path}/{name_part}.fasta"
        ret: list[tuple[int, int, str, str]] = get_position_in_fasta(fasta_file_path,
                              [i[0] - sequence_len_before_gene -1 for i in gene_boundary],
                              [i[0] - 1 for i in gene_boundary],
                              ret_homer_motifs)
        write_to_output_folder(output_folder, name_part, ret)


def make_homer2_output(fasta_folder_path: str,
                       gff3_folder_path: str,
                       homer2_output_folder_path: str,
                       sequence_len_before_gene: int
                       ):
    for file in os.listdir(gff3_folder_path):
        name_part, _ = os.path.splitext(file)
        out_path: str = f"{homer2_output_folder_path}/{name_part}_promoter.bed"
        extract_promoters_to_bed(f"{gff3_folder_path}/{file}", out_path, sequence_len_before_gene)  # all gene positive intervals
        # -----
        fasta_path: str = f"{fasta_folder_path}/{name_part}.fasta"
        promoter_fasta_path: str = f"{homer2_output_folder_path}/{name_part}_promoter.fasta"
        get_fasta_from_bed(fasta_path, out_path, promoter_fasta_path)
        output_genome: str = f"{homer2_output_folder_path}/{name_part}.genome"
        create_genome_file(fasta_path, output_genome)
    pass

