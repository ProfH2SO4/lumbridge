import os
import pybedtools
import subprocess
from Bio import SeqIO

from .common import (
    parse_gff3_file,
    create_folder,
    extract_promoters_to_bed,
    get_fasta_from_bed,
    create_genome_file,
)


__all__ = ["annotate_homer2_motifs", "make_homer2_output"]


def is_value_lower(value_to_check: float, upper_bound: float) -> bool:
    return value_to_check <= upper_bound


def get_motif_files_with_tolerance(
    motif_path_folder: str, max_tolerated_ambiguity: int = 0, p_threshold: float = 0.05
) -> list[tuple[str, str]]:
    file_ending: str = ".motif"
    # Define ambiguous bases
    dna_bases: set = {"A", "T", "C", "G"}
    motif_file_paths: list[tuple[str, str]] = []
    for file in os.listdir(motif_path_folder):
        if file.endswith(file_ending):
            file_path = os.path.join(motif_path_folder, file)
            with open(file_path, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        sequence_part: str = line.split("\t", 1)[0]
                        sequence_part = sequence_part[1:]
                        p_part: str = line.rsplit(",", 1)[-1].strip()
                        if not is_value_lower(float(p_part[2:]), p_threshold):
                            continue
                        # Count ambiguous bases in the motif
                        count = sum(base not in dna_bases for base in sequence_part)
                        if count <= max_tolerated_ambiguity:
                            motif_file_paths.append((sequence_part, p_part[2:]))
    return motif_file_paths


def get_position_in_fasta(
    fasta_file_path: str,
    start_intervals: list[int],
    end_intervals: list[int],
    sequences: list[tuple[str, str]],
) -> list[tuple[int, int, str, str]]:
    ret: list[tuple[int, int, str, str]] = []
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        for start_interval, end_interval in zip(start_intervals, end_intervals):
            subsequence = str(record.seq[start_interval - 1 : end_interval])
            for seq in sequences:
                position = subsequence.find(seq[0])
                if position != -1:
                    ret.append(
                        (
                            start_interval + (position - 1),
                            start_interval + (position - 1) + len(seq[0]) - 1,
                            seq[0],
                            seq[1],
                        )
                    )
    return ret


def write_to_output_folder(
    output_folder: str,
    file_body: str,
    lines: list[tuple[int, int, str, str]],
):
    homer2_output_path: str = f"{output_folder}/homer2_annotation"
    create_folder(homer2_output_path)
    with open(f"{homer2_output_path}/{file_body}.txt", "w") as file:
        file.write("start\tend\tp_value\tseq\n")
        for interval in lines:
            file.write(f"{interval[0]}\t{interval[1]}\t{interval[3]}\t{interval[2]}\n")


def annotate_homer2_motifs(
    fasta_folder_path: str,
    gff3_folder_path: str,
    homer2_folder_path: str,
    sequence_len_before_gene_start: int,
    p_threshold: float,
    output_folder: str,
) -> None:
    """
    Homer2 gives me motifs of promotors in desired sequences, but it does annotate them.
    :param fasta_folder_path:
    :param gff3_folder_path:
    :param homer2_folder_path:
    :param sequence_len_before_gene_start:
    :param p_threshold:
    :param output_folder:
    :return:
    """
    homer2_folder_post_fix: str = "_output_dir"
    # /motifResults/knownResults
    for file in os.listdir(gff3_folder_path):
        gene_boundary: list[tuple[int, int]] = parse_gff3_file(
            f"{gff3_folder_path}/{file}"
        )  # all gene intervals
        name_part, _ = os.path.splitext(file)
        # open homer2
        homer_folder_path: str = (
            f"{homer2_folder_path}/{name_part}{homer2_folder_post_fix}"
        )
        homer_known_motif_path: str = f"{homer_folder_path}/knownResults"
        ret_homer_motifs: list[tuple[str, str]] = get_motif_files_with_tolerance(
            homer_known_motif_path, p_threshold=p_threshold
        )
        # open fasta
        fasta_file_path: str = f"{fasta_folder_path}/{name_part}.fasta"
        ret: list[tuple[int, int, str, str]] = get_position_in_fasta(
            fasta_file_path,
            [i[0] - sequence_len_before_gene_start - 1 for i in gene_boundary],
            [i[0] - 1 for i in gene_boundary],
            ret_homer_motifs,
        )
        write_to_output_folder(output_folder, name_part, ret)


def create_background_file(
    fasta_path: str,
    genome_file: str,
    promoter_bed_path: str,
    homer2_output_background_path: str,
    homer2_output_folder: str,
    sequence_len_before_gene_start: int,
    seq_len_after_inclusive_gene_start: int,
    file_name: str,
) -> None:
    # Create random intervals and save to a file
    len_of_seq: int = (
        sequence_len_before_gene_start + seq_len_after_inclusive_gene_start
    )
    random_intervals = pybedtools.BedTool().random(l=len_of_seq, n=10000, g=genome_file)
    random_intervals.saveas(f"{homer2_output_folder}/{file_name}_random_intervals.bed")

    # Subtract promoters from random intervals and save to a file
    promoters = pybedtools.BedTool(promoter_bed_path)
    background_intervals = random_intervals.subtract(promoters)
    background_intervals.saveas(
        f"{homer2_output_folder}/{file_name}_background_intervals.bed"
    )

    # Get fasta sequences for the background intervals
    background_intervals.sequence(
        fi=fasta_path, fo=homer2_output_background_path, name=True
    )


def run_find_motifs(
    promoters_fasta: str, output_directory: str, background_fasta: str, cpu_cores: int
):
    # Set the PATH environment variable
    homer_bin_path = "/home/matej/homer2/bin"
    current_env = os.environ.copy()
    current_env["PATH"] = homer_bin_path + os.pathsep + current_env["PATH"]

    cmd = [
        "findMotifs.pl",
        promoters_fasta,
        "fasta",
        output_directory,
        "-len",
        "6,8,10",
        "-p",
        str(cpu_cores),
        "-fasta",
        background_fasta,
    ]
    subprocess.run(cmd, env=current_env)


def make_homer2_output(
    fasta_folder_path: str,
    gff3_folder_path: str,
    homer2_output_folder_path: str,
    sequence_len_before_gene_start: int,
    sequence_len_after_gene_start: int,
    cpu_cores: int,
) -> None:
    """

    :param fasta_folder_path:
    :param gff3_folder_path:
    :param homer2_output_folder_path:
    :param sequence_len_before_gene_start: if 200, and gene_start at pos 1000 => [1000- 200, 1000)
    :param sequence_len_after_gene_start: if 200, and gene_start at pos 1000 => [1000, 1000 + 200)
    :param cpu_cores:
    :return:
    """
    create_folder(homer2_output_folder_path)
    for file in os.listdir(gff3_folder_path):
        name_part, _ = os.path.splitext(file)
        out_path: str = f"{homer2_output_folder_path}/{name_part}_promoter.bed"
        fasta_path: str = f"{fasta_folder_path}/{name_part}.fasta"
        promoter_fasta_path: str = (
            f"{homer2_output_folder_path}/{name_part}_promoter.fasta"
        )
        output_genome: str = f"{homer2_output_folder_path}/{name_part}.genome"
        background_fasta_path: str = (
            f"{homer2_output_folder_path}/{name_part}_background.fasta"
        )
        output_dir_path: str = f"{homer2_output_folder_path}/{name_part}_output_dir"

        extract_promoters_to_bed(
            f"{gff3_folder_path}/{file}",
            out_path,
            sequence_len_before_gene_start,
            sequence_len_after_gene_start,
        )  # all gene positive intervals
        get_fasta_from_bed(fasta_path, out_path, promoter_fasta_path)
        create_genome_file(fasta_path, output_genome)
        create_background_file(
            fasta_path,
            output_genome,
            out_path,
            background_fasta_path,
            homer2_output_folder_path,
            sequence_len_before_gene_start,
            sequence_len_after_gene_start,
            name_part,
        )

        run_find_motifs(
            promoter_fasta_path,
            output_dir_path,
            background_fasta_path,
            cpu_cores=cpu_cores,
        )
