#!/usr/bin/env python
import sys
import os
import concurrent.futures
import multiprocessing
import argparse
import logging

args_parser = argparse.ArgumentParser(
    prog="SELEXderep",
    description="SELEXderep combines FASTA files acquired from a SELEX process. "
                "The FASTA files should be preprocessed (filtered, merged and primers cut). "
                "SELEXderep writes to STDOUT or writes to FASTA files with the primerless sequences as identifier "
                "and the corresponding read counts per SELEX round.\n"
                "\n"
                "Example for sequence \'CATCATCAT\' for 5 SELEX rounds. "
                "Primers \'GG\' are attached on both sides. \n"
                ">CATCATCAT 0-2-9-17-90\n"
                "GGCATCATCATGG\n"
                "\n"
                "Author: Ulrich Aschl (ulrich.aschl@tuwien.ac.at)",
    formatter_class=argparse.RawTextHelpFormatter
)
args_parser.add_argument("SELEX_Rounds", metavar="SELEX_Round_FASTA", type=str, nargs='+', help=
"FASTA files from SELEX process. Files have to be in sequential order")

args_parser.add_argument("-o", "--out-fasta", type=str, default=None, help=
"Output fasta file location. Default: stdout")
args_parser.add_argument("--force", action="store_true", default=False, help=
"Overwrite output files if they already exists")

args_parser.add_argument("-c", "--out-csv", type=str, default=None, help=
"Output csv file location. Default: no csv output")

args_parser.add_argument("-t", "--threads", type=int, default=None, help=
"Define how many threads should be used. Uses all cores by default.")

# args_parser.add_argument("--input-untrimmed", action="store_true", default=False, help=
# "Specifies if input FASTA files have their SELEX primers still attached.")

args_parser.add_argument("--add-primers", action="store_true", default=False, help=
"Add forward and reverse primers to sequences. If used, please also specify primer sequences.")
args_parser.add_argument("--primer-forward", type=str, default="")
args_parser.add_argument("--primer-reverse", type=str, default="")


def get_round_name(fasta_file_path):
    return os.path.splitext(os.path.basename(fasta_file_path))[0]


def read_file(fasta_file_path, selex_dict: dict):
    round_name = get_round_name(fasta_file_path)
    fasta_file = open(fasta_file_path, "r")

    for seq_id in fasta_file:
        seq = fasta_file.readline().rstrip('\n')

        # dicts are thread safe
        selex_dict[seq] = selex_dict.get(seq, {})
        selex_dict[seq][round_name] = selex_dict[seq].get(round_name, 0) + 1


def read_files(fasta_files, threads, sequences: dict()):
    # Reading of sequences
    thread_executor = concurrent.futures.ThreadPoolExecutor(max_workers=threads)
    futures = []
    for fasta_file in fasta_files:
        futures.append(thread_executor.submit(read_file, fasta_file, sequences))
    return futures


def print_fasta_file(out_file, round_names, sequences, p5, p3):
    # Print counted sequences to out_file
    for sequence, round_counts in sequences.items():
        str_round_counts = '-'.join(str(round_counts.get(round_name, 0)) for round_name in round_names)

        # WRITE TO FILE:
        out_file.write(">{} {}\n".format(sequence, str_round_counts))
        out_file.write("{}{}{}\n".format(p5, sequence, p3))

    out_file.flush()


def print_csv_file(out_file, round_names, sequences, p5, p3):
    # CSV Header
    out_file.write("id\tseq\tp5_seq_p3\t{}\n".format('\t'.join(round_names)))

    id = 1
    for sequence, round_counts in sequences.items():
        str_round_counts = '\t'.join(str(round_counts.get(round_name, 0)) for round_name in round_names)

        # WRITE TO FILE:
        out_file.write("{}\t{}\t{}{}{}\t{}\n".format(
            str(id),
            sequence,
            p5, sequence, p3,
            str_round_counts)
        )
        id = id + 1

    out_file.flush()


def main():
    args = args_parser.parse_args()

    # ======= Parameter verification =======
    threads = args.threads
    if threads is None:
        try:
            threads = os.cpu_count()
        except NotImplementedError:
            threads = 1

    out_fasta = args.out_fasta
    if (out_fasta is not None) and (not args.force) and os.path.exists(out_fasta):
        print("Output fasta file {} already exists. Add --force to overwrite the output file".format(
            out_fasta))
        print("Halting.")
        return 1

    out_csv = args.out_csv
    if (out_csv is not None) and (not args.force) and os.path.exists(out_csv):
        print("Did not start! Output csv file {} already exists. Add --force to overwrite the output file".format(
            out_csv))
        print("Halting.")
        return 1

    add_primers = args.add_primers
    if add_primers:
        p5 = args.primer_forward
        p3 = args.primer_reverse
        if p5 == "" or p3 == "":
            print("Did not start! If flag \'--add-primer\' is set, please specify forward and reverse primers using "
                  "\'--primer-forward\' and \'--primer-reverse\'")
            return 1
    else:
        p5 = ""
        p3 = ""

    fasta_files = args.SELEX_Rounds
    # check if fasta files exist
    missing_file_flag = False
    for fasta_file in fasta_files:
        if not os.path.exists(fasta_file):
            print("Input file {} does not seem to be at the specified path.".format(fasta_file))
            missing_file_flag = True
    if missing_file_flag:
        print("Halting.")
        return 1

    # ======= Read fasta files =======
    sequences = {}
    thread_futures = read_files(fasta_files, threads, sequences)
    # wait to finish reading
    concurrent.futures.wait(thread_futures, return_when=concurrent.futures.ALL_COMPLETED)

    # ======= Print fasta file =======
    round_names = [get_round_name(fasta_file) for fasta_file in fasta_files]
    if out_fasta is not None:
        out_file = open(out_fasta, "w")
        print_fasta_file(out_file, round_names, sequences, p5, p3)
        out_file.close()
    else:
        print_fasta_file(sys.stdout, round_names, sequences, p5, p3)

    # ======= Print csv file =======
    if out_csv is not None:
        out_file = open(out_csv, "w")
        print_csv_file(out_file, round_names, sequences, p5, p3)
        out_file.close()


main()
