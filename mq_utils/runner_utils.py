#!/bin/python3

import os
import sys
import argparse
import time
import psutil
import pandas as pd
import csv
from pathlib import Path

start_time = time.time()


def save_parameters(script_desc: str, arguments):
    """
    Save command line options into local variables.

    :return: values assigned to input arguments
    """
    descr = "\n############################################################################\n"
    descr += "################# MaxQuantHandler - %(prog)s ##################\n"
    descr += script_desc
    descr += "\n############################################################################\n"
    descr += "\nusage: python3 %(prog)s [required arguments] [optional arguments]\n"
    epilo = _get_epilog(script_name=os.path.basename(sys.argv[0]))
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawTextHelpFormatter, epilog=epilo,
                                     usage=argparse.SUPPRESS, add_help=False)
    required_args = parser.add_argument_group("required arguments")
    if 'qf' in arguments:
        required_mut = required_args.add_mutually_exclusive_group(required=True)
        required_mut.add_argument('-q', '--maxquant_file', type=str, help='MaxQuant file', default=None)
        required_mut.add_argument('-s', '--single_file', type=str, help='Single file', default=None)
    if 'q' in arguments:
        required_args.add_argument('-q', '--maxquant_file', type=str, help='MaxQuant file', required=True)
    if 'f_req' in arguments:
        required_args.add_argument('-f', '--fasta_file', type=str, help='Fasta file', required=True)
    if 'or_req' in arguments:
        required_args.add_argument('-or', '--organism', choices=["human", "mouse", "rat", "rabbit"], type=str,
                                   required=True, help='Specify organism the ids should match to.')
    if 'tor_req' in arguments:
        required_args.add_argument('-or', '--organism', choices=["human", "mouse", "rat", "rabbit"], type=str,
                                   required=True, help='Specify organism the ids are mapped to.')
        required_args.add_argument('-tor', '--tar_organism', choices=["human", "mouse", "rat", "rabbit"], type=str,
                                   required=True, help='Specify organism from which orthologs should be mapped.')
    if 'm' in arguments:
        required_args.add_argument('-m', '--mode',
                                   choices=['all', 'fasta', 'uniprot', 'uniprot_one', 'uniprot_primary'],
                                   type=str, required=True, help='Mode of refilling. See below for more infos.')
    if 'i' in arguments:
        required_args.add_argument('-i', '--in_type', choices=['protein', 'gene'], required=True,
                                   help='Define what type should be the source.')
    optional_args = parser.add_argument_group("optional arguments")
    if 'c' in arguments:
        optional_args.add_argument('-pc', '--protein_column', type=str, default=None,
                                   help='Name of column with protein IDs [Default=None]')
        optional_args.add_argument('-gc', '--gene_column', type=str, default=None,
                                   help='Name of column with gene names [Default=None]')
    if 'f' in arguments:
        optional_args.add_argument('-f', '--fasta_file', type=str, help='Fasta file', default=None)
    if 'l' in arguments:
        optional_args.add_argument('-l', '--fill', action='store_false', default=True,
                                   help='Use this flag, if filled values should be skipped. [Default=True]')
    if 'a' in arguments:
        optional_args.add_argument('-a', '--action', type=str, default="delete", choices=['keep', 'delete'],
                                   help='What to do, if IDs cell is empty after filtering. '
                                        'Keep empty cell or delete it.')
    if 'or' in arguments:
        optional_args.add_argument('-or', '--organism', choices=["human", "mouse", "rat", "rabbit"], type=str,
                                   default=None, help='Specify organism the ids should match to.')
    if 'r' in arguments:
        optional_args.add_argument('-r', '--reviewed', action='store_true', default=False,
                                   help='Bool to indicate if only reviewed protein IDs should be kept.')
    if 'd' in arguments:
        optional_args.add_argument('-d', '--decoy', action='store_true', default=False,
                                   help='Set flag if protein ids from decoy fasta (REV__, CON__) should be kept.')
    if 'o' in arguments:
        optional_args.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory. [Default=./]')
    optional_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    args = parser.parse_args()
    if 'qf' in arguments:
        if args.maxquant_file is not None:
            args.data = pd.read_table(args.maxquant_file, sep=find_delimiter(args.maxquant_file)).fillna("")
            args.file_name = Path(args.maxquant_file).stem
        else:
            args.data = pd.read_table(args.single_file).fillna("")
            args.file_name = Path(args.single_file).stem
    if 'q' in arguments:
        args.data = pd.read_table(args.maxquant_file, sep=find_delimiter(args.maxquant_file)).fillna("")
        args.file_name = Path(args.maxquant_file).stem
    return args


def _get_epilog(script_name):
    epilog = ""
    if script_name == 'remap_genenames.py':
        epilog += "\n----------------------------------------------------------------------------\n"
        epilog += "\nsupported modes\n"
        epilog += "  all\t\t\tUse primarly fasta infos and additionally uniprot infos.\n"
        epilog += "  fasta\t\t\tUse information extracted from fasta headers.\n"
        epilog += "  uniprot\t\tUse mapping information from uniprot and use all gene names.\n"
        epilog += "  uniprot_primary\tUse mapping information from uniprot and only all primary gene names.\n"
        epilog += "  uniprot_one\t\tUse mapping information from uniprot and only use most frequent single gene name.\n"
    epilog += "\n############################################################################\n"
    return epilog


def print_current_usage(text):
    memory_usage = '{0:.2f}'.format(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    time_usage = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print('[{}|{}MB] '.format(time_usage, memory_usage) + text)


def find_delimiter(filename):
    sniffer = csv.Sniffer()
    with open(filename) as fp:
        delimiter = sniffer.sniff(fp.readline()).delimiter
    return delimiter
