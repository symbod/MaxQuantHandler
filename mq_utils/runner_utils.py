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
    if 'd' in arguments:
        required_args.add_argument('-d', '--data', type=str, help='Data file', default=None, required=True)
    if 'f_req' in arguments:
        required_args.add_argument('-f', '--fasta_file', type=str, help='Fasta file', required=True)
    if 'or_req' in arguments:
        required_args.add_argument('-or', '--organism', choices=["human", "mouse", "rat", "rabbit"], type=str,
                                   required=True, help='Specify organism the ids should match to.')
    if 'tor_req' in arguments:
        required_args.add_argument('-tor', '--tar_organism', choices=["human", "mouse", "rat", "rabbit"], type=str,
                                   required=True, help='Specify organism from which orthologs should be mapped.')
    if 'pc_req' in arguments:
        required_args.add_argument('-pc', '--protein_column', type=str,
                                   help='Name of column with protein IDs.', required=True)
    if 'gc_req' in arguments:
        required_args.add_argument('-gc', '--gene_column', type=str, default=None,
                                   help='Name of column with gene names.', required = True)
    if 'm' in arguments:
        required_args.add_argument('-m', '--mode',
                                   choices=['all', 'fasta', 'uniprot', 'uniprot_one', 'uniprot_primary'],
                                   type=str, required=True, help='Mode of refilling. See below for more infos.')
    if 'rm' in arguments:
        required_args.add_argument( '-m', '--mode',
                                    choices=[ 'ensembl', 'mygeneinfo', 'HGNC', 'enrichment'],
                                    type=str, required=True, help='Mode of reducing. See below for more infos.' )
    if 'i' in arguments:
        required_args.add_argument('-i', '--in_type', choices=['protein', 'gene'], required=True,
                                   help='Define what type should be the source.')

    optional_args = parser.add_argument_group("optional arguments")
    if 'pc' in arguments:
        optional_args.add_argument('-pc', '--protein_column', type=str, default=None,
                                   help='Name of column with protein IDs [Default=None]')
    if 'gc' in arguments:
        optional_args.add_argument('-gc', '--gene_column', type=str, default=None,
                                   help='Name of column with gene names [Default=None]')

    if 'ke' in arguments:
        optional_args.add_argument('-ke', '--keep_empty', action='store_true', default=True,
                                   help = "Bool to indicate whether empty rows should be kept. [Default=True]")
    if 'f' in arguments:
        optional_args.add_argument('-f', '--fasta_file', type=str, help='Fasta file', default=None)
    if 'l' in arguments:
        optional_args.add_argument('-l', '--fill', action='store_false', default=True,
                                   help='Use this flag, if filled values should be skipped. [Default=True]')
    if 'a' in arguments:
        optional_args.add_argument('-a', '--action', type=str, default="delete", choices=['keep', 'delete'],
                                   help='What to do, if IDs cell is empty after filtering. '
                                        'Keep empty cell or delete it.')
    if 'hm' in arguments:
        optional_args.add_argument('-hm', '--hgnc_mode', type=str, choices=["mostfrequent", "all"],
                                   default="mostfrequent",
                                   help="What to do if reduce_mode is HGNC. Take most frequent gene name or all. "
                                        "[Default=mostfrequent]")
    if 'rc' in arguments:
        optional_args.add_argument('-rc', '--res_column', type=str, default=None,
                                   help='Name of output column. If None, input column will be edited. [Default = None]')
    if 'or' in arguments:
        optional_args.add_argument('-or', '--organism', choices=["human", "mouse", "rat", "rabbit"], type=str,
                                   default=None, help='Specify organism the ids should match to.')
    if 'r' in arguments:
        optional_args.add_argument('-r', '--reviewed', action='store_true', default=False,
                                   help='Bool to indicate if only reviewed protein IDs should be kept.')
    if 'rv' in arguments:
        optional_args.add_argument('-rv', '--rev_con', action='store_true', default=False,
                                   help='Set flag if decoy and contaminants IDs (REV__, CON__) should be kept.')
    if 'o' in arguments:
        optional_args.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory. [Default=./]')
    optional_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    args = parser.parse_args()
    if 'd' in arguments:
        args.data = pd.read_table(args.data).fillna("")
        args.file_name = Path(args.single_file).stem
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
    if script_name == "reduce_genenames.py":
        epilog += "\n----------------------------------------------------------------------------\n"
        epilog += "\nsupported modes\n"
        epilog += "  ensembl\t\tUse gProfiler to reduce gene names to those have an Ensembl ID.\n"
        epilog += "  mygeneinfo\t\tUse mygeneinfo database to reduce gene names to those having an entry in mygeneinfo.\n"
        epilog += "  HGNC\t\tUse HGNC database to reduce gene names to those having an entry in HGNC (only for human).\n"
        epilog += "  enrichment\tUse gProfiler to reduce gene names to those having a functional annotation.\n"
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
