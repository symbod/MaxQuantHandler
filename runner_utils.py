#!/bin/python3

import os
import sys
import argparse
import time
import psutil

start_time = time.time()
organisms = {"human": "Homo sapiens (Human)",
             "rat": "Rattus norvegicus (Rat)"}

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
    if 'f_req' in arguments:
        required_args.add_argument('-f', '--fasta_file', type=str, help='Fasta file', required=True)
    if 'r_req' in arguments:
        required_args.add_argument('-r', '--organism', choices=organisms.keys(), type=str, required=True,
                                   help='Specify organism the ids should match to.')
    if 'q' in arguments:
        required_args.add_argument('-q', '--maxquant_file', type=str, help='MaxQuant file', required=True)
    if 'm' in arguments:
        required_args.add_argument('-m', '--mode', choices=['all', 'fasta','uniprot','uniprot_one'], type=str,
                                   required=True, help='Mode of refilling. See below for more infos.')
    if 'i' in arguments:
        required_args.add_argument('-i', '--in_type', choices=['proteinID', 'genename'], required=True,
                                   help='Define what type should be the source.')
    optional_args = parser.add_argument_group("optional arguments")
    if 'f' in arguments:
        required_args.add_argument('-f', '--fasta_file', type=str, help='Fasta file', default=None)
    if 'l' in arguments:
        optional_args.add_argument('-l', '--fill', action='store_true', default=False,
                                   help='Use this flag, if only missing values should be filled.')
    if 'r' in arguments:
        optional_args.add_argument('-r', '--organism', choices=organisms.keys(), type=str, default=None,
                                   help='Specify organism the ids should match to.')
    if 'o' in arguments:
        optional_args.add_argument('-o', '--out_dir', type=str, default='./', help='Output directory. [Default=./]')
    optional_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    args = parser.parse_args()
    return args

def _get_epilog(script_name):
    epilog = ""
    if script_name == 'remap_genenames.py':
        epilog += "\n----------------------------------------------------------------------------\n"
        epilog += "\nsupported modes\n"
        epilog += "  all\t\t\tUse primarly fasta infos and additionally uniprot infos.\n"
        epilog += "  fasta\t\t\tUse information extracted from fasta headers.\n"
        epilog += "  uniprot\t\tUse mapping information from uniprot and use all gene names.\n"
        epilog += "  uniprot_one\t\tUse mapping information from uniprot and only use most frequent single gene name.\n"
    epilog += "\n############################################################################\n"
    return epilog


def print_current_usage(text):
    memory_usage = '{0:.2f}'.format(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2)
    time_usage = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print('[{}|{}MB] '.format(time_usage, memory_usage) + text)
