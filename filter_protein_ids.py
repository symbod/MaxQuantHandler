#!/usr/bin/python3

import csv
import pandas as pd
import runner_utils as ru
import uniprot_handler as uh
from pathlib import Path


def filter_protein_ids(mq_file, organism=None, decoy=False):
    """
    Filter protein ids in MaxQuant file by chosen organism.

    :param mq_file: MaxQuant file
    :param organism: Organism to map to
    :return: Filtered MaxQuant file as dataframe
    """
    handler = uh.UniprotHandler()
    max_quant = pd.read_table(mq_file, sep=" ").fillna("")
    max_quant['Protein IDs'] = max_quant['Protein IDs'].apply(
        lambda x: handler.get_filtered_ids(ids=x.split(";"), organism=organism, decoy=decoy))
    handler.save_mappings()
    max_quant = max_quant[max_quant['Protein IDs'] != ""]
    return max_quant


if __name__ == "__main__":
    description = "               Filter protein ids by organism and/or decoy names."
    parameters = ru.save_parameters(script_desc=description, arguments=('q', 'r', 'd', 'o'))
    df = filter_protein_ids(mq_file=parameters.maxquant_file, organism=parameters.organism, decoy=parameters.decoy)
    df.to_csv(parameters.out_dir + Path(parameters.maxquant_file).stem + "_filtered.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
