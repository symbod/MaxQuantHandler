#!/usr/bin/python3

import pandas as pd
import runner_utils as ru
import uniprot_handler as uh
from pathlib import Path


def filter_ids(mq_file, organism):
    """
    Filter protein ids in MaxQuant file by chosen organism.

    :param mq_file: MaxQuant file
    :param organism: Organism to map to
    :return: Filtered MaxQuant file as dataframe
    """
    handler = uh.Uniprot_Handler()
    max_quant = pd.read_table(mq_file, sep=" ").fillna("")
    max_quant['Protein IDs'] = max_quant['Protein IDs'].apply(
        lambda x: handler.get_filtered_ids(ids=x.split(";"), organism=organism))
    handler.save_mappings()
    return max_quant


if __name__ == "__main__":
    description = "                   Filter protein ids by organism."
    parameters = ru.save_parameters(script_desc=description, arguments=('q', 'r_req', 'o'))
    df = filter_ids(mq_file=parameters.maxquant_file, organism=parameters.organism)
    df.to_csv(parameters.out_dir + Path(parameters.maxquant_file).stem + "_filtered.txt", header=True)
