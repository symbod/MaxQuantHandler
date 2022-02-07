#!/usr/bin/python3

import pandas as pd
import runner_utils as ru
import uniprot_handler as uh


def get_uniprot_mappings(mq_file, organism=None):
    """
    Get uniprot mapping to protein ids in MaxQuant file. Optionally only mapped to chosen organism.

    :param mq_file: MaxQuant file
    :param organism: Organism to map to
    :return: Uniprot mapping to protein ids as dataframe
    """
    handler = uh.Uniprot_Handler()
    max_quant = pd.read_table(mq_file, sep=" ").fillna("")
    mappings = pd.DataFrame(columns=['Gene names', 'Gene names  (primary )', 'Status', 'Organism', 'Protein ID'])
    for index, row in max_quant.iterrows():
        mappings = pd.concat(
            [mappings, handler.get_mapping(ids=row['Protein IDs'].split(";"), organism=organism)])
    handler.save_mappings()
    return mappings


if __name__ == "__main__":
    description = "                   Get uniprot mapping to protein ids optionally by organism."
    parameters = ru.save_parameters(script_desc=description, arguments=('q', 'r', 'o'))
    df = get_uniprot_mappings(mq_file=parameters.maxquant_file, organism=parameters.organism)
    df.to_csv(parameters.out_dir + "uniprot_mapping.txt", header=True)
