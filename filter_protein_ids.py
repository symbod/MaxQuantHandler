#!/usr/bin/python3

import csv
import pandas as pd
import runner_utils as ru
import uniprot_handler as uh


def filter_protein_ids(data: pd.DataFrame, organism: str = None, decoy: bool = False, action: str = "delete",
                       reviewed: bool = True):
    """
    Filter protein ids in given data by chosen organism.

    :param data: MaxQuant data or data with one single column
    :param organism: Organism to map to
    :param decoy: bool to indicate if decoy IDs (REV_, ) should be kept
    :param action: What to do, if IDs cell is empty after filtering. Keep empty cell, delete it
    or fill it based on gene name (if it exists).
    :param reviewed: Set True if during action=fill only reviewed protein IDs should be taken
    :return: Filtered MaxQuant file as dataframe
    """
    id_column = "Protein IDs"
    handler = uh.UniprotHandler()
    # Get all existing mappings in one batch
    handler.get_mapping(ids=";".join(data[id_column]).split(";"),
                        in_type="proteinID", organism=organism)

    # filter row wise
    data[id_column] = data[id_column].apply(
        lambda x: handler.get_filtered_ids(ids=x.split(";"), organism=organism, decoy=decoy))
    handler.save_mappings()
    # keep fill or remove
    if action == "delete":
        data = data[data[id_column] != ""]  # remove
    elif action == "fill":
        if 'Gene names' in data.columns:
            data[id_column] = data.apply(lambda row:
                                         handler.get_ids_from_gene(genenames=row['Gene names'].split(";"),
                                                                   organism=organism, reviewed=reviewed)
                                         if row[id_column] == "" else row[id_column], axis=1)
        data = data[data[id_column] != ""]  # remove no matchable
    return data


if __name__ == "__main__":
    description = "               Filter protein ids by organism and/or decoy names."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'or_req', 'r', 'd', 'a', 'o'))
    df = filter_protein_ids(data=parameters.data, organism=parameters.organism, decoy=parameters.decoy,
                            action=parameters.action, reviewed=parameters.reviewed)
    df.to_csv(parameters.out_dir + parameters.file_name + "_filtered.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
