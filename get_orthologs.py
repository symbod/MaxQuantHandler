#!/usr/bin/python3

import csv
import pandas as pd
import runner_utils as ru
import mapping_handler as mh


def get_orthologs(data: pd.DataFrame, organism: str, tar_organism: str):
    """
    Filter protein ids in given data by chosen organism.

    :param data: MaxQuant data or data with one single column
    :param organism: Organism of the input ids
    :param tar_organism: Organism to map to

    :return: MaxQuant file as dataframe with ortholog ids
    """
    handler = mh.MappingHandler()

    # Get all existing mappings in one batch
    handler.get_mapping(ids=";".join(data['Gene names']).split(";"),
                        in_type="orthologs", organism=organism, tar_organism=tar_organism)
    handler.save_mappings()

    # ==== If Input was single column file ====
    if len(data.columns) == 1:
        data['Ortholog Gene names'] = data['Gene names'].apply(
            lambda x: handler.get_orthologs(ids=x.split(";"), organism=organism, tar_organism=tar_organism))

    # ==== If Input was MaxQuant file ====
    else:
        data['Gene names'] = data['Gene names'].apply(
            lambda x: handler.get_orthologs(ids=x.split(";"), organism=organism, tar_organism=tar_organism))

    return data


if __name__ == "__main__":
    description = "               Get ortholog gene names."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'tor_req', 'r', 'd', 'a', 'o'))
    df = get_orthologs(data=parameters.data, organism=parameters.organism, tar_organism=parameters.torganism)
    df.to_csv(parameters.out_dir + parameters.file_name + "_ortholog.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
