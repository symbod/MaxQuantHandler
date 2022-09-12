#!/usr/bin/python3

import csv
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru


def get_orthologs(data: pd.DataFrame, gene_column: str, organism: str, tar_organism: str):
    """
    Map gene names of origin organism to orthologs of target organism.

    :param data: MaxQuant data or data with one single column
    :param gene_column: Name of column in gene names
    :param organism: Organism of the input ids
    :param tar_organism: Organism to map to

    :return: Data as dataframe with ortholog ids
    """
    handler = mh.MappingHandler(mapping_dir="mappings/")

    # Get all existing mappings in one batch
    handler.get_mapping(ids=";".join(data[gene_column]).split(";"),
                        in_type="orthologs", organism=organism, tar_organism=tar_organism)
    handler.save_mappings(mapping_dir="mappings/")

    # ==== If Input was single column file ====
    if len(data.columns) == 1:
        data['Ortholog Gene names'] = data[gene_column].apply(
            lambda x: handler.get_orthologs(ids=x.split(";"), organism=organism, tar_organism=tar_organism))

    # ==== If Input was MaxQuant file ====
    else:
        data[gene_column] = data[gene_column].apply(
            lambda x: handler.get_orthologs(ids=x.split(";"), organism=organism, tar_organism=tar_organism))
    print(data[data[gene_column] != ""])
    return data


if __name__ == "__main__":
    description = "                       Get ortholog gene names."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'tor_req', 'c', 'o'))
    df = get_orthologs(data=parameters.data, organism=parameters.organism, tar_organism=parameters.tar_organism,
                       gene_column=parameters.gene_column)
    df.to_csv(parameters.out_dir + parameters.file_name + "_ortholog.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
