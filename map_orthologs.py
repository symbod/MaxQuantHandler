#!/usr/bin/python3

import csv
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru
from mq_utils.logger import get_ortholog_genenames_logging


def get_orthologs(data: pd.DataFrame, gene_column: str, organism: str, tar_organism: str,
                  res_column: str = None, return_log: bool = True):
    """
    Map gene names of origin organism to orthologs of target organism.

    :param data: Dataframe containing a column with gene names
    :param gene_column: Column name with gene names
    :param organism: Organism of the input ids
    :param tar_organism: Organism to map to
    :param res_column: Set column name for ortholog results. If None, the gene_column will be overridden.
    :param return_log: Set True if log dataframes should be returned

    :return: Data as dataframe with ortholog ids
    """

    data_copy = data.copy(deep=True)

    handler = mh.MappingHandler(mapping_dir="mappings/")
    # ==== Get all existing mappings in one batch ====
    handler.get_mapping(ids=";".join(data_copy[gene_column]).split(";"),
                        in_type="orthologs", organism=organism, tar_organism=tar_organism)
    handler.save_mappings(mapping_dir="mappings/")

    ortholog_gene_names= data_copy[gene_column].apply(
        lambda x: handler.get_orthologs(ids=x.split(";"), organism=organism, tar_organism=tar_organism))

    # ==== Logging ====
    log_dict = dict()
    if return_log:
        log_dict = get_ortholog_genenames_logging(data_copy[gene_column], ortholog_gene_names)

    # ==== If target column depending if res_column is set ====
    column = res_column if res_column is not None else gene_column

    # ==== Set Reduced Gene Names To DataFrame ====
    data_copy[column] = ortholog_gene_names

    return data_copy, log_dict


if __name__ == "__main__":
    description = "                       Get ortholog gene names."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'tor_req', 'c', 'o'))
    df = get_orthologs(data=parameters.data, organism=parameters.organism, tar_organism=parameters.tar_organism,
                       gene_column=parameters.gene_column)
    df.to_csv(parameters.out_dir + parameters.file_name + "_ortholog.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
