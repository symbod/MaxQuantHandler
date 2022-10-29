#!/usr/bin/python3

import csv
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru
from mq_utils.logger import get_filter_ids_logging

def filter_protein_ids(data: pd.DataFrame, id_column: str = "Protein IDs", organism: str = None,
                       decoy: bool = False, action: str = "delete", gene_column:str = "Gene names",
                       reviewed: bool = True, inplace: bool = True, return_log: bool = True):
    """
    Filter protein ids in given data by chosen organism.

    :param data: MaxQuant data or data with one single column
    :param id_column: Column name with protein IDs
    :param organism: Organism to map to
    :param decoy: Bool to indicate if decoy IDs (REV_, ) should be kept
    :param action: What to do, if IDs cell is empty after filtering. Keep empty cell, delete it
    or fill it based on gene name (if it exists).
    :param gene_column: if action fill is chosen, this indicated the column name with the gene names
    :param reviewed: Set True if during action=fill only reviewed protein IDs should be taken
    :param inplace: Set True if filtering should be applied on original dataframe (False for returning a new dataframe)
    :param return_log: Set True if a log dataframe should be returned
    :return: Filtered data as dataframe
    """

    # ==== Check if user want to work on input or create a deep copy ====
    if not inplace:
        data = data.copy( deep=True )

    # ==== Get all existing mappings in one batch ====
    handler = mh.MappingHandler(mapping_dir="mappings/")
    handler.get_mapping(ids=";".join(data[id_column]).split(";"),
                        in_type="protein", organism=organism)

    # ==== Filter row wise ====
    filtered_ids = data[id_column].apply(
        lambda x: handler.get_filtered_ids(ids=x.split(";"), in_type="protein", organism=organism, decoy=decoy))

    # ==== Logging ====
    if return_log:
        log_df = get_filter_ids_logging( data[ id_column ], filtered_ids )

    # ==== Set filtered ids to dataframe ====
    data[ id_column ] = filtered_ids

    # ==== Save current mapping to files
    handler.save_mappings(mapping_dir="mappings/")

    # ==== Keep fill or remove
    if action == "delete":
        data = data[data[id_column] != ""]  # remove
    elif action == "fill":
        if gene_column in data.columns:
            data[id_column] = data.apply(lambda row:
                                         handler.get_ids_from_gene(genenames=row[gene_column].split(";"),
                                                                   organism=organism, reviewed=reviewed)
                                         if row[id_column] == "" else row[id_column], axis=1)
        data = data[data[id_column] != ""]  # remove no matchable

    # ==== Return data and logging dataframe if requested
    if return_log:
        return data, log_df
    else:
        return data


def filter_gene_names(data: pd.DataFrame, id_column: str = "Gene names", organism: str = None,
                      decoy: bool = False, action: str = "delete", protein_column: str = "Protein IDs",
                      reviewed: bool = True):
    """
    Filter gene names in given data by chosen organism.

    :param data: MaxQuant data or data with one single column
    :param id_column: Column name with gene names
    :param organism: Organism to map to
    :param decoy: Bool to indicate if decoy IDs (REV_, ) should be kept
    :param action: What to do, if IDs cell is empty after filtering. Keep empty cell, delete it
    or fill it based on protein IDs (if it exists).
    :param protein_column: if action fill is chosen, this indicated the column name with the protein IDs
    :param reviewed: Set True if during action=fill only reviewed protein IDs should be taken
    :return: Filtered data as dataframe
    """
    handler = mh.MappingHandler(mapping_dir="mappings/")
    # Get all existing mappings in one batch
    handler.get_mapping(ids=";".join(data[id_column]).split(";"),
                        in_type="gene", organism=organism)
    # TODO
    raise Exception("Not yet implemented")
    # # filter row wise
    # data[id_column] = data[id_column].apply(
    #     lambda x: handler.get_filtered_ids(ids=x.split(";"), in_type="gene", organism=organism, decoy=decoy))
    # handler.save_mappings(mapping_dir="mappings/")
    # # keep fill or remove
    # if action == "delete":
    #     data = data[data[id_column] != ""]  # remove
    # elif action == "fill":
    #     if 'Protein IDs' in data.columns:
    #         data[id_column] = data.apply(lambda row:
    #                                      handler.get_ids_from_gene(genenames=row['Protein IDs'].split(";"),
    #                                                                organism=organism, reviewed=reviewed)
    #                                      if row[id_column] == "" else row[id_column], axis=1)
    #     data = data[data[id_column] != ""]  # remove no matchable
    # return data


if __name__ == "__main__":
    description = "        Filter proteins or genes by organism and/or decoy names."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'or_req', 'c', 'r', 'd', 'i', 'a', 'o'))
    if parameters.in_type == "protein":
        df = filter_protein_ids(data=parameters.data, organism=parameters.organism, decoy=parameters.decoy,
                                id_column=parameters.protein_column, gene_column=parameters.gene_column,
                                action=parameters.action, reviewed=parameters.reviewed)
    elif parameters.in_type == "gene":
        df = filter_gene_names(data=parameters.data, organism=parameters.organism, decoy=parameters.decoy,
                               id_column=parameters.gene_column, protein_column=parameters.protein_column,
                               action=parameters.action, reviewed=parameters.reviewed)
    else:
        df = parameters.data
    df.to_csv(parameters.out_dir + parameters.file_name + "_filtered.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
