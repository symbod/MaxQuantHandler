#!/usr/bin/python3

import csv
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru
from mq_utils.logger import get_filter_ids_logging
from pathlib import Path


def filter_protein_ids(data: pd.DataFrame, protein_column: str, organism: str = None,
                       rev_con: bool = False, keep_empty: bool = True,
                       reviewed: bool = True, res_column: str = None):
    """
    Filter protein ids in given data by chosen organism.

    :param data: Dataframe containing a column with protein IDs
    :param protein_column: Column name with protein IDs
    :param organism: Organism to map to
    :param rev_con: Bool to indicate if decoy IDs (REV_) and contaminants (CON_) should be kept
    :param keep_empty: Set True if empty rows should be kept
    :param reviewed: Set True if only reviewed protein IDs should be kept
    :param res_column: Set column name for remap genenames results. If None, the gene_column will be overridden.
    :return: Filtered data as dataframe
    """
    data_copy = data.copy(deep=True)
    data_copy = data_copy.fillna("")
    data_copy[protein_column] = data_copy[protein_column].astype("string")

    handler = mh.MappingHandler()
    # ==== Get all existing mappings in one batch ====
    handler.get_mapping(ids=";".join(data_copy[protein_column]).split(";"),
                        in_type="protein", organism=organism)

    # ==== Filter row wise ====
    filtered_ids = data_copy[protein_column].apply(
        lambda x: get_filtered_ids(ids=x.split(";"), handler=handler, organism=organism, rev_con=rev_con,
                                   reviewed=reviewed))

    # ==== Logging ====
    log_dict = get_filter_ids_logging(original=data_copy[protein_column], filtered=filtered_ids, handler=handler,
                                      organism=organism)

    # ==== If target column depending if res_column is set ====
    column = res_column if res_column is not None else protein_column

    # ==== Set filtered ids to dataframe ====
    data_copy[column] = filtered_ids

    # ==== Save current mapping to files
    # handler.save_mappings(mapping_dir="mappings/")

    # ==== Remove rows with empty protein IDs ====
    if keep_empty is False:
        data_copy = data_copy[data_copy[column] != ""]  # remove

    return data_copy, log_dict


def get_filtered_ids(ids, handler: mh.MappingHandler, organism: str = None, rev_con: bool = False,
                     reviewed: bool = False) -> str:
    """
    Filter given set of protein ids based on organism, decoy/contaminants and/or review status.

    :param ids: Set of protein IDs
    :param handler: MappingHandler object
    :param organism: Organism the IDs should belong to
    :param rev_con: Bool to indicate if decoy and contaminant IDs should be kept
    :param reviewed: Bool to indicate if only reviewed IDs should be kept
    :return: filtered IDs combined into a string
    """
    # ==== Get mapping on protein IDs ====
    mapping = handler.get_mapping(ids=ids, in_type="protein", organism=organism, ignore_missing=True)
    if mapping.empty:
        return ""

    # ==== Keep or remove decoy/contaminants IDs based on flag ====
    keep = set([x for x in ids if x.startswith(("REV", "CON"))]) if rev_con else set()

    # ==== Keep only reviewed IDs based on flag ====
    if reviewed:
        mapping = mapping[mapping['Reviewed'] == "reviewed"]

    # ==== Combine mapped IDs with kept or left decoy IDs ====
    prot_ids = set(mapping['Protein ID']).union(keep)
    return ';'.join(prot_ids)


if __name__ == "__main__":
    description = "        Filter proteins by organism and/or decoy/contaminants names in data file."
    parameters = ru.save_parameters(script_desc=description,
                                    arguments=('d', 'pc_req', 'or', 'rv', 'ke', 'r', 'rc', 'o'))
    df, log = filter_protein_ids(data=parameters.data, protein_column=parameters.protein_column,
                                 organism=parameters.organism, rev_con=parameters.rev_con,
                                 keep_empty=parameters.keep_empty, reviewed=parameters.reviewed,
                                 res_column=parameters.res_column)
    df.to_csv(parameters.out_dir + Path(parameters.file_name).stem + "_filtered.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")

    log["Overview_Log"].to_csv(
        parameters.out_dir + Path(parameters.file_name).stem + "_filtered_overview_log.txt",
        header=True, index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")
    log["Detailed_Log"].to_csv(
        parameters.out_dir + Path(parameters.file_name).stem + "_filtered_detailed_log.txt",
        header=True, index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")
