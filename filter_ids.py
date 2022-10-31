#!/usr/bin/python3

import csv
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru
from mq_utils.logger import get_filter_ids_logging


def filter_protein_ids(data: pd.DataFrame, id_column: str = "Protein IDs", organism: str = None,
                       decoy: bool = False, action: str = "delete",
                       reviewed: bool = True, return_log: bool = True):
    """
    Filter protein ids in given data by chosen organism.

    :param data: Dataframe containing a column with protein IDs
    :param id_column: Column name with protein IDs
    :param organism: Organism to map to
    :param decoy: Bool to indicate if decoy IDs (REV_, ) should be kept
    :param action: What to do, if IDs cell is empty after filtering. Keep empty cell or delete it.
    :param reviewed: Set True if only reviewed protein IDs should be kept
    :param return_log: Set True if a log dataframe should be returned
    :return: Filtered data as dataframe
    """

    data = data.copy(deep=True)

    # ==== Get all existing mappings in one batch ====
    handler = mh.MappingHandler(mapping_dir="mappings/")
    handler.get_mapping(ids=";".join(data[id_column]).split(";"),
                        in_type="protein", organism=organism)

    # ==== Filter row wise ====
    filtered_ids = data[id_column].apply(
        lambda x: handler.get_filtered_ids(ids=x.split(";"), organism=organism, decoy=decoy, reviewed=reviewed))

    # ==== Logging ====
    if return_log:
        log_dict = get_filter_ids_logging(original=data[id_column], filtered=filtered_ids, handler=handler, organism=organism)

    # ==== Set filtered ids to dataframe ====
    data[id_column] = filtered_ids

    # ==== Save current mapping to files
    handler.save_mappings(mapping_dir="mappings/")

    # ==== Remove empty rows if flag is set
    if action == "delete":
        data = data[data[id_column] != ""]  # remove

    if return_log:
        return data, log_dict
    else:
        return data


if __name__ == "__main__":
    description = "        Filter proteins by organism and/or decoy names."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'or_req', 'c', 'r', 'd', 'a', 'o'))
    df, log = filter_protein_ids(data=parameters.data, organism=parameters.organism, decoy=parameters.decoy,
                                 id_column=parameters.protein_column,
                                 action=parameters.action, reviewed=parameters.reviewed)
    df.to_csv(parameters.out_dir + parameters.file_name + "_filtered.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
