#!/usr/bin/python3

import pandas as pd
import itertools
from mq_utils.logger import get_reduced_genenames_logging
from mq_utils import mapping_handler as mh, runner_utils as ru
from pathlib import Path
import csv


def reduce_genenames(data: pd.DataFrame, gene_column: str, mode:str, organism: str,
                     res_column: str = None, keep_empty: bool = True, HGNC_mode: str = "mostfrequent"):
    """
    Reduce gene names in data file based on chosen mode.

    :param data: Dataframe containing a column with gene names
    :param mode: Mode on how to reduce gene names
    :param gene_column: Column name with gene names
    :param res_column: Set column name for reduced results. If None, the gene_column will be overridden.
    :param keep_empty: Set True if rows with no gene names should be kept
    :param organism: Organism to map to
    :param HGNC_mode: Mode on how to select the gene names in HGNC (mostfrequent, all)
    :param return_log: Set True if log dataframes should be returned

    :return: Reduced data as dataframe
    """
    data_copy = data.copy(deep=True)
    data_copy[gene_column] = data_copy[gene_column].astype("string")

    handler = mh.MappingHandler(mapping_dir="mappings/")
    # ==== Preload info for all IDs ====
    handler.get_mapping(ids=";".join(data_copy[gene_column]).split(";"),
                        in_type="reduced_genes", organism=organism, reduction_mode=mode)

    # ==== If gene_column not in data ====
    if gene_column not in data_copy.columns:
        raise Exception("Gene Column Not in Data Column!")

    # ==== If Mode is HGNC check if organism is Human ====
    if (mode == "HGNC") and (organism != "human"):
        raise Exception("HGNC Database only for Human Genes!")

    # ==== Reduce Gene Names ====
    reduced_gene_names = data_copy[gene_column].apply(
        lambda row: get_reduced_genenames(ids=row.split(";"), handler=handler,
                                          reduction_mode=mode, HGNC_mode=HGNC_mode, organism=organism))

    # ==== Logging ====
    log_dict = get_reduced_genenames_logging(data_copy[gene_column], reduced_gene_names, handler, organism, mode)

    # ==== If target column depending if res_column is set ====
    column = res_column if res_column is not None else gene_column

    # ==== Set Reduced Gene Names To DataFrame ====
    data_copy[column] = reduced_gene_names

    # ==== Remove Rows with Empty Gene Names ====
    if keep_empty is False:
        data_copy = data_copy[data_copy[column] != ""]  # remove

    # ==== Save Current Mappings To Files ====
    handler.save_mappings(mapping_dir="mappings/")

    return data_copy, log_dict


def get_reduced_genenames(ids, handler, organism=None, reduction_mode="ensembl", HGNC_mode="mostfrequent"):
    mapping = handler.get_mapping(ids=ids, in_type="reduced_genes", organism=organism, reduction_mode=reduction_mode)
    if mapping.empty:
        return ""
    else:
        # check orga
        if reduction_mode == "HGNC":
            # separate case because we have two modes (mostfrequent and all)
            HGNC_output = mapping[mapping["Mode"] == reduction_mode].dropna()
            HGNC_list = list(itertools.chain.from_iterable(list(HGNC_output["Reduced Gene Name"])))
            if HGNC_mode == "mostfrequent":
                reduced_genenames = list(pd.Series(HGNC_list).mode())
            elif HGNC_mode == "all":
                reduced_genenames = HGNC_list
        else:
            # remove None
            reduced_genenames = list(mapping[-mapping["Reduced Gene Name"].isin(["None", None])][
                                         "Reduced Gene Name"])  # "None" and None because Ensembl returns for example "None"
        return ";".join(list(set(reduced_genenames)))


if __name__ == "__main__":
    description = "                  Reduce gene names in data file."
    parameters = ru.save_parameters(script_desc=description,
                                    arguments=('d', 'gc_req', 'rm', 'or_req', 'rc', 'ke', 'hm', 'o'))
    res, log = reduce_genenames(data=parameters.data, gene_column=parameters.gene_column, mode=parameters.reduce_mode,
                                organism=parameters.organism, res_column = parameters.res_column,
                                keep_empty=parameters.keep_empty, HGNC_mode= parameters.hgnc_mode)
    res.to_csv(parameters.out_dir + Path(parameters.file_name).stem + "_reduced.txt", header=True,
               index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")

    log["Overview_Log"].to_csv(parameters.out_dir + Path(parameters.file_name).stem + "_reduced_overview_log.txt",
                               header = True, index = False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")
    log["Detailed_Log" ].to_csv(
        parameters.out_dir + Path( parameters.file_name ).stem + "_reduced_detailed_log.txt",
        header=True, index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" " )
