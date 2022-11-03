#!/usr/bin/python3

import pandas as pd
import itertools
from mq_utils import mapping_handler as mh
from mq_utils.logger import get_reduced_genenames_logging


def reduce_genenames(data: pd.DataFrame, mode, gene_column: str, organism: str,
                     res_column: str = None, keep_empty: bool = False, HGNC_mode: str = "mostfrequent",
                     return_log: bool = True):
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
    log_dict = dict()
    if return_log:
        log_dict = get_reduced_genenames_logging(data_copy[gene_column], reduced_gene_names)

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
