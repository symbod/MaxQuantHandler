#!/usr/bin/python3

import csv
from pathlib import Path
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru
from mq_utils.logger import get_ortholog_genenames_logging


def map_orthologs(data: pd.DataFrame, gene_column: str, organism: str, tar_organism: str,
                  keep_empty: bool = True, res_column: str = None):
    """
    Map gene names of origin organism to orthologs of target organism.

    :param data: Dataframe containing a column with gene names
    :param gene_column: Column name with gene names
    :param organism: Organism of the input ids
    :param tar_organism: Organism to map to
    :param keep_empty: Set True if empty rows should be kept
    :param res_column: Set column name for ortholog results. If None, the gene_column will be overridden.
    :return: Data as dataframe with ortholog ids
    """
    data_copy = data.copy(deep=True)
    data_copy[gene_column] = data_copy[gene_column].astype("string")

    handler = mh.MappingHandler(mapping_dir="mappings/")
    # ==== Get all existing mappings in one batch ====
    handler.get_mapping(ids=";".join(data_copy[gene_column]).split(";"),
                        in_type="orthologs", organism=organism, tar_organism=tar_organism)
    handler.save_mappings(mapping_dir="mappings/")

    ortholog_gene_names = data_copy[gene_column].apply(
        lambda x: get_orthologs(ids=x.split(";"), handler=handler, organism=organism, tar_organism=tar_organism))

    # ==== Logging ====
    log_dict = get_ortholog_genenames_logging(original=data_copy[gene_column], orthologs=ortholog_gene_names,
                                                  handler=handler, organism=organism, tar_organism=tar_organism)

    # ==== If target column depending if res_column is set ====
    column = res_column if res_column is not None else gene_column

    # ==== Set ortholog gene names to dataframe ====
    data_copy[column] = ortholog_gene_names

    # ==== Remove rows with empty ortholog gene names ====
    if keep_empty is False:
        data_copy = data_copy[data_copy[column] != ""]  # remove

    return data_copy, log_dict


def get_orthologs(ids, handler, organism: str, tar_organism: str):
    """
    Get orthologs of genes from one organism to another.

    :param ids: Set of gene names
    :param handler: Handler for mappings
    :param organism: Organism of the input ids
    :param tar_organism: Organism to map to

    :return:
    """
    mapping = handler.get_mapping(ids=ids, in_type="orthologs", organism=organism,
                                  tar_organism=tar_organism, ignore_missing=True)
    if mapping.empty:
        return ""
    else:
        orthologs = {x for x in mapping['target_symbol'] if pd.notna(x)}
        return ';'.join(orthologs)


if __name__ == "__main__":
    description = "                       Map ortholog gene names in data file."
    parameters = ru.save_parameters(script_desc=description,
                                    arguments=('d','gc_req','or_req','tor_req','ke','rc','o'))
    df, log = map_orthologs(data=parameters.data, gene_column=parameters.gene_column,organism=parameters.organism,
                            tar_organism=parameters.tar_organism, keep_empty=parameters.keep_empty,
                            res_column=parameters.res_column)

    df.to_csv(parameters.out_dir + Path( parameters.file_name ).stem + "_ortholog.txt", header=True, index=False,
              quoting=csv.QUOTE_NONNUMERIC, sep=" ")
    log[ "Overview_Log" ].to_csv(
        parameters.out_dir + Path( parameters.file_name ).stem + "_ortholog_overview_log.txt",
        header=True, index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" " )
    log[ "Detailed_Log" ].to_csv(
        parameters.out_dir + Path( parameters.file_name ).stem + "_ortholog_detailed_log.txt",
        header=True, index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" " )

