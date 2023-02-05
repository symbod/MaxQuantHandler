#!/usr/bin/python3

import csv
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru
from mq_utils.logger import get_remapped_genenames_logging
from fasta_grepper import grep_header_info
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')


def remap_genenames(data: pd.DataFrame, mode: str, protein_column: str, gene_column: str = None,
                    skip_filled: bool = False, organism: str = None, fasta: str = None, keep_empty: bool = True,
                    res_column: str = None):
    """
    Remap gene names in data file based on chosen mode.

    :param data: Dataframe containing a column with protein IDs and optionally gene names
    :param mode: Mode on how to map gene names
    :param protein_column: Column name with protein IDs
    :param gene_column: Column name with gene names
    :param skip_filled: Set True if rows with already filled gene names should be ignored
    :param organism: Organism to map to
    :param fasta: Fasta file
    :param keep_empty: Set True if empty rows should be kept
    :param res_column: Set column name for remap genenames results. If None, the gene_column will be overridden.
    :return: Remapped data as dataframe
    """
    data_copy = data.copy(deep=True)
    data_copy = data_copy.fillna("")
    data_copy[protein_column] = data_copy[protein_column].astype("string")

    handler = mh.MappingHandler()
    # ==== Preload info for all IDs ====
    handler.get_mapping(ids=";".join(data_copy[protein_column]).split(";"),
                        in_type="protein", organism=organism)

    # ==== If gene_column does not exist in given dataframe ====
    temp_gene_column = "Remapped Gene Names" if gene_column is None else gene_column
    data_copy[temp_gene_column] = "" if temp_gene_column not in data_copy.columns else data_copy[temp_gene_column].fillna("")

    # ==== Get fasta mapping ====
    if fasta is not None and mode in ['all', 'fasta']:
        fasta_mapping = grep_header_info(fasta=parameters.mapping_file)
        remapped_gene_names = data_copy.apply(
            lambda row: get_fasta_mapping(ids=row[protein_column].split(";"), genename=row[temp_gene_column],
                                          mapping=fasta_mapping, skip_filled=skip_filled), axis=1)
        skip_filled = True
    else:
        remapped_gene_names = data_copy[temp_gene_column]

    # ==== Get uniprot mappings ====
    if mode != 'fasta':
        data_copy["temp"] = remapped_gene_names  # copy into df for apply function
        remapped_gene_names = data_copy.apply(
            lambda row: get_uniprot_mapping(ids=row[protein_column].split(";"), genename=row["temp"],
                                            mode=mode, organism=organism, handler=handler,
                                            skip_filled=skip_filled), axis=1)
        del data_copy["temp"]

    # ==== Logging ====
    log_dict = get_remapped_genenames_logging(original=data_copy[temp_gene_column], remapped=remapped_gene_names,
                                              protein_ids=";".join(data_copy[protein_column]).split(";"),
                                              handler=handler, organism=organism)

    # ==== If target column depending if res_column is set ====
    column = res_column if res_column is not None else temp_gene_column

    # ==== Set remapped gene names to dataFrame ====
    data_copy[column] = remapped_gene_names

    # Remove temporary gene column
    if gene_column is None and res_column is not None:
        data_copy = data_copy.drop(temp_gene_column, axis=1)

    # ==== Remove rows with empty gene names ====
    if keep_empty is False:
        data_copy = data_copy[data_copy[column] != ""]  # remove

    return data_copy, log_dict


def get_fasta_mapping(ids, genename, mapping=None, skip_filled=False):
    """
    Get gene names from fasta file for empty entries or all if skip_filles is set to false.

    :param ids: List of protein ids
    :param genename: Mapped gene name
    :param mapping: Mapping from fasta file
    :param skip_filled: Set True if skip mapping when genename is not empty
    :return: Gene name
    """
    if genename == "" or not skip_filled:
        symbols = set(mapping[mapping["uniprot"].isin(ids)]["symbol"].dropna())
        return ";".join(list(symbols))
    else:
        return genename


def get_uniprot_mapping(ids, genename, mode, handler, organism=None, skip_filled=False):
    """
    Get gene names from uniprot for empty entries or all if skip_filles is set to false.

    :param ids: List of protein ids
    :param genename: Mapped gene name
    :param mode: Mode on how to map gene names
    :param handler: Handler for uniprot mappings
    :param organism: Organism to map to
    :param skip_filled: Set True if skip mapping when genename is not empty
    :return: Gene name
    """
    if genename == "" or not skip_filled:
        if mode == "uniprot_one":
            return get_single_genename(ids=ids, organism=organism, handler=handler)
        elif mode == "uniprot_primary":
            return get_primary_genenames(ids=ids, organism=organism, handler=handler)
        else:
            return get_all_genenames(ids=ids, organism=organism, handler=handler)
    else:
        return genename


def get_single_genename(ids, handler: mh.MappingHandler, organism=None):
    """
    Get most frequent gene name from uniprot.

    :param ids: List of protein IDs
    :param organism: Organism to map to
    :param handler: Handler for uniprot mappings
    :return: Single gene name
    """
    mapping = handler.get_mapping(ids=ids, in_type="protein", organism=organism)
    # ==== Get primary gene name first if only one for all ====
    prim_keys = set(mapping['Gene Names (primary)'].fillna(""))
    if len(prim_keys) == 1:
        return ";".join(prim_keys)
    # ==== Most frequent out of all ====
    mapping['Gene Names'] = mapping['Gene Names'].fillna("")
    gene_names = mapping['Gene Names'].apply(lambda x: set(x.split(";")))
    lst = [x for z in gene_names for x in z if x != ""]
    # ==== Return most frequent ====
    return max(set(lst), key=lst.count)


def get_primary_genenames(ids, handler: mh.MappingHandler, organism=None):
    """
    Get only primary gene names from uniprot.

    :param ids: List of protein IDs
    :param organism: Organism to map to
    :param handler: Handler for uniprot mappings
    :return: Primary gene names
    """
    mapping = handler.get_mapping(ids=ids, in_type="protein", organism=organism)
    if mapping.empty:
        return ""
    else:
        genenames = {x for x in mapping['Gene Names (primary)'].str.split(";").explode() if pd.notna(x)}  # set()
        return ';'.join(genenames)


def get_all_genenames(ids, handler: mh.MappingHandler, organism=None):
    """
    Get all gene names from uniprot.

    :param ids: List of protein IDs
    :param organism: Organism to map to
    :param handler: Handler for uniprot mappings
    :return: All gene names
    """
    mapping = handler.get_mapping(ids=ids, in_type="protein", organism=organism)
    if mapping.empty:
        return ""
    else:
        mapping['Gene Names'] = mapping['Gene Names'].fillna("")
        gene_names_series = mapping['Gene Names'].apply(lambda x: set(x.split(";")))
        genenames = set([x for y in gene_names_series for x in y if x != ""])
        return ';'.join(genenames)


if __name__ == "__main__":
    description = "                  Re-mapp gene names in data file."
    parameters = ru.save_parameters(script_desc=description,
                                    arguments=('d', 'm', 'pc_req', 'gc', 'l', 'or', 'f', 'ke', 'rc', 'o'))
    res, log = remap_genenames(data=parameters.data, mode=parameters.mode, protein_column=parameters.protein_column,
                               gene_column=parameters.gene_column, skip_filled=parameters.fill,
                               organism=parameters.organism, fasta=parameters.fasta_file,
                               keep_empty=parameters.keep_empty,
                               res_column=parameters.res_column)

    res.to_csv(parameters.out_dir + Path(parameters.file_name).stem + "_remapped.txt", header=True,
               index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")

    log["Overview_Log"].to_csv(parameters.out_dir + Path(parameters.file_name).stem + "_remapped_overview_log.txt",
                               header=True, index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")
    log["Detailed_Log"].to_csv(
        parameters.out_dir + Path(parameters.file_name).stem + "_remapped_detailed_log.txt",
        header=True, index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")
