#!/usr/bin/python3

import csv
import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru
from fasta_grepper import grep_header_info
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

full_mapping = pd.DataFrame(columns=['Gene Names', 'Gene Names (primary)', 'Reviewed', 'Organism', 'Protein ID'])


def remap_genenames(data, mode, protein_column: str, gene_column:str = "Gene names",
                    skip_filled=False, organism=None, fasta=None):
    """
    Remap gene names in MaxQuant file based on chosen mode.

    :param data: MaxQuant data
    :param mode: Mode on how to map gene names
    :param protein_column:
    :param gene_column:
    :param skip_filled: Set True if rows with already filled gene names should be ignored
    :param organism: Organism to map to
    :param fasta: Fasta file
    :return: Remapped MayQuant file as dataframe
    """
    handler = mh.MappingHandler(mapping_dir="mappings/")
    # ==== Preload info for all IDs ====
    handler.get_mapping(ids=";".join(data[protein_column]).split(";"),
                        in_type="protein", organism=organism)

    # ==== If Input was single column file ====
    if gene_column not in data.columns:
        data[gene_column] = ""
    # ==== Get fasta mapping ====
    if fasta is not None and mode in ['all', 'fasta']:
        fasta_mapping = grep_header_info(fasta=parameters.mapping_file)
        data[gene_column] = data.apply(
            lambda row: run_fasta_mapping(ids=row[protein_column].split(";"), genename=row[gene_column],
                                          mapping=fasta_mapping, skip_filled=skip_filled), axis=1)
        skip_filled = True

    # ==== Get uniprot mappings ====
    if mode != 'fasta':
        data[gene_column] = data.apply(
            lambda row: run_uniprot_mapping(ids=row[protein_column].split(";"), genename=row[gene_column],
                                            mode=mode, organism=organism, handler=handler,
                                            skip_filled=skip_filled), axis=1)
    handler.save_mappings(mapping_dir="mappings/")
    return data


def run_fasta_mapping(ids, genename, mapping=None, skip_filled=False):
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
        return ";". join(list(symbols))
    else:
        return genename


def run_uniprot_mapping(ids, genename, mode, handler, organism=None, skip_filled=False):
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
            return handler.get_primary_genenames(ids=ids, organism=organism)
        else:
            return handler.get_all_genenames(ids=ids, organism=organism)
    else:
        return genename


def get_single_genename(ids, organism=None, handler: mh.MappingHandler = mh.MappingHandler(mapping_dir="mappings/")):
    """
    Get gene name from uniprot.

    :param ids: List of protein ids
    :param organism: Organism to map to
    :param handler: Handler for uniprot mappings
    :return: Single gene name
    """
    df = handler.get_mapping(ids=ids, in_type="protein", organism=organism)
    # ==== Get primary gene name first ====
    prim_keys = set(df['Gene Names (primary)'].fillna(""))
    if len(prim_keys) == 1:
        return ";".join(prim_keys)
    # ==== Check all gene names ====
    df['Gene Names'] = df['Gene Names'].fillna("").str.upper()
    gene_names = df['Gene Names'].apply(mh.series_to_set)
    lst = [x for z in gene_names for x in z if x != ""]
    # ==== Return most frequent ====
    return max(set(lst), key=lst.count)


if __name__ == "__main__":
    description = "                  Re-mapp gene names in max quant file."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'f', 'c', 'or', 'l', 'm', 'o'))
    res = remap_genenames(data=parameters.data, mode=parameters.mode, skip_filled=parameters.fill,
                          protein_column=parameters.protein_column, gene_column=parameters.gene_column,
                          organism=parameters.organism, fasta=parameters.fasta_file)
    res.to_csv(parameters.out_dir + Path(parameters.file_name).stem + "_remapped.txt", header=True,
               index=False, quoting=csv.QUOTE_NONNUMERIC, sep=" ")
