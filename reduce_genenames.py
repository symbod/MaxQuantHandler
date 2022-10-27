#!/usr/bin/python3

from gprofiler import GProfiler
import mygene
import pandas as pd

from mq_utils.HGNC_mapping import get_HGNC_mapping
from mq_utils import mapping_handler as mh

def reduce_genenames(data, mode, gene_column: str = "Gene names",
                     keep_empty=False, organism=None, HGNC_mode="mostfrequent"):
    """
    Reduce gene names in MaxQuant file based on chosen mode.

    :param data: MaxQuant data
    :param mode: Mode on how to reduce gene names
    :param gene_column:
    :param keep_empty: Set True if rows with no gene names should be kept
    :param organism: Organism to map to
    :param HGNC_mode: Mode on how to select the gene names in HGNC (mostfrequent, all)
    :return: Remapped MaxQuant file as dataframe
    """

    handler = mh.MappingHandler(mapping_dir="mappings/")

    # ==== Preload info for all IDs ====
    handler.get_mapping(ids=";".join(data[gene_column]).split(";"),
                                     in_type="reduced_genes", organism=organism, reduction_mode=mode)

    # ==== If Input was single column file ====
    if gene_column not in data.columns:
        print("Gene Column Not in Data Column!")
        raise SystemExit


    # ==== If Mode == HGNC check if organism = Human
    if (mode == "HGNC") and (organism != "human"):
        print("HGNC Database only for Human Genes!")
        raise SystemExit

    # === If Mode == Ensembl or enrichment --> check that organism is not None
    if (mode in ["ensembl", "enrichment"]) and (organism is None):
        print("For ensembl and enrichment modes, an organism is required!")
        raise SystemExit

    # === Reduce Gene Names ====
    data[gene_column] = data[gene_column].apply(lambda row: handler.get_reduced_genenames(
                                            ids = row.split(";"),
                                            reduction_mode = mode, HGNC_mode = HGNC_mode,
                                            organism = organism)
                                            )

    # === Remove Rows with Empty Gene Names ====
    if keep_empty is False:
        data = data.drop(data[data[gene_column] == ""].index)

    handler.save_mappings(mapping_dir="mappings/")
    return data
