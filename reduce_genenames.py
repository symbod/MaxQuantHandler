#!/usr/bin/python3


from mq_utils import mapping_handler as mh
from mq_utils.logger import get_reduced_genenames_logging

def reduce_genenames(data, mode, gene_column: str = "Gene names",
                     keep_empty=False, organism=None, HGNC_mode="mostfrequent",
                     inplace=True, return_log=True):
    """
    Reduce gene names in MaxQuant file based on chosen mode.

    :param data: MaxQuant data
    :param mode: Mode on how to reduce gene names
    :param gene_column:
    :param keep_empty: Set True if rows with no gene names should be kept
    :param organism: Organism to map to
    :param HGNC_mode: Mode on how to select the gene names in HGNC (mostfrequent, all)
    :param inplace: Set True if reduction should be applied on original dataframe (False for returning a new dataframe)
    :param return_log: Set True if a log dataframe should be returned
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

    if organism == None:
        print("Organism is required!")
        raise SystemExit


    # ==== If Mode == HGNC check if organism = Human
    if (mode == "HGNC") and (organism != "human"):
        print("HGNC Database only for Human Genes!")
        raise SystemExit


    # === Reduce Gene Names ====
    reduced_gene_names = data[gene_column].apply(lambda row: handler.get_reduced_genenames(
                                            ids = row.split(";"),
                                            reduction_mode = mode, HGNC_mode = HGNC_mode,
                                            organism = organism)
                                            )

    # === Logging ===
    log_df = get_reduced_genenames_logging(data[gene_column], reduced_gene_names)

    # === Remove Rows with Empty Gene Names ====
    if keep_empty is False:
        data = data.drop( data[ data[ gene_column ] == "" ].index )

    handler.save_mappings( mapping_dir="mappings/" )

    # === Set Reduced Gene Names To DataFrame
    if inplace:
        data[gene_column] = reduced_gene_names
        if return_log:
            return data, log_df
        else:
            return data
    else:
        reduced_data = data.copy(deep=True)
        reduced_data[gene_column] = reduced_gene_names
        if return_log:
            return reduced_data, log_df
        return reduced_data


