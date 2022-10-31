#!/usr/bin/python3


from mq_utils import mapping_handler as mh
from mq_utils.logger import get_reduced_genenames_logging


def reduce_genenames(data, mode, gene_column: str = "Gene names",
                     keep_empty: bool = False, organism: str = None, HGNC_mode: str = "mostfrequent",
                     return_log: bool = True):
    """
    Reduce gene names in MaxQuant file based on chosen mode.

    :param data: MaxQuant data
    :param mode: Mode on how to reduce gene names
    :param gene_column:
    :param keep_empty: Set True if rows with no gene names should be kept
    :param organism: Organism to map to
    :param HGNC_mode: Mode on how to select the gene names in HGNC (mostfrequent, all)
    :param return_log: Set True if a log dataframe should be returned
    :return: Remapped MaxQuant file as dataframe
    """

    data = data.copy( deep=True )

    # ==== Preload info for all IDs ====
    handler = mh.MappingHandler(mapping_dir="mappings/")
    handler.get_mapping(ids=";".join(data[gene_column]).split(";"),
                                     in_type="reduced_genes", organism=organism, reduction_mode=mode)

    # ==== If Input was single column file ====
    if gene_column not in data.columns:
        print("Gene Column Not in Data Column!")
        raise SystemExit

    # ==== If Organism is not set ====
    if organism == None:
        print("Organism is required!")
        raise SystemExit


    # ==== If Mode == HGNC check if organism = Human
    if (mode == "HGNC") and (organism != "human"):
        print("HGNC Database only for Human Genes!")
        raise SystemExit


    # ==== Reduce Gene Names ====
    reduced_gene_names = data[gene_column].apply(lambda row: handler.get_reduced_genenames(
                                            ids = row.split(";"),
                                            reduction_mode = mode, HGNC_mode = HGNC_mode,
                                            organism = organism)
                                            )

    # ==== Logging ====
    if return_log:
        log_dict = get_reduced_genenames_logging(data[gene_column], reduced_gene_names)

    # ==== Remove Rows with Empty Gene Names ====
    if keep_empty is False:
        data = data.drop( data[ data[ gene_column ] == "" ].index )

    # ==== Save Current Mappings To Files ====
    handler.save_mappings( mapping_dir="mappings/" )

    # ==== Set Reduced Gene Names To DataFrame ====
    data[ gene_column ] = reduced_gene_names

    if return_log:
        return data, log_dict
    else:
        return data


