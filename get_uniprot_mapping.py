#!/usr/bin/python3

import pandas as pd
import runner_utils as ru
import uniprot_handler as uh


def get_uniprot_mappings(mq_file, in_type, organism=None):
    """
    Get uniprot mapping to protein ids or gene names in MaxQuant file. Optionally only mapped to chosen organism.

    :param mq_file: MaxQuant file
    :param in_type: Type which column should be taken as a reference
    :param organism: Organism to map to
    :return: Uniprot mapping to protein ids or gene names as dataframe
    """
    handler = uh.UniprotHandler()
    max_quant = pd.read_table(mq_file, sep=" ").fillna("")
    if in_type == "proteinID":
        mappings = pd.DataFrame(columns=['Gene names', 'Gene names  (primary )', 'Status', 'Organism', 'Protein ID'])
    else:  # if in_type == "genename":
        mappings = pd.DataFrame(columns=['Protein ID', 'Status', 'Organism', 'Gene names'])
    for index, row in max_quant.iterrows():
        if in_type == "proteinID":
            mappings = pd.concat(
                [mappings, handler.get_mapping(ids=row['Protein IDs'].split(";"), in_type=in_type, organism=organism)])
        elif in_type == "genename":
            def agg_func(x):
                return ';'.join(set(x))

            mappings = pd.concat(
                [mappings, handler.get_mapping(ids=row['Gene names'].split(";"), in_type=in_type, organism=organism)])
            mappings= mappings.groupby('Gene name').agg({'Protein ID': agg_func, 'Status': agg_func,
                                                         'Organism': agg_func}).reset_index()
    handler.save_mappings()
    return mappings


if __name__ == "__main__":
    description = "  Get uniprot mapping to protein ids or gene names optionally by organism."
    parameters = ru.save_parameters(script_desc=description, arguments=('q', 'r_req', 'i', 'o'))
    df = get_uniprot_mappings(mq_file=parameters.maxquant_file, in_type=parameters.in_type,
                              organism=parameters.organism)
    df.to_csv(parameters.out_dir + "uniprot_" + parameters.in_type + "_mapping.csv", header=True, index=False)
