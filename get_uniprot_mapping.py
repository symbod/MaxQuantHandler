#!/usr/bin/python3

import pandas as pd
from mq_utils import mapping_handler as mh, runner_utils as ru


def get_uniprot_mappings(data: pd.DataFrame, in_type: str, organism=None):
    """
    Get uniprot mapping to protein ids or gene names in given data. Optionally only mapped to chosen organism.

    :param data: MaxQuant data or data with one single column
    :param in_type: Type which column should be taken as a reference
    :param organism: Organism to map to
    :return: Uniprot mapping to protein ids or gene names as dataframe
    """
    handler = mh.MappingHandler(mapping_dir="mappings/")
    if in_type == "protein":
        mappings = handler.get_mapping(ids=";".join(data['Protein IDs']).split(";"),
                                       in_type=in_type, organism=organism)
    else:  # in_type == "gene":
        print("Gene names option deferred.")
        def agg_func(x):
            return ';'.join(set(x))

        mappings = handler.get_mapping(ids=";".join(data['Gene names']).split(";"),
                                       in_type=in_type, organism=organism)
        mappings = mappings.groupby(['Gene name', 'Status']).agg({'Protein ID': agg_func,
                                                                  'Organism': agg_func}).reset_index()
    handler.save_mappings()
    return mappings


if __name__ == "__main__":
    description = "  Get uniprot mapping to protein ids or gene names optionally by organism."
    parameters = ru.save_parameters(script_desc=description, arguments=('qf', 'or_req', 'i', 'o'))
    df = get_uniprot_mappings(data=parameters.data, in_type=parameters.in_type,
                              organism=parameters.organism)
    df.to_csv(parameters.out_dir + parameters.file_name + "_uniprot_" + parameters.in_type + "_mapping.csv",
              header=True, index=False)
