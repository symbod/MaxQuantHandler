#!/usr/bin/python3

import re
import pandas as pd
import runner_utils as ru


def grep_header_info(fasta: str) -> pd.DataFrame:
    """
    Grep information from headers in given fasta file.

    :param fasta: Fasta file
    :return: Grepped information as dataframe with columns := 'uniprot','description','symbol','organism','ox','pe','sv'
    """
    infos = list()
    with open(fasta) as file:
        for line in file:
            if line.startswith('>'):
                # get gene name
                matches = re.search("GN=(.*)\sPE", line)
                gene_name = matches.group(1) if matches else ""
                gene_string = "GN=.*\s" if matches else ""
                # get full info
                matches = re.search("\|(.*)\|\S*\s(.*)\sOS=(.*)\sOX=(.*)\s"+gene_string+"PE=(.*)\sSV=(.*)", line)
                infos.append([matches.group(1), matches.group(2), gene_name, matches.group(3),
                              matches.group(4), matches.group(5), matches.group(6)])
    return(pd.DataFrame(infos, columns=['uniprot','description','symbol','organism','ox','pe','sv']))


if __name__ == "__main__":
    description = "                  Grep protein info from fasta file."
    parameters = ru.save_parameters(script_desc=description, arguments=('f_req','o'))
    df = grep_header_info(fasta=parameters.fasta_file)
    df.to_csv(parameters.out_dir + "grepped_info.csv", header=True, index=False)