#!/usr/bin/python3

from os.path import dirname, abspath, join
import pandas as pd
from pathlib import Path
from gprofiler import GProfiler
import requests
from .HGNC_mapping import get_HGNC_mapping
import mygene
import numpy as np

here = dirname(dirname(abspath(__file__)))  # needed to make files reachable in python package


class MappingHandler:
    full_protein_mapping = pd.DataFrame(columns=['Gene Names', 'Gene Names (primary)', 'Reviewed', 'Organism',
                                                 'Protein ID'])
    full_ortholog_mapping = pd.DataFrame(columns=['source_symbol', 'source_organism', 'ensg', 'ortholog_ensg',
                                                  'target_symbol', 'target_organism', 'description'])
    full_reduced_gene_mapping = pd.DataFrame(columns=["Gene Name", "Reduced Gene Name", "Organism", "Mode"])

    def __init__(self, mapping_dir):
        mapping_dir = join(here, mapping_dir)
        if Path(mapping_dir + "protein_to_genenames.csv").exists():
            self.full_protein_mapping = pd.read_csv(mapping_dir + "protein_to_genenames.csv")
        if Path(mapping_dir + "genenames_to_orthologs.csv").exists():
            self.full_ortholog_mapping = pd.read_csv(mapping_dir + "genenames_to_orthologs.csv")
        if Path(mapping_dir + "genenames_to_reduced_genenames.csv").exists():
            self.full_reduced_gene_mapping = pd.read_csv(mapping_dir + "genenames_to_reduced_genenames.csv",
                                                         na_values=None)

    # === UniProt Mapping ====
    def get_uniprot_mapping(self, ids, organism: str = None):
        """
        Get UniProt mapping for UniProt protein IDs. Optionally filter for organism.

        :param ids: Set of protein IDs
        :param organism: Organism to map to
        :return: dataframe with mapping to each mappable ID
        """
        organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)", "mouse": "Mus musculus (Mouse)",
                     "rabbit": "Oryctolagus cuniculus (Rabbit)"}
        url = 'https://rest.uniprot.org/uniprotkb/accessions'
        mapping = pd.DataFrame()
        # ==== Get mappings in <= 500 IDs chunks ====
        for i in range(0, len(ids), 500):
            # ==== Remove contaminated and reverse mapped IDs ====
            ids_chunk = [x for x in ids[i:i + 500] if not x.startswith(("REV", "CON"))]
            params = {'format': 'tsv',
                      'accessions': ",".join(ids_chunk),
                      'fields': 'gene_names,gene_primary,reviewed,organism_name,accession'}
            f = requests.get(url=url, params=params)
            # ==== If at least one ID doesn't exist ====
            if f.status_code == 400:
                mapping_chunk = pd.DataFrame()
                # ==== Check each id separately ====
                for id in ids_chunk:
                    params["accessions"] = id
                    f = requests.get(url=url, params=params)
                    if f.status_code != 400:
                        mapping_chunk2 = pd.read_csv(f.url, sep="\t")
                        mapping_chunk = mapping_chunk2 if mapping_chunk.empty \
                            else pd.concat([mapping_chunk, mapping_chunk2])
            # ==== All IDs were mapped ====
            else:
                mapping_chunk = pd.read_csv(f.url, sep="\t")
            # ==== Combine to one mapping dataframe ====
            mapping = mapping_chunk if mapping.empty else pd.concat([mapping, mapping_chunk])
        # ==== Changes inside final mapping dataframe ====
        if not mapping.empty:
            mapping.columns = [*mapping.columns[:-1], 'Protein ID']  # change name of last column
            # ==== Clean accidental white spaces from UniProt ====
            mapping['Gene Names'] = mapping['Gene Names'].str.replace(r';?\s+;?', ';', regex=True)  # change separation to ;
            mapping['Gene Names (primary)'] = mapping['Gene Names (primary)'].str.replace(r'\s+', '', regex=True)
            mapping['Protein ID'] = mapping['Protein ID'].str.replace(r'\s+', '', regex=True)
            # ==== Split and explode to create one row for each ID ====
            mapping['Protein ID'] = mapping['Protein ID'].apply(lambda x: x.split(","))
            mapping = mapping.explode('Protein ID')
            # ==== Save to global mapping ====
            self.full_protein_mapping = pd.concat([self.full_protein_mapping, mapping])
            # ==== Filter for organism if given ====
            if organism is not None:
                mapping = mapping[mapping['Organism'] == organisms[organism]]
        return mapping

    # === Ortholog Mapping ====
    def get_ortholog_mapping(self, ids, organism, tar_organism):
        """
        Get ortholog mapping from source to target organism using gProfiler.

        :param ids: Set of gene names
        :param organism: Organism of the input ids
        :param tar_organism: Organism to map to
        :return: dataframe with mapped ortholog pairs
        """
        organisms = {"human": "hsapiens", "mouse": "mmusculus", "rat": "rnorvegicus", "rabbit": "ocuniculus"}
        gp = GProfiler(return_dataframe=True)
        mapping = gp.orth(organism=organisms[organism], query=ids, target=organisms[tar_organism])
        mapping = mapping.fillna('')
        # ==== Subset and rename columns ====
        mapping = mapping[['incoming', 'converted', 'ortholog_ensg', 'name', 'description']]
        mapping.columns = ['source_symbol', 'ensg', 'ortholog_ensg', 'target_symbol', 'description']
        # ==== Save organism info ====
        mapping.insert(loc=1, column='source_organism', value=organism)
        mapping.insert(loc=5, column='target_organism', value=tar_organism)
        # ==== Save to global mapping ====
        self.full_ortholog_mapping = pd.concat([self.full_ortholog_mapping, mapping])
        return mapping

    # === Reduced Mapping ====
    def get_reduced_mapping(self, ids, organism, reduction_mode="ensembl"):
        if reduction_mode == "ensembl":
            mapping = self.get_ensembl_reduction(ids, organism)  # organism must be set
        elif reduction_mode == "HGNC":
            mapping = self.get_HGNC_reduction(ids)  # organism must be human
        elif reduction_mode == "enrichment":
            mapping = self.get_enrichment_reduction(ids, organism)  # organism must be set
        elif reduction_mode == "mygeneinfo":
            mapping = self.get_mygeneinfo_reduction(ids)  # add for all organisms directly
        if mapping is not None:
            mapping["Mode"] = reduction_mode
            self.full_reduced_gene_mapping = pd.concat([self.full_reduced_gene_mapping, mapping])
            return mapping
        else:
            return None

    def get_ensembl_reduction(self, ids, organism):  # organism required
        organisms = {"human": "hsapiens", "mouse": "mmusculus", "rat": "rnorvegicus", "rabbit": "ocuniculus"}
        gp = GProfiler(return_dataframe=True)
        gp_df = gp.convert(organism=organisms[organism], query=ids, target_namespace="ENSG")
        if len(gp_df) == 0:
            return None
        else:

            # Case 1: in ensembl id and name is None --> remove them
            gp_df = gp_df.drop( gp_df[ (gp_df[ "name" ] == "None") & (gp_df[ "converted" ] == "None") ].index )

            # Case 2: remove entries that have already an entry in gp_df with a corresponding name
            gp_df = gp_df.drop( gp_df[ (gp_df[ "n_converted" ] > 1) & (gp_df[ "name" ] == "None") ].index )

            # Case 3: in name is also the ensembl id saved --> save in new name the incoming name else take name
            gp_df[ "new_name" ] = np.where( gp_df[ "converted" ] == gp_df[ "name" ], gp_df[ "incoming" ],
                                            gp_df[ "name" ] )

            # Case 4: for entries with ensembl id but without name --> take input name
            gp_df[ "new_name" ] = np.where( (gp_df[ "converted" ] != "None") & (gp_df[ "name" ] == "None"),
                                            gp_df[ "incoming" ], gp_df[ "new_name" ] )

            mapping = gp_df[ [ "incoming", "new_name" ] ]

            # Case 5: get entries that have more than 1 entry in gProfiler
            duplicates = mapping[mapping.duplicated("incoming", keep=False)]

            if len(duplicates) > 0:
                # remove them from initial df
                mapping = mapping.drop(duplicates.index)
                # take input names of duplicates as new names
                mapping_chunk = pd.DataFrame({"incoming": duplicates["incoming"].unique(), "new_name": duplicates["incoming"].unique()})
                mapping = pd.concat([mapping, mapping_chunk])

            # if after removing nothing remains --> return empty dataframe
            if len(mapping) > 0:
                mapping.columns = ["Gene Name", "Reduced Gene Name"]
                mapping.loc[:, "Organism"] = organism
                return mapping
            else:
                return None

    def get_HGNC_reduction(self, ids):  # human organism required
        mapping_dict = {}
        for id in ids:
            alias_df = get_HGNC_mapping(id, "alias_symbol")
            symbol_df = get_HGNC_mapping(id, "symbol")
            if alias_df is None and symbol_df is None:
                # no information in HGNC for this id
                mapping_dict[id] = None
            else:
                if alias_df is None:
                    result_df = symbol_df
                elif symbol_df is None:
                    result_df = alias_df
                else:
                    result_df = pd.concat([alias_df, symbol_df])
                mapping_dict[id] = list(result_df["Symbol"])

        mapping = pd.DataFrame({"Gene Name": k, "Reduced Gene Name": v} for k, v in mapping_dict.items())
        mapping["Organism"] = "human"
        return mapping

    def get_enrichment_reduction(self, ids, organism):
        organisms = {"human": "hsapiens", "mouse": "mmusculus", "rat": "rnorvegicus", "rabbit": "ocuniculus"}
        gp = GProfiler(return_dataframe=True)
        # no_evidences = False returns intersections column with ids that match the specific annotation
        gp_df = gp.profile(organism=organisms[organism], query=ids, no_evidences=False)

        if len(gp_df) == 0:
            reduced_ids = [None for id in ids]
        else:
            enrichment_ids = gp_df["intersections"].explode().to_list()
            reduced_ids = [id if id in enrichment_ids else None for id in ids]
        mapping = pd.DataFrame({"Gene Name": ids, "Reduced Gene Name": reduced_ids})
        mapping["Organism"] = organism
        return mapping

    def get_mygeneinfo_reduction(self, ids):
        tax_ids = {"human": 9606, "mouse": 10090, "rat": 10116, "rabbit": 9986}
        inv_tax_ids = {tax_id: organism for organism, tax_id in tax_ids.items()}

        mg = mygene.MyGeneInfo()
        mg_output = mg.querymany(ids, scopes="symbol", fields="symbol,taxid", as_dataframe=True, returnall=True)
        mg_df = mg_output["out"]
        mg_df = mg_df[["symbol", "taxid"]]

        # only get the species that we support
        mg_df = mg_df[mg_df["taxid"].isin(inv_tax_ids.keys())]
        mg_df = mg_df.replace({"taxid": inv_tax_ids})

        mg_df["query"] = mg_df.index.values

        # complete missing information
        # TODO --> do it nicer
        for tax_id, organism in inv_tax_ids.items():
            for id in ids:
                if not ((mg_df["query"] == id) & (mg_df["taxid"] == organism)).any():
                    row = pd.DataFrame({"query": [id], "symbol": ["None"], "taxid": [organism]})
                    mg_df = mg_df.append(row, ignore_index=True)

        mapping = pd.DataFrame(
            {"Gene Name": mg_df["query"], "Reduced Gene Name": mg_df["symbol"], "Organism": mg_df["taxid"]})
        mapping = mapping.replace({np.nan: "None"})
        return mapping

    def get_mapping(self, ids, in_type, organism=None, tar_organism=None, ignore_missing=False,
                    reduction_mode="ensembl"):
        """
        Load prefetched mappings to set of IDs and add missing entries.

        :param ids: Set of either protein IDs or gene names
        :param in_type: Type of needed mapping [protein, orthologs, reduced_genes]
        :param organism: Organism the input IDs (should) belong to
        :param tar_organism: (Orthologs mode) Target organism to find the orthologs of
        :param ignore_missing: Bool indicating if not previously fetched IDs should be ignored
        :param reduction_mode: Mode of how to reduce the gene names
        :return: Dataframe with mapping
        """
        # ===== get precalculated =====
        df, missing = self.get_preloaded(in_list=ids, in_type=in_type, organism=organism, tar_organism=tar_organism,
                                         reduction_mode=reduction_mode)
        # ===== get missing =====
        if len(missing) > 0 and not ignore_missing:
            # ==== Filter protein IDs ====
            if in_type == "protein":
                df2 = self.get_uniprot_mapping(ids=missing, organism=organism)
                if df2 is not None:
                    df = pd.concat([df, df2])
                if organism is not None:
                    df = df[df['Organism'] == organism]
            # ==== Map orthologs ====
            if in_type == "orthologs":
                df2 = self.get_ortholog_mapping(ids=missing, organism=organism, tar_organism=tar_organism)
                if df2 is not None:
                    df = pd.concat([df, df2])
            # ==== Reduce gene names ====
            if in_type == "reduced_genes":
                df2 = self.get_reduced_mapping(ids=missing, organism=organism, reduction_mode=reduction_mode)
                if df2 is not None:
                    df = pd.concat([df, df2])
        return df

    # === Check existing mapping entries and return missing ones ====
    def get_preloaded(self, in_list: list, in_type: str, organism=None, tar_organism=None, reduction_mode="ensembl"):
        # ==== Filter protein IDs ====
        if in_type == "protein":
            organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)",
                         "mouse": "Mus musculus (Mouse)", "rabbit": "Oryctolagus cuniculus (Rabbit)"}
            cur_mapping = self.full_protein_mapping[self.full_protein_mapping["Protein ID"].isin(in_list)]
            if organism is not None:
                cur_mapping = cur_mapping[cur_mapping['Organism'] == organisms[organism]]
            return cur_mapping, list(set(in_list) - set(self.full_protein_mapping["Protein ID"]))
        # ==== Map orthologs ====
        elif in_type == "orthologs":
            cur_mapping = self.full_ortholog_mapping[self.full_ortholog_mapping["source_symbol"].isin(in_list)]
            cur_mapping = cur_mapping[cur_mapping['source_organism'] == organism]
            cur_mapping = cur_mapping[cur_mapping['target_organism'] == tar_organism]
            return cur_mapping, list(set(in_list) - set(cur_mapping["source_symbol"]))
        # ==== Reduce gene names ====
        elif in_type == "reduced_genes":
            cur_mapping = self.full_reduced_gene_mapping[self.full_reduced_gene_mapping["Gene Name"].isin(in_list)]
            cur_mapping = cur_mapping[cur_mapping['Organism'] == organism]
            cur_mapping = cur_mapping[cur_mapping['Mode'] == reduction_mode]
            return cur_mapping, list(set(in_list) - set(cur_mapping["Gene Name"]))
        else:
            return None

    # === Save new mappings to files ====
    def save_mappings(self, mapping_dir):
        mapping_dir = join(here, mapping_dir)
        self.full_protein_mapping.to_csv(mapping_dir + "protein_to_genenames.csv", index=False)
        self.full_ortholog_mapping.to_csv(mapping_dir + "genenames_to_orthologs.csv", index=False)
        self.full_reduced_gene_mapping.to_csv(mapping_dir + "genenames_to_reduced_genenames.csv", index=False)
