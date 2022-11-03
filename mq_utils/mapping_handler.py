#!/usr/bin/python3

from os.path import dirname, abspath, join
import pandas as pd
from pathlib import Path
from gprofiler import GProfiler
import requests
from HGNC_mapping import get_HGNC_mapping
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

    # === Uniprot Mapping ====
    def get_uniprot_mapping(self, ids, organism=None):
        organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)", "mouse": "Mus musculus (Mouse)",
                     "rabbit": "Oryctolagus cuniculus (Rabbit)"}
        url = 'https://rest.uniprot.org/uniprotkb/accessions'
        mapping = pd.DataFrame()
        for i in range(0, len(ids), 500):
            ids_chunk = ids[i:i + 500]
            params = {
                'format': 'tsv',
                'accessions': ",".join([x for x in ids_chunk if not x.startswith(("REV", "CON"))]),
                'fields': 'gene_names,gene_primary,reviewed,organism_name,accession'}
            f = requests.get(url=url, params=params)
            mapping_chunk = pd.read_csv(f.url, sep="\t")
            if organism is not None:
                mapping_chunk = mapping_chunk[mapping_chunk['Organism'] == organisms[organism]]
            if mapping.empty:
                mapping = mapping_chunk
            else:
                mapping = pd.concat([mapping, mapping_chunk])
        mapping.columns = [*mapping.columns[:-1], 'Protein ID']
        mapping['Gene Names'] = mapping['Gene Names'].str.replace(' ', ';')
        mapping['Protein ID'] = mapping['Protein ID'].apply(lambda x: x.split(","))
        mapping = mapping.explode('Protein ID')
        self.full_protein_mapping = pd.concat([self.full_protein_mapping, mapping])
        return mapping

    # === Ortholog Mapping ====
    def get_ortholog_mapping(self, ids, organism, tar_organism):
        organisms = {"human": "hsapiens", "mouse": "mmusculus", "rat": "rnorvegicus", "rabbit": "ocuniculus"}
        gp = GProfiler(return_dataframe=True)
        mapping = gp.orth(organism=organisms[organism], query=ids, target=organisms[tar_organism])
        mapping = mapping[['incoming', 'converted', 'ortholog_ensg', 'name', 'description']]
        mapping.columns = ['source_symbol', 'ensg', 'ortholog_ensg', 'target_symbol', 'description']
        # save organism info
        mapping.insert(loc=1, column='source_organism', value=organism)
        mapping.insert(loc=5, column='target_organism', value=tar_organism)
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
        mapping["Mode"] = reduction_mode
        self.full_reduced_gene_mapping = pd.concat([self.full_reduced_gene_mapping, mapping])
        return mapping

    def get_ensembl_reduction(self, ids, organism):  # organism required
        organisms = {"human": "hsapiens", "mouse": "mmusculus", "rat": "rnorvegicus", "rabbit": "ocuniculus"}
        gp = GProfiler(return_dataframe=True)
        gp_df = gp.convert(organism=organisms[organism], query=ids, target_namespace="ENSG")
        if len(gp_df) == 0:
            return None
        else:
            # check if ensembl id == name --> then take the incoming name
            gp_df["new_name"] = np.where(gp_df["converted"] == gp_df["name"], gp_df["incoming"], gp_df["name"])
            # replace Nones
            gp_df["new_name"] = np.where(gp_df["name"].isin(["None", None]), gp_df["name"], gp_df["new_name"])

            mapping = gp_df[["incoming", "new_name"]]
            mapping.columns = ["Gene Name", "Reduced Gene Name"]
            mapping.loc[:, "Organism"] = organism
            return mapping

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
        if in_type == "protein":
            organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)",
                         "mouse": "Mus musculus (Mouse)", "rabbit": "Oryctolagus cuniculus (Rabbit)"}
            cur_mapping = self.full_protein_mapping[self.full_protein_mapping["Protein ID"].isin(in_list)]
            if organism is not None:
                cur_mapping = cur_mapping[cur_mapping['Organism'] == organisms[organism]]
            return cur_mapping, list(set(in_list) - set(self.full_protein_mapping["Protein ID"]))
        elif in_type == "orthologs":
            cur_mapping = self.full_ortholog_mapping[self.full_ortholog_mapping["source_symbol"].isin(in_list)]
            cur_mapping = cur_mapping[cur_mapping['source_organism'] == organism]
            cur_mapping = cur_mapping[cur_mapping['target_organism'] == tar_organism]
            return cur_mapping, list(set(in_list) - set(cur_mapping["source_symbol"]))
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
