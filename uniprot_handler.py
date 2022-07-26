#!/usr/bin/python3

import pandas as pd
from pathlib import Path
import requests

class UniprotHandler:

    full_proteinID_mapping = pd.DataFrame(columns=['Gene Names', 'Gene Names  (primary )', 'Reviewed', 'Organism', 'Protein ID'])
    full_genenames_mapping = pd.DataFrame(columns=['Protein ID', 'Status', 'Organism', 'Gene Name'])

    def __init__(self):
        if Path("protein_to_genenames.csv").exists():
            self.full_proteinID_mapping = pd.read_csv("protein_to_genenames.csv")
        if Path("genenames_to_protein.csv").exists():
            self.full_genenames_mapping = pd.read_csv("genenames_to_protein.csv")

    def get_uniprot_mapping(self, ids, in_type, organism=None):
        organisms = {'Homo sapiens (Human)': '9606', 'Mus musculus (Mouse)': '10090', 'Rattus norvegicus (Rat)': '10116'}
        setup = {'proteinID': {'fields': 'gene_names,gene_primary,reviewed,organism_name,accession'},
                 'genename': {'fields': 'id,reviewed,organism'}}
        if in_type == "proteinID":
            mapping = self.get_uniprot_protein_mapping(ids=ids, organism=organism)

        else: # in_type == "genename":
            # mapping.columns = ['Protein ID', *mapping.columns[1:-1], 'Gene name']
            # mapping['Gene name'] = mapping['Gene name'].apply(lambda x: x.split(" "))
            # mapping = mapping.explode('Gene name')
            # mapping = mapping[mapping['Gene name'].isin(ids)]
            # self.full_genenames_mapping = pd.concat([self.full_genenames_mapping, mapping])
            mapping = pd.DataFrame()
        return mapping

    def get_uniprot_protein_mapping(self, ids, organism=None):
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
                mapping_chunk = mapping_chunk[mapping_chunk['Organism'] == organism]
            if mapping.empty:
                mapping = mapping_chunk
            else:
                mapping = pd.concat([mapping, mapping_chunk])
        mapping.columns = [*mapping.columns[:-1], 'Protein ID']
        mapping['Gene Names'] = mapping['Gene Names'].str.replace(' ', ';')
        mapping['Protein ID'] = mapping['Protein ID'].apply(lambda x: x.split(","))
        mapping = mapping.explode('Protein ID')
        self.full_proteinID_mapping = pd.concat([self.full_proteinID_mapping, mapping])
        return mapping

    def get_mapping(self, ids, in_type, organism=None, ignore_missing=False):
        # ===== get precalculated =====
        df, missing = self.get_preloaded(in_list=ids, in_type=in_type, organism=organism)
        # ===== get missing =====
        if len(missing) > 0 and not ignore_missing:
            df2 = self.get_uniprot_mapping(ids=missing, in_type=in_type, organism=organism)
            if df2 is not None:
                df = pd.concat([df, df2])
        if organism is not None:
            df = df[df['Organism'] == organism]
        return df

    def get_primary_genenames(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, in_type="proteinID", organism=organism)
        if mapping.empty:
            return ""
        else:
            genenames = {x for x in mapping['Gene names  (primary )'] if pd.notna(x)}  # set()
            return ';'.join(genenames)

    def get_all_genenames(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, in_type="proteinID", organism=organism)
        if mapping.empty:
            return ""
        else:
            mapping['Gene names'] = mapping['Gene names'].fillna("").str.upper()
            gene_names_series = mapping['Gene names'].apply(series_to_set)
            genenames = set([x for y in gene_names_series for x in y if x != ""])
            return ';'.join(genenames)

    def get_filtered_ids(self, ids, organism=None, decoy=False):
        if decoy:
            keep = set([x for x in ids if x.startswith(("REV", "CON"))])
        else:
            keep = set()
        mapping = self.get_mapping(ids=ids, in_type="proteinID", organism=organism, ignore_missing=True)
        if mapping.empty:
            return ""
        else:
            prot_ids = set(mapping['Protein ID']).union(keep)
            return ';'.join(prot_ids)

    def get_ids_from_gene(self, genenames, organism=None, reviewed=True):
        mapping = self.get_mapping(ids=genenames, in_type="genename", organism=organism)
        if mapping.empty:
            return ""
        else:
            if reviewed:
                mapping = mapping[mapping["Status"] == "reviewed"]
            prot_ids = set(mapping['Protein ID'])
            return ';'.join(prot_ids)

    def get_preloaded(self, in_list: list, in_type: str, organism=None):
        if in_type == "proteinID":
            cur_mapping = self.full_proteinID_mapping[self.full_proteinID_mapping["Protein ID"].isin(in_list)]
            if organism is not None:
                cur_mapping = cur_mapping[cur_mapping['Organism'] == organism]
            return cur_mapping, list(set(in_list) - set(self.full_proteinID_mapping["Protein ID"]))
        elif in_type == "genename":
            cur_mapping = self.full_genenames_mapping[self.full_genenames_mapping["Gene name"].isin(in_list)]
            if organism is not None:
                cur_mapping = cur_mapping[cur_mapping['Organism'] == organism]
            return cur_mapping, list(set(in_list) - set(self.full_genenames_mapping["Gene name"]))
        else:
            return None

    def save_mappings(self):
        self.full_proteinID_mapping.to_csv("protein_to_genenames.csv", index=False)
        self.full_genenames_mapping.to_csv("genenames_to_protein.csv", index=False)


def series_to_set(x):
    return set(x.split(";"))
