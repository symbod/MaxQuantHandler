#!/usr/bin/python3

import pandas as pd
import urllib.parse
import urllib.request
import io
from pathlib import Path


class UniprotHandler:

    full_proteinID_mapping = pd.DataFrame(
        columns=['Gene names', 'Gene names  (primary )', 'Status', 'Organism', 'Protein ID'])
    full_genenames_mapping = pd.DataFrame(columns=['Protein ID', 'Status', 'Organism', 'Gene names'])

    def __init__(self):
        if Path("protein_to_genenames.csv").exists():
            self.full_proteinID_mapping = pd.read_csv("protein_to_genenames.csv")
        if Path("genenames_to_protein.csv").exists():
            self.full_genenames_mapping = pd.read_csv("genenames_to_protein.csv")

    def get_uniprot_mapping(self, ids, in_type, organism=None):
        setup = {'proteinID': {'from': 'ACC+ID', 'columns': 'genes,genes(PREFERRED),reviewed,organism'},
                 'genename': {'from': 'Genenames', 'columns': 'id,reviewed,organism'}}
        url = 'https://www.uniprot.org/uploadlists/'
        params = {
            'from': setup[in_type]['from'],
            'to': 'ACC',
            'format': 'tab',
            'query': " ".join(ids),
            'columns': setup[in_type]['columns']}

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        mapping = pd.read_csv(io.StringIO(response.decode('utf-8')), sep="\t")
        if organism is not None:
            mapping = mapping[mapping[organism] == mapping[organism]]
        if in_type == "proteinID":
            mapping.columns = [*mapping.columns[:-1], 'Protein ID']
            self.full_proteinID_mapping = pd.concat([self.full_proteinID_mapping, mapping])
        elif in_type == "genename":
            mapping.columns = [*mapping.columns[:-1], 'Gene name']
            self.full_genenames_mapping = pd.concat([self.full_genenames_mapping, mapping])
        return mapping

    def get_mapping(self, ids, in_type, organism=None):
        # ===== get precalculated =====
        df, missing = self.get_preloaded(in_list=ids, in_type=in_type)
        # ===== get missing =====
        if len(missing) > 0:
            df2 = self.get_uniprot_mapping(ids=missing, in_type=in_type, organism=organism)
            if df2 is not None:
                df = pd.concat([df, df2])
        return df

    def get_primary_genenames(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, in_type="proteinID", organism=organism)
        if mapping.empty:
            return ""
        else:
            genenames = set(mapping['Gene names  (primary )'])
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

    def get_filtered_ids(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, in_type="proteinID", organism=organism)
        if mapping.empty:
            return ""
        else:
            prot_ids = set(mapping['Protein ID'])
            return ';'.join(prot_ids)

    def get_preloaded(self, in_list: list, in_type: str):
        if in_type == "proteinID":
            cur_mapping = self.full_proteinID_mapping[self.full_proteinID_mapping["Protein ID"].isin(in_list)]
            return cur_mapping, list(set(in_list) - set(self.full_proteinID_mapping["Protein ID"]))
        elif in_type == "genename":
            cur_mapping = self.full_genenames_mapping[self.full_genenames_mapping["Gene name"].isin(in_list)]
            return cur_mapping, list(set(in_list) - set(self.full_genenames_mapping["Gene name"]))
        else:
            return None

    def save_mappings(self):
        self.full_proteinID_mapping.to_csv("protein_to_genenames.csv", index=False)
        self.full_genenames_mapping.to_csv("genenames_to_protein.csv", index=False)


def series_to_set(x):
    return set(x.split(" "))
