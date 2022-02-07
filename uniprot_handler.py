#!/usr/bin/python3

import pandas as pd
import urllib.parse
import urllib.request
import io
from pathlib import Path

class Uniprot_Handler:

    full_mapping = pd.DataFrame(columns=['Gene names', 'Gene names  (primary )', 'Status', 'Organism', 'Protein ID'])

    def __init__(self):
        if Path("protein_to_genenames.csv").exists():
            self.full_mapping = pd.read_csv("protein_to_genenames.csv")

    def get_uniprot_mapping(self, ids, organism=None):
        url = 'https://www.uniprot.org/uploadlists/'
        params = {
                'from': "ACC+ID",
                'to': 'ACC',
                'format': 'tab',
                'query': " ".join(ids),
                'columns': 'genes,genes(PREFERRED),reviewed,organism'}

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        if len(response.decode('utf-8')) == 0:
            return None
        mapping = pd.read_csv(io.StringIO(response.decode('utf-8')), sep="\t")
        mapping.columns = [*mapping.columns[:-1], 'Protein ID']
        if organism is not None:
            mapping = mapping[mapping[organism] == mapping[organism]]
        self.full_mapping = pd.concat([self.full_mapping, mapping])
        return mapping

    def get_mapping(self, ids, organism=None):
        # ===== get precalculated =====
        df, missing = self.get_preloaded(ids)
        # ===== get missing =====
        if len(missing) > 0:
            df2 = self.get_uniprot_mapping(ids=missing, organism=organism)
            if df2 is not None:
                df = pd.concat([df, df2])
        return df

    def get_primary_genenames(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, organism=organism)
        if mapping.empty:
            return ""
        else:
            genenames = set(mapping['Gene names  (primary )'])
            return ';'.join(genenames)

    def get_all_genenames(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, organism=organism)
        if mapping.empty:
            return ""
        else:
            mapping['Gene names'] = mapping['Gene names'].fillna("").str.upper()
            gene_names_series = mapping['Gene names'].apply(series_to_set)
            genenames = set([x for y in gene_names_series for x in y if x != ""])
            return ';'.join(genenames)

    def get_filtered_ids(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, organism=organism)
        if mapping.empty:
            return ""
        else:
            prot_ids = set(mapping['Protein ID'])
            return ';'.join(prot_ids)

    def get_preloaded(self, in_list:list):
        cur_mapping = self.full_mapping[self.full_mapping["Protein ID"].isin(in_list)]
        return cur_mapping, list(set(in_list)-set(self.full_mapping["Protein ID"]))

    def save_mappings(self):
        self.full_mapping.to_csv("protein_to_genenames.csv", index=False)


def series_to_set(x):
    return set(x.split(" "))