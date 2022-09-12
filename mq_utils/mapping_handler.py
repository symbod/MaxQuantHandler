#!/usr/bin/python3

from os.path import dirname, abspath, join
import pandas as pd
from pathlib import Path
from gprofiler import GProfiler
import requests

here=dirname(dirname(abspath(__file__)))

class MappingHandler:

    full_protein_mapping = pd.DataFrame(columns=['Gene Names', 'Gene Names (primary)', 'Reviewed', 'Organism',
                                                 'Protein ID'])
    full_gene_mapping = pd.DataFrame(columns=['Protein ID', 'Status', 'Organism', 'Gene Name'])
    full_ortholog_mapping = pd.DataFrame(columns=['source_symbol', 'source_organism', 'ensg', 'ortholog_ensg',
                                                  'target_symbol', 'target_organism', 'description'])

    def __init__(self, mapping_dir):
        mapping_dir = join(here, mapping_dir)
        if Path(mapping_dir+"protein_to_genenames.csv").exists():
            self.full_protein_mapping = pd.read_csv(mapping_dir+"protein_to_genenames.csv")
        if Path(mapping_dir+"genenames_to_protein.csv").exists():
            self.full_gene_mapping = pd.read_csv(mapping_dir+"genenames_to_protein.csv")
        if Path(mapping_dir+"genenames_to_orthologs.csv").exists():
            self.full_ortholog_mapping = pd.read_csv(mapping_dir+"genenames_to_orthologs.csv")

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

    def get_mapping(self, ids, in_type, organism=None, tar_organism=None, ignore_missing=False):
        # ===== get precalculated =====
        df, missing = self.get_preloaded(in_list=ids, in_type=in_type, organism=organism, tar_organism=tar_organism)
        # ===== get missing =====
        if len(missing) > 0 and not ignore_missing:
            if in_type == "protein":
                df2 = self.get_uniprot_mapping(ids=missing, organism=organism)
                if df2 is not None:
                    df = pd.concat([df, df2])
                if organism is not None:
                    df = df[df['Organism'] == organism]
            if in_type == "gene":
                # TODO
                df = pd.DataFrame()
            if in_type == "orthologs":
                df2 = self.get_ortholog_mapping(ids=missing, organism=organism, tar_organism=tar_organism)
                if df2 is not None:
                    df = pd.concat([df, df2])
        return df

    def get_orthologs(self, ids, organism: str, tar_organism: str):
        mapping = self.get_mapping(ids=ids, in_type="orthologs", organism=organism,
                                   tar_organism=tar_organism, ignore_missing=True)
        if mapping.empty:
            return ""
        else:
            orthologs = {x for x in mapping['target_symbol'] if pd.notna(x)}
            return ';'.join(orthologs)

    def get_primary_genenames(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, in_type="protein", organism=organism)
        if mapping.empty:
            return ""
        else:
            genenames = {x for x in mapping['Gene Names (primary)'] if pd.notna(x)}  # set()
            return ';'.join(genenames)

    def get_all_genenames(self, ids, organism=None):
        mapping = self.get_mapping(ids=ids, in_type="protein", organism=organism)
        if mapping.empty:
            return ""
        else:
            mapping['Gene names'] = mapping['Gene names'].fillna("").str.upper()
            gene_names_series = mapping['Gene names'].apply(series_to_set)
            genenames = set([x for y in gene_names_series for x in y if x != ""])
            return ';'.join(genenames)

    def get_filtered_ids(self, ids, in_type, organism=None, decoy=False):
        mapping = self.get_mapping(ids=ids, in_type=in_type, organism=organism, ignore_missing=True)
        if mapping.empty:
            return ""
        if in_type == "protein":
            if decoy:
                keep = set([x for x in ids if x.startswith(("REV", "CON"))])
            else:
                keep = set()
            prot_ids = set(mapping['Protein ID']).union(keep)
            return ';'.join(prot_ids)
        elif in_type == "gene":
            gene_names = set(mapping['Gene name'])
            return ';'.join(gene_names)
        else:
            return set()

    def get_ids_from_gene(self, genenames, organism=None, reviewed=True):
        mapping = self.get_mapping(ids=genenames, in_type="gene", organism=organism)
        if mapping.empty:
            return ""
        else:
            if reviewed:
                mapping = mapping[mapping["Status"] == "reviewed"]
            protein_ids = set(mapping['Protein ID'])
            return ';'.join(protein_ids)

    def get_preloaded(self, in_list: list, in_type: str, organism=None, tar_organism=None):
        if in_type == "protein":
            organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)",
                         "mouse": "Mus musculus (Mouse)"}
            cur_mapping = self.full_protein_mapping[self.full_protein_mapping["Protein ID"].isin(in_list)]
            if organism is not None:
                cur_mapping = cur_mapping[cur_mapping['Organism'] == organisms[organism]]
            return cur_mapping, list(set(in_list) - set(self.full_protein_mapping["Protein ID"]))
        elif in_type == "gene":
            organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)",
                         "mouse": "Mus musculus (Mouse)"}
            cur_mapping = self.full_gene_mapping[self.full_gene_mapping["Gene name"].isin(in_list)]
            if organism is not None:
                cur_mapping = cur_mapping[cur_mapping['Organism'] == organisms[organism]]
            return cur_mapping, list(set(in_list) - set(self.full_gene_mapping["Gene name"]))
        elif in_type == "orthologs":
            cur_mapping = self.full_ortholog_mapping[self.full_ortholog_mapping["source_symbol"].isin(in_list)]
            cur_mapping = cur_mapping[cur_mapping['source_organism'] == organism]
            cur_mapping = cur_mapping[cur_mapping['target_organism'] == tar_organism]
            return cur_mapping, list(set(in_list) - set(self.full_ortholog_mapping["source_symbol"]))
        else:
            return None

    def save_mappings(self, mapping_dir):
        mapping_dir = join(here, mapping_dir)
        self.full_protein_mapping.to_csv(mapping_dir+"protein_to_genenames.csv", index=False)
        self.full_gene_mapping.to_csv(mapping_dir+"genenames_to_protein.csv", index=False)
        self.full_ortholog_mapping.to_csv(mapping_dir+"genenames_to_orthologs.csv", index=False)


def series_to_set(x):
    return set(x.split(";"))
