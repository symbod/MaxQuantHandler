#!/usr/bin/python3

import os
import requests
import pandas as pd
from collections import Counter
from upsetplot import plot
from matplotlib import pyplot as plt


def count_intersection(data: dict, threshold=1):
    full_list = []
    [full_list.extend(data[key]) for key in data]
    intersections = pd.DataFrame.from_dict(Counter(full_list), orient='index').reset_index()
    intersections = intersections.sort_values(by=[0], ascending=False).rename(columns={"index": "ID", 0: "Count"})
    return intersections[intersections["Count"] >= threshold]


def plot_intersections(data: dict, out_dir, file_type="png"):
    names = list()
    all_elements = set()
    data_sets = dict()
    for key in data:
        names.append(key)
        data_sets[key] = set(data[key])
        all_elements = all_elements.union(data_sets[key])
    df = pd.DataFrame([[element in data_sets[key] for key in data_sets] for element in all_elements], columns=names)
    df_up = df.groupby(names).size()

    fig = plt.figure(figsize=(12, 6), dpi=80)
    plot(df_up, fig=fig, element_size=None, show_counts=True, orientation='horizontal')
    plt.suptitle('Overview of Intersections')
    plt.show()
    fig.savefig(os.path.join(out_dir, f"overview_intersections.{file_type}"), bbox_inches='tight')


def inspect_for_drugs(genes: list):
    url = 'https://api.drugst.one/create_network'
    myobj = {"network": {'nodes': [{"id": gene, "group": "gene"} for gene in genes]}}
    result = requests.post(url, json=myobj)
    return "https://drugst.one?id="+result.json()


def load_multi_files(files: [list, str], columns: [list, str]):
    # Read column names
    if type(columns) is str:
        column_names = pd.read_csv(columns, names=['columns', 'file'])
    elif type(columns) is list:
        column_names = pd.DataFrame(columns, columns=['columns']).reset_index()
    else:
        column_names = pd.DataFrame()

    # Read files
    data = dict()
    if type(files) is list:
        for index, file in enumerate(files):

            if "file" in column_names.columns:
                col = column_names[column_names['file'] == os.path.basename(file)]["columns"].iloc[0]
            elif "index" in column_names.columns:
                col = column_names[column_names['index'] == index]["columns"].iloc[0]
            else:
                col = None

            if type(file) is pd.DataFrame:
                if col is not None:
                    data[str(index)] = file[col].to_list()
                else:
                    data[str(index)] = file.iloc[:, 0].to_list()
            elif type(file) is list:
                data[str(index)] = file
            elif os.path.isfile(file):
                if col is not None:
                    data[os.path.basename(file)] = pd.read_csv(file)[col].to_list()
                else:
                    data[os.path.basename(file)] = pd.read_csv(file).iloc[:, 0].to_list()
    elif type(files) is str:
        if os.path.isfile(files):
            if "file" in column_names.columns:
                col = column_names[column_names['file'] == os.path.basename(files)]["columns"].iloc[0]
                data[os.path.basename(files)] = pd.read_csv(files)[col].to_list()
            elif "index" in column_names.columns:
                col = column_names[column_names['index'] == 0]["columns"].iloc[0]
                data[os.path.basename(files)] = pd.read_csv(files)[col].to_list()
            else:
                data[os.path.basename(files)] = pd.read_csv(files).iloc[:, 0].to_list()
        elif os.path.isdir(files):
            if "file" in column_names.columns:
                for file in column_names["file"]:
                    col = column_names[column_names['file'] == os.path.basename(file)]["columns"].iloc[0]
                    data[os.path.basename(file)] = pd.read_csv(os.path.join(files, file))[col].to_list()
            else:
                for index, file in enumerate(os.listdir(files)):
                    if "index" in column_names.columns:
                        col = column_names[column_names['index'] == index]["columns"].iloc[0]
                        data[os.path.basename(file)] = pd.read_csv(os.path.join(files, file))[col].to_list()
                    else:
                        data[os.path.basename(file)] = pd.read_csv(os.path.join(files, file)).iloc[:, 0].to_list()

    # split and explode
    for data_key in data:
        seen = set()
        data[data_key] = [x for x in ";".join(data[data_key]).split(";") if x not in seen and not seen.add(x)]

    return data
