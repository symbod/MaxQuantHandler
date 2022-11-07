#!/usr/bin/python3

import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


# ==== Get Log Plots ====
def create_overview_plot(logging, out_dir, file_type="png"):
    # ==== Prepare dataframe ====
    logging_names = [x for x in logging.columns if x.startswith(("Nr"))]
    df = pd.melt(logging[logging_names])
    # ==== Create plot ====
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = sns.boxplot(data=df, x="variable", y="value")
    ax.set(title="Distribution of the number of entries per line",
           ylabel="Number of entries per line", xlabel="")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"overview_log.{file_type}"), bbox_inches='tight')


def create_filter_detailed_plot(logging, organism, reviewed, decoy, out_dir, file_type="png"):
    organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)",
                 "mouse": "Mus musculus (Mouse)", "rabbit": "Oryctolagus cuniculus (Rabbit)"}
    # ==== Prepare dataframe ====
    df_dict = dict()
    if decoy is False:
        df_dict["Decoys"] = len(logging[logging["Organism"] == "Decoy"].index)
    if reviewed:
        df_dict["Unreviewed"] = len(
            logging[(logging["Reviewed"] == "unreviewed") & (logging["Organism"] == organisms[organism])].index)
    df_dict["Not found IDs"] = len(logging[logging["Organism"] == "Not found"].index)
    df_dict["Wrong organism"] = len(
        logging[~logging["Organism"].isin(["Not found", "Decoy", organisms[organism]])].index)
    df = pd.melt(pd.DataFrame(df_dict, index=[0]))
    # ==== Create plot ====
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = sns.barplot(data=df, x="variable", y="value")
    ax.set(title="Number of entries removed by cause",
           ylabel="Number of entries", xlabel="")
    ax.bar_label(ax.containers[0])
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"filter_detailed_log.{file_type}"), bbox_inches='tight')
