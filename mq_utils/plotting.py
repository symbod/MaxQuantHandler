#!/usr/bin/python3

import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


# ==== Get Log Plots ====
def create_overview_plot(logging, out_dir, file_type="png"):
    """
    Overview plot for logging data of proharmed's functions

    :param logging: overview logging Dataframe that has been returned by one of the four methods of proharmed
    :param out_dir: (optional) output directory to save the plots
    :param file_type: (optional) file type for the plots (png, pdf, jpg, ...)
    :return:
    """
    # ==== Prepare dataframe ====
    logging_names = [x for x in logging.columns if x.startswith(("Nr"))]
    # ==== Create box plot ====
    df = pd.melt(logging[logging_names])
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = sns.boxplot(data=df, x="variable", y="value")
    ax.set(title="Distribution of the number of entries per line",
           ylabel="Number of entries per line", xlabel="")
    plt.xticks(rotation=45)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"overview_log_box.{file_type}"), bbox_inches='tight')
    # ==== Create bar plot ====
    df = logging[logging_names].sum().to_frame().reset_index(level=0)
    # ==== Create plot ====
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = sns.barplot(data=df, x="index", y=0)
    ax.set(title="Number of entries added and removed",
           ylabel="Number of entries", xlabel="")
    ax.bar_label(ax.containers[0])
    plt.xticks(rotation=45)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"overview_log_bar.{file_type}"), bbox_inches='tight')


def create_filter_detailed_plot(logging, organism, out_dir, file_type="png"):
    """
    Detailed plot for logging data of proharmed's filter_protein_ids method

    :param logging: detailed logging Dataframe that has been returned by the filtered_protein_ids method of proharmed
    :param organism: organism that has been specified in the filter_protein_ids method of proharmed
    :param out_dir: (optional) output directory to save the plot
    :param file_type: (optional) file type for the plot
    :return:
    """
    organisms = {"human": "Homo sapiens (Human)", "rat": "Rattus norvegicus (Rat)",
                 "mouse": "Mus musculus (Mouse)", "rabbit": "Oryctolagus cuniculus (Rabbit)"}
    # ==== Prepare dataframe ====
    df_dict = dict()
    if logging.empty:
        df_dict["Decoys"] = 0
        df_dict["Unreviewed"] = 0
        df_dict["Not found IDs"] = 0
        df_dict["Wrong organism"] = 0
    else:
        df_dict["Decoys"] = len(logging[logging["Organism"] == "Decoy"].index)
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
    fig.savefig(os.path.join(out_dir, f"detailed_log.{file_type}"), bbox_inches='tight')


def create_reduced_detailed_plot(logging, out_dir, file_type="png"):
    """
    Detailed plot for logging data of proharmed's reduce_genenames method

    :param logging: detailed logging Dataframe that has been returned by the reduce_genenames method of proharmed
    :param out_dir: (optional) output directory to save the plot
    :param file_type: (optional) file type for the plot
    :return:
    """
    # ==== Prepare dataframe ====
    df_dict = dict()
    if logging.empty:
        df_dict["Not found IDs"] = 0
        df_dict["Different Name"] = 0
    else:
        df_dict["Not found IDs"] = len(logging[logging["Reduced Gene Name"] == "Not found"].index)
        df_dict["Different Name"] = len(
            logging[logging["Reduced Gene Name"] != "Not found"].index)
    df = pd.melt(pd.DataFrame(df_dict, index=[0]))
    # ==== Create plot ====
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = sns.barplot(data=df, x="variable", y="value")
    ax.set(title="Number of entries removed by cause",
           ylabel="Number of entries", xlabel="")
    ax.bar_label(ax.containers[0])
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"detailed_log.{file_type}"), bbox_inches='tight')


def create_ortholog_detailed_plot(logging, organism, out_dir, file_type="png"):
    """
    Detailed plot for logging data of proharmed's map_orthologs method

    :param logging: detailed logging Dataframe that has been returned by the map_orthologs method of proharmed
    :param organism: organism that has been specified in the map_orthologs method of proharmed
    :param out_dir: (optional) output directory to save the plot
    :param file_type: (optional) file type for the plot    :return:
    """
    # ==== Prepare dataframe ====
    df_dict = dict()
    if logging.empty:
        df_dict["Not found target names"] = 0
        df_dict["Not found target IDs"] = 0
        df_dict["Not found source IDs"] = 0
        df_dict["Wrong organism"] = 0
    else:
        df_dict["Not found target names"] = len( logging[(logging["ortholog_ensg"] != "") & (logging["target_symbol"] == "")].index)
        df_dict["Not found target IDs"] = len( logging[(logging["ortholog_ensg"] == "") & (logging["target_symbol"] == "")].index)
        df_dict["Not found source IDs"] = len(logging[logging["ortholog_ensg"] == "Not found"].index)
        df_dict["Wrong organism"] = len(
            logging[~logging["source_organism"].isin(["Not found", organism])].index)
    df = pd.melt(pd.DataFrame(df_dict, index=[0]))
    # ==== Create plot ====
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = sns.barplot(data=df, x="variable", y="value")
    ax.set(title="Number of entries removed by cause",
           ylabel="Number of entries", xlabel="")
    ax.bar_label(ax.containers[0])
    plt.xticks(rotation=45)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"detailed_log.{file_type}"), bbox_inches='tight')
