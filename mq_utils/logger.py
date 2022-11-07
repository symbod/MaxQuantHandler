#!/usr/bin/python3

import pandas as pd


# ==== Logging DataFrame For Filtering Ids ====
def get_filter_ids_logging(original, filtered, handler, organism):
    # ==== create dataframe with for each row original ids, filtered ids, nr ids, etc. ====
    log_df = pd.DataFrame({"IDs": original.str.split(";"), "Filtered IDs": filtered.str.split(";")})
    log_df["Removed IDs"] = log_df.apply(
        lambda row: list(set(row["IDs"]).difference(set(row["Filtered IDs"]))), axis=1)
    log_df["Nr IDs"] = log_df["IDs"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Filtered IDs"] = log_df["Filtered IDs"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Removed IDs"] = log_df["Nr IDs"] - log_df["Nr Filtered IDs"]

    # ==== for removed ids --> find out the cause and create df ====
    if sum(log_df["Nr Removed IDs"]) == 0:
        return {"Overview_Log": log_df, "Detailed_Log": pd.DataFrame()}
    else:
        removed_ids = [x for x in log_df["Removed IDs"].explode() if str(x) != 'nan']
        decoy = [x for x in removed_ids if x.startswith(("REV", "CON"))]  # save decoys
        removed_ids_nondecoy = [x for x in removed_ids if (x not in decoy)]
        # ==== Get Information for removed IDs ====
        df, missing = handler.get_preloaded(in_list=removed_ids_nondecoy, in_type="protein", organism=organism)
        df = pd.concat([df, pd.DataFrame({"Protein ID": missing})], ignore_index=True).fillna("Not found")
        df = pd.concat([df, pd.DataFrame({"Protein ID": decoy})], ignore_index=True).fillna("Decoy")
        return {"Overview_Log": log_df, "Detailed_Log": df}


# ==== Logging DataFrame For Remapped Gene Names ====
def get_remapped_genenames_logging(original, remapped, protein_ids, handler, organism):
    # ==== create dataframe with for each row original names, remapped names, nr names, etc. ====
    log_df = pd.DataFrame({"Gene Names": original.str.split(";"), "Remapped Gene Names": remapped.str.split(";")})
    log_df["Added Gene Names"] = log_df.apply(
        lambda row: list(set(row["Remapped Gene Names"]).difference(set(row["Gene Names"]))), axis=1)
    log_df["Removed Gene Names"] = log_df.apply(
        lambda row: list(set(row["Gene Names"]).difference(set(row["Remapped Gene Names"]))), axis=1)
    log_df["Nr Gene Names"] = log_df["Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Remapped Gene Names"] = log_df["Remapped Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Added Gene Names"] = log_df["Added Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Removed Gene Names"] = log_df["Removed Gene Names"].apply(lambda x: len(list(filter(None, x))))

    # ==== for added names --> find out the cause and create df ====
    df, missing = handler.get_preloaded(in_list=protein_ids, in_type="protein", organism=organism)
    df = df[df["Gene Names (primary)"].isin(log_df["Added Gene Names"].explode().to_list())]
    return {"Overview_Log": log_df, "Detailed_Log": df}


# ==== Logging DataFrame For Remapped Gene Names ====
def get_ortholog_genenames_logging(original, orthologs, handler, organism, tar_organism):
    # ==== Get Information for all original names ====
    original_names = [x for x in ";".join(original).split(";") if str(x) != 'nan']
    df, missing = handler.get_preloaded(in_list=original_names, in_type="orthologs",
                                        organism=organism, tar_organism=tar_organism)
    removed_gene_names = set(original_names).difference(set(df["source_symbol"]))

    # ==== create dataframe with for each row original names, ortholog names, nr names, etc. ====
    log_df = pd.DataFrame({"Gene Names": original.str.split(";"), "Ortholog Gene Names": orthologs.str.split(";")})
    log_df["Removed Gene Names"] = log_df.apply(
        lambda row: list(set(row["Gene Names"]).intersection(removed_gene_names)), axis=1)
    log_df["Nr Gene Names"] = log_df["Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Ortholog Gene Names"] = log_df["Ortholog Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Removed Gene Names"] = log_df["Removed Gene Names"].apply(lambda x: len(list(filter(None, x))))

    # ==== for removed ids --> find out the cause and create df ====
    if len(removed_gene_names) == 0:
        return {"Overview_Log": log_df, "Detailed_Log": pd.DataFrame()}
    else:
        # ==== Get Information for removed names ====
        df = df[df["target_symbol"] == ""]  # no ortholog found
        df = pd.concat([df, pd.DataFrame({"source_symbol": missing})], ignore_index=True).fillna("Not found")
        return {"Overview_Log": log_df, "Detailed_Log": df}


# ==== Logging DataFrame For Reducing Gene Names ====
def get_reduced_genenames_logging(original, reduced, handler, organism, mode):
    # ==== Get Information for all original names ====
    df, missing = handler.get_preloaded(in_list=original, in_type="reduced_genes",
                                        organism=organism, reduction_mode=mode)
    # ==== create dataframe with for each row original names, reduced names, nr names, etc. ====
    log_df = pd.DataFrame({"Gene Names": original.str.split(";"), "Reduced Gene Names": reduced.str.split(";")})
    log_df["Added Gene Names"] = log_df.apply(
        lambda row: list(set(row["Reduced Gene Names"]).difference(set(row["Gene Names"]))), axis=1)
    log_df["Removed Gene Names"] = log_df.apply(
        lambda row: list(set(row["Gene Names"]).difference(set(row["Reduced Gene Names"]))), axis=1)
    log_df["Nr Gene Names"] = log_df["Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Reduced Gene Names"] = log_df["Reduced Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Added Gene Names"] = log_df["Added Gene Names"].apply(lambda x: len(list(filter(None, x))))
    log_df["Nr Removed Gene Names"] = log_df["Removed Gene Names"].apply(lambda x: len(list(filter(None, x))))

    # ==== for reduced ids --> find out the cause and create df ====
    df = pd.concat([df, pd.DataFrame({"Gene Name": missing})], ignore_index=True).fillna("Not found")
    df = df[df["Reduced Gene Name"] == "Not found"]
    return {"Overview_Log": log_df, "Detailed_Log": df}
