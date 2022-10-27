#!/usr/bin/python3

import pandas as pd

def get_reduced_genenames_logging(original, reduced):
    log_df = pd.DataFrame({"Gene Names": original.str.split(";"), "Reduced Gene Names": reduced.str.split(";")}) # easy to compute the next steps
    log_df["Removed Gene Names"] = log_df.apply(lambda row: list(set(row["Gene Names"]).symmetric_difference(set(row["Reduced Gene Names"]))), axis=1)
    log_df["Nr Gene Names"] = log_df["Gene Names"].str.len().mask(log_df["Gene Names"].eq(""),0)
    log_df["Nr Reduced Gene Names"] = log_df["Reduced Gene Names"].str.len().mask(log_df["Reduced Gene Names"].eq(""),0)
    log_df["Nr Removed Gene Names"] = log_df["Nr Gene Names"] - log_df["Nr Reduced Gene Names"]
    return log_df
