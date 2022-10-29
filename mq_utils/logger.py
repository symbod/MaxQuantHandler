#!/usr/bin/python3

import pandas as pd


# === Logging DataFrame For Reducing Gene Names ===
def get_reduced_genenames_logging(original, reduced):
    log_df = pd.DataFrame({"Gene Names": original.str.split(";"), "Reduced Gene Names": reduced.str.split(";")}) # easy to compute the next steps
    log_df["Removed Gene Names"] = log_df.apply(lambda row: list(set(row["Gene Names"]).symmetric_difference(set(row["Reduced Gene Names"]))), axis=1)
    log_df["Nr Gene Names"] = log_df["Gene Names"].str.len().mask(log_df["Gene Names"].eq(""),0)
    log_df["Nr Reduced Gene Names"] = log_df["Reduced Gene Names"].str.len().mask(log_df["Reduced Gene Names"].eq(""),0)
    log_df["Nr Removed Gene Names"] = log_df["Nr Gene Names"] - log_df["Nr Reduced Gene Names"]
    return log_df

# === Logging DataFrame For Filtering Ids ===
def get_filter_ids_logging(original, filtered):
    log_df = pd.DataFrame({"IDs":original.str.split(";"), "Filtered IDs": filtered.str.split(";")})
    log_df[ "Removed IDs" ] = log_df.apply(
        lambda row: list( set( row[ "IDs" ] ).symmetric_difference( set( row[ "Filtered IDs" ] ) ) ),
        axis=1 )
    log_df[ "Nr IDs" ] = log_df[ "IDs" ].str.len().mask( log_df[ "IDs" ].eq( "" ), 0 )
    log_df[ "Nr Filtered IDs" ] = log_df[ "Filtered IDs" ].str.len().mask(
        log_df[ "Filtered IDs" ].eq( "" ), 0 )
    log_df[ "Nr Removed IDs" ] = log_df[ "Nr IDs" ] - log_df[ "Nr Filtered IDs" ]
    return log_df
