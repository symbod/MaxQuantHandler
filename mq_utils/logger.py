#!/usr/bin/python3

import pandas as pd


# === Logging DataFrame For Reducing Gene Names ===
def get_reduced_genenames_logging(original, reduced):
    log_df = pd.DataFrame({"Gene Names": original.str.split(";"), "Reduced Gene Names": reduced.str.split(";")}) # easy to compute the next steps
    log_df["Removed Gene Names"] = log_df.apply(lambda row: list(set(row["Gene Names"]).symmetric_difference(set(row["Reduced Gene Names"]))), axis=1)
    log_df["Nr Gene Names"] = log_df["Gene Names"].str.len().mask(log_df["Gene Names"].eq(""),0)
    log_df["Nr Reduced Gene Names"] = log_df["Reduced Gene Names"].str.len().mask(log_df["Reduced Gene Names"].eq(""),0)
    log_df["Nr Removed Gene Names"] = log_df["Nr Gene Names"] - log_df["Nr Reduced Gene Names"]
    return {"Overview_Log":log_df}

# === Logging DataFrame For Filtering Ids ===
def get_filter_ids_logging(original, filtered, handler, organism):
    # create dataframe with for each row original ids, filtered ids, nr ids, etc.
    log_df = pd.DataFrame({"IDs":original.str.split(";"), "Filtered IDs": filtered.str.split(";")})
    log_df[ "Removed IDs" ] = log_df.apply(
        lambda row: list( set( row[ "IDs" ] ).symmetric_difference( set( row[ "Filtered IDs" ] ) ) ),
        axis=1 )
    log_df[ "Nr IDs" ] = log_df[ "IDs" ].str.len().mask( log_df[ "IDs" ].eq( "" ), 0 )
    log_df[ "Nr Filtered IDs" ] = log_df[ "Filtered IDs" ].str.len().mask(
        log_df[ "Filtered IDs" ].eq( "" ), 0 )
    log_df[ "Nr Removed IDs" ] = log_df[ "Nr IDs" ] - log_df[ "Nr Filtered IDs" ]

    # for removed ids --> find out the cause and create df
    if sum(log_df["Nr Removed IDs"])==0:
        return {"Overview_Log": log_df}
    else:
        removed_ids = log_df["Removed IDs"].explode().to_list()
        df, missing = handler.get_preloaded( in_list=removed_ids, in_type="protein", organism=organism) # missing should be empty (because has already been added in the steps before)
        return {"Overview_Log":log_df, "Removed_Log": df}