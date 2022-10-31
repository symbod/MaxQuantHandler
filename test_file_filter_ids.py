import pandas as pd
from filter_ids import filter_protein_ids


data = pd.read_csv("/Users/lisiarend/Desktop/Hiwi/proteomics_analysis/data/TMT_data/rat_plasma_diabetic/proteinGroups_filtered_remapped.txt", sep=" ", header=0)
data = data[["Protein IDs"]]
data = data[10:20]



organism = "rat" # Specify organism the ids should match to
in_type = "protein" #  Define what type should be the source
protein_column = "Protein IDs" # Name of column with protein IDs
gene_column = "Gene names" # Name of column with gene names
action = "delete" # What to do, if IDs cell is empty after filtering. Keep empty cell, delete it or fill it based on gene name.
reviewed = True # Bool to indicate if newly retrieved protein IDs should be reduced to reviewed ones
decoy = False # Bool to indicate if protein ids from decoy fasta (REV__, CON__) should be kept


filtered_data, logging = filter_protein_ids(data = data, id_column = protein_column, organism = organism,
                                      decoy = decoy, action = action, gene_column = gene_column,
                                      reviewed = reviewed, return_log=True)
print(logging)


