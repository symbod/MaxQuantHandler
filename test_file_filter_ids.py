import pandas as pd
from filter_ids import filter_protein_ids

# data = pd.read_csv("/Users/lisiarend/Desktop/Hiwi/proteomics_analysis/data/TMT_data/rat_plasma_diabetic/proteinGroups_filtered_remapped.txt", sep=" ", header=0)
# data = data[["Protein IDs"]]
# data = data[10:100]

# data = pd.DataFrame({"Protein IDs" : ["A9UMVB", "M0QXN2;H7C441"]})
data = pd.DataFrame({"Protein IDs": ["A9UMVB", "A8MTI9;Q96A49"]})

organism = "human"  # Specify organism the ids should match to
in_type = "protein"  # Define what type should be the source
protein_column = "Protein IDs"  # Name of column with protein IDs
gene_column = "Gene names"  # Name of column with gene names
keep_empty = True  # Bool to indicate if empty rows should be kept
reviewed = False  # Bool to indicate if newly retrieved protein IDs should be reduced to reviewed ones
decoy = False  # Bool to indicate if protein ids from decoy fasta (REV__, CON__) should be kept

filtered_data, logging = filter_protein_ids(data=data, protein_column=protein_column, organism=organism,
                                            decoy=decoy, keep_empty=keep_empty,
                                            reviewed=reviewed, return_log=True)
print(filtered_data)
