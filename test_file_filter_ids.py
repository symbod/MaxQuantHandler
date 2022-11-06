import pandas as pd
from filter_ids import filter_protein_ids
from mq_utils.plotting import create_overview_plot, create_filter_detailed_plot

# data = pd.read_csv("/Users/lisiarend/Desktop/Hiwi/proteomics_analysis/data/TMT_data/rat_plasma_diabetic/proteinGroups_filtered_remapped.txt", sep=" ", header=0)
data = pd.read_csv("in/MQ_Uniprots.txt", sep=",", header=0)
#data = data[["Protein IDs"]]
#data = data[100:200]

# data = pd.DataFrame({"Protein IDs" : ["CON_A9UMVB", "M0QXN2;H7C441"]})
# data = pd.DataFrame({"Protein IDs": ["Q71DI3"]})

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

# create_overview_plot(logging=logging["Overview_Log"], out_dir="out/")
# create_filter_detailed_plot(logging=logging["Detailed_Log"], organism=organism, reviewed=reviewed, decoy=decoy,
#                             out_dir="out/")
print(logging["Detailed_Log"])