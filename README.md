# MaxQuantHandler
Multiple simple small scripts to help fill missing information in MaxQuant output using fasta files and/or retrieved uniprot mappings.
## Setup
1. Clone this repository
```shell
git clone git@github.com:MaxQuantHandler.git
```

2. Make sure you install python 3
## Main Scripts
### [filter_ids.py](filter_ids.py)
Sometimes MaxQuant maps protein ids of different organisms. Run this script to filter ids out which do not match your given organism.
```
############################################################################
################# MaxQuantHandler - filter_ids.py ##################
        Filter proteins or genes by organism and/or decoy names.
############################################################################

usage: python3 filter_ids.py [required arguments] [optional arguments]

required arguments:
  -q MAXQUANT_FILE, --maxquant_file MAXQUANT_FILE
                        MaxQuant file
  -s SINGLE_FILE, --single_file SINGLE_FILE
                        Single file
  -or {human,mouse,rat,rabbit}, --organism {human,mouse,rat,rabbit}
                        Specify organism the ids should match to.
  -i {protein,gene}, --in_type {protein,gene}
                        Define what type should be the source.

optional arguments:
  -pc PROTEIN_COLUMN, --protein_column PROTEIN_COLUMN
                        Name of column with protein IDs [Default="Protein IDs"]
  -gc GENE_COLUMN, --gene_column GENE_COLUMN
                        Name of column with gene names [Default="Gene names"]
  -a {keep,delete,fill}, --action {keep,delete,fill}
                        What to do, if IDs cell is empty after filtering. Keep empty cell, delete it or fill it based on gene name.
  -r, --reviewed        Bool to indicate if newly retrieved protein IDs should be reduced to reviewed ones.
  -d, --decoy           Set flag if protein ids from decoy fasta (REV__, CON__) should be kept.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################

```
### [remap_genenames.py](remap_genenames.py) 
Remap gene names in MaxQuant file based on fasta file and/or uniprot mappings with multiple other options.
```
############################################################################
################# MaxQuantHandler - remap_genenames.py ##################
                  Re-mapp gene names in max quant file.
############################################################################

usage: python3 remap_genenames.py [required arguments] [optional arguments]

required arguments:
  -q MAXQUANT_FILE, --maxquant_file MAXQUANT_FILE
                        MaxQuant file
  -s SINGLE_FILE, --single_file SINGLE_FILE
                        Single file
  -m {all,fasta,uniprot,uniprot_one,uniprot_primary}, --mode {all,fasta,uniprot,uniprot_one,uniprot_primary}
                        Mode of refilling. See below for more infos.

optional arguments:
  -pc PROTEIN_COLUMN, --protein_column PROTEIN_COLUMN
                        Name of column with protein IDs [Default="Protein IDs"]
  -gc GENE_COLUMN, --gene_column GENE_COLUMN
                        Name of column with gene names [Default="Gene names"]
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        Fasta file
  -l, --fill            Use this flag, if filled values should be skipped. [Default=True]
  -or {human,mouse,rat,rabbit}, --organism {human,mouse,rat,rabbit}
                        Specify organism the ids should match to.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

----------------------------------------------------------------------------

supported modes
  all			Use primarly fasta infos and additionally uniprot infos.
  fasta			Use information extracted from fasta headers.
  uniprot		Use mapping information from uniprot and use all gene names.
  uniprot_primary	Use mapping information from uniprot and only all primary gene names.
  uniprot_one		Use mapping information from uniprot and only use most frequent single gene name.

############################################################################

```
## [map_orthologs.py](map_orthologs.py) 
```
############################################################################
################# MaxQuantHandler - map_orthologs.py ##################
                       Get ortholog gene names.
############################################################################

usage: python3 map_orthologs.py [required arguments] [optional arguments]

required arguments:
  -q MAXQUANT_FILE, --maxquant_file MAXQUANT_FILE
                        MaxQuant file
  -s SINGLE_FILE, --single_file SINGLE_FILE
                        Single file
  -or {human,mouse,rat,rabbit}, --organism {human,mouse,rat,rabbit}
                        Specify organism the ids are mapped to.
  -tor {human,mouse,rat,rabbit}, --tar_organism {human,mouse,rat,rabbit}
                        Specify organism from which orthologs should be mapped.

optional arguments:
  -pc PROTEIN_COLUMN, --protein_column PROTEIN_COLUMN
                        Name of column with protein IDs [Default="Protein IDs"]
  -gc GENE_COLUMN, --gene_column GENE_COLUMN
                        Name of column with gene names [Default="Gene names"]
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################


```

## Smaller Scripts
### [fasta_grepper.py](fasta_grepper.py)
Get informations stored in the headers of the given fasta file and saving result to grepped_info.csv
```
############################################################################
################### MaxQuantHandler - fasta_grepper.py #####################
                   Grep protein info from fasta file.
############################################################################

usage: python3 fasta_grepper.py [required arguments] [optional arguments]

required arguments:
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        Fasta file

optional arguments:
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################
```

### [get_uniprot_mapping.py](get_uniprot_mapping.py) 
Get a file with uniprot mappings to either protein ids or gene names in the file. Optionally those can be filtered by given organism, which is highly recommended for gene names option.
```
############################################################################
################# MaxQuantHandler - get_uniprot_mapping.py ##################
  Get uniprot mapping to protein ids or gene names optionally by organism.
############################################################################

usage: python3 get_uniprot_mapping.py [required arguments] [optional arguments]

required arguments:
  -q MAXQUANT_FILE, --maxquant_file MAXQUANT_FILE
                        MaxQuant file
  -i {proteinID,genename}, --in_type {proteinID,genename}
                        Define what type should be the source.

optional arguments:
  -r {human,rat}, --organism {human,rat}
                        Specify organism the ids should match to.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################

```