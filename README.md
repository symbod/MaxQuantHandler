# MaxQuantHandler
Multiple simple small scripts to help fill missing information in MaxQuant output using fasta files and/or retrieved uniprot mappings.
## Setup
1. Clone this repository
```
git clone git@github.com:MaxQuantHandler.git
```
2. Make sure you install python 3
## Smaller Scripts
### fasta_grepper
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
### filter_ids
Sometimes MaxQuant maps protein ids of different organisms. Run this script to filter ids out which do not match your given organism.
```
############################################################################
#################### MaxQuantHandler - filter_ids.py #######################
                   Filter protein ids by organism.
############################################################################

usage: python3 filter_ids.py [required arguments] [optional arguments]

required arguments:
  -r {human,rat}, --organism {human,rat}
                        Specify organism the ids should match to.
  -q MAXQUANT_FILE, --maxquant_file MAXQUANT_FILE
                        MaxQuant file

optional arguments:
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################
```
### get_uniprot_mapping 
Get a file uniprot mappings to protein ids in the file. Optionally those can be filtered by given organism.
```
############################################################################
################### MaxQuantHandler - get_uniprot_mapping.py ###############
                   Get uniprot mapping to protein ids optionally by organism.
############################################################################

usage: python3 get_uniprot_mapping.py [required arguments] [optional arguments]

required arguments:
  -q MAXQUANT_FILE, --maxquant_file MAXQUANT_FILE
                        MaxQuant file

optional arguments:
  -r {human,rat}, --organism {human,rat}
                        Specify organism the ids should match to.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################
```
## Main Scripts
### remap_genenames
Remap gene names in MaxQuant file based on fasta file and/or uniprot mappings with multiple other options.
```
############################################################################
#################### MaxQuantHandler - remap_genenames.py ##################
                   Re-mapp gene names in max quant file.
############################################################################

usage: python3 remap_genenames.py [required arguments] [optional arguments]

required arguments:
  -q MAXQUANT_FILE, --maxquant_file MAXQUANT_FILE
                        MaxQuant file
  -m {all,fasta,uniprot,uniprot_one}, --mode {all,fasta,uniprot,uniprot_one}
                        Mode of refilling. See below for more infos.

optional arguments:
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        Fasta file
  -l, --fill            Use this flag, if only missing values should be filled.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

----------------------------------------------------------------------------

supported modes
  all			Use primarly fasta infos and additionally uniprot infos.
  fasta			Use information extracted from fasta headers.
  uniprot		Use mapping information from uniprot and use all gene names.
  uniprot_one		Use mapping information from uniprot and only use most frequent single gene name.

############################################################################
```
