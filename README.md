# MaxQuantHandler

##Introduction
<span style="color:red">**TODO**</span>

This repository comprises four main functionalities:
- filter protein IDs
- remap gene names
- reduce gene names
- map orthologs

## Installation

You have multiple options to use the methods of this repo. You can either use the 
[mqhandler python package](https://pypi.org/project/mqhandler/) or use the [R package](https://github.com/symbod/MaxQuantHandler-R). 
Or you can simply clone this github repository.

If you want to clone this repository:

```shell
git clone git@github.com:MaxQuantHandler.git
```

and make sure you install python 3.

By cloning this github repo, you can either simply import the functions in a python script or even call the 
functions in a terminal/console. 


## Filter Protein IDs ([filter_ids.py](filter_ids.py))
For a protein assignment using MaxQuant, Fasta files are required. Since MaxQuant can also be used to run several data collectively, 
it can also happen that results are provided with protein IDs of several organisms.

This method makes it possible to check the protein IDs for their organism by directly accessing the Uniprot database, and to 
remove incorrectly assigned IDs. Additionally, decoy (REV_) and contaminants (CON_) IDs and/or unreviewed protein IDs can be removed.

One might be interested to know how many IDs were filtered out, in total and per row. Therefore, with this call, you can generate 2 data frames that display this information as a table.

In addition to the information as a table, it can also be displayed directly as plots with a simple call.

This function requires multiple arguments:
- **mandatory arguments:**
  - data: pandas dataframe
  - protein_column: name of colum with protein IDs
- **optional arguments:**
  - organism: specify organism the ids should match to
  - rev_con: bool to indicate if protein IDs from decoy (REV__) ad contaminants (CON__) should be kept
  - reviewed: bool to indicate if newly retrieved protein IDs should be reduced to reviewed ones
  - keep_empty: bool to indicate if empty ID cells should be kept or deleted
  - res_column: name of column for filter protein IDs results, if None, the protein_column will be overridden
  
  
To call the function in a python script:
```{python}
from mqhandler import filter_ids
data = pd.read_csv(<file>)
filtered_data, logging = filter_protein_ids(data=data, protein_column=protein_column, organism=organism,
                                            decoy=decoy, keep_empty=keep_empty,
                                            reviewed=reviewed, return_log=True)
```

If you want to run the function from console/terminal:

```
############################################################################
################# MaxQuantHandler - filter_ids.py ##################
        Filter proteins by organism and/or decoy/contaminants names in data file.
############################################################################

usage: python3 filter_ids.py [required arguments] [optional arguments]

required arguments:
  -d DATA, --data DATA  Data file
  -pc PROTEIN_COLUMN, --protein_column PROTEIN_COLUMN
                        Name of column with protein IDs.

optional arguments:
  -ke, --keep_empty     Bool to indicate whether empty rows should be kept. [Default=True]
  -rc RES_COLUMN, --res_column RES_COLUMN
                        Name of output column. If None, input column will be edited. [Default = None]
  -or {human,mouse,rat,rabbit}, --organism {human,mouse,rat,rabbit}
                        Specify organism the ids should match to.
  -r, --reviewed        Bool to indicate if only reviewed protein IDs should be kept.
  -rv, --rev_con        Set flag if decoy and contaminants IDs (REV__, CON__) should be kept.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################

```

## Remap Gene Names ([remap_genenames.py](remap_genenames.py))

Besides protein IDs, gene names are also taken out of the respective Fasta files and mapped. These are needed for easier naming in plots and in analytical procedures such as enrichment analysis. Unfortunately, Fasta files are not always complete in terms of gene names.

This method makes it possible to retrieve the assigned gene names based on the protein IDs with direct access to the Uniprot database and to fill the empty entries in the user file or even replace existing entries. There are multiple possible modes for which names should be taken.

Here, too, it is possible to subsequently obtain information on how many gene names were found for how many rows. This can also be displayed as a plot with a simple call.

In this tutorial, we will call the remap gene names function on the data that has already been processed using the filter IDs method.


This function requires multiple arguments:
- **mandatory arguments:**
  - data: pandas dataframe
  - mode: mode of refilling, see below for more infos
  - protein_column: name of colum with protein IDs
- **optional arguments:**
  - gene_colum: name of column with gene names
  - skip_filled: bool to indicate if already filled gene names should be skipped
  - organism: specify organism the ids should match to
  - fasta: path of fasta file when mode all or fasta
  - keep_empty: bool to ind
  - res_column: name of column for remap gene names results, if None, the gene_column will be overridden
  
There are five **modes of refilling**:
- all: use primarily fasta infos and additionally uniprot infos
- fasta: use information extracted from fasta headers
- uniprot: use mapping information from uniprot and use all gene names
- uniprot_primary: use mapping information from uniprot and only all primary gene names
- uniprot_one: use mapping information from uniprot and only use most frequent single gene name

To call the function in a python script:

```{python}
from mqhandler import remap_genenames
data = pd.read_csv(<file>)
remapped_data, logging = remap_genenames(data=data, mode=mode, protein_column=protein_column, 
gene_column=gene_column, skip_filled=skip_filled, organism=organism, fasta=fasta, keep_empty=keep_empty, 
res_column = res_column)
```

If you want to run the function from console/terminal:

```
############################################################################
################# MaxQuantHandler - remap_genenames.py ##################
                  Re-mapp gene names in data file.
############################################################################

usage: python3 remap_genenames.py [required arguments] [optional arguments]

required arguments:
  -d DATA, --data DATA  Data file
  -pc PROTEIN_COLUMN, --protein_column PROTEIN_COLUMN
                        Name of column with protein IDs.
  -m {all,fasta,uniprot,uniprot_one,uniprot_primary}, --mode {all,fasta,uniprot,uniprot_one,uniprot_primary}
                        Mode of refilling. See below for more infos.

optional arguments:
  -gc GENE_COLUMN, --gene_column GENE_COLUMN
                        Name of column with gene names [Default=None]
  -ke, --keep_empty     Bool to indicate whether empty rows should be kept. [Default=True]
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        Fasta file
  -l, --fill            Use this flag, if filled values should be skipped. [Default=True]
  -rc RES_COLUMN, --res_column RES_COLUMN
                        Name of output column. If None, input column will be edited. [Default = None]
  -or {human,mouse,rat,rabbit}, --organism {human,mouse,rat,rabbit}
                        Specify organism the ids should match to.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

----------------------------------------------------------------------------

supported modes
  all                   Use primarly fasta infos and additionally uniprot infos.
  fasta                 Use information extracted from fasta headers.
  uniprot               Use mapping information from uniprot and use all gene names.
  uniprot_primary       Use mapping information from uniprot and only all primary gene names.
  uniprot_one           Use mapping information from uniprot and only use most frequent single gene name.

############################################################################
```

## Reduce Gene Names ([reduce_genenames.py](reduce_genenames.py))

<span style="color:red">**TODO --> write description**</span>

This function requires multiple arguments:
- **mandatory arguments:**
  - data: pandas dataframe
  - mode: mode of reduction, see below for more infos
  - gene_column: name of colum with gene names
  - organism: specify organism the IDs match to
- **optional arguments:**
  - keep_empty: bool to ind
  - res_column: name of column for remap gene names results, if None, the gene_column will be overridden
  - HGNC_mode: mode on how to reduce the gene names using HGNC (mostfrequent or all)
  
 
There are four **modes of reduction**:
 - ensembl: use gProfiler to reduce gene names to those having an Ensembl ID
 - HGNC: use HGNC database to reduce gene names to those having an entry in HGNC (only for human)
 - mygeneinfo: Use mygeneinfo database to reduce gene names to those having an entry in mygeneinfo
 - enrichment: Use gProfiler to reduce gene names to those having a functional annotation
 
```{python}
from mqhandler import reduce_genenames
data = pd.read_csv(<file>)
reduced_data, logging = reduce_genenames(data=data, mode=mode, gene_column=gene_column, 
organism=organism, keep_empty=keep_empty, res_column = res_column, HGNC_mode = HGNC_mode)
```
 
If you want to run the function from console/terminal:

```
############################################################################
################# MaxQuantHandler - reduce_genenames.py ##################
                  Reduce gene names in data file.
############################################################################

usage: python3 reduce_genenames.py [required arguments] [optional arguments]

required arguments:
  -d DATA, --data DATA  Data file
  -or {human,mouse,rat,rabbit}, --organism {human,mouse,rat,rabbit}
                        Specify organism the ids should match to.
  -gc GENE_COLUMN, --gene_column GENE_COLUMN
                        Name of column with gene names.
  -m {ensembl,mygeneinfo,HGNC,enrichment}, --mode {ensembl,mygeneinfo,HGNC,enrichment}
                        Mode of reducing. See below for more infos.

optional arguments:
  -ke, --keep_empty     Bool to indicate whether empty rows should be kept. [Default=True]
  -hm {mostfrequent,all}, --hgnc_mode {mostfrequent,all}
                        What to do if reduce_mode is HGNC. Take most frequent gene name or all. [Default=mostfrequent]
  -rc RES_COLUMN, --res_column RES_COLUMN
                        Name of output column. If None, input column will be edited. [Default = None]
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

----------------------------------------------------------------------------

supported modes
  ensembl               Use gProfiler to reduce gene names to those have an Ensembl ID.
  mygeneinfo            Use mygeneinfo database to reduce gene names to those having an entry in mygeneinfo.
  HGNC                  Use HGNC database to reduce gene names to those having an entry in HGNC (only for human).
  enrichment            Use gProfiler to reduce gene names to those having a functional annotation.

############################################################################

```



## Map Orthologs ([map_orthologs.py](map_orthologs.py))

Suppose you want to compare data between organisms, for example if you want to do a review across several species, you come across a known problem. Gene names differ between species, making it necessary to map all IDs to a selected organism through an ortholog mapping.

Using the commonly used gProfiler, this method simply maps the gene names from the current organism to the target organism.

Unfortunately, depending on the original and target organism, there are more or less cases where no orthologous gene could be found. For a simplified overview of how many cases this was the case, this method can be used to obtain this information.

As with the previous tasks, the log information can be displayed in plots.


This function requires multiple arguments:
- **mandatory arguments:**
  - data: pandas dataframe
  - gene_column: name of colum with gene names
  - organism: specify organism the IDs match to
  - tar_organism: specify organism the IDs should me mapped to
- **optional arguments:**
  - keep_empty: bool to ind
  - res_column: name of column for remap gene names results, if None, the gene_column will be overridden
  
```{python}
from mqhandler import map_orthologs
data = pd.read_csv(<file>)
map_orthologs_data, logging = map_orthologs(data=data, genee_column=gene_column, organism=organism, 
tar_organism=tar_organism, keep_empty=keep_empty, res_column = res_column)
```
  

If you want to run the function from console/terminal:

```
############################################################################
################# MaxQuantHandler - map_orthologs.py ##################
                       Map ortholog gene names in data file.
############################################################################

usage: python3 map_orthologs.py [required arguments] [optional arguments]

required arguments:
  -d DATA, --data DATA  Data file
  -or {human,mouse,rat,rabbit}, --organism {human,mouse,rat,rabbit}
                        Specify organism the ids should match to.
  -tor {human,mouse,rat,rabbit}, --tar_organism {human,mouse,rat,rabbit}
                        Specify organism from which orthologs should be mapped.
  -gc GENE_COLUMN, --gene_column GENE_COLUMN
                        Name of column with gene names.

optional arguments:
  -ke, --keep_empty     Bool to indicate whether empty rows should be kept. [Default=True]
  -rc RES_COLUMN, --res_column RES_COLUMN
                        Name of output column. If None, input column will be edited. [Default = None]
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -h, --help            show this help message and exit

############################################################################

```