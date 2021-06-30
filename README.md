# Phylo Distribution

Allow to create a figure similar to the figure 4 of Denise et al 2020 PLoS Biology

## Requirement

Python 3.8

### Libraries

pandas==1.2.0  
numpy==1.18.3  
matplotlib==3.3.3  
biopython==1.76  
seaborn==0.11.1  

## Input format

### Report 

This is a tabular file, the separator is a tabulation (\t). This is a file that contains all the hits that you have for the set of genomes that you use as database.

Need to have a column with an identifier for the genome (chromosome or plasmid) [`Replicon_name`] and a column that contains the annotation that you want to see the distribution [`Your_column_name`]

#### Exemple

Replicon_name	Gene  
ACAC001.D.00001.C001	serC  
ACAC001.D.00001.C001	pdxP  
ACAC001.D.00001.P002	yggS  
ACAC001.D.00001.P002	DXS  
BAMY001.D.00002.C001	pdxP  
BAMY001.D.00002.C001	pdxS  
BAMY001.D.00002.C001	pdxT  
BAMY001.D.00003.C001	pdxT  
BAMY001.D.00003.C001	pdxS  

### INFO_TAB

This is a tabular file, the separator is a tabulation (\t). This is the full description of all the species of the database that you used.

The information needed is the same identifier for the genome (chromosome or plasmid) as the Report file [`Replicon_name`], an identifier for the species to know which replicons belong to the same species [`Species_Id`] and the full lineage of the species [`Lineage`]

#### Exemple

Replicon_name	Species_Id	Lineage  
ACAC001.D.00001.C001	ACAC001.D.00001	Bacteria;Actinobacteria;Propionibacteriales;Propionibacteriaceae;Acidipropionibacterium  
ACAC001.D.00001.P002	ACAC001.D.00001	Bacteria;Actinobacteria;Propionibacteriales;Propionibacteriaceae;Acidipropionibacterium  
BAMY001.D.00002.C001	BAMY001.D.00002	Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus cereus group  
BAMY001.D.00003.C001	BAMY001.D.00003	Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus cereus group  


### Tree file 

A tree of the phyla you want to plot the distribution against in newick format. The phyla need to be present in the lineage of the genomes present in your database.

#### Exemple 

((Phylum_A,Phylum_B),(Phylum_C,Phylum_D),Phylum_E)

## Help

More help with the option -h/--help 

```
usage: phylo_distribution.py [-h] -r <REPORT> -i <INFO_TAB> -c <NAME OF THE COLUMNS OF INTEREST> -t <FILE_TREE>
                             [-o <FOLDER>] [-order <GENE_NAME/SYSTEM_NAME> [<GENE_NAME/SYSTEM_NAME> ...]]
                             [-size <SIZE>] [-main_rowspan <INT>] [-tree_colspan <INT>] [-heatmap_colspan <INT>]
                             [-barplot_rowspan <INT>] [-base_span <INT>] [-height <INT>] [-width <INT>] [-p]
                             [-color <COLOR>]

Phylegentic Dictribution Figure Creation

optional arguments:
  -h, --help            show this help message and exit

General input dataset options:
  -r <REPORT>, --report <REPORT>
                        Tabulated File that contain the detection of the genes/systems with the columns
                        Replicon_name (name of the replicon [for exemple each chromosome of a species is a
                        different replicon]) and the columns that you want to see the distribution
  -i <INFO_TAB>, --infotab <INFO_TAB>
                        Name of the tabulated file that contain the lineage information for each replicon. Need
                        to have a columns name Replicon_name (same as in the report), Species_Id (could be the
                        name of the species, but need to be unique for each genome), Lineage (full lineage of the
                        genome)
  -c <NAME OF THE COLUMNS OF INTEREST>, --columns_distribution <NAME OF THE COLUMNS OF INTEREST>
                        The name of the column in the report that will be use to do the distribution. The column
                        contains the annotation for the sequence
  -t <FILE_TREE>, --tree_file <FILE_TREE>
                        The name of tree file of the phyla in newick format
  -o <FOLDER>, --output <FOLDER>
                        Name of the figure where to put the database (default: $PWD/phylo_distribution.pdf)
  -order <GENE_NAME/SYSTEM_NAME> [<GENE_NAME/SYSTEM_NAME> ...], --order_columns <GENE_NAME/SYSTEM_NAME> [<GENE_NAME/SYSTEM_NAME> ...]
                        Ordered list of the genes/systems you want in the heatmap (default: The order of
                        appearance in the report)

Figure options:
  -size <SIZE>, --size_font <SIZE>
                        Font size in points or as a string (e.g., 'large') (default: 8pt)
  -main_rowspan <INT>, --main_rowspan <INT>
                        Number of row to span the main figure (tree and heatmap) (default: 5)
  -tree_colspan <INT>, --tree_colspan <INT>
                        Number of columns to span the tree (default: 2)
  -heatmap_colspan <INT>, --heatmap_colspan <INT>
                        Number of columns to span the heatmap (default: 8)
  -barplot_rowspan <INT>, --barplot_rowspan <INT>
                        Number of rows to span the barplot under the heatmap (default: 1)
  -base_span <INT>, --base_span <INT>
                        Number of rows/columns to span that are not cited above (defaul: 1)
  -height <INT>, --height <INT>
                        Total height of the figure in point (default: 10)
  -width <INT>, --width <INT>
                        Total width of the figure in point (default: 12)
  -color <COLOR>, --color_heatmap <COLOR>
                        Color of the gradient in the heatmap (default: red)
```
