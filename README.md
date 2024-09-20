# Data Load For Computational Bio & Pharma - Working Group for Rare Disease

This script loads the data used by the mocomakers meetup [Computational Bio & Pharma - Working Group for Rare Disease](https://www.meetup.com/mocomakers/events/301844982/).

The data files should be placed in the `data` directory. The sources are:
- `approved_compounds.txt` - approved compounds for difference diseases - text copied from [NIH/NCI web site](https://www.cancer.gov/about-cancer/treatment/types/targeted-therapies/approved-drug-list)
- `family_target_ligand_all_gtopdb.csv` - info on compounds - downloaded from [hugging face](https://huggingface.co/datasets/dwb2023/gtopdb-target-family-interaction-detail)    
- `OmicsSomaticMutationsMatrixDamaging.csv` - the mutation file that holds medication and target pivot table indicating level of mutation - for access to the input data files, please contact Matt Zamora via `matt@mocomakers.com`.  
- `secondary-screen-dose-response-curve-parameters.csv` - the response file that has medication effect on disease and target - for access to the input data files, please contact Matt Zamora via `matt@mocomakers.com`.

The output of this python script will upload the 2 data tables into a sqlite database called `rare_disease.db`. Please note that the database calculation may need about 2TB free disk on the hard drive - SSD preferred.

It will also create a tables that merge and process the data from both tables and outputs:
- `response` - the response data file loaded.
- `mutation_matrix_full` - the mutation file loaded.
- `select_genes2` - a list of select genes that have a mutation of level 2.
- `select_mutation_matrix2` -  a subset of the mutation matrix containing only genes that have a mutation of level 2.
- `mutated2` - a join of the response file and the subset of the mutation matrix containing only select genes.  
- `aggregated_screen_id_gene_tissue_moa_2` - an aggregation of `mutated2` that is grouped by screen_id, gene, tissue, moa and calculates some statistics of the response results.   
- `select_aggregated_screen_id_gene_tissue_moa_2_0` - a subset of `aggregated_screen_id_gene_tissue_moa_2` where mutation value is 0.
- `select_aggregated_screen_id_gene_tissue_moa_2_1` - a subset of `aggregated_screen_id_gene_tissue_moa_2` where mutation value is 1.
- `select_aggregated_screen_id_gene_tissue_moa_2_2` - a subset of `aggregated_screen_id_gene_tissue_moa_2` where mutation value is 2.
- `diff_aggregated_screen_id_gene_tissue_moa_name_family_types_family_names_approved_2_approved` -  a joined table containing the difference between `select_aggregated_screen_id_gene_tissue_moa_2_2` and `select_aggregated_screen_id_gene_tissue_moa_2_0` for each group of screen_id, gene, tissue, moa, name, family_types, family_names, and approved.    

The differences in `diff_aggregated_screen_id_gene_tissue_moa_2` can be used for further analysis.

LICENSE
-------
<a rel="license" href="http://creativecommons.org/publicdomain/zero/1.0/"> <img src="https://licensebuttons.net/p/zero/1.0/88x31.png" style="border-style: none;" alt="CC0" />  </a>

To the extent possible under law, Jacob Barhak has waived all copyright and 
related or neighboring rights to Data Load For Computational Bio & Pharma - Working Group for Rare Disease
This work is published from: Israel.
