# What it does

Normalization takes a tabular table of raw read counts (it's better to use a filtered matrix), each column corresponding to a sample count, and returns a normalized table of counts.  
 
You can normalized your data :  
 - Using R package scater : Log Normalization to either size factors or library size 
 - Using R package Seurat : Log2 normalization to the median of UMI counts.

NB : To be able to use size factor with scater you need to have at least 100 cells in your dataset.

# Generation of output test files 

The different outputs in `test-data/` have been generated thanks to the Rscript normalization_singleCells.R and the input file `test-data/counts.tab`.

> Rscript normalization.R -f test-data/counts.tab -c TRUE -s tab-m scater --min_mean 0.1 --size_factor FALSE -o test-data/counts_normalized_scater.tab

> Rscript normalization.R -f test-data/counts.tab -c TRUE -s tab-m Seurat -o test-data/counts_normalized_Seurat.tab
