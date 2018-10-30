#What it does

 Filter_cells_and_genes takes a tabular table of raw read counts, each column corresponding to a sample count, and returns (1) a PDF with QC (quality controls) plots and (2) a filtered table of counts.  
 
 You can filter out low quality cells from your raw expression table :  
 - Based on expression level of mitochondrial genes :  
     - Using R package scater: remove cells that detect mitochondrial genes 3 times higher than the median of all cells.  
    - Using R package Seurat: if cells expressed more than 20% of mitochondrial genes, they are removed.  
 - Using R package scater : search cells with low library size and few detected genes based on the deviation from the median of all cells  
 - Using R package Seurat : Remove all cells that detected at least n genes.  
 
 Also you can filter out low quality genes from your raw expression table :  
 - Using R package scater :  
     - Remove all genes that aren't detected in at least n cells  
    - Filter out genes based on its average count (which is adjust for library size between cells). It is more stringent.  
 - Using R package Seurat : Remove all genes that aren't detected in at least n cells.  

#Generation of output test files 

The different outputs in `test-data/` have been generated thanks to the Rscript filter_cells.R and either of two inputs files `test-data/counts.tab` or `test-data/counts_with_mito.tab`

> Rscript filter_cells.R -f test-data/counts.tab -c TRUE -s tab -m FALSE -p scater -g low.abundances --min_cells 1 --output_matrix test-data/counts_scater_low.abundances_filtered.tab --output_pdf test-data/counts_scater_low.abundances_filter_plots.pdf

> Rscript filter_cells.R -f test-data/counts.tab -c TRUE -s tab -m FALSE -p scater -g min.cells --min_cells 3 --output_matrix test-data/counts_scater_min.cells_filtered.tab --output_pdf test-data/counts_scater_min.cells_filter_plots.pdf

> Rscript filter_cells.R -f test-data/counts.tab -c TRUE -s tab -m FALSE -p Seurat --min_cells 3 --min_genes 5 --output_matrix test-data/counts_Seurat_filtered.tab --output_pdf test-data/counts_Seurat_filter_plots.pdf

> Rscript filter_cells.R -f test-data/counts_with_mito.tab -c TRUE -s tab -m TRUE --header_mito MT- -p Seurat --min_cells 3 --min_genes 5 --output_matrix test-data/counts_with_mito_Seurat_filtered.tab --output_pdf test-data/counts_with_mito_Seurat_filter_plots.pdf

> Rscript filter_cells.R -f test-data/counts_with_mito.tab -c TRUE -s tab -m TRUE --header_mito MT- -p scater -g low.abundances --min_cells 1 --output_matrix test-data/counts_with_mito_scater_low.abundances_filtered.tab --output_pdf test-data/counts_with_mito_scater_low.abundances_filter_plots.pdf

> Rscript filter_cells.R -f test-data/counts_with_mito.tab -c TRUE -s tab -m TRUE --header_mito MT- -p scater -g min.cells --min_cells 3 --output_matrix test-data/counts_with_mito_scater_min.cells_filtered.tab --output_pdf test-data/counts_with_mito_scater_min.cells_filter_plots.pdf