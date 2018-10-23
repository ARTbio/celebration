## Filtering methods on single cell data
# required packages : scater, SingleCellExperiment, Seurat

data = "./test-data/counts_with_mito.tab" #input file
data.sep = "\t"                           #input file separator
data.header = TRUE                        #Does your input file have a header (column names) or not ?
method = "scater"                         #what method do you want to use ? scater, Seurat
mito = TRUE                               # Is there mitochondrial genes in dataset ? : TRUE or FALSE
if(mito == TRUE) header.mito = "MT-"      #How to recognize mitochondrial genes in gene set. 
if(method == "Seurat"){
  min.cells = 3                           #Keep all genes that at least detect in n cells
  min.genes = 10                          #Keep all cells that detect at least n genes
}
if(method == "scater") {
  gene_filter = "low.abundances"          #what method do you want to use for filter genes ? low.abundances or min.cells
  min.cells = 3                           #Keep all genes that at least detect in n cells, threshold for gene filtering
}

#Import dataset
data.counts = read.table(
  data,
  header = data.header,
  stringsAsFactors = F,
  sep = data.sep,
  check.names = FALSE,
  row.names = 1
)


pdf(file = paste(tools::file_path_sans_ext(data), ifelse(method == "scater", paste(method, gene_filter, sep = "_"),method), "filter_plots.pdf", sep = "_"), paper = "a4")
if(method == "scater") {
  sce <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(data.counts)))
  
  if (mito == TRUE) {
    #retrieve mitochondrial genes in the dataset
    mito.genes = grep(header.mito, rownames(sce))
    if (length(mito.genes) == 0)
      stop(paste("No genes match mitochondrial pattern :", header.mito))
    
    #calculate QC metrics
    sce <-
      scater::calculateQCMetrics(sce, feature_controls = list(Mito = mito.genes))
    
    #search which cells had a high level of expressed mitochondrial genes
    high.mito <-
      scater::isOutlier(sce$pct_counts_Mito, nmads = 3, type = "higher")
    
    #remove those cells
    sce <- sce[,!high.mito]
    
    #some QC plots
    hist(
      sce$total_counts,
      breaks = 20,
      col = "grey80",
      xlab = "Log-total UMI count"
    )
    hist(
      sce$log10_total_features,
      breaks = 20,
      col = "grey80",
      xlab = "Log-total number of expressed features"
    )
    hist(
      sce$pct_counts_Mito,
      breaks = 20,
      col = "grey80",
      xlab = "Proportion of counts in mitochondrial genes"
    )
    #scater::plotExprsFreqVsMean(sce, feature_controls = rowData(sce)$is_feature_control_Mito) #Error due to feature controls in rowData, resolved in version with R>=3.5
    #Inspecting the most highly expressed genes
    scater::plotQC(sce, type = "highest-expression", n = ifelse(nrow(sce) < 50, nrow(sce), 50))
    
  } else{
    #calculate QC metrics
    sce = scater::calculateQCMetrics(sce)
    
    #search which cells had a low counts
    libsize.drop <-
      scater::isOutlier(sce$total_counts,
                        nmads = 3,
                        type = "lower",
                        log = TRUE)
    
    #search for cells with few detected genes
    feature.drop <-
      scater::isOutlier(
        sce$total_features,
        nmads = 3,
        type = "lower",
        log = TRUE
      )
    
    #remove low quality cells
    sce <- sce[,!(libsize.drop | feature.drop)]
    print(data.frame(
      ByLibSize = sum(libsize.drop),
      ByFeature = sum(feature.drop),
      Remaining = ncol(sce)
    ))
    
    #Some QC plots
    hist(
      sce$total_counts,
      breaks = 20,
      col = "grey80",
      xlab = "Log-total UMI count"
    )
    hist(
      sce$log10_total_features,
      breaks = 20,
      col = "grey80",
      xlab = "Log-total number of expressed features"
    )
    #Verify that the frequency of expression (i.e., number of cells with non-zero expression) and the mean are positively correlated
    scater::plotQC(sce, type = "exprs-freq-vs-mean")
    #Inspecting the most highly expressed genes
    scater::plotQC(sce, type = "highest-expression", n = ifelse(nrow(sce) < 50, nrow(sce), 50))
  }
  if (gene_filter == "min.cells") {
    numcells <- scater::nexprs(sce, byrow = TRUE)
    #Filter genes detected in less than n cells
    numcells2 <- numcells >= min.cells
    sce <- sce[numcells2, ]
  }
  if (gene_filter == "low.abundances") {
    ave.counts <- scater::calcAverage(sce)
    
    num.cells <- scater::nexprs(sce, byrow = TRUE)
    smoothScatter(
      log10(ave.counts),
      num.cells,
      ylab = "Number of cells",
      xlab = expression(Log[10] ~ "average count")
    )
    
    to.keep <- num.cells > min.cells
    sce <- sce[to.keep, ]
    print(summary(to.keep))
  }
}
if (method == "Seurat") {
  sce = Seurat::CreateSeuratObject(raw.data = data.counts,
                                   min.cells = min.cells,
                                   min.genes = min.genes)
  if (mito == TRUE) {
    #retrieve mitochondrial genes in the dataset
    mito.genes = grep(header.mito, rownames(sce@raw.data))
    if (length(mito.genes) == 0)
      stop(paste("No genes match mitochondrial pattern :", header.mito))
    
    #calculate QC metrics
    percent.mito <-
      Matrix::colSums(sce@raw.data[mito.genes,]) / Matrix::colSums(sce@raw.data)
    sce <-
      AddMetaData(object = sce,
                  metadata = percent.mito,
                  col.name = "percent.mito")
    
    #Filter low quality cells
    sce <-
      FilterCells(sce, subset.names = "percent.mito", high.thresholds = 0.2)
    
    #QC plot after filtering
    print(Seurat::VlnPlot(
      object = sce,
      features.plot = c("nGene", "nUMI", "percent.mito"),
      nCol = 3
    ))
    
    
  } else{
    print(Seurat::VlnPlot(
      object = sce,
      features.plot = c("nGene", "nUMI"),
      nCol = 2
    ))
  }
  
  #Some QC plots
  hist(
    sce@meta.data$nUMI,
    breaks = 20,
    col = "grey80",
    main = "Number of UMI/Cells"
  )
  hist(
    sce@meta.data$nGene,
    breaks = 20,
    col = "grey80",
    main = "Number of genes detected/Cells"
  )
  
}

dev.off()

save(sce,
     file = paste(
       tools::file_path_sans_ext(data),
       ifelse(method == "scater", paste(method, gene_filter, sep = "_"), method),
       "filter_sce.rds",
       sep = "_"
     ))


if (method == "scater") {
  filtered_data_counts = sce@assays$data$counts
} else filtered_data_counts = sce@data

write.table(
  filtered_data_counts,
  file = paste(
    tools::file_path_sans_ext(data),
    ifelse(method == "scater", paste(method, gene_filter, sep = "_"), method),
    "filtered.tab",
    sep = "_"
  ),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)
