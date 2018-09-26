## Filtering methods on single cell data
# required packages : scater, SingleCellExperiment

data = "/home/lea/Darmanis_data/GBM_raw_gene_counts.csv"
data.sep = " "
data.header = TRUE
method = "scater" #what method do you want to use ? scater, Seurat
mito = FALSE # Is there mitochondrial genes in dataset ? : TRUE or FALSE
if(mito == TRUE) header.mito = "" #How to recognize mitochondrial genes in gene set. 

#Import dataset
data.counts = read.table(
  data,
  header = data.header,
  stringsAsFactors = F,
  sep = data.sep,
  check.names = FALSE
)


pdf(file = "../Filter_cells_plots.pdf", paper = "a4")
if(method == "scater") {
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(data.counts)))
  
  if (mito == TRUE) {
    #retrieve mitochondrial genes in the dataset
    mito.genes = grep(header.mito, rownames(sce))
    
    #calculate QC metrics
    sce <-
      scater::calculateQCMetrics(sce, feature_controls = list(Mito = mito.genes))
    
    #search which cells had a high level of expressed mitochondrial genes
    high.mito <-
      scater::isOutlier(sce$pct_counts_Mito, nmads = 3, type = "higher")
    
    #remove those cells
    sce <- sce[, !high.mito]
    
    #some QC plots
    hist(sce$total_counts, breaks=20, col="grey80", xlab="Log-total UMI count")
    hist(sce$log10_total_features, breaks=20, col="grey80", xlab="Log-total number of expressed features")
    hist(sce$pct_counts_Mito, breaks=20,col="grey80", xlab="Proportion of counts in mitochondrial genes")
    
  }else{
    #calculate QC metrics
    sce = scater::calculateQCMetrics(sce)
    
    #search which cells had a low counts
    libsize.drop <-
      scater::isOutlier(sce$total_counts,
                nmads = 3,
                type = "lower",
                log = TRUE)
    
    #search low abundance genes
    feature.drop <-
      scater::isOutlier(
        sce$total_features,
        nmads = 3,
        type = "lower",
        log = TRUE
      )
    
    #remove low quality cells and genes
    sce <- sce[, !(libsize.drop | feature.drop)]
    print(data.frame(
      ByLibSize = sum(libsize.drop),
      ByFeature = sum(feature.drop),
      Remaining = ncol(sce)
    ))
    
    #Some QC plots
    hist(sce$total_counts, breaks=20, col="grey80", xlab="Log-total UMI count")
    hist(sce$log10_total_features, breaks=20, col="grey80", xlab="Log-total number of expressed features")
    #Inspecting the most highly expressed genes
    plotQC(sce, type = "highest-expression", n=50)
    }
}

dev.off()


