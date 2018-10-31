## Normalisation script for scRNAseq data

## Setup R error handling to go to stderr
options(
  show.error.messages = F,
  error = function (){
    cat(geterrmessage(), file = stderr())
    q("no", 1, F)
  }
)

option_list <- list(
  optparse::make_option(c("-f", "--file"),
                        type = "character",
                        help = "Path to the input file"),
  optparse::make_option(c("-c", "--colnames"),
                        type = "logical", help =
                          "Consider first line as header"),
  optparse::make_option(c("-s", "--sep"),
                        type = "character",
                        help = "Input file separator"),
  optparse::make_option(c("-m", "--method"),
                        type = "character",
                        help = "Name of the method/package used for the
                        normalization : scater or Seurat"),
  optparse::make_option("--min_mean",
                        type = "character",
                        help = "Minimum average count of genes to be used
                        for normalization with scater"),
  optparse::make_option("--size_factor",
                        type = "logical",
                        help = "Using size factors based on deconvolution
                        from cell pools (TRUE) or considering library sizes
                        as size Factors (FALSE) ? NB: Mandatory if there is
                        less than 100 cells in dataset, you must choose library
                        sizes"),
  optparse::make_option(c("-o", "--output_matrix"),
                        type = "character",
                        help = "Path to the filtered matrix")
)

parser <- optparse::OptionParser(usage = "%prog [options] file",
                                 option_list = option_list)
args <- optparse::parse_args(parser)

## Import dataset
data.counts <- read.table(
  args$file,
  header = args$colnames,
  stringsAsFactors = F,
  sep = ifelse(args$sep == "tab", "\t", args$sep),
  check.names = FALSE,
  row.names = 1
)

## Checking parameters
if (ncol(data.counts) < 100 & args$size_factor == TRUE){
  stop("You need at least 100 cells in your dataset
       to use the deconvolution method.")
}

if (args$method == "scater"){
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(data.counts)))
  if (args$size_factor == TRUE){
    sce <- scran::computeSumFactors(
      sce,
      min.mean = args$min_mean,
      sf.out = FALSE,
      get.spikes = FALSE,
      assay.type = "counts"
    )
  }
  sce <- scater::normalize(sce)
  normalized_data_counts <- SingleCellExperiment::logcounts(sce)
}

if (args$method == "Seurat"){
  sce <- Seurat::CreateSeuratObject(raw.data = data.counts)
  sce <- Seurat::NormalizeData(sce,
                               scale.factor = median(sce@meta.data[, "nUMI"]))
  normalized_data_counts <- as.matrix(sce@data)
}

write.table(
  normalized_data_counts,
  file = args$output_matrix,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)