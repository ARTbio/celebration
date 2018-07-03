## functions

Down_Sample_Counts <- function (expr_vector, prob) 
    {
        down_sample <- function(x, prob.=prob) {
            return(unlist(lapply(x, function(y) {
                rbinom(1, y, prob.)
            })))
        }
        down_sampled_df <- apply(as.data.frame(expr_vector), 2, down_sample)
        down_sampled_df = data.frame(down_sampled_df)
        row.names(down_sampled_df) = rownames(expr_vector)
        return(down_sampled_df)
    }

prob.range = c(0.001, 0.002, 0.004, 0.008, 0.01, 0.02, 0.04, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)


Downsampling.curve <- function(exprV, prob.range){
    res_matrix <- data.frame(matrix(ncol = length(prob.range), nrow = dim(res_matrix)[1]))
    i = 0
    for (prob in prob.range) {
        i = i + 1
        res_matrix[, i] = Down_Sample_Counts(exprV, prob)
        colnames(res_matrix) <- prob.range
        }
    return(res_matrix)
    }
    
Compute_features_matrix <- function(rawcounts, prob.range){
    Features.matrix = data.frame(matrix(ncol = dim(rawcounts)[2], nrow = length(prob.range)))
    i = 0
    for (label in colnames(rawcounts)) {
        i = i + 1
        Features.matrix[, i] = colSums(Downsampling.curve(rawcounts[, i], prob.range) != 0)
        # or >= 5, etc
        colnames(Features.matrix) = colnames(rawcounts)
        if (i == 1) {
            plot(prob.range, Features.matrix[, i], type="l")
            } else {
            lines(prob.range, Features.matrix[, i])
        }
    }
    return(Features.matrix)
}

Visualise_features_matrix <- function(features_matrix, prob.range, Quality_Treshold=20000) {
    i = 0
    max = max(tail(features_matrix, 1))
    Vmaxs = unlist(features_matrix[2,]/.002)
    colors = ifelse(Vmaxs<Quality_Treshold, "red", "black")
    for (label in features_matrix) {
        i = i + 1
        if (i == 1) {
            plot(prob.range, features_matrix[, i], type="l", ylim=c(0,max), col=colors[i])
            } else {
            lines(prob.range, features_matrix[, i], col=colors[i])
        }
    }
    return(Vmaxs)
}

## main

CountArray = read.delim("raw_counts.tab", header=T)
rownames(CountArray) = CountArray[,1]
CountArray = CountArray[,-1]
mymatrix = Compute_features_matrix(CountArray, prob.range)
vmaxs = Visualise_features_matrix(mymatrix, prob.range)