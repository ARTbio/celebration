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


Downsampling.curve <- function (exprV, prob.range) {
    res_matrix <- data.frame(matrix(ncol = length(prob.range), nrow = dim(res_matrix)[1]))
    i = 0
    for (prob in prob.range) {
        i = i + 1
        res_matrix[, i] = Down_Sample_Counts(exprV, prob)
        colnames(res_matrix) <- prob.range
        }
#    plot(prob.range, colSums(res_matrix != 0))
#    plot(prob.range, colSums(res_matrix > 4))
#    return(colSums(res_matrix != 0))
    return(res_matrix)
    }
    
Downsampling.dataframe = data.frame(row.names = prob.range)
for (label in colnames(CountArray)) {
    newdownsampling = Downsampling.curve(CountArray[,label], prob.range)
    colnames(newdownsampling) = label
    Downsampling.dataframe = cbind(Downsampling.dataframe, newdownsampling)
    }