### R script for quantile normalization 
### quantile normalization (individual gene with mean of 0 and sd of 1)

suppressMessages(library("argparse"))


parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Given sample read count files, return normalized read count data in one combined file.")
    parser$add_argument("-norm_count", "--norm_count", help = "library-size normalized read count in tab-separated file. ")
    parser$add_argument("-O", "--output", help = "output prefix. ")
    parser$parse_args()
}

###
argv <- parse_arg()
norm_count <- argv$norm_count
out <- argv$output

norm_df <- read.table(norm_count, sep = "\t", header = T, row.names = 1)
norm_mat <- as.matrix(norm_df)

exp_mat <- t(apply(norm_mat, 1, rank, ties.method = "average"))

quantile_norm <- qnorm(exp_mat/(ncol(exp_mat)+1))
write.table(quantile_norm, file = paste0(out, ".txt"), sep = "\t", quote = F)
