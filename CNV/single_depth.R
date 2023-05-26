# Updated on Oct-04-2021
# Developed by Meiyuan Ji (for questions contact: meiyuan.ji@utah.edu)
# This program estimates the parameters (mean and standard devidation) for the distribution of single-copy read depth via normal mixture models with one component. 

# update on Feb-1-2023
# remove the disperse parameters 

suppressMessages(library("argparse"))
suppressMessages(library("stats4"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Use the median read-depth value as representative observation for each gene region. Observed values was fitted using normal mixture model given one component. ")
    parser$add_argument("statdepth", help = "mean and median of read-depth of individual gene. ")
    parser$add_argument("coverage", help = "estimated single read-depth (integer). ")
    parser$add_argument("-O", "--outdir", help = "output directory. ")
    parser$parse_args()
}

options(warn=-1)
argv <- parse_arg()
rd <- argv$statdepth
cov <- as.numeric(argv$coverage)
outdir <- argv$outdir
# 
rd_df <- read.table(rd, sep = "\t", header = T)
pdf(file.path(outdir, 'histogram.pdf'), width = 6, height = 6)
low_cov <- cov*0.1
high_cov <- cov*20
sub_rd <- rd_df[rd_df$read_depth > low_cov & rd_df$read_depth < high_cov, ]
# 
min_dp <- min(sub_rd$read_depth)
max_dp <- max(sub_rd$read_depth)
rd_hist_freq <- hist(sub_rd$read_depth, breaks = (max_dp - min_dp)%/%1, freq = T, xlab = "read depth", ylab = "frequency", main = NULL) 
# 
hist_count <- rd_hist_freq$counts
hist_mids <- rd_hist_freq$mids
# 
hist_freq <- data.frame(bar_mid = hist_mids, bar_count = hist_count)
total_count <- sum(hist_freq$bar_count)
hist_freq$bar_prop <- hist_freq$bar_count/total_count 
# 
mode_mid <- hist_freq[hist_freq$bar_count == max(hist_freq$bar_count), ]$bar_mid # find mode value
abline(v = mode_mid, col = "red")
write.table(mode_mid, file = file.path(outdir, "single_cov.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
# 
dev.off()
