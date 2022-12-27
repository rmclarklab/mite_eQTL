###########################################################
### This script is for library-size normalization 
### using DESeq2. Accept read count data file with 
### column information. 
###########################################################

suppressMessages(library("DESeq2"))
suppressMessages(library("argparse"))
suppressMessages(library("stringr"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Given sample read count files, return normalized read count data in one combined file.")
    parser$add_argument("-count", "--count", help = "raw read count in tab-separated file. ")
    parser$add_argument("-O", "--output", help = "output prefix. ")
    parser$parse_args()
}

###
argv <- parse_arg()
count_file <- argv$count
out <- argv$output

### generate sample file based on read count data in column
count_df <- read.table(count_file, sep = "\t", header = T, row.names = 1)
col_info <- data.frame(samples = colnames(count_df))

### use DESeq2 to normalizing read count 
dds <- DESeqDataSetFromMatrix(countData = count_df, colData = col_info, design = ~ 1) # no need to give experimental design when only for normalization 
dds <- estimateSizeFactors(dds)
norm_cts <- counts(dds, normalized = T)

# to transform the count data to the log2 scale which minimizes differences between samples for rows with sample counts, and normalizes with respect to library size (rlog). 
# rld <- rlog(dds)
# norm_cts <- as.data.frame(assay(rld))

write.table(norm_cts, file = paste0(out, ".txt"), sep = "\t", quote = F, row.names = T)
