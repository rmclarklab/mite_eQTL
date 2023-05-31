
suppressMessages(library("argparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("plyr"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "Usage: Rscript block2bin.R -genodir <block_directory> -chrLen <chrlen.txt> -SNP <SNP_pos.txt> ")
    parser$add_argument("-genotype", "--genotype", help = 'all samples genotype information. ')
    parser$add_argument("-chrLen", "--chrLen", help = "chromosome length file. ")
    parser$add_argument("-out", "--output", help = "output name")
    parser$parse_args()
}

args <- parse_arg()
genotype <- args$genotype
chrlen <- args$chrLen
out <- args$output

# all sample genotype block information (only for chromosome)
df <- read.table(genotype, sep = "\t", header = T)
df <- df[grepl("chr", df$chromosome), ] # only chromosome
df <- df[df$left != "-", ] # only breakpoint position

### chromosome length file
len_df <- read.table(chrlen, sep = "\t", header = T)

### assign bin and get the representative SNP
bin_start_vec <- c()
bin_end_vec <- c()
chr_vec <- c()

chrs <- unique(df$chromosome)
for (i in chrs) {
    print(i)
    bin_start <- 1
    chr_len <- len_df[len_df$chromosome == i, 'length']
    sub_df <- df[df$chromosome == i, ]
    # sub_snp <- snp_df[snp_df$chromosome == i, ]
    sub_df <- sub_df[order(sub_df$end), ]
    min_mid <- min(sub_df$end)
    max_mid <- max(sub_df$end)
    rownames(sub_df) <- seq(1, nrow(sub_df))
    sub_df$left <- as.numeric(sub_df$left)
    sub_df$right <- as.numeric(sub_df$right)
    sub_df$end <- as.numeric(sub_df$end)
    ###
    l0 <- 0
    r0 <- 0
    for (n in seq(1, nrow(sub_df))) {
        l <- sub_df[n, "left"]
        r <- sub_df[n, "right"]
        e <- sub_df[n, "end"]
        gap <- r-l
        overlaplr <- intersect(seq(l0, r0), seq(l, r))
        if (length(overlaplr) > 1) {
            print(paste0("overlap with ", i, " on ", l0, "-", r0))
        } else {
            bin_start_vec <- c(bin_start_vec, bin_start)
            bin_end_vec <- c(bin_end_vec, e)
            chr_vec <- c(chr_vec, i)
            l0 <- l
            r0 <- r
            bin_start <- e+1
        }
    }
    bin_start_vec <- c(bin_start_vec, bin_start)
    bin_end_vec <- c(bin_end_vec, chr_len)
    chr_vec <- c(chr_vec, i)
}


bin_df <- data.frame(chromosome = chr_vec, bin_start = bin_start_vec, bin_end = bin_end_vec)
write.table(bin_df, file = paste0(out, ".txt"), sep = "\t", quote = F, row.names = F)
