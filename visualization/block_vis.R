suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("ggbio"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Usage: Rscript ")
    parser$add_argument("-geno", "--geno", help = 'Input genotype block for visualization. ')
    #parser$add_argument("-chromosome", "--chromosome", nargs = "+", help = "assign the chromosomes for ploting crossover events. ")
    parser$add_argument("-O", "--output", default = "default", help = "the output name")
    parser$parse_args()
}

# arguments
args <- parse_arg()
geno <- args$geno
#chr <- args$chromosome
out <- args$output

### assign output file name
if (out == "default") {
    gen <- basename(geno)
    pre <- strsplit(gen, ".txt")[[1]]
} else {
    pre <- out
}

### plot genotype data
df_genotype <- read.table(geno, sep = '\t', header = T)
win.geno <- df_genotype[df_genotype$chromosome %in% c('chromosome_1', 'chromosome_2', 'chromosome_3'), c("chromosome", "start", "end", "genotype")]
### transform to GRange object
win.range <- makeGRangesFromDataFrame(win.geno, keep.extra.columns=TRUE)
# autoplot function: a generic function to visualize various data object
pdf(file = paste0(pre, ".pdf"))
win.range$genotype <- factor(win.range$genotype, levels = c('homo', "hete"))
autoplot(win.range, layout = "karyogram", aes(color = genotype, fill = genotype)) + scale_color_manual(values = c("#48D1CC", "#ADFF2F")) + scale_fill_manual(values = c("#48D1CC", "#ADFF2F"))+ theme_bw() + theme(axis.text.y = element_blank(), legend.position = "none") + ylab(NULL)
dev.off()
# #48D1CC == similar to blue 'homozygous'
# #ADFF2F == similar to yellow 'heterzygous'
