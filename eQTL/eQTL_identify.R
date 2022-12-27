
suppressMessages(library("argparse"))
suppressMessages(library("MatrixEQTL"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "Usage: Rscript eQTL_identify.R -genotype <genotype.txt> -expression <expression.txt> -O <output> ")
    parser$add_argument("-genotype", "--genotype", help = 'genotype data for all samples in a tab-separated file. ')
    parser$add_argument("-expression", "--expression", help = "gene expression data for all samples (recommend: quantile normalization performed on gene expression). ")
    parser$add_argument("-p", "--pvalue", default = 0.01, help = "p cutoff value (default: 0.01). ")
    parser$add_argument("-O", "--output", help = "output file prefix. ")
    parser$parse_args()
}

###
args <- parse_arg()
geno <- args$genotype
expre <- args$expression
pv <- as.numeric(args$pvalue)
out <- args$output
### make sure samples are matched between gneotype and expression file
genotype <- read.table(geno, sep = "\t", header = T, row.names = 1)
expression <- read.table(expre, sep = "\t", header = T, row.names = 1)
geno_col <- colnames(genotype) 
exp_col <- colnames(expression)

if (identical(geno_col,exp_col)) {
    run <- T
    print("genotype and expression files all match!")
} else {
    run <- F
    print("genotype and expression files do NOT match!")
}

###
if (run == T) {
    snp_mat <- as.matrix(genotype)
    snp_new <- SlicedData$new()
    snp_new$fileSliceSize = 2000
    snp_new$CreateFromMatrix(snp_mat)
    # phenotype data (quantile normalized)
    gene_mat <- as.matrix(expression)
    gene_new <- SlicedData$new()
    gene_new$fileSliceSize <- 2000
    gene_new$CreateFromMatrix(gene_mat)
    ### 
    # Set pvOutputThreshold > 0, pvOutputThreshold.cis = 0: eQTL analysis without using gene/SNP locations 
    # Set pvOutputThreshold = 0, pvOutputThreshold.cis > 0: eQTL analysis for local gene-SNP pairs only
    # Set pvOutputThreshold > 0, pvOutputThreshold.cis > 0: eQTL analysis with separate p-value thresholds for local and distant eQTLs.
    tetur_ANOVA <- Matrix_eQTL_main(snps = snp_new, gene = gene_new, pvOutputThreshold = pv, pvOutputThreshold.cis = 0, useModel = modelANOVA, pvalue.hist = F, noFDRsaveMemory = F,  output_file_name = paste0(out, ".anova.txt"))
    tetur_LINEAR <- Matrix_eQTL_main(snps = snp_new, gene = gene_new, pvOutputThreshold = pv, pvOutputThreshold.cis = 0, useModel = modelLINEAR, pvalue.hist = F, noFDRsaveMemory = F, output_file_name = paste0(out, ".linear.txt"))
} 
