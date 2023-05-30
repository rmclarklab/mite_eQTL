
suppressMessages(library("reshape2"))
suppressMessages(library("qtl"))
suppressMessages(library("tidyr"))

# take genotype_bin.txt with marker on row and sample on column (0 for homozygous, 1 for heterzygous)
genotype <- read.table("data/eQTL_genotype_all.txt", sep = "\t", header = T, row.names = 1)
genotype <- ifelse(genotype == 0, "AA", "AB")
genotype <- as.data.frame(genotype)
samples <- colnames(genotype)

# re-format genotype_bin file as a csv format
genotype.t <- t(as.matrix(genotype))
genotype.t <- as.data.frame(genotype.t)
write.csv(genotype.t, file = "data/genotype_bin.csv", quote = F, row.names = T)

# read genotype file using read.cross
mapthis <- read.cross("csv", file = "genotype_bin.csv", genotypes = c("AA", "AB"), alleles = c("A", "B"), estimate.map = F, crosstype = 'bc')
# summary(mapthis)
# estimate recombination fraction between all pairs of genetic markers
rf <- est.rf(mapthis)

# pull out recombination fraction as calculated by est.rf
rf.rf <- pull.rf(rf, what = "rf")
# pull out LOD scores as calculated by est.rf
rf.lod <- pull.rf(rf, what = "lod")
# plot(rf.rf[1,], rf.lod[1,], xlab = "rec frac", ylab = "LOD score")

prf <- as.data.frame(rf.rf)
plod <- as.data.frame(rf.lod)

prf['marker1'] <- rownames(prf)
rf.melt <- melt(prf, id = c("marker1"), na.rm = T)
colnames(rf.melt) <- c("marker1", "marker2", "rf")

plod['marker1'] <- rownames(plod)
lod.melt <- melt(plod, id = c("marker1"), na.rm = T)
colnames(lod.melt) <- c("marker1", "marker2", "LOD")

marker_asso <- merge(rf.melt, lod.melt, by = c("marker1", "marker2"))
write.table(marker_asso, file = "marker_asso.txt", sep = "\t", quote = F, row.names = F)

lg <- formLinkageGroups(rf, max.rf = 0.4, min.lod = 3)
# 2 markers will be placed in the same linkage groups if they have estimated recombination fractioin < max.rf and LOD score > min.lod
write.table(lg, file = "linkage_group.txt", sep = "\t", quote = F)

lng <- read.table("linkage_group.txt", sep = "\t", header = T, row.names = 1)
lng2 <- lng[lng$LG == 2, ]
lng1 <- lng[lng$LG == 1, ]
lng3 <- lng[lng$LG == 3, ]
write.table(lng1, file = "lng1.txt", sep = "\t", row.names = T, quote = F)
write.table(lng2, file = "lng2.txt", sep = "\t", row.names = T, quote = F)
write.table(lng3, file = "lng3.txt", sep = "\t", row.names = T, quote = F)

