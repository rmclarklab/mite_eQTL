# Expression quantitative trait locus (eQTL) mapping in <i>Tetranychus urticae</i> (a generalist spider mite herbivore) and associated analyses

This is a repository for performing and analyzing an eQTL mapping experiment with <i>Tetranychus urticae</i>, the two-spotted spider mite, and follow up studies.

## Authorship

This site, and the respective pipeline documentation and code, has been assembled by Meiyuan Ji, a graduate student in Richard Michael Clark's laboratory at the University of Utah. This site supports a manuscript under review, and will be linked to the resulting publication when it is published. The work is a collaboration with Thomas Van Leeuwen's laboratory at Ghent University. <br> <br>
For questions regarding to the use of scripts and data analysis steps, contact <i>meiyuan.ji@utah.edu</i> or <i>clark@biology.utah.edu</i>. 

## Experimental design
   
We used a pesticide susceptible strain ROS-ITi (diploid mother, ♀, **S**) and more resistant strain MR-VPi (haploid father, ♂, **R**) as inbred parent strains for the generation of isogenic F3 full sibling families from which RNA was collected for RNA-seq for use in eQTL mapping (the resulting set of F3 isogenic families constituted the eQTL mapping population).

## Table of contents

- [DNA-seq for variants calling](#DNA-seq-for-variants-calling)
- [Map RNA-seq against the reference genome](#Map-RNA-seq-against-the-reference-genome)
- [Genotype call for eQTL mapping populations based on RNA-seq alignment](#Genotype-call-for-eQTL-mapping-populations-based-on-RNA-seq-alignment)
- [Update GFF3 file for the reference genome](#Update-GFF3-file-for-the-reference-genome)
- [Gene expression level quantification](#Gene-expression-level-quantification)
- [Association analysis between genotype and gene expression](#Association-analysis-between-genotype-and-gene-expression)
- [Gene copy number variation estimation](#Gene-copy-number-variation-estimation)

## Programs used / Dependencies

- python3+ (packages: pysam v0.15.3; biopython v1.76; pandas v0.25; numpy v1.21; mpi4py v3.0) and upper
- BWA v0.7.17-r1188
- STAR v2.7.3a
- GATK v4.2 and upper
- samtools v1.15 and upper
- htseq-count v2.0.1 and upper
- R v4.1.3 (packages: DESeq2 v1.34; MatrixEQTL v2.3; R/qtl v1.46) and upper

[NOTE]To enable parallel processing, python model mpi4py need to be installed. 

## DNA-seq for variants calling
For variant calling see folder "VCF". <br> 
To call variants for the inbred **R** and **S** strains, we mapped illumina DNA-seq against the three-chromosome reference genome (London strain, see [Wybouw, Kosterlitz, et al., 2019](https://academic.oup.com/genetics/article/211/4/1409/5931522)). <br>
For more about GATK best practices for variants calling, see [here](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). <br>
1. First, prepare the index for the genome fasta file;
```bash
# make directory for bwa index files
mkdir bwa_index
# change working directory to the folder
cd bwa_index
# generate index files using bwa index command
bwa index Tetranychus_urticae_2017.11.21.fasta
```
2. Then, map DNA-seq of **R** and **S** onto the reference fasta genome using BWA;
```bash
# make directory for bwa mapping
mkdir bwa_map
# change working directory to the mapping folder
cd bwa_map
# run bwa mapping for ROS-ITi sample (paired-end DNA sequences)
bwa mem -t 20 -R "@RG\tID:20190412_8\tSM:ROS-ITi\tPL:Illumina\tLB:ROS-ITi" bwa_index/Tetranychus_urticae_2017.11.21.fasta r1.fastq.gz r2.fastq.gz | samtools view -Su - | samtools sort -@ 20 - -o ROS-ITi.BWA.bam
# run bwa mapping for MR-VPi sample (paired-end DNA sequences)
bwa mem -t 20 -R "@RG\tID:20190312\tSM:MR-VPi\tPL:Illumina\tLB:MR-VPi" bwa_index/Tetranychus_urticae_2017.11.21.fasta r1.fastq.gz r2.fastq.gz | samtools view -Su - | samtools sort -@ 20 - -o MR-VPi.BWA.bam
```
3. Mark duplicated reads (PCR duplicates) using picard MarkDuplicate;
```bash
# mark duplicated reads in the BWA generated BAM files
picard MarkDuplicates I=ROS-ITi.BWA.bam O=ROS-IT_duplicate.bam M=ROS-IT_metrics.txt && samtools index ROS-IT_duplicate.bam
picard MarkDuplicates I=MR-VPi.BWA.bam O=MR-VP_duplicate.bam M=MR-VP_metrics.txt && samtools index MR-VP_duplicate.bam
# left align indels (optional)
gatk LeftAlignIndels -R Tetranychus_urticae_2017.11.21.fasta -I ROS-IT_duplicate.bam -O ROS-IT_leftalign.bam
gatk LeftAlignIndels -R Tetranychus_urticae_2017.11.21.fasta -I MR-VP_duplicate.bam -O MR-VP_leftalign.bam
```
4. Run the Best practice of GATK pipeline for variant calling;
```bash
# run gatk HaplotypeCaller to generate g.vcf files
gatk HaplotypeCaller -R <ref> -I ROS-IT_leftalign.bam -ERC GVCF -O ROS-IT.g.vcf.gz
gatk HaplotypeCaller -R <ref> -I MR-VP_leftalign.bam -ERC GVCF -O MR-VP.g.vcf.gz
# Merge GVCFs from the two samples
gatk GenomicsDBImport --genomicsdb-workspace-path <path> -R <ref> -sequence-dictionary <dict> --sample-name-map <samples>
# Perform joint genotype calling
gatk GenotypeGVCFs -R <ref> -V gendb://your_database -O <output.vcf.gz>
# sort the resulting VCF and generate its index file
bcftools sort -o <sorted.vcf.gz> -O z input.vcf
bcftools index -t sorted.vcf.gz
```
5. Select variants from the unfiltered vcf file;
```bash
# select SNPs
gatk SelectVariants -R <ref> -V input.vcf.gz -select-type-to-include SNP -O SNP.vcf.gz
```
6. Filter out SNPs based on the following criteria (MQ > 40, QD > 2, SOR < 3, and different between the two **R** and **S** strains), and output a new vcf file with variants passed the critieria.
```bash
vcf_filter_2sample.py -vcf SNP.vcf.gz -O filtered.snp.vcf.gz
```
For Variants filtering, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) also for hard-filtering. <br>

7. To retrieve variants information from vcf file into a tab-separated file, which will be used for following eQTL analysis. 

```bash
variant_retrieve.py -vcf filtered.snp.vcf.gz -O variant_RS
```
Use the output file named "variant_RS.allele.txt" for further analysis. 

## Map RNA-seq against the reference genome
Align the Illumina RNA read data for the eQTL mapping population to the three-chromosome reference genome (that is, the same genome used for DNA-seq read alignment). <br>
For RNA-seq read alignment we used the RNA-seq aligner [STAR](https://github.com/alexdobin/STAR) 

1. Generate indices for the genome fasta file.

```bash
# STAR index generation
STAR --runMode genomeGenerate --runThreadN 30 --genomeDir STAR_index --genomeFastaFiles Tetranychus_urticae_2017.11.21.fasta --genomeSAindexNbases 12
```

2. Map RNA-seq onto the reference genome using the index folder. 

```bash
# STAR mapping, sort and index BAM alignment file
STAR --genomeDir STAR_index --runThreadN 20 --readFilesIn <r1.fastq.gz> <r2.fastq.gz> --twopassMode Basic --sjdbOverhang 99 --outFileNamePrefix <sample_name>. --readFilesCommand zcat --alignIntronMax 30000 --outSAMtype BAM Unsorted && samtools sort <sample_name>.Aligned.out.bam -o <sample_name>_sorted.bam -@ 8 && samtools index <sample_name>_sorted.bam 
```

## Genotype calling for eQTL mapping populations based on RNA-seq alignment
See folder "genotype". <br>
We developed a customized pipeline to genotype the F3 isogenic sibling families using the RNA-seq read alignments and SNP data obtained from the DNA-seq read alignments and variant calls for the **R** and **S** parents (F0 generation). 
Inputs:
  - BAM RNA-seq data file from each F3 isogenic family (the input must be coordinate sorted and have a bai index file);
  - SNP information (tab-separated) for the two F0 strains (the fixed differences between strains, see final output of [DNA-seq for variants calling](#DNA-seq-for-variants-calling)).

1. Count allele-specific reads at SNP sites for each F3 BAM file

```bash
# run SNP_allele_count.py to count allele-specifc reads on the SNP sites
# this is a multiple-core processing program, adjust core usage via "-n"
mpiexec -n 10 SNP_allele_count.py -V variant_RS.allele.txt -bam <sample.bam> -O <sample>_allele_count
```
To call genotype at SNP site, see script record in ```allele2genotype.py```. <br>
After running for all F3 BAM files, place all of them in the same folder (raw_count). <br>
The output genotype file for each file with name "<sample>_allele_count.geno.txt". <br>

2. Call genotype blocks for each F3 family based on allele-specific read counts at SNP sites.

You need to set up the chromosomes of interested for genotype block assignment, and also provide chromosome length information in a tab-separated file. See example input under the "data" folder (chrlen.txt).

```bash
# run genotype_block.py to call genotype blocks for each F3 isogenic family
genotype_block.py -chrLen chrlen.txt -count <sample>_allele_count.geno.txt -O <sample>_genotype_block
```

3 (optional) Visulization of genotype blocks on chromosome level

```bash
Rscript block_vis.R -geno sample_genotype_block.txt
```

  Example output (in PDF): <br>
<img width="300" alt="Screen Shot 2022-05-01 at 3 15 45 PM" src="https://user-images.githubusercontent.com/63678158/166165078-eaeace45-abfc-48ca-9301-684e0f670db4.png">

## Update GFF3 file for the reference genome
See folder "GFF_update" <br>
  To integrate existing annotated gene information onto the <i>T. urticae</i> three-chromosome reference genome (see above), we added gene models from [Orcae database](https://bioinformatics.psb.ugent.be/gdb/tetranychus/) (version of 01252019) to create an updated GFF3 file for the current study. <br>
  To transfer gene models from the prior <i>T. urticae</i> genome that was not assembled to chromosomes, to the three-chromosome genome, see script ```GFF_record.py``` under GFF_update folder. <br>
  To combine the added gene models to the current GFF3 file and sort it using [gff3sort.pl](https://github.com/billzt/gff3sort). <br>
  Transform from GFF3 to GTF format using script ```gff2gtf.py``` under GFF_update folder. <br>
  To compress and add index for GFF/GTF, see below (the output was converted to the GTF file used for subsequence analyses, see below):

```bash
# sort gff using gff3sort.pl
gff3sort.pl input.gff > output.gff
# compress gff
bgzip output.gff
# add index for compressed gff
tabix -p gff output.gff.gz
```

## Gene expression level quantification
See folder "normalization" <br> 
1. Using the updated GTF, run htseq-count on the RNA-seq alignment BAM files for the eQTL mapping population and output the read counts per gene.

```bash
# htseq-count command line (adjust number of CPU "-n" based on sample BAMs number, per BAM per CPU)
htseq-count -r pos -s reverse -t exon -i gene_id --nonunique none --with-header -n 3 -c sample1-3.tsv sample1.bam sample2.bam sample3.bam $GTF 
```
Notice that the new version htseq-count (v2.0) can process multiple BAM files in parallel.

2. From the raw read-count data obtained by running htseq-count (see above), performed library-size normalization using DESeq2 (function estimateSizeFactors). <br>

```bash 
# merge htseq-count output of all samples into one file, -O for output name
Rscript DESeq2_norm.R -count all_sample.txt -O all_sample_normalizedbyDESeq2
```
3. Perform quantile normalization on the library-size normalized read data to alleviate the effects from outliers and alleviate any systematic inflation (see also [here](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/faq.html)). Briefly, the normalization is on individual gene level, with gene expression across all samples is fit to a normal distribution while preserving the relative rankings.

```bash
# use the code quantile_norm.R 
Rscript quantile_norm.R -norm_count all_sample_normalizedbyDESeq2.txt -O all_sample_quantile.txt
```

## Association analysis between genotype and gene expression for eQTL mapping
See folder "eQTL" <br>
  For each F3 isogenic sibling family, the genotype blocks and gene expression data are now available (see above). Association analysis can then be performed.

1. First, to alleviate computational pressure for association tests that run for each combination of SNP genotype and individual gene, we assigned genotype bins based on the overlap of genotype blocks among F3 isogenic populations (for logic please see also [Ranjan et. al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5074602/)). The output of this step are the genotype bins (genetic markers) to be used directly for eQTL mapping.

```bash
# run block2bin.R to assign genotype bins based on overlapped blocks across all F3 isogenic sibling families
# for this step, you also need to provide chromosome length information (chrlen.txt) and the SNP position file (SNP_loc.txt)
Rscript block2bin.R -genodir sample_genotype_block/ -chrLen chrlen.txt -SNP SNP_loc.txt
```
All samples genotype files should be stored under the "sample_genotype_block" folder. 

2. Perform genotype-expression association analysis using [MatrixeQTL](https://github.com/andreyshabalin/MatrixEQTL). <br>
  For instructions on preparing input files for MatrixeQTL, see the respective tutorial [here](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html).

```bash 
# command for using MatrixeQTL for association test (without taking gene/SNP location information)
# both ANOVA and LINEAR models are used for the association analysis, and the output files are *.anova.txt and *.linear.txt (for our study, we only subsequently used the output of the linear model)
Rscript eQTL_identify.R -genotype <genotype.txt> -expression <expression.txt> -O <output>
```

3. Significant associations can arise from linkage disequilibruim (LD). To alleviate its effect, we take the bin genotype data to reconstruct linkage groups using [R/qtl](https://rqtl.org/download/). <br>

  For linkage group construction, see R file ```marker_association.R``` that has the needed commandlines.

4. Based on the linkage group information, we then parsed significant associations for each bin and its target gene. The most significant bin for a given linkage group was chosen as the location for the given eQTL. 

```bash
# use the script parse_eQTL.py (support multiple-core running, adjust core usage in "-n")
mpiexec -n 5 eQTL_parse.py -eQTL MatrixeQTL_output -assoc marker_association.txt -O parsed_eQTL
```
Note that this process is memory consuming. 

## Gene copy number variation estimation 
See folder "CNV" <br>
1. Count coverage depth on gene coding regions (default stepsize 1 bp). <br>
```bash
python gene_coverage.py [ref] [gtf] [bam] -O [out] 
```
required arguments for this script
ref: reference genome in fasta file <br>
gtf: gtf file for the reference genome <br>
bam: bam of reads aligned to reference genome <br>
out: output folder
- A new "out" folder will be created (if not existing) and all output files will be written under the folder. <br>
- File named "pos_depth.txt" (coverage at gene coding positions) will be generated under the \[out\] folder. <br>

2. Report the single-copy coverage depth for the BAM file.  <br>
```bash
Rscript single_depth.R [out]/pos_depth.txt [cov_est] -O [out] 
```
Arguments from last step:
cov_est: the estimated coverage for the BAM file <br>
out: output folder (it should be the same as in "Step 1") <br>
- File named "single_cov.txt" (single copy coverage) will be generated under the \[out\] folder <br>
- File named "histogram.pdf" (coverage distribution) will be generated under the \[out\] folder <br>

3. Estimate gene CNV based on gene coding region coverage depth <br>
```bash
python gene_CNV.py [out]/pos_depth.txt [out]/single_cov.txt -O [out] 
```
Arguments from last step:
out: output folder (it should be the same as in "Step 1" and "Step 2") <br>
- File named "gene_cnv.txt" will be generated under the \[out\] folder. <br>

