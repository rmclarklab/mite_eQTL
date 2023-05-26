""" Updated on Oct-04-2021. 
Program developed by Meiyuan Ji (for questions please contact: meiyuan.ji@utah.edu).
This program takes genome fasta file, gene annotation (GTF) and DNA-read alignment (BAM) files, and measuring the read-coverage on gene region with assigned step size. 

The purpose of this program is for providing a big-picture of read-coverage on gene region, which provide cues for estimating single-copy coverage.

Update on Jan-27-2023.
Mapping quality control with parameters setting (MQ, ignore_orphans, flag_filter)
Update on Jan-30-2023.
Combine bam_parse and gene_stat running in order;
go through gene exon region for read-coverage calculation (reduce running time);
Update on Jan-31-2023.
Not write temporary file in the folder. Instead, store data in a dataframe and collect all from all running cores. 
update on Feb-01-2023.
Do not collect the median and mean coverage for each gene, use the raw site coverage on gene CDS region. 
"""

import pandas as pd
import argparse
import os
import pysam
from statistics import mean
from statistics import median
from Bio import SeqIO
from mpi4py import MPI
import statistics
import random 
import time

def endtime(start):
    end = time.time()
    t = end-start
    if t < 60:
        print('{:.2f} seconds elapsed'.format(t))
    elif t < 3600:
        print('{:.2f} minutes elapsed'.format(t/60))
    else:
        print('{:.2f} hours elapsed'.format(t/3600))
    return end

def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, description="Count read-coverage across gene coding region.")
    parser.add_argument("ref", help = "reference FASTA file. ")
    parser.add_argument("gtf", help = "GTF file of the reference genome (indexed and sorted). ")
    parser.add_argument("bam", help = "BAM alignment file of the reference genome (indexed and sorted). ")
    parser.add_argument("-contig", "--contig", help = "assign contig of interest in a table. ", required = False)
    parser.add_argument("-step", "--stepsize", default = 1, type = int, help = "step size for coverage in BAM (default: 1 bp)")
    parser.add_argument("-MQ", "--mapping_quality", default = 0, type = int, help = "minimum mapping quality required (default: 0). ")
    parser.add_argument("-ignore_orphans", "--ignore_orphans", action = "store_true", help = "ignore orphan reads (read not properly paired). ")
    parser.add_argument("-flag_filter", "--flag_filter", default = 256, type = int, help = "ignore reads where any of the bits in the flag are set (default: 256). ")
    parser.add_argument("-O", "--outdir", help = "output directory. ")
    args = parser.parse_args()
    return args

def GTF_CDS_pos(args):
    '''Takes gtf file and write each gene CDS position per row in a new data frame (format: gene|contig|start|end). '''
    gtf_handle = pysam.TabixFile(args.gtf, parser = pysam.asGTF())
    if args.contig == None:
        geneinfo = [(gene.gene_id.split(";")[0].strip("\""), gene.contig, gene.start, gene.end) for gene in gtf_handle.fetch() if gene.feature == "CDS"]
    else:
        tigdf = pd.read_table(args.contig, header = None)
        tigs = tigdf[0].tolist()
        geneinfo = [(gene.gene_id.split(";")[0].strip("\""), gene.contig, gene.start, gene.end) for c in tigs for gene in gtf_handle.fetch(c) if gene.feature == "CDS"]
    #
    gene_pos = pd.DataFrame(geneinfo, columns = ["gene", "contig", "start", "end"])
    gene_pos = gene_pos.sort_values(by = ["gene"])
    return gene_pos

def bam_parse(args, pos_df, seq_dict, bam_handle):
    '''Count read-coverage of given position. '''
    # pos_df is gene CDS position
    step = args.stepsize
    outdir = args.outdir
    mq = args.mapping_quality
    #
    gene_set = sorted(set(pos_df["gene"]))
    g_df_list = list()
    for g in gene_set:
        g_pos = pos_df[pos_df["gene"] == g]
        chr = list(set(g_pos["contig"]))[0]
        pos_cov = {p:len({bam_read.alignment.query_name for read_column in bam_handle.pileup(row.contig, p, p+1, min_mapping_quality = mq, ignore_orphans = args.ignore_orphans, flag_filter = args.flag_filter, truncate = True) for bam_read in read_column.pileups}) for row in g_pos.itertuples() for p in range(row.start, row.end, step) if seq_dict[chr][p] != "N"}
        g_df = pd.DataFrame({"gene":g, "chromosome":chr, "position":list(pos_cov.keys()), "read_depth":list(pos_cov.values())})
        g_df_list.append(g_df)
    #
    df = pd.concat(g_df_list, axis = 0)
    return df

def main():
    '''Split work based on gene number and thread number, run for read-coverage measuring via reading BAM file.'''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0:
        start = time.time()
    else:
        start = None
    start = comm.bcast(start, root = 0)
    # all threads parsing bam file
    args = args_parser()
    bam = args.bam
    outdir = args.outdir
    bam_handle = pysam.AlignmentFile(bam, "rb")
    # all threads parsing reference fasta 
    ref = args.ref
    seq = SeqIO.parse(ref, "fasta")
    seq_dict = SeqIO.to_dict(seq)
    # all threads parsing position dataframe
    pos_df = GTF_CDS_pos(args) # gene CDS region 
    # split work based on gene number
    if rank == 0:
        if outdir not in os.listdir():
            os.mkdir(outdir)
        worker_tasks = {w:[] for w in range(size)} # a dictionary
        w_idx = 0
        gid_set = sorted(set(pos_df["gene"]))
        print("#### On " + os.path.basename(args.bam) + " with ignore_orphans=" + str(args.ignore_orphans))
        print("#### See output " + args.outdir)
        for i in gid_set:
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    # running coverage screen on BAM file
    i = worker_tasks[rank] # i stores gene names
    subpos = pos_df[pos_df.gene.isin(i)]
    #
    print("Process BAM on core ", rank)
    sub_bam_count = bam_parse(args, subpos, seq_dict, bam_handle) # read-coverage of exon positions by parsing BAM file
    bam_count_gather = comm.gather(sub_bam_count, root = 0)
    if rank == 0:
        print("Gather bam count data from all processors ......")
        bam_count = pd.concat(bam_count_gather, axis = 0)
        bam_count.to_csv(os.path.join(outdir, "pos_depth.txt"), sep = "\t", index = False)
        

################
#### Run it ####
################

if __name__ == "__main__":
    main()
