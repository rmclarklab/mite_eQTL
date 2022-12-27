#!/usr/bin/env python

import pysam
import pandas as pd
import argparse
import os
from mpi4py import MPI
import time 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Allele-specifi read count at retrieved SNP positions. \n \
        Usage: mpiexec -n <cores> genotype_allele.py -variant <var> -bam <bam> -O <output> ")
    parser.add_argument("-V", "--variant", help = "variant table with genotype alleles with their position in 1-basis. ")
    parser.add_argument("-bam", "--bam", help = "the sorted bam with index file in the same directory. ")
    parser.add_argument("-O", "--output", help = "Output name of the read count at allele positions. ")
    args = parser.parse_args()
    return args

def genotype_call(args, chr, pos, alt1, alt2):
    '''Read count at each allele position. '''
    bam_handle = pysam.AlignmentFile(args.bam, "rb")
    ###
    alt1_read = set()
    alt2_read = set()
    other_read = set()
    for read_column in bam_handle.pileup(chr, pos-1, pos): 
        # pysam is for 0-based position
        if read_column.reference_pos == pos-1:
            for bam_read in read_column.pileups:
                read_name = bam_read.alignment.query_name
                if (("NH", 1) in bam_read.alignment.tags and bam_read.query_position != None):
                    ###
                    allele = bam_read.alignment.query_sequence[bam_read.query_position]
                    if allele == alt1:
                        alt1_read.add(read_name)
                    elif allele == alt2:
                        alt2_read.add(read_name)
                    else:
                        other_read.add(read_name)
    alt1_count = len(alt1_read)
    alt2_count = len(alt2_read)
    other = len(other_read)
    return alt1_count, alt2_count, other

def check_rm(args, size):
    output = args.output
    if output + ".txt" in os.listdir():
        for s in range(size):
            os.remove(output + "_temp" + str(s) + ".txt")

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

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    ###
    args = args_parser()
    var = args.variant
    bam = args.bam
    output = args.output
    ### split work
    if rank == 0:
        start = time.time()
        var_df = pd.read_table(var, sep = "\t", header = 0)
        var_num = len(var_df)
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for i in range(var_num):
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
        print("split work!")
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0) # broadcast work division
    ###
    outdf = open(output + "_temp" + str(rank) + ".txt", "w")
    outdf.write("chromosome\tposition\talt1\talt2\talt1_count\talt2_count\tothers\n")
    var_df = pd.read_table(var, sep = "\t", header = 0)
    rank_df = var_df.loc[worker_tasks[rank]]
    for row in rank_df.itertuples():
        chr = row.chromosome
        pos = row.position
        alt1 = row.alt1
        alt2 = row.alt2
        alt1_count, alt2_count, other = genotype_call(args, chr, pos, alt1, alt2)
        ### for each variant site, write read count
        outdf.write(f"{chr}\t{pos}\t{alt1}\t{alt2}\t{alt1_count}\t{alt2_count}\t{other}\n")
    outdf.close()
    print("Run with {} cores, finished on processor {}".format(size, rank))
    df = pd.read_table(output + "_temp" + str(rank) + ".txt", header = 0, sep = "\t")
    ### gather data
    df = comm.gather(df, root = 0) # return list object
    if rank == 0:
        df = pd.concat(df, axis = 0)
        df = df.sort_values(by = ["chromosome", "position"])
        df.to_csv(output + ".txt", sep = "\t", index = False)
        print("Write genotyping data successfully! ")
        endtime(start)
        check_rm(args, size)

##################
##### Run it #####
##################

if __name__ == "__main__":
    main()

