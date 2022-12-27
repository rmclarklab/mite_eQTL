#!/usr/bin/env python

"""

This script was developed as a part of eQTL identification pipeline. 
For significant eQTLs associated with one single target gene, they are often closely located, which is arised from linkage disequilibruim (LD). To maintain eQTLs in approximate linkage equilibruim and remove noise eQTLs in LD, I developed this script. 
Briefly, for all eQTLs associated with one single gene, if they are in a linkage group, only the one with the most significant p-value was kept, all others are thought from LD effects. 

For any questions, please directly email: 
Meiyuan Ji (meiyuan.ji@utah.edu)
last update: 10/01/2021

"""


import pandas as pd
import argparse
from mpi4py import MPI
import time 
import os

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=" \n \
        Example: mpiexec -n 8 eQTL_parse.py -eQTL <eQTL> -assoc <mark_assoc> ")
    parser.add_argument("-eQTL", "--eQTL", help = 'eQTL result file from MatrixeQTL. ')
    parser.add_argument("-assoc", "--association", help = "association file with recombination fraction and LOD values. ")
    parser.add_argument("-FDR", "--FDR", default = 0.01, type = float, help = "FDR cutoff (default: 0.01). ")
    parser.add_argument("-O", "--output", required = None, help = "Prefix for the output file. ")
    args = parser.parse_args()
    return args

def min_sig(trans_gene):
    snp_sig = sorted(list(trans_gene[trans_gene["FDR"] == min(trans_gene['FDR'])]['SNP']))
    return snp_sig

def check_trans(snp_sig, asso, trans_gene):
    if len(snp_sig) > 1:
        for s1 in snp_sig:
            for s2 in snp_sig:
                if s1 < s2:
                    rf = float(asso[(asso['marker1'] == s1) & (asso["marker2"] == s2)]['rf'])
                    lod = float(asso[(asso['marker1'] == s1) & (asso["marker2"] == s2)]['LOD'])
                    if rf < 0.4 or lod > 3:
                        print("link: " + s1 + "-" + s2)
                    else:
                        print("not link" + s1 + "-" + s2)
                    trans_gene = remove_link(s1, asso, trans_gene)
                    trans_gene = remove_link(s2, asso, trans_gene)
    else:
        s = snp_sig[0]
        trans_gene = remove_link(s, asso, trans_gene)
    return trans_gene

def remove_link(s, asso, trans_gene):
    snplink = list(asso[(asso['marker1'] == s) & ((asso['rf'] < 0.4) | (asso['LOD'] > 3))]["marker2"])
    snplink.append(s)
    trans_gene = trans_gene[~trans_gene['SNP'].isin(snplink)]
    return trans_gene

def iter_check(trans_gene, asso):
    trans_snp = list()
    while True:
        if len(trans_gene) == 0:
            break
        else:
            snp_sig = min_sig(trans_gene)
            trans_snp.append(trans_gene[trans_gene['SNP'].isin(snp_sig)])
            # delete any SNPs which are significant because LD
            trans_gene = check_trans(snp_sig, asso, trans_gene) 
    return trans_snp

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

def check_rm(output, size):
    if output + ".txt" in os.listdir():
        for s in range(size):
            os.remove(str(s) + "_rank.txt")

def main():
    ### parallel
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    ### 
    args = args_parser()
    eQTL = args.eQTL
    assoc = args.association
    out = args.output
    if out == None:
        out = os.path.basename(eQTL)
        out = ".".join(out.split(".")[:-1]) + ".parsed"
    p = args.FDR
    eQTL = pd.read_table(eQTL, sep = "\t", header = 0)
    asso = pd.read_table(assoc, sep = "\t", header = 0)
    eQTL_filt = eQTL[(eQTL["FDR"] < p) & (eQTL['SNP'].str.contains("chr"))]
    # eQTL_other = eQTL[(eQTL["FDR"] < p) & ~(eQTL["SNP"].str.contains("chr"))]
    # print(eQTL_other)
    unigene = list(sorted(set(eQTL_filt['gene'])))
    gene_num = len(unigene)
    if rank == 0:
        print("Total " + str(gene_num) + " genes.")
    ### split work
    if rank == 0:
        start = time.time()
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for i in unigene:
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
        print("split work ......")
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    ### parallel run
    g_trans_unlink = list()
    for g in worker_tasks[rank]:
        eQTL_gene = eQTL_filt[eQTL_filt["gene"] == g]
        trans_snp = iter_check(eQTL_gene, asso)
        g_clean = pd.concat(trans_snp, axis = 0)
        g_trans_unlink.append(g_clean)
    g_all = pd.concat(g_trans_unlink, axis = 0)
    g_all.to_csv(str(rank) + "_rank.txt", sep = "\t", index = False)
    ### collect all dataframe
    df = comm.gather(g_all, root = 0) # return list object
    if rank == 0:
        df = pd.concat(df, axis = 0)
        df = df.sort_values(by = ["gene", "SNP"])
        df.to_csv(out + ".txt", sep = "\t", index = False)
        print("Write trans-eQTL for genes successfully! ")
        endtime(start)
        check_rm(out, size)

##################
##### Run it #####
##################

if __name__ == "__main__":
    main()
