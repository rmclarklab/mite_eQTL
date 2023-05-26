
import pandas as pd 
import numpy as np 
import argparse
from collections import Counter
import os 
from mpi4py import MPI

def args_parser():
    parser = parser=argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, description="Based on read-coverage on gene exon region, to get gene CNV estimation.")
    parser.add_argument("RD", help = "Read-depth of gene exon region. ")
    parser.add_argument("single", help = "single copy coverage. ")
    parser.add_argument("-min_depth", "--min_depth", default = 0.25, help = "minimum depth to be counted at least one copy (default: 0.25). ")
    parser.add_argument("-representative_prop", "--representative_prop", default = 0.80, help = "minimum proportion to be indicated as representative value for gene (default: 0.80).")
    parser.add_argument("-noise_prop", "--noise_prop", default = 0.05, help = "site CNV counts as noisy when proportion less than the value (default: 0.05). ")
    parser.add_argument("-O", "--outdir", help = "the output directory")
    parser.parse_args()
    args = parser.parse_args()
    return args

def CNV_pos(rd_df, single_value, min_dp):
    '''normalized read depth based on single-copy coverage value. '''
    rd_df["depth_norm"] = rd_df["read_depth"]/single_value
    # CNV of 0 is unique in that the round cutoff is 0.3
    rd_df["CNV"] = np.where(rd_df["depth_norm"] < min_dp, 0, np.where((rd_df["depth_norm"] >= min_dp) & (rd_df["depth_norm"] < 1), 1, round(rd_df["depth_norm"])))
    rd_df["CNV"] = rd_df.CNV.astype("int32")
    return rd_df

def CNV_gene(rd_df, representative_prop, noise_prop):
    gid = sorted(set(rd_df["gene"]))
    cnv_type_list = list()
    cnv_freq_list = list()
    cnv_final_list = list()
    represent_gid = list()
    for g in gid:
        g_df = rd_df[rd_df["gene"] == g]
        gid_cnv = Counter(g_df["CNV"])
        site_g = sum(gid_cnv.values())
        gid_cnv_correct = {i[0]:i[1] for i in gid_cnv.items() if i[1]/site_g > noise_prop}
        site_correct = sum(gid_cnv_correct.values())
        # use the unfiltered CNV values
        cnv_type_list.append(",".join([str(i[0]) for i in gid_cnv.items()]))
        cnv_freq_list.append(",".join([str(i[1]) for i in gid_cnv.items()]))
        # use the filtered CNV values
        cnv_value = round(sum([i[0]*i[1] for i in gid_cnv_correct.items()])/site_correct, 2)
        for i in gid_cnv_correct.items():
            if i[1]/site_correct > representative_prop:
                cnv_value = int(i[0])
                represent_gid.append(g)
                break
        cnv_final_list.append(cnv_value)
    gid_df = pd.DataFrame({"gene":gid, "CNV list":cnv_type_list, "Freq list":cnv_freq_list, "CNV":cnv_final_list})
    gid_df["label"] = np.where(gid_df["gene"].isin(represent_gid), "rep", "non")
    return gid_df

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # all threads parsing bam file
    args = args_parser()
    minp = args.representative_prop
    noisep = args.noise_prop
    min_dp = args.min_depth
    rd_df = pd.read_table(args.RD, sep = "\t", header = 0)
    single_df = pd.read_table(args.single, sep = "\t", header = None)
    single_cov = single_df[0].item()
    if rank == 0:
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        gid_set = sorted(set(rd_df["gene"]))
        for i in gid_set:
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    i = worker_tasks[rank]
    rand_df = rd_df[rd_df["gene"].isin(i)].copy()
    cnv_pos = CNV_pos(rand_df, single_cov, min_dp)
    cnv_gene = CNV_gene(cnv_pos, minp, noisep)
    #
    cnv_gene_all = comm.gather(cnv_gene, root = 0)
    if rank == 0:
        cnv_out = pd.concat(cnv_gene_all, axis = 0)
        cnv_out.to_csv(os.path.join(args.outdir, "gene_cnv.txt"), sep = "\t", index = False)

################
#### Run it ####
################

if __name__ == "__main__":
    main()
