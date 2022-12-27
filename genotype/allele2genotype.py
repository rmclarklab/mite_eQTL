import os 
import pandas as pd 
from mpi4py import MPI
import numpy as np 

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

new_snp = "/home/clarklab/Desktop/Tetranychus_urticae/eQTL/count/allele_specific/SNP_STAR/gatk4.2_new/"

new_fs = sorted([f for f in os.listdir(new_snp) if f.endswith(".txt")])

if rank == 0:
    dir_num = len(new_fs)
    worker_tasks = {w:[] for w in range(size)}
    w_idx = 0
    for i in range(dir_num):
        worker_tasks[w_idx].append(new_fs[i])
        w_idx = (w_idx + 1) % size
    print("split work successfully!")
else:
    worker_tasks = None
worker_tasks = comm.bcast(worker_tasks, root = 0)

log_rank = open(str(rank)+"_rank.txt", "w")
log_rank.write("file\ttotal\tinformative\tproportion\n")

for f in worker_tasks[rank]:
    pre = f.split(".txt")[0]
    new_df = pd.read_table(os.path.join(new_snp, f), sep = "\t", header = 0)
    total_var = len(new_df)
    new_df["total_count"] = new_df['alt1_count'] + new_df['alt2_count'] + new_df['others']
    new_df["total_allele"] = new_df['alt1_count'] + new_df['alt2_count']
    ### informative SNP
    count_df = new_df[((new_df["others"] <= 5) | (new_df["others"]/new_df["total_count"]*100 < 1)) & ((new_df["alt1_count"] >= 5) | (new_df["alt2_count"] >= 5))].copy()
    informative_var = len(count_df)
    prop = round(informative_var/total_var*100, 2)
    log_rank.write(f"{f}\t{total_var}\t{informative_var}\t{prop}\n")
    ### assign genotype
    count_df["genotype"] = np.where((count_df['alt2_count'] < 8) & (count_df['alt2_count']/count_df['total_allele']*100 < 5), 0, 1)
    count_df.to_csv(pre + ".geno.txt", sep = "\t", index = False)

log_rank.close()

comm.barrier()

if rank == 0:
    fs = sorted([f for f in os.listdir() if f.endswith("_rank.txt")])
    df_list = list()
    for f in fs:
        df_log = pd.read_table(f, sep = "\t", header = 0)
        df_list.append(df_log)
    df_all = pd.concat(df_list, axis = 0)
    df_all.to_csv("variants.log", sep = "\t", index = False)

