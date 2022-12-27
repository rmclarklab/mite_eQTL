#!/usr/bin/env python

import pandas as pd
import os
import argparse
import time
import shutil


def args_parser():
    '''parser the arguments from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=". \n\
        Usage: ")
    parser.add_argument("-genotype", "--genotype", help = "genotype block data for all samples. ")
    parser.add_argument("-bin", '--bin', help = 'the bin file. ')
    parser.add_argument("-O", "--output", help = "the output file. ")
    args = parser.parse_args()
    return args

def collect_genotype(genotype):
    '''Collect genotype block data and build nested dictionary for all samples.'''
    geno_df = pd.read_table(genotype, sep = "\t", header = 0)
    sample_set = sorted(set(geno_df["sample"]))
    geno_dict = dict()
    for s in sample_set:
        geno_dict[s] = dict()
        sub_geno = geno_df[geno_df["sample"] == s]
        chr_set = sorted(set(sub_geno["chromosome"]))
        for chr in chr_set:
            geno_dict[s][chr] = dict()
            chr_sub = sub_geno[sub_geno['chromosome'] == chr]
            for row in chr_sub.itertuples():
                rstart = row.start
                rend = row.end
                rgeno = str(row.genotype)
                geno_dict[s][chr][range(rstart, rend+1)] = rgeno
    return geno_dict

def snp_genotype(bin, geno_dict):
    bin_df = pd.read_table(bin, sep = '\t', header = 0)
    bin_df["position"] = (bin_df["bin_start"] + bin_df["bin_end"])//2
    for g in geno_dict:
        bin_df[g] = "-"
        sample_chrs = [sample_chr for sample_chr in geno_dict[g]]
        for sample_chr in sample_chrs:
            bin_sub = bin_df[bin_df["chromosome"] == sample_chr]
            for row in bin_sub.itertuples():
                pos = row.position
                grange = [r for r in geno_dict[g][sample_chr] if pos in r]
                sample_genotype = geno_dict[g][sample_chr][grange[0]]
                bin_df.loc[(bin_df["chromosome"] == sample_chr) & (bin_df["position"] == pos), g] = sample_genotype
    return bin_df

def main():
    args = args_parser()
    genotype = args.genotype
    bin = args.bin
    geno_dict = collect_genotype(genotype)
    snp_genotype(bin, geno_dict)

##################
##### Run it #####
##################

if __name__ == "__main__":
    main()