#! /usr/bin/env python
import pysam
import argparse
from Bio import SeqIO
import subprocess
import os 
import pandas as pd
import random 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Filter variants and rewrite the passed the variants into a new vcf file. \n\
    Usage: vcf_pass.py -vcf <snp.vcf> -R <reference.fa> -O <vcf_output>")
    parser.add_argument("-vcf", "--vcf", help="filtered clean vcf file to subset. ")
    parser.add_argument("-tab", "--table", required = False, help = "table with id of variant set for information retrieve")
    parser.add_argument("-O", "--output", help="output name for the new vcf file with filtered variants. ")
    parser.add_argument("-window", "--window", required = False, type =  int, help = "window size to pick variant data in approximate linkage equilibruim. ")
    parser.add_argument("-if_SNPseq", "--if_SNPseq", action = "store_true", help = "if write SNP sequence. ")
    args=parser.parse_args()
    return args

def get_snpid(args):
    tab = args.table
    vcf = args.vcf
    output = args.output
    ### 
    tab_df = pd.read_table(tab, sep = "\t", header = 0)
    pos_id = list(tab_df.iloc[:, 0]) # 1-based (from R)
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    header = hd.strip().split("\n")[-1].strip("#")
    ### write to new file
    fw = open(output + ".txt", "w")
    fw.write(f"{header}\n")
    ###
    vcf_handle = pysam.TabixFile(vcf, parser = pysam.asVCF())
    i = 0
    vcf_dict = dict()
    for v in vcf_handle.fetch():
        i += 1
        vcf_dict[i] = v
    ###
    for p in pos_id:
        v = vcf_dict[p]
        fw.write(str(v) + "\n")
    fw.close()

def vcf_sample(vcf):
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    header = hd.strip().split("\n")[-1]
    samples_list = header.split("FORMAT\t")[-1].split("\t")
    return samples_list

def vcf_chrom_len(vcf):
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    tig_len = {i.split("ID=")[-1].split(",")[0]:int(i.split("length=")[-1].split(",")[0]) for i in hd.split("\n") if i.startswith("##contig=")}
    return tig_len

def get_vcf_real_chrom(vcf_handle):
    tigs = sorted({v.contig for v in vcf_handle.fetch()})
    return tigs

def random_pick_vcf(args, tig_len):
    vcf = args.vcf
    win = args.window
    output = args.output
    vcf_handle = pysam.TabixFile(vcf, parser = pysam.asVCF())
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    header = hd.strip().split("\n")[-1].strip("#")
    vcf_win = open(output + ".txt", "w")
    vcf_win.write(f"{header}\n")
    tig_space = get_vcf_real_chrom(vcf_handle)
    for t in tig_len:
        if tig_len[t] > win and t in tig_space:
            for w in range(0, tig_len[t], win):
                vw_dict = dict()
                i = 0
                for v in vcf_handle.fetch(t, w, w+win):
                    vw_dict[i] = v
                    i += 1
                if i > 0:
                    v_pick = random.sample(range(i), 1)
                    v_pick_sub = vw_dict[v_pick[0]]
                    vcf_win.write(str(v_pick_sub) + "\n")
    vcf_win.close()

def all_vcf(args):
    vcf = args.vcf
    output = args.output
    vcf_handle = pysam.TabixFile(vcf, parser = pysam.asVCF())
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    header = hd.strip().split("\n")[-1].strip("#")
    vcf_variant = open(output + ".txt", "w")
    vcf_variant.write(f"{header}\n")
    for v in vcf_handle.fetch():
        vcf_variant.write(str(v) + "\n")
    vcf_variant.close()

def parse_geno(output, samples_list):
    genotype_tab = pd.read_table(output + ".txt", sep = "\t", header = 0)
    for s in samples_list:
        genotype_tab[s] = genotype_tab[s].str.split(":", expand = True)[0]
    # genotype_tab.to_csv(output + ".GT.txt", sep = "\t", index = False)
    return genotype_tab

def parse_allele(parsed_genotype, samples_list, output):
    # genotype_df = pd.read_table(output + ".GT.txt", sep = "\t", header = 0)
    genotype_df = parsed_genotype
    col_list = list(genotype_df.columns.values)
    col_index = [col_list.index(s)+1 for s in samples_list]
    fw = open(output + ".allele.txt", "w")
    sample_join = "\t".join(samples_list)
    fw.write(f"CHROM\tPOS\t{sample_join}\n")
    for row in genotype_df.itertuples():
        rchr = row.CHROM
        rpos = row.POS
        rref = row.REF
        ralt = row.ALT.split(",")
        fw.write(f"{rchr}\t{rpos}\t")
        sample_allele = list()
        ###
        for i in col_index:
            r0 = row[i][0]
            r1 = row[i][-1]
            if r0 == r1 and r0 != ".":
                i_index = int(r0)
                if i_index == 0:
                    i_allele = rref
                else:
                    i_allele = ralt[i_index - 1]
            else:
                i_allele = "N"
            sample_allele.append(i_allele)
        ###
        sample_allele = "\t".join(sample_allele)
        fw.write(f"{sample_allele}\n")
    fw.close()

def snp_pick(output, samples_list):
    allele_df = pd.read_table(output + ".allele.txt", sep = "\t", header = 0, index_col = ["CHROM", "POS"])
    sample_num = len(samples_list)
    df_out_new = list()
    for row in allele_df.itertuples():
        for s in range(1, sample_num+1):
            row_s = row[s]
            if len(row_s) > 1:
                df_out_new.append(row.Index)
                break
    allele_sub = allele_df[~allele_df.index.isin(df_out_new)]
    allele_sub.to_csv(output + ".allele.new.txt", sep = "\t", index = True)

def sequence_get(output, samples_list):
    allele_df = pd.read_table(output + ".allele.new.txt", sep = "\t", header = 0, index_col = ["CHROM", "POS"])
    fw = open(output + ".fasta", "w")
    for s in samples_list:
        s_seq = "".join(list(allele_df[s]))
        s_seq = s_seq.replace("*", "-")
        fw.write(f">{s}\n{s_seq}\n")
    fw.close()

def sequence_align(output):
    muscle_command = f"muscle -align {output}.fasta -output {output}.afa"
    subprocess(muscle_command, shell = True)

def fasttree(output):
    command = f"FastTreeMP -gtr -gamma < {output}.afa > {output}.tree"
    subprocess.call(command, shell = True)

def main():
    args = args_parser()
    output = args.output
    vcf = args.vcf
    tig_len = vcf_chrom_len(vcf)
    samples_list = vcf_sample(vcf)
    if args.table != None:
        print("### Write variants by given table ......")
        get_snpid(args)
    if args.window != None:
        print("### Write variants across window ......")
        random_pick_vcf(args, tig_len)
    if args.table == None and args.window == None:
        print("### Write all variants in VCF ......")
        all_vcf(args)
    parsed_genotype = parse_geno(output, samples_list)
    #
    parse_allele(parsed_genotype, samples_list, output)
    if args.window != None:
        snp_pick(output, samples_list)
    if args.if_SNPseq == True:
        sequence_get(output, samples_list)
    # sequence_align(output)
    # fasttree(output)

############################
######### Run it ###########
############################

if __name__ == "__main__":
    main()


# snpset_LD0.4