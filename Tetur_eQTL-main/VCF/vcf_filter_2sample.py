#! /usr/bin/env python

"""
This script is for read and filter variants in vcf_snp file,
and variants passed our filter will be written into a new vcf file.
NOTE: The filter condition as GT == "1/1" or "1|1" and MQ > 40.0.
By Meiyuan Ji, meiyuan.ji@utah.edu
Last update 11/19/2019
"""

import pysam
import argparse
from Bio import SeqIO
import subprocess
import os 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Filter variants and rewrite the passed the variants into a new vcf file. \n\
    Usage: vcf_pass.py -vcf <snp.vcf> -R <reference.fa> -O <vcf_output>")
    parser.add_argument("-vcf", "--vcf", help = "unfiltered vcf file. ")
    parser.add_argument("-map_quality", "--mapping_quality", type = float, default = 40.0, help = "apping quality for variant calling. ")
    parser.add_argument("-qualdepth", "--qualdepth", default = 2, type = float, help = "quality of depth")
    parser.add_argument("-strandodds", "--strandodds", default = 3, type = float, help = "strand odds parameter. ")
    # parser.add_argument("-genotype_cutoff", "--genotype_cutoff", type = float, default = 0.8, help = "cutoff for the percentage of available genotype across all population (default:0.8).")
    parser.add_argument("-O", "--output", help="output name for the new vcf file with filtered variants. ")
    args=parser.parse_args()
    return args

def new_vcf(args):
    # pysam.VariantFile read vcf file and write header into new file
    vcf = args.vcf
    mq = args.mapping_quality
    sor = args.strandodds
    qd = args.qualdepth
    # cutoff = args.genotype_cutoff
    ## write header line for new vcf
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    output = args.output
    vcf_w = open(output + ".vcf", "w")
    vcf_w.write(hd)
    print("###### Write vcf header ######")
    sample_num = len(hd.strip().split("\n")[-1].split("FORMAT")[-1].strip().split("\t"))
    # cutoff_sample = sample_num*cutoff # reverse to missing data cutoff
    ###
    vcf_handle = pysam.TabixFile(vcf, parser = pysam.asVCF())
    for v in vcf_handle.fetch():
        v_boolean = False
        if ("MQ" in v.info) and ("SOR" in v.info) and ("QD" in v.info):
            MQ = float(v.info.split("MQ=")[-1].split(";")[0])
            SOR = float(v.info.split("SOR=")[-1].split(";")[0])
            QD = float(v.info.split("QD=")[-1].split(";")[0])
            if MQ >= mq and SOR <= sor and QD >= qd:
                sample1_GT = v[0].split(":")[0]
                sample2_GT = v[1].split(":")[0]
                if ("." not in sample1_GT) and ("." not in sample2_GT) and (sample1_GT[0] == sample1_GT[-1]) and (sample2_GT[0] == sample2_GT[-1]) and (sample1_GT[0] != sample2_GT[0]):
                    v_boolean = True
        ###
        if v_boolean == True:
            v.filter = "PASS"
            vcf_w.write(str(v) + "\n")
    vcf_r.close()
    vcf_w.close()

def sort_index_vcf(args):
    output = args.output
    sort_command = f"bcftools sort -o {output}.vcf.gz -O z {output}.vcf"
    index_command = f"bcftools index -t {output}.vcf.gz"
    subprocess.call(sort_command, shell = True)
    subprocess.call(index_command, shell = True)

def main():
    args = args_parser()
    new_vcf(args)
    sort_index_vcf(args)

############################
######### Run it ###########
############################

if __name__ == "__main__":
    main()
