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
    parser.add_argument("-vcf", "--vcf", help="vcf file with raw snp")
    parser.add_argument("-R", "--reference", help="reference genome in fasta file")
    parser.add_argument("-O", "--output", help="output name for the new vcf file with filtered variants. ")
    args=parser.parse_args()
    return args

def chrom_set(args):
    '''chromosome/contig with snp'''
    vcf = args.vcf
    pvf=pysam.VariantFile(vcf)
    uniq_chr=set()
    for vf_fetch in pvf.fetch():
        uniq_chr.add(vf_fetch.chrom)
    print("###### Build seq id container with variants successfully! ######")
    return sorted(uniq_chr)

def fasta_dict(args):
    '''parser reference genome and build dictionary'''
    ref_fasta = args.reference
    seq_parse=SeqIO.parse(ref_fasta, "fasta")
    seq_dict=SeqIO.to_dict(seq_parse)
    seq_dict_mut={}
    for entry in seq_dict:
        seq_dict_mut[seq_dict[entry].id] = str(seq_dict[entry].seq)
    print("###### Building reference genome dictionary successfully! ######")
    return seq_dict_mut

def new_vcf(args):
    # pysam.VariantFile read vcf file and write header into new file
    vcf = args.vcf
    vcf_r = pysam.VariantFile(vcf, "r")
    hd = str(vcf_r.header)
    output = args.output
    vcf_w = open(output + ".vcf", "w")
    vcf_w.write(hd)
    print("###### Write vcf header ######")
    # chromosome set
    chrom = chrom_set(args)
    seq_dict = fasta_dict(args)
    for cr in chrom:
        seq_length = len(seq_dict[cr])
        # pysam.TabixFile read vcf file
        vcf_handle = pysam.TabixFile(vcf, parser = pysam.asVCF())
        vcf_ps = vcf_handle.fetch(cr, 0, seq_length)
        for v in vcf_ps:
            # here, v is in class 'pysam.libctabixproxies.VCFProxy'
            contig = str(v.contig)
            # NOTE: here, the position is 0-based, if I want to keep the vcf as the input 1-based then change the pos-value
            pos = str(v.pos+1)
            id = str(v.id)
            ref = str(v.ref)
            alt = str(v.alt)
            quality = str(v.qual)
            # filter = v.filter
            info = str(v.info)
            # info as string
            MQ = v.info.split("MQ=")[-1].split(";")[0]
            MQ = float(MQ)
            format = str(v.format)
            genotype = str(v[0])
            # genotype as string
            GT = genotype.split(":")[0]
            # print(contig, pos, id, ref, alt, quality, info, MQ, genotype, GT)
            if ((GT == "1|1" or GT == "1/1") and MQ >= 40.0):
                vcf_w.write(contig + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + quality + "\tPASS\t" + info + "\t" + format + "\t" + genotype + "\n")
    vcf_r.close()
    vcf_w.close()

def zip_vcf(args):
    output = args.output
    print("###### compressing vcf file ...... ")
    zip_command = f"bgzip -c {output}.vcf > {output}.vcf.gz"
    subprocess.call(zip_command, shell = True)

def index_vcf(args):
    output = args.output
    print("###### index vcf.gz file ...... ")
    index_command = f"tabix -p vcf {output}.vcf.gz"
    subprocess.call(index_command, shell = True)

def rm_unzip(args):
    output = args.output
    gzvcf_file = output + ".vcf.gz"
    if gzvcf_file in os.listdir():
        rm_command = f"rm {output}.vcf"
        subprocess.call(rm_command, shell = True)

def main():
    args = args_parser()
    new_vcf(args)
    zip_vcf(args)
    index_vcf(args)
#    rm_unzip(args)

############################
######### Run it ###########
############################

if __name__ == "__main__":
    main()
