#!/usr/bin/env python

import pandas as pd
import argparse
import numpy as np
import os
import random 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description=" \n \
        Usage: genotype_block.py -chrLen <chrlen> -count <cnt.txt>")
    parser.add_argument("-count", "--count", help = "allele-specific read count with genotype data. ")
    parser.add_argument("-chrLen", '--chrLen', help = "chromosome length. ")
    parser.add_argument("-extension_SNP", "--extension_SNP", type = int, default = 10, help = "extension SNPs for possible breakpoint check over the shift genotypes barrier. \n Recommend number of > 5. (default: 10). ")
    parser.add_argument("-extension_range", "--extension_range", type = int, default = 250000, help = "extension range for possible breakpoint check over the shift genotypes barrier. \n Recommend number of > 100000. (default: 250000 bp). ")
    parser.add_argument("-coverage_control", "--coverage_control", action = "store_true", help = "Control coverage depth for genotype information. ")
    parser.add_argument("-quality_coverage", "--quality_coverage", type = int, required = False, default = 25, help = "read-depth coverage as a quality site for genotype desicion (default: 25). ")
    parser.add_argument("-dilute", "--dilute_SNP", action = "store_true", help = "Dilute SNP information to prevent over-density SNP overtaking genotype call. ")
    parser.add_argument("-dilute_step", "--dilute_step", type = int, required = False, default = 100, help = "Step site for diluting over-density SNPs (default: 100 bp). ")
    parser.add_argument("-O", "--output", required = False, default = "block", help = "Prefix for the output file. ")
    args = parser.parse_args()
    return args

def initial_correction(subcnt):
    # the start point genotype of a chromosome should be double check and is correct, so that start from a correct start
    # the first 10 rows are picked
    start_index = subcnt.index.min()
    end_index = start_index + 10
    checkcnt = subcnt[subcnt.index.isin(range(start_index, end_index))]
    checktype = sorted(set(checkcnt.genotype))
    if len(checktype) == 1:
        t1 = checktype[0]
        return (t1, [])
    else:
        check_dict = {c:list(checkcnt.genotype).count(c) for c in checktype}
        ## to get the correct start genotype
        c0 = check_dict[checktype[0]]
        for c in checktype:
            if check_dict[c] > c0:
                c0 = check_dict[c]
        t1 = [g for g in check_dict if check_dict[g] == c0][0]
        drop_list = checkcnt[checkcnt["genotype"] != t1].index.tolist()
        return (t1, drop_list)

def locate_recombination(args, cnt, r_list, chrlen_dict, out):
    locate_f = open(out+".txt", "w") # chromosome genotype start end info1 info2 total noise
    locate_f.write("chromosome\tgenotype\tstart\tend\tleft\tright\ttotal\tnoise\n")
    for r in r_list:
        print(f"Processing {r} ......")
        subcnt = cnt[cnt['chromosome'] == r].copy()
        subcnt = subcnt.sort_values(by = ["position"])
        subcnt.index = list(range(0, len(subcnt)))
        # correct the starting point of genotype
        t1, drop_list = initial_correction(subcnt)
        if len(drop_list) > 0:
            subcnt.drop(drop_list, axis = 0, inplace = True)
        ###
        rlen = chrlen_dict[r]
        # noise = []
        sr = 1 # block start 
        pos1 = subcnt.position.min()
        end_pos = subcnt.position.max()
        nis = 0 # total noisy SNP number
        tsnp = 0 # total SNP number
        ###
        for row in subcnt.itertuples():
            tsnp += 1
            ti = t1 # store last row genotype
            posi = pos1 # store last row position
            t1 = row.genotype
            pos1 = row.position # information for current row
            if (t1 != ti): # potential break
                if break_bh(args, subcnt, pos1, t1): # real breakpoint
                    bp = (posi + pos1)//2 # midpoint as breakpoint
                    print(f"({posi}-{pos1}) transition on {r}")
                    locate_f.write(f"{r}\t{ti}\t{sr}\t{bp}\t{posi}\t{pos1}\t{tsnp}\t{nis}\n")
                    sr = bp + 1 # new start of a genotype block
                    nis = 0
                    tsnp = 0
                else: # false breakpoint
                    nis += 1
                    subcnt.loc[(subcnt.position == pos1), 'genotype'] = ti # correct genotype from last row 
                    t1 = ti
                    # noise.append(pos1) # add the position as noise position
            if row.position == end_pos: # capture the last row
                locate_f.write(f"{r}\t{t1}\t{sr}\t{rlen}\t-\t-\t{tsnp}\t{nis}\n") # write the corrected genotype into dataframe 
    locate_f.close()
    locate_df = pd.read_table(out+".txt", sep = "\t", header = 0)
    os.remove(out+".txt")
    return locate_df
    
def refine_recombination(locate_df, count_all):
    index_list = locate_df.index.values.tolist()
    for i in index_list:
        rchr = locate_df.loc[i].chromosome
        rleft = locate_df.loc[i].left
        rright = locate_df.loc[i].right
        count_locate = count_all[(count_all["chromosome"] == rchr) & (count_all.position.isin(range(int(rleft), int(rright)+1)))].copy()
        gleft = count_locate[count_locate.position == int(rleft)].genotype.values[0]
        gright = count_locate[count_locate.position == int(rright)].genotype.values[0]
        count_left = count_locate[count_locate.genotype == gleft].copy()
        count_right = count_locate[count_locate.genotype == gright].copy()
        if len(count_left) > 1:
            left_block = consecutiveSNP_genotype(count_left)
            rleft = left_block.position.max()
        if len(count_right) > 1:
            right_block = consecutiveSNP_genotype(count_right)
            rright = right_block.position.min()
        # print(rleft, rright, type(rleft), type(rright))
        if int(rleft) < int(rright): 
            # only if rleft and rright in the correct order, to update the dataframe information
            locate_df.loc[locate_df.index == i, "left"] = rleft
            locate_df.loc[locate_df.index == i, "right"] = rright
    return locate_df # updated with more allele signal when it is available

def consecutiveSNP_genotype(count_sub):
    pos_order = np.array(count_sub.index.tolist())
    pos_diff = [0] + np.diff(pos_order).tolist()
    area_bk = [idx for idx,val in enumerate(pos_diff) if val > 2]
    # count_sub.index = list(range(0, len(count_sub)))
    if len(area_bk) > 0:
        n = 0
        df_list = list()
        for a in area_bk:
            sub_cnt = count_sub.iloc[n:a]
            n = a
            df_list.append(sub_cnt)
        sub_cnt = count_sub.iloc[n:]
        df_list.append(sub_cnt)
        # return the longest part of dataframe genotype area
        df_longest = max(df_list, key = len)
        return df_longest
    else:
        return count_sub

def break_bh(args, subcnt, pos, t1):
    snp_num = args.extension_SNP
    range_size = args.extension_range
    short_break_check = False
    ### short distance check (the SNP positioin and its following 5 number of SNPs in every 100 bp)
    pos_index = subcnt[subcnt["position"] == pos].index.values[0]
    if subcnt.index.max() - pos_index < snp_num:
        return False
    else:
        shortcnt = subcnt[subcnt["position"].index.isin(range(pos_index, pos_index+snp_num))]
        if list(shortcnt.genotype).count(t1)/len(shortcnt)*100 <= 60:
            return False
        else:
            short_break_check = True
    ### long distance check (the SNP position: following genotype of some distance and some number of SNPs)
    if subcnt.position.max() - pos > range_size*0.2 and short_break_check:
        geno_cnt = subcnt[subcnt['position'].isin(range(pos, pos+range_size))]
        snp_row = len(geno_cnt)
        if snp_row < int(snp_num*2.5): # control factor for regions when SNP density is scarce
            # print("pick more SNPs")
            geno_cnt = subcnt[subcnt["position"].index.isin(range(pos_index, pos_index + int(snp_num*2.5)))]
        t_perc = list(geno_cnt.genotype).count(t1)/len(geno_cnt)*100
        if (t_perc >= 80): 
            return True
    else:
        return False

def generate(start, stop, step):
    for p in range(start, stop, step):
        yield (p, p + step)

def dilute_SNP(count_all, dilute_step):
    count_chr = sorted(set(count_all["chromosome"]))
    count_dilute = list()
    for c in count_chr:
        sample_pos = list()
        count_sub = count_all[count_all["chromosome"] == c]
        all_pos = np.array(count_sub.position.tolist())
        min_pos = count_sub.position.min()
        max_pos = count_sub.position.max()
        for range1, range2 in generate(min_pos, max_pos, dilute_step):
            sub_pos = all_pos[(all_pos>=range1)*(all_pos<range2)].tolist()
            if len(sub_pos) > 0:
                sample_pos.append(random.sample(sub_pos, 1)[0])
        sample_chr = count_sub[count_sub.position.isin(sample_pos)]
        count_dilute.append(sample_chr)
    count_dilute_all = pd.concat(count_dilute, axis = 0)
    return count_dilute_all

def main():
    args = args_parser()
    out = args.output
    # chromosome length directionary 
    chrLen = args.chrLen
    chrlen = pd.read_table(chrLen, sep = "\t", header = 0, index_col = 0)
    chrlen_dict = chrlen.to_dict()["length"]
    # 
    out = args.output
    count = args.count
    count_all = pd.read_table(count, sep = '\t', header = 0, dtype = {"genotype": "str"})
    if args.dilute_SNP:
        print("Dilute over-density SNPs ......")
        count_all = dilute_SNP(count_all, args.dilute_step)
        print("Write new count file ......")
        count_all.to_csv(out+".scarse.txt", sep = "\t", index = False)
    if args.coverage_control:
        print(f"Read depth of {args.quality_coverage} is controling for genotype information ......")
        coverage = args.quality_coverage
        count_all = count_all[count_all["total_allele"] > coverage]
    # chr_quality = sorted(set(count_quality["chromosome"]))
    # print(chr_quality)
    # r1_recombination = locate_recombination(args, count_quality, chr_quality, chrlen_dict, out+".temp1")
    # print(r1_recombination)
    # refine_df = r1_recombination[r1_recombination["left"] != "-"].copy()
    # non_refine_df = r1_recombination[r1_recombination["left"] == "-"].copy()
    # refine_new = refine_recombination(refine_df, count_all)
    # assigned_df = pd.concat([refine_new, non_refine_df], axis = 0)
    # print(assigned_df)
    # unassigned_chr = set(count_all["chromosome"]).difference(set(assigned_df["chromosome"]))
    chrs = sorted(set(count_all['chromosome']))
    r_recombination = locate_recombination(args, count_all, chrs, chrlen_dict, out+".temp1")
    # 
    # all_assigned =  pd.concat([assigned_df, r2_recombination], axis = 0)
    all_assigned = r_recombination.sort_values(by = ["chromosome", "start"])
    all_assigned.to_csv(out+".txt", sep = "\t", index = False)

##################
##### Run it #####
##################

if __name__ == "__main__":
    main()