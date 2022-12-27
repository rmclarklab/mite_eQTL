from Bio import SeqIO
import pandas as pd 

fa = "Orcae_2019.protein.fa"
seq = SeqIO.parse(fa, "fasta")
seq_dict = SeqIO.to_dict(seq)

good_protein = list()
bad_protein = list()
for entry in seq:
    id = entry.id
    seq = str(entry.seq)
    if len(seq) > 0 and seq.startswith("M") and seq.endswith("*") and seq.strip("*").find("*") == -1:
        good_protein.append(id)
    else:
        bad_protein.append(id)
    

####
gff3_orcae = pd.read_table("Orcae_01252019/scaffold.gff3", comment = "#", header = None)
gff3_clark = pd.read_table("Clark_2019/Tetranychus_urticae_2018.06.01.gff3", comment = "#", header = None)

def gff_geneid(gff):
    gff = gff.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
    gene_gff = gff[gff["feature"] == "gene"]
    geneid_list = list(gene_gff.attributes.str.split(";", expand = True)[[0]][0])
    return geneid_list

orcae_geneid = gff_geneid(gff3_orcae)
clark_geneid = gff_geneid(gff3_clark)
orcae_geneid = [o.split("=")[-1] for o in orcae_geneid]
clark_geneid = [c.split(":")[-1] for c in clark_geneid]

Orcae_gdf = pd.DataFrame({"gene" : orcae_geneid})
clark_gdf = pd.DataFrame({"gene" : clark_geneid})

Orcae_gdf.to_csv("data/Orcae_gid.txt", sep = "\t", index = False)
clark_gdf.to_csv("data/Clark_gid.txt", sep = "\t", index = False)

####
orcae_gid = pd.read_table("data/Orcae_gid.txt", sep = "\t", header = 0)
clark_gid = pd.read_table("data/Clark_gid.txt", sep = "\t", header = 0)

orcae_gid = list(orcae_gid["gene"])
clark_gid = list(clark_gid["gene"])

overlapped_gid = list(set(orcae_gid).intersection(set(clark_gid)))

overlapped_gdf = pd.DataFrame({"gene" : overlapped_gid})
overlapped_gdf.to_csv("data/overlapped_gid.txt", sep = "\t", index = False)

####
overlapped_gid = pd.read_table("data/overlapped_gid.txt", sep = "\t", header = 0)
overlapped_gid = list(overlapped_gid["gene"])

fa1 = "Clark_2019/Tetur_2018.protein.fa"
fa2 = "Orcae_01252019/Orcae_2019.protein.fa"

seq1 = SeqIO.parse(fa1, "fasta")
seq2 = SeqIO.parse(fa2, "fasta")
seq1_dict = SeqIO.to_dict(seq1)
seq2_dict = SeqIO.to_dict(seq2)

agree_gid = list()
disagree_gid = list()
for g in overlapped_gid:
    seq1_str = str(seq1_dict[g].seq)
    seq2_str = str(seq2_dict[g].seq)
    if seq1_str == seq2_str:
        agree_gid.append(g)
    else:
        disagree_gid.append(g)

disagree_df = pd.DataFrame({"gene" : disagree_gid})
disagree_df.to_csv("data/disagree_gid.txt", sep = "\t", index = False)

###
disagree_df = pd.read_table("data/disagree_gid.txt", sep = "\t", header = 0)

fa1 = "Clark_2019/Tetur_2018.protein.fa"
fa2 = "Orcae_01252019/Orcae_2019.protein.fa"

seq1 = SeqIO.parse(fa1, "fasta")
seq2 = SeqIO.parse(fa2, "fasta")
seq1_dict = SeqIO.to_dict(seq1)
seq2_dict = SeqIO.to_dict(seq2)

fw = open("data/disagree_gene.fa", "w")

for row in disagree_df.itertuples():
    gene = row.gene
    seq1_str = str(seq1_dict[gene].seq)
    seq2_str = str(seq2_dict[gene].seq)
    fw.write(f">{gene}_clark\n{seq1_str}\n")
    fw.write(f">{gene}_orcae\n{seq2_str}\n")

fw.close()

###
orcae_gid = pd.read_table("data/Orcae_gid.txt", sep = "\t", header = 0)
clark_gid = pd.read_table("data/Clark_gid.txt", sep = "\t", header = 0)

orcae_gid = list(orcae_gid["gene"])
clark_gid = list(clark_gid["gene"])

orcae_unique = list(set(orcae_gid).difference(set(clark_gid)))
orcae_unique_gdf = pd.DataFrame({"gene" : orcae_unique})
clark_unique = list(set(clark_gid).difference(set(orcae_gid)))
clark_unique_gdf = pd.DataFrame({"gene" : clark_unique})
orcae_unique_gdf.to_csv("data/Orcae_unique.txt", sep = "\t", index = False)
clark_unique_gdf.to_csv("data/Clark_unique.txt", sep = "\t", index = False)

fa1 = "Clark_2019/Tetur_2018.protein.fa"
fa2 = "Orcae_01252019/Orcae_2019.protein.fa"

seq1 = SeqIO.parse(fa1, "fasta")
seq2 = SeqIO.parse(fa2, "fasta")
seq1_dict = SeqIO.to_dict(seq1)
seq2_dict = SeqIO.to_dict(seq2)

clark_fw = open("data/clark_unique.fa", "w")
for c in clark_unique:
    clark_seq = str(seq1_dict[c].seq)
    clark_fw.write(f">{c}\n{clark_seq}\n")

clark_fw.close()

orcae_fw = open("data/orcae_unique.fa", "w")
for o in orcae_unique:
    orcae_seq = str(seq2_dict[o].seq)
    orcae_fw.write(f">{o}\n{orcae_seq}\n")

orcae_fw.close()

###
fa2 = "Orcae_01252019/Orcae_2019.protein.fa"
seq2 = SeqIO.parse(fa2, "fasta")
seq2_dict = SeqIO.to_dict(seq2)

fa2_2 = "Orcae_01252019/mRNA_tetur_active_pep_20190125.tfa"
seq2_2 = SeqIO.parse(fa2_2, "fasta")
seq2_2_dict = SeqIO.to_dict(seq2_2)

seq_translate = seq2_dict.keys()
seq_db = seq2_2_dict.keys()

translated_unique = set(seq_translate).difference(set(seq_db)) # 15
db_unique = set(seq_db).difference(set(seq_translate)) # 3
both_on = list(set(seq_db).intersection(set(seq_translate)))

agree_gid1 = list()
disagree_gid1 = list()
for b in both_on:
    seqb1 = str(seq2_dict[b].seq)
    seqb2 = str(seq2_2_dict[b].seq)
    if seqb1 == seqb2:
        agree_gid1.append(b)
    else:
        disagree_gid1.append(b)

len(agree_gid1)
len(disagree_gid1)

####

orcae_unique_translate = SeqIO.parse("data/orcae_unique.fa", "fasta")
orcae_prot = SeqIO.parse("Orcae_01252019/mRNA_tetur_active_pep_20190125.tfa", "fasta")
orcae_unique_translate_dict = SeqIO.to_dict(orcae_unique_translate)
orcae_prot_dict = SeqIO.to_dict(orcae_prot)

agree_gid2 = list()
disagree_gid2 = list()
weired_gid = list()
for entry in orcae_unique_translate_dict:
    unique_id1 = orcae_unique_translate_dict[entry].id
    if unique_id1 in orcae_prot_dict:
        unique_seq1 = str(orcae_unique_translate_dict[entry].seq)
        unique_seq2 = str(orcae_prot_dict[unique_id1].seq)
        if unique_seq1 == unique_seq2:
            agree_gid2.append(unique_id1) # 28
        else:
            disagree_gid2.append(unique_id1) # 5
    else:
        weired_gid.append(unique_id1) # 12 genes in GFF but not in protein tfa

def write_gene(gene, out):
    fw = open(out + ".txt", "w")
    for g in gene:
        fw.write(f"{g}\n")
    fw.close()



#### 45 unique genes in the Orcae database
orcae_df = pd.read_table("data/Orcae_unique.txt", sep = "\t", header = 0)
orcae_unique = list(orcae_df["gene"])
snoeck = pd.read_table("Snoeck_2019_NHR_from_Marilou.txt", sep = "\t", header = None, comment = "#")
snoeck_nhr = list(snoeck[0])

snoeck_miss = set(snoeck_nhr).intersection(set(orcae_unique))

orcae_other = set(orcae_unique).difference(set(snoeck_nhr))

### 12 genes not in protein tfa file
noprotein = pd.read_table("data/weired_12w45.txt", sep = "\t", header = None)
noprotein = list(noprotein[0])
translate_prot = SeqIO.parse("Orcae_01252019/Orcae_2019.protein.fa", "fasta")
translate_dict = SeqIO.to_dict(translate_prot)

fw2 = open("weired_protein.fa", "w")
for n in noprotein:
    seq_str = str(translate_dict[n].seq)
    fw2.write(f">{n}\n{seq_str}\n")

fw2.close()

#### If Simon's gene good in Clark GFF
# fa1 = "Clark_2019/Tetur_2018.protein.fa"
fa1 = "Orcae_01252019/Orcae_2019.protein.fa"
seq1 = SeqIO.parse(fa1, "fasta")
seq1_dict = SeqIO.to_dict(seq1)

fa3 = "Orcae_01252019/mRNA_tetur_active_pep_20190125.tfa"
seq3 = SeqIO.parse(fa3, "fasta")
seq3_dict = SeqIO.to_dict(seq3)

snoeck = pd.read_table("Snoeck_2019_NHR_from_Marilou.txt", sep = "\t", header = None, comment = "#")
snoeck_nhr = list(snoeck[0])

ks = [k for k in seq1_dict.keys() if k in snoeck_nhr] # 42

mis_nhr = list()
i = 0
for k in ks:
    seq1_str = str(seq1_dict[k].seq)
    seq3_str = str(seq3_dict[k].seq)
    if seq1_str == seq3_str:
        i += 1
    else:
        mis_nhr.append(k)

write_gene(mis_nhr, "orcae_gff_mis_nhr")



### This is a record for GFF reconstruction from Tetranychus_urticae_2018.06.01.gff3

# First, pick NHRs and make sure all genes (except tetur11g05480) can be translated perfectly
gene_info = pd.read_table("data/Tetur_gene_info_20190125.txt", sep = "\t", header = 0)
nhr_info = gene_info[(gene_info.protein_domains.str.contains("nuclear receptor", case = False)) | (gene_info.protein_domains.str.contains("hormone receptor", case = False)) | (gene_info.description.str.contains("nuclear receptor", case = False)) | (gene_info.description.str.contains("hormone receptor", case = False))]
nhr_info[["Gene_ID"]].to_csv("NHR_id.txt", sep = "\t", index = False)

# temp (make sure all simon's NHRs96 in nhr_info)
snoeck = pd.read_table("data/Snoeck_2019_NHR_from_Marilou.txt", sep = "\t", header = None, comment = "#")
snoeck_nhr = list(snoeck[0])
set(snoeck_nhr).difference(set(nhr_info["Gene_ID"])) # YES!

# Then, make sure all selected NHRs can be correctly translated in Orcae GFF3
gff3_orcae = pd.read_table("Orcae_01252019/scaffold.gff3", comment = "#", header = None)

def gff3_subset_select(gff, geneid):
    gff_df = gff.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
    gff_df = gff_df[gff_df["feature"] != "sequence_assembly"]
    gene_list = list()
    for row in gff_df.itertuples():
        if row.feature == "gene":
            gid = row.attributes.split(";")[0].split("ID=")[-1]
        else:
            gid = row.attributes.split("Parent=")[-1].split(";")[0].split(".")[0]
        gene_list.append(gid)
    gff_df["gene"] = gene_list
    gene_df_new = gff_df[gff_df.gene.isin(geneid)]
    return gene_df_new

nhr_gff3 = gff3_subset_select(gff3_orcae, list(nhr_info["Gene_ID"]))
nhr_gff3 = nhr_gff3.drop("gene", axis = 1)
nhr_gff3.to_csv("NHR.gff3", sep = "\t", index = False, header = False)

### run translation in shell
python code/cds_translate.py -gff3 NHR.gff3 -fasta Orcae_01252019/Tetranychus_urticae.main_genome_200909.scaffolds.fasta -O NHR_translate

# Then, report genes not good as it is in tfa
tfa = "Orcae_01252019/mRNA_tetur_active_pep_20190125.tfa"
tfa_seq = SeqIO.parse(tfa, "fasta")
tfa_dict = SeqIO.to_dict(tfa_seq)

nhrfa = "NHR_translate.protein.fa"
nhr_seq = SeqIO.parse(nhrfa, "fasta")
nhr_dict = SeqIO.to_dict(nhr_seq)

for entry in nhr_dict:
    translate_prot = str(nhr_dict[entry].seq)
    tfa_prot = str(tfa_dict[entry.split(".")[0]].seq)
    if translate_prot == tfa_prot:
        pass
    else:
        print(entry)

# tetur11g05480.1 ***
# tetur30g02190.1 :)
# tetur46g90320.1 :)
# tetur46g90340.1 :)
# tetur87g90401.1 (need to check in IGV)

# Then, drop all NHR genes in Clark gff (version 1)
gff3_clark = pd.read_table("Clark_2019/Tetranychus_urticae_2018.06.01.gff3", comment = "#", header = None)
nhr = pd.read_table("NHR_id.txt", sep = "\t", header = 0)

def gff3_subset_delete(gff, geneid):
    gff_df = gff.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
    gff_df = gff_df[gff_df["feature"] != "assembly"]
    gene_list = list()
    for row in gff_df.itertuples():
        if row.feature == "gene":
            gid = row.attributes.split(";")[0].split("ID=gene:")[-1]
        else:
            gid = row.attributes.split("Parent=")[-1].split(":")[1].split(";")[0].split(".")[0]
        gene_list.append(gid)
    gff_df["gene"] = gene_list
    gene_df_new = gff_df[~gff_df.gene.isin(geneid)]
    return gene_df_new

clark_v1 = gff3_subset_delete(gff3_clark, list(nhr["Gene_ID"]))
clark_v1 = clark_v1.drop("gene", axis = 1)
clark_v1.to_csv("Tetranychus_urticae_v1.gff3", sep = "\t", index = False, header = False)

# Then, put tetur11g05480 coordinates from Marilou to gff as version 2
### see Tetranychus_urticae_v2.gff3

# read Orcae GFF3 and Clark GFF3 for unique geneid in each dataset
gff3_orcae = pd.read_table("Orcae_01252019/scaffold.gff3", comment = "#", header = None)
gff3_clark = pd.read_table("Clark_2019/Tetranychus_urticae_v2.gff3", comment = "#", header = None)

def gff_geneid(gff):
    gff = gff.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
    gene_gff = gff[gff["feature"] == "gene"]
    geneid_list = list(gene_gff.attributes.str.split(";", expand = True)[[0]][0])
    return geneid_list

orcae_geneid = gff_geneid(gff3_orcae)
clark_geneid = gff_geneid(gff3_clark)

orcae_geneid = [o.split("=")[-1] for o in orcae_geneid]
clark_geneid = [c.split(":")[-1] for c in clark_geneid]

orcae_unique = sorted(set(orcae_geneid).difference(set(clark_geneid))) # 151
clark_unique = sorted(set(clark_geneid).difference(set(orcae_geneid))) # 18

write_gene(orcae_unique, "orcae_shift")
write_gene(clark_unique, "clark_unique")

# Pull out subset of Orcae GFF3 that are ready for shift to clark GFF
gff3_orcae = pd.read_table("Orcae_01252019/scaffold.gff3", comment = "#", header = None)
orcae2shift = pd.read_table("orcae_shift.txt", sep = "\t", header = None)

gff3_orcae2shift = gff3_subset_select(gff3_orcae, list(orcae2shift[0]))
gff3_orcae2shift = gff3_orcae2shift.drop("gene", axis = 1)

gff3_orcae2shift.to_csv("data/Orcae2add.gff3", sep = "\t", index = False)

# esatablish scaffold to chromosome alignment relationship (critical step)
coord = pd.read_table("data/scaffold2chr_coordinate.txt", sep = "\t", header = 0)

range_dict = dict()
for row in coord.itertuples():
    scaff_post = row.scaffold_break
    scaff_direction = row.direction
    if scaff_direction == "+":
        range_dict[scaff_post] = [range(row.scaffold_start, row.scaffold_end, 1), range(row.start, row.end, 1)]
    elif scaff_direction == "-":
        range_dict[scaff_post] = [range(row.scaffold_end, row.scaffold_start, -1), range(row.start, row.end, 1)]
    else:
        print("strand not indicated. ")

# shift the location from scaffold to chromosome based on coordinates association 
gff3_orcae_use = pd.read_table("data/Orcae2add.gff3", sep = "\t", header = 0)

row_list = list()
strand_dict = {"+":"-", "-":"+", ".":"."}
for index, row in gff3_orcae_use.iterrows():
    rtig = row.tig
    rstart = row.start
    rend = row.end
    rstrand = row.strand
    scaf_search_sub = coord[(coord["scaffold_before_break"] == rtig) & (coord["scaffold_start"] < rstart) & (coord["scaffold_end"] > rend)]
    if len(scaf_search_sub) > 0:
        scaf_search = scaf_search_sub["scaffold_break"].item()
        chr_search = scaf_search_sub["chromosome"].item()
        sstart = range_dict[scaf_search][0].index(rstart)
        cstart = list(range_dict[scaf_search][1])[sstart]
        send = range_dict[scaf_search][0].index(rend)
        cend = list(range_dict[scaf_search][1])[send]
        row.tig = "chromosome_" + str(chr_search)
        row.start = cstart
        row.end = cend
        if row.start > row.end:
            row.strand = strand_dict[rstrand]
            row.start = cend
            row.end = cstart
    row_dict = row.to_dict()
    row_list.append(row_dict)

df_new = pd.DataFrame(row_list)
df_new.to_csv("data/Orcae2add_chr.gff3", sep = "\t", index = False)

# drop all CDS: in Orcae2add_chr.gff3 

# format chromosome scale gff3 
gff_new = pd.read_table("data/Orcae2add_chr.gff3", sep = "\t", header = 0)
parent_dict = {"CDS":"mRNA", "exon":"mRNA", "mRNA":"gene", "three_prime_UTR":"mRNA", "five_prime_UTR":"mRNA"}

gene_list_record = list()
gff3_format = open("data/Orcae2add_chr_format.gff3", "w")

for index, row in gff_new.iterrows():
    rfeature = row.feature
    rattributes = row.attributes
    if rfeature != "gene":
        id = rattributes.split("ID=")[-1].split(";")[0]
        parent = rattributes.split("Parent=")[-1].split(";")[0]
        id_label = "ID=" + rfeature +  ":" + id + ";"
        name_label = "Parent=" + parent_dict[rfeature] + ":" + parent + ";"
        parent_left = ";".join(rattributes.split("Parent=")[-1].split(";")[1:])
        new_attributes = id_label + name_label + parent_left
    else:
        id = rattributes.split("ID=")[-1].split(";")[0]
        id_left = ";".join(rattributes.split("ID=")[-1].split(";")[1:])
        id_label = "ID=" + rfeature + ":" + id + ";"
        new_attributes = id_label + id_left
        #if id not in gene_list_record:
        #    gene_list_record.append(id)
        #    gff3_format.write("###\n")
    row.attributes = new_attributes
    row_l = [str(r) for r in row.tolist()]
    row_w = "\t".join(row_l)
    gff3_format.write(f"{row_w}\n")

gff3_format.close()

# format Tetranychus_urticae_v2.gff3 for mRNA number label (.1 suffix)
gff_v2 = pd.read_table("Clark_2019/Tetranychus_urticae_v2.gff3", sep = "\t", header = None)
gff_v2 = gff_v2.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
gff_2_gene = sorted(set(gff_v2[gff_v2["feature"] == "gene"]["attributes"].str.split(";", expand = True)[0].str.split(":", expand = True)[1]))
gff_v2_mRNA = gff_v2[gff_v2["feature"] == "mRNA"]
for gg in gff_2_gene:
    gff_v2_mRNA_sub = gff_v2_mRNA[gff_v2_mRNA["attributes"].str.contains(gg)]
    if len(gff_v2_mRNA_sub) == 1:
        pass
    else:
        print(gg + " more than one mRNA") # examine all gene has and only has one mRNA (transcript) YES!

row_list = list()
for index, row in gff_v2.iterrows():
    attr = row.attributes
    attr = ";".join([a + ".1" if "mRNA:" in a else a for a in attr.split(";")])
    attr = ";".join([a.split(".")[0] + ".1." + a.split(".")[-1] if ("CDS:" in a) or ("exon:" in a) or ("five_prime_UTR:" in a) or ("three_prime_UTR:" in a) else a for a in attr.split(";")])
    row.attributes = attr
    row_dict = row.to_dict()
    row_list.append(row_dict)

df_new = pd.DataFrame(row_list)
df_new.to_csv("Clark_2019/Tetranychus_urticae_v3.gff3", sep = "\t", index = False, header = False)

# combine Tetranychus_urticae_v2.gff3 and Orcae2add_chr_format.gff3
### run: cat Clark_2019/Tetranychus_urticae_v3.gff3 data/Orcae2add_chr_format.gff3 > Tetranychus_urticae_v4.gff3

# translate gff3 to check the quality 
### run: mpiexec -n 6 python code/cds_translate.py -fasta Clark_2019/Tetranychus_urticae_2017.11.21.fasta -gff3 Tetranychus_urticae_v4.gff3 -O v4_translate

import re 
import pandas as pd 

gff3_v3 = pd.read_table("Clark_2019/Tetranychus_urticae_v4.gff3", sep = "\t", header = None)
gff3_v3 = gff3_v3.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
row_list = list()
for index, row in gff3_v3.iterrows():
    attr = row.attributes
    attr_dict = {a.split(":")[0]:a.split(":")[-1] for a in attr.split(";")}
    attr_dict_new = dict()
    for a in attr_dict:
        if "mRNA" in a:
            attr_dict_new[a] = re.sub("\.1$", "", attr_dict[a])
        else:
            attr_dict_new[a] = attr_dict[a]
    attr_list = [a + ":" + attr_dict_new[a] for a in attr_dict_new]
    attr_new = ";".join(attr_list)
    row.attributes = attr_new
    row_dict = row.to_dict()
    row_list.append(row_dict)

df_new = pd.DataFrame(row_list)
df_new.to_csv("Tetranychus_urticae_v5.gff3", sep = "\t", index = False, header = False)

# add header line to Tetranychus_urticae_v4.gff3
### run: cat header.gff3 Tetranychus_urticae_v5.gff3 > Tetranychus_urticae_v6.gff3

# sort gff3 (on Adele)
### run: perl ~/Downloads/gff3sort-master/gff3sort.pl Tetranychus_urticae_v6.gff3 > Tetranychus_urticae_v7.gff3 

# fix tetur14g02500 overlaped CDS error in Tetranychus_urticae_v6.gff3 
# save as Tetranychus_urticae_v7.gff3 

# provide GTF format (change the gff3 and gtf file name accordingly) for Tetranychus_urticae_v7.gtf
### run: python gff2gtf.py

# add header for GTF 
### run: cat header.gtf Tetranychus_urticae_v7.gtf > Tetranychus_urticae_v8.gtf

# save the final version as 2022 GFF
### run: cp Tetranychus_urticae_v8.gtf Tetranychus_urticae_2022.gtf
