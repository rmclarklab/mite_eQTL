import pandas as pd 
import csv 

gff3 = "Clark_2019/Tetranychus_urticae_v7.gff3"
gtf = "Clark_2019/Tetranychus_urticae_v7.gtf"
gff3 = pd.read_table(gff3, sep = "\t", header = None, comment = "#")
gff3 = gff3.rename(columns = {0:"tig", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"score", 6:"strand", 7:"phase", 8:"attributes"})
gff3_clean = gff3[gff3["feature"] != "assembly"]

### gene:mRNA dictionary
gff3_mRNA = gff3[gff3["feature"] == "mRNA"]
gene_mRNA_dict = dict()
for row in gff3_mRNA.itertuples():
    row_att = row.attributes
    row_mRNA = row_att.split(";")[0].split(":")[-1]
    row_gene = row_att.split(";")[1].split(":")[-1]
    gene_mRNA_dict[row_mRNA] = row_gene

### row by row to update
row_list = list()
for index, row in gff3_clean.iterrows():
    att = row.attributes
    if row.feature == "gene":
        gid = att.split("gene:")[-1].split(";")[0]
        glength = att.split("length=")[-1]
        row.attributes = f'gene_id \"{gid}\";length \"{glength}\"'
    elif row.feature == "mRNA":
        mid = att.split("mRNA:")[-1].split(";")[0]
        gid = gene_mRNA_dict[mid]
        mlen = att.split("length=")[-1]
        row.feature = "transcript"
        if "nb_exon" in att:
            exnb = att.split("nb_exon=")[-1].split(";")[0]
            row.attributes = f"gene_id \"{gid}\";transcript_id \"{mid}\";exon_number \"{exnb}\";length \"{mlen}\""
        else:
            row.attributes = f"gene_id \"{gid}\";transcript_id \"{mid}\";length \"{mlen}\""
    else:
        mid = att.split("mRNA:")[-1].split(";")[0]
        gid = gene_mRNA_dict[mid]
        keyid = att.split("ID=")[-1].split(";")[0].split(mid+".")[-1]
        onto = att.split("Ontology_term=")[-1]
        row.attributes = f"gene_id \"{gid}\";transcript_id \"{mid}\";exon_label \"{keyid}\";Ontology_term \"{onto}\""
    new_row = row.to_frame().T
    row_list.append(new_row)

row_df = pd.concat(row_list, axis = 0)
row_df.to_csv(gtf, sep = "\t", index = False, quoting=csv.QUOTE_NONE, header = False)
