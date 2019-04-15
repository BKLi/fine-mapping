# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:03:13 2019

@author: libingkun

Can be run block by block
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
trait = "BP"
cell_types = ["astrocyte", "motor", "cortical", "hippocampal"]


## overall pre-processing
header = "index rsid chromosome position allele1 allele2 maf beta se z prob log10bf group corr_group prob_group log10bf_group mean sd mean_incl sd_incl"
header = header.split(" ")
all_snps = pd.read_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\snp\{}_all.snp'.format(trait, trait, trait), delim_whitespace=True, names=header)
sig_snps = all_snps[(all_snps["prob"] > 0.001) & (all_snps["maf"] > 0.001)] \
.drop_duplicates(subset=['rsid'],keep='first')

# dropped missed snps stes
## prepare chipseeker input -- **NEXT STEP: liftover**
for_cs = sig_snps[['chromosome', 'position', "rsid"]]
# add in finemap-non-eligible snps --- AD_noneli.set
notelg = pd.read_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_noteli.set'.format(trait, trait, trait), names=['chromosome', 'position', "rsid"], sep="\t")
for_cs = pd.concat([for_cs, notelg],ignore_index=True).drop_duplicates(subset=['rsid'],keep='first')
for_cs["end"] = for_cs['position'].apply(lambda x: x+1)
identifier = ["{}_{}".format(trait, i) for i in range(0, for_cs.shape[0])]
for_cs.loc[:,"ID"] = identifier
for_cs = for_cs.rename(index=str, columns={"position":"start", "chromosome":"chr"})
for_cs = for_cs[["chr","start","end","rsid","ID"]]
# input for liftover --> chipseeker --> motifbreakR
for_cs.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}_all.snps.chipseeker.bed'.format(trait), sep="\t", index=False, header=True)
## next step -- Chioseeker

## gtf preprocessing
gtf = pd.read_csv(r"C:\Users\libin\UCSF\gene_annotation\gencode.v19.annotation.gtf", sep="\t", comment="#",
                         names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
gtf_genes = gtf[gtf["type"] == "gene"]
gtf_genes["gene_id"] = gtf_genes["annotation"].str.extract(r'gene_id "(.+?)";')
gtf_genes["gene_name"] = gtf_genes["annotation"].str.extract(r'gene_name "(.+?)";')
gtf_genes["gene_type"] = gtf_genes["annotation"].str.extract(r'gene_type "(.+?)";')

gtf_genes_pos = gtf_genes[gtf_genes["strand"] == "+"]
gtf_genes_neg = gtf_genes[gtf_genes["strand"] == "-"]

gtf_genes_pos_TSS = gtf_genes_pos[["chr","start","gene_id", "gene_name"]]
gtf_genes_neg_TSS = gtf_genes_neg[["chr","end","gene_id", "gene_name"]]
gtf_genes_pos_TSS.to_csv(r"C:\Users\libin\UCSF\gene_annotation\gtf_pos_tss.bed", sep="\t", index=False, header=True)
gtf_genes_neg_TSS.to_csv(r"C:\Users\libin\UCSF\gene_annotation\gtf_neg_tss.bed", sep="\t", index=False, header=True)


## prepare bedtools intersect input (atac-peak and snps)
for cell_type in cell_types:
    peak_file = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}.atac.seq.narrowPeak'.format(cell_type), delim_whitespace=True, 
                        names=["chr","start","end","strand","null","score","fold_change","-log10pvalue","-log10qvalue","summitPosition"])
    peak_identifier=["{}_{}".format(cell_type, i) for i in range(0, peak_file.shape[0])]
    peak_file["peakID"] = peak_identifier
    peak_file = peak_file[['chr',
                     'start',
                     'end',
                     '-log10pvalue',
                     'peakID']]
    # peak_file['-log10pvalue'] = peak_file['-log10pvalue'].apply("float64")
    peak_file.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}.atac.peak.bed'.format(cell_type), sep="\t", index=False, header=False)

hg19_SNP_file = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_all.snps.chipseeker.hg19.bed'.format(trait, trait), sep="\t")
# hg19_SNP_file = hg19_SNP_file[['seqnames', 'start', 'end','rsid', 'ID']]
hg19_SNP_file.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_all.snps.bedtools.hg19.bed'.format(trait, trait), sep="\t", index=False, header=False)
## next step -- bedtools intersect


### chipseeker output data cleaning
cs_out = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_snp.chipseek.output.bed'.format(trait), sep="\t")
cs_out_anno_logical_full = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_anno_full'.format(trait, trait), sep="\t", names=["five_UTR", "three_UTR", "intron", "exon"])
# remove reformat annotation
# cs_out = cs_out[~cs_out['annotation'].str.contains("Downstream")]
cs_out.loc[cs_out['annotation'].str.contains("Intergenic"), "annotation"] = "Intergenic"
cs_out.loc[cs_out['annotation'].str.contains("Downstream"), "annotation"] = "Intergenic"
cs_out.loc[cs_out['annotation'].str.contains("Intron"), "annotation"] = "Intron"
cs_out.loc[cs_out['annotation'].str.contains("Exon"), "annotation"] = "Exon"
# annotate promoter : +/- 500bp TSS
cs_out["promoter"] = False
cs_out["promoter"] = np.where(abs(cs_out['distanceToTSS']) < 500, True, cs_out["promoter"])
cs_out = pd.concat([cs_out, cs_out_anno_logical_full], axis=1)
cs_out["annotation"] = np.where((cs_out["promoter"] == True) & (cs_out["annotation"] != "Exon"), "promoter", cs_out["annotation"])
cs_out["annotation"] = np.where((((cs_out["five_UTR"]==True) | (cs_out["three_UTR"]==True)) & (cs_out["promoter"] != True) & (cs_out["annotation"] != "Exon")), "UTR", cs_out["annotation"])
anno_all_logical = cs_out[
['promoter',
 "five_UTR",
 "three_UTR",
 'intron',
 'exon']]
# change true/false to 1/0
anno_all_logical[["five_UTR", "three_UTR", 'intron', 'exon']] = anno_all_logical[["five_UTR", "three_UTR", 'intron', 'exon']].astype(int)
anno_all_logical.to_csv(r'C:\Users\libin\R_projects\variant_ann\AD_snp.chipseek.logical.ann.input', sep="\t", index=False)
value_counts_df = cs_out['annotation'].value_counts().to_frame().reset_index().rename(columns={"index":"annotation","annotation":"count"})
sns.barplot(x="annotation", y="count", data=value_counts_df, palette="Pastel1")


## processing intersect results
# peak intersection
peak_annotations = []
for cell_type in cell_types:
    atac_peak_intersect = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\AD_{}_peaks.intersect'.format(cell_type), sep="\t", names=\
                                  ["snp_chr","snp_start","snp_end","rsid","snpID","{}_peak_chr".format(cell_type),"{}_peak_start".format(cell_type),"{}_peak_end".format(cell_type),"{}_peak_pvalue".format(cell_type),"{}_peakID".format(cell_type)])

# merge on snp
    atac_peak_intersect_merge_snp = (atac_peak_intersect \
                                     .groupby(['snpID','rsid','snp_chr'], as_index=False) \
                                     .agg(",".join))
    peak_annotations.append(atac_peak_intersect_merge_snp)
# reformat cs_out for merging
    cs_out = cs_out.rename(columns={"seqnames":"snp_chr", "ID":"snpID"}) 
    anno_peak = reduce(lambda left,right: pd.merge(left, right, on=["snp_chr","snpID","rsid"], how="outer").fillna("0"), peak_annotations)
    anno_cs_peak = pd.merge(cs_out, anno_peak,on=["snp_chr","snpID","rsid"], how="outer").fillna("0")
# interaction intersect
interaction_annotations = []
for cell_type in cell_types:
    interactions_intersect_lh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\AD_{}_lh.interaction'.format(cell_type), sep="\t", names=\
                                      ["snp_chr","snp_start","snp_end","rsid","snpID","{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"score","{}_lh_interactionID".format(cell_type)])
    interactions_intersect_rh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\AD_{}_rh.interaction'.format(cell_type), sep="\t", names=\
                                      ["snp_chr","snp_start","snp_end","rsid","snpID","{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"score","{}_rh_interactionID".format(cell_type)])

    # merge on snp
    interactions_intersect_lh = (interactions_intersect_lh \
    .groupby(['snpID','rsid','snp_chr'], as_index=False) \
    .agg(",".join))                           
    interactions_intersect_rh = (interactions_intersect_rh \
    .groupby(['snpID','rsid','snp_chr'], as_index=False) \
    .agg(",".join))   
    interaction_annotations.append(interactions_intersect_lh)
    interaction_annotations.append(interactions_intersect_rh)                             
# merge interaction annotations
anno_interaction = reduce(lambda left,right: pd.merge(left, right, on=["snp_chr","snpID","rsid"], how="outer").fillna("0"), interaction_annotations)
anno_cs_peak_interaction = pd.merge(anno_cs_peak, anno_interaction, on=["snp_chr","snpID","rsid"], how="outer").fillna("0")
# final cleaning & format before downstream filtering analysis








# ----------------------- deprecated codes below ---------------------------
# plt.pie(cs_out['annotation'].value_counts(), labels=cs_out['annotation'].value_counts().index.tolist(),
# labeldistance=0.45, colors=["red","yellow","green","blue"])
# my_circle=plt.Circle((0,0), 0.7, color='white')
# p=plt.gcf()
# p.gca().add_artist(my_circle)
# plt.show()

## prepare homer input -- now use chipseeker instead
# for_homer = sig_snps[['chromosome', 'position']]
# for_homer["end"] = for_homer['position'].apply(lambda x: x+1)
# identifier = ["{}_{}".format(trait,i) for i in range(0, for_homer.shape[0])]
# for_homer.loc[:,"ID"] = identifier
# for_homer.loc[:, "empty"] = "."
# for_homer.loc[: , "strand"] = "+"
# for_homer.to_csv(r"C:\Users\libin\UCSF\FINEMAP\AD_set\AD_snp\snp\AD_all.snp.homer.bed", sep="\t", header=False, index=False)
