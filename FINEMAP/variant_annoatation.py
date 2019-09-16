# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:03:13 2019

@author: bingkunli

Can be run block to block
"""

import sys
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
import itertools
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

trait = "ASD"
cell_types = ["astrocyte", "motor", "cortical", "hippocampal"]

gtf = True
prepare_peak = False
FUMA = False
prepare_motifbreakr_input = False
prepare_promoter = False
annotate_interaction_promoter = False



### ------------------- the following codes are only needed once for each session/disease ----------------------
## gtf preprocessing
# needed for preparing FUMA input
if gtf:
    gtf = pd.read_csv(r"C:\Users\libin\UCSF\gene_annotation\gencode.v19.annotation.gtf", sep="\t", comment="#",
                         names=["chr", "source", "type", "start", "end", "score", "strand", "genomic_phase", "annotation"])
    gtf_genes = gtf[gtf["type"] == "gene"]
    gtf_genes["gene_id"] = gtf_genes["annotation"].str.extract(r'gene_id "(.+?)";')
    gtf_genes["gene_name"] = gtf_genes["annotation"].str.extract(r'gene_name "(.+?)";')
    gtf_genes["gene_type"] = gtf_genes["annotation"].str.extract(r'gene_type "(.+?)";')

if prepare_promoter:
    gtf_genes_pos = gtf_genes[gtf_genes["strand"] == "+"]
    gtf_genes_neg = gtf_genes[gtf_genes["strand"] == "-"]
    gtf_genes_pos_TSS = gtf_genes_pos[["chr","start","gene_id", "gene_name", "gene_type"]]
    gtf_genes_neg_TSS = gtf_genes_neg[["chr","end","gene_id", "gene_name", "gene_type"]]
    
    gtf_genes_pos_TSS["pro_start"] = gtf_genes_pos_TSS["start"].apply(lambda x:x-500)
    gtf_genes_pos_TSS["pro_end"] = gtf_genes_pos_TSS["start"].apply(lambda x:x+500)
    gtf_genes_pos_TSS = gtf_genes_pos_TSS[['chr', 'pro_start', 'pro_end', 'gene_id', 'gene_name', 'gene_type']]
    
    gtf_genes_neg_TSS["pro_start"] = gtf_genes_neg_TSS["end"].apply(lambda x:x-500)
    gtf_genes_neg_TSS["pro_end"] = gtf_genes_neg_TSS["end"].apply(lambda x:x+500)
    gtf_genes_neg_TSS = gtf_genes_neg_TSS[['chr', 'pro_start', 'pro_end', 'gene_id', 'gene_name', 'gene_type']]
   
    gtf_promoter_all = pd.concat([gtf_genes_pos_TSS, gtf_genes_neg_TSS], ignore_index=True)
    gtf_promoter_all.to_csv(r'C:\Users\libin\R_projects\variant_ann\promoters_all.bed', sep="\t", index=False, header=False)
    # next step : intersect promoter with interactions
        
## prepare bedtools intersect input (atac-peak and snps)
if prepare_peak:
    for cell_type in cell_types:
        peak_file = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}.atac.seq.narrowPeak'.format(cell_type), 
                                delim_whitespace=True, 
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

# add in promoter annotation
if annotate_interaction_promoter:
    annotate_interaction_promoter_all = pd.DataFrame()
    for cell_type in cell_types:
        promoter_interaction_lh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_lh.promoter'.format(cell_type), sep="\t",
                                                    names = ["{}_lh_pro_chr".format(cell_type),"{}_lh_pro_start".format(cell_type),"{}_lh_pro_end".format(cell_type),
                                                             "{}_lh_gene_id".format(cell_type),"{}_lh_gene_name".format(cell_type),"{}_lh_gene_type".format(cell_type),
                                                             "{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)])
        promoter_interaction_rh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_rh.promoter'.format(cell_type), sep="\t",
                                                    names = ["{}_rh_pro_chr".format(cell_type),"{}_rh_pro_start".format(cell_type),"{}_rh_pro_end".format(cell_type),
                                                             "{}_rh_gene_id".format(cell_type),"{}_rh_gene_name".format(cell_type),"{}_rh_gene_type".format(cell_type),
                                                             "{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)])
    
        # merge on interaction ID
        promoter_interaction_lh = (promoter_interaction_lh \
        .groupby(["{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)], as_index=False) \
        .agg(",".join))                           
        promoter_interaction_rh = (promoter_interaction_rh \
        .groupby(["{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)], as_index=False) \
        .agg(",".join))

        annotate_interaction_promoter_all = pd.concat([annotate_interaction_promoter_all, promoter_interaction_lh, promoter_interaction_rh], axis=1, ignore_index=False, sort=False).fillna(0)

### --------------------------------- block ends here ----------------------------------


## overall pre-processing
# read in original results of fine-mapping 
header = "index rsid chromosome position allele1 allele2 maf beta se z prob log10bf group corr_group prob_group log10bf_group mean sd mean_incl sd_incl"
header = header.split(" ")
all_snps = pd.read_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\snp\{}_all.snp'.format(trait, trait, trait), delim_whitespace=True, names=header)
sig_snps = all_snps[(all_snps["prob"] > 0.001) & (all_snps["maf"] > 0.001)].drop_duplicates(subset=['rsid'],keep='first')


## prepare chipseeker input -- **NEXT STEP: liftover**
for_cs = sig_snps[['chromosome', 'position', "rsid"]]
# add in finemap-non-eligible snps --- AD_noneli.set
# excluded missed snps stes
notelg = pd.read_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_noteli.set'.format(trait, trait, trait), names=['chromosome', 'position', "rsid"], sep="\t")
for_cs = pd.concat([for_cs, notelg],ignore_index=True).drop_duplicates(subset=['rsid'],keep='first')
for_cs["end"] = for_cs['position'].apply(lambda x: x+1)
identifier = ["{}_{}".format(trait, i) for i in range(0, for_cs.shape[0])]
for_cs.loc[:,"ID"] = identifier
for_cs = for_cs.rename(index=str, columns={"position":"start", "chromosome":"chr"})
for_cs = for_cs[["chr","start","end","rsid","ID"]]
# input for liftover --> chipseeker --> motifbreakR
for_cs.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_all.snps.chipseeker.bed'.format(trait, trait), sep="\t", index=False, header=True)
## next step --liftover ->  chipseeker


### ----------------- hg19 needed for analyses below -----------------------
if prepare_peak:
    # prepare SNP file for intersecting
    hg19_SNP_file = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_all.snps.chipseeker.hg19.bed'.format(trait, trait), sep="\t")
    # hg19_SNP_file = hg19_SNP_file[['seqnames', 'start', 'end','rsid', 'ID']]
    hg19_SNP_file.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_all.snps.bedtools.hg19.bed'.format(trait, trait), sep="\t", index=False, header=False)
## ! next step -- bedtools intersect


### chipseeker output data cleaning
cs_out = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_snp.chipseek.output.bed'.format(trait, trait), sep="\t")
cs_out_anno_logical_full = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_anno_full'.format(trait, trait), 
                                       sep="\t", names=["five_UTR", "three_UTR", "intron", "exon"])
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
cs_out["annotation"] = np.where((((cs_out["five_UTR"]==True) | (cs_out["three_UTR"]==True)) & \
                                  (cs_out["promoter"] != True) & (cs_out["annotation"] != "Exon")), "UTR", cs_out["annotation"])
cs_out["intergenic"] = False
cs_out["intergenic"] = np.where((cs_out['annotation'] == "Intergenic"), True, cs_out["intergenic"])

anno_all_logical = cs_out[
['promoter',
 "five_UTR",
 "three_UTR",
 'intron',
 'exon',
 'intergenic']]
# change true/false to 1/0
anno_all_logical[['promoter', "five_UTR", "three_UTR", 'intron', 'exon', 'intergenic']] = anno_all_logical[['promoter', "five_UTR", "three_UTR", 'intron', 'exon', 'intergenic']].astype(int)
# write out input for plotting upset plot
anno_all_logical.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_snp.chipseek.logical.ann.input'.format(trait, trait), sep="\t", index=False)


# prepare input for FUMA
if FUMA:
    cs_out_sig = pd.merge(cs_out, sig_snps, on=["rsid"], how="inner")
    cs_out_sig = cs_out_sig[cs_out_sig["prob"] > 0.1]
    cs_out_gene = cs_out_sig[cs_out_sig["annotation"] == "Exon"].rename(columns={"geneId":"gene_id"})
    cs_out_gene = pd.merge(cs_out_gene, gtf_genes, on=["gene_id"], how="inner")
    cs_out_gene_protein_coding = cs_out_gene[cs_out_gene['gene_type']=="protein_coding"]
    cs_out_gene_protein_coding_geneName = cs_out_gene_protein_coding[['gene_name']].drop_duplicates(keep='first')
    cs_out_gene_protein_coding_geneName.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_snp.codingGene.forFUMA.bed'.format(trait, trait), sep="\t", index=False, header=False)
    cs_out_gene_all = cs_out_gene[['gene_name']].drop_duplicates(keep='first')
    cs_out_gene_all.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_snp.allGene.forFUMA.bed'.format(trait, trait), sep="\t", index=False, header=False)


### ------------------------------------ motif analysis -------------------------------------
# prepare motifbreakR input
if prepare_motifbreakr_input:
    cs_out_sig = pd.merge(cs_out, sig_snps, on=["rsid"], how="inner")
    cs_out_sig = cs_out_sig[cs_out_sig["prob"] > 0.1]
    cs_out_nonCoding = cs_out_sig[cs_out_sig["annotation"] != "Exon"].rename(columns={"geneId":"gene_id"})
    for_mb = cs_out_nonCoding[['seqnames','start','end',"rsid"]]
    for_mb["real_start"] = for_mb['start'].apply(lambda x:x-1)
    for_mb["score"] = "0"
    for_mb["strand"] = "+"
    for_mb = for_mb[['seqnames', 'real_start', 'start', 'rsid', 'score', 'strand']].sort_values(by=['seqnames', 'start'])
    for_mb.to_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_all.snps.motifbreakR.bed'.format(trait, trait), sep="\t", index=False, header=False)
    for_mb_snps = for_mb[["rsid"]]
    for_mb_snps.to_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_all.snps.motifbreakR.rsid'.format(trait, trait), sep="\t", index=False, header=False)
    ## sampling random non-coding variants as control
    number_of_sig_nonCoding_var = for_mb.shape[0]
    # cs_out_random_backgroud = pd.merge(cs_out, sig_snps, on=["rsid"], how="inner")
    cs_out_random_backgroud = cs_out[cs_out["annotation"] != "Exon"].rename(columns={"geneId":"gene_id"})
    cs_out_random_backgroud = cs_out_random_backgroud.sample(n=number_of_sig_nonCoding_var, random_state=1)
    for_mb_control = cs_out_random_backgroud[['seqnames','start','end',"rsid"]]
    for_mb_control["real_start"] = for_mb_control['start'].apply(lambda x:x-1)
    for_mb_control["score"] = "0"
    for_mb_control["strand"] = "+"
    for_mb_control = for_mb_control[['seqnames', 'real_start', 'start', 'rsid', 'score', 'strand']].sort_values(by=['seqnames', 'start'])
    for_mb_control.to_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_all.snps.motifbreakR.control.bed'.format(trait, trait), sep="\t", index=False, header=False)
    for_mb_snps_control = for_mb_control[["rsid"]]
    for_mb_snps_control.to_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_all.snps.motifbreakR.control.rsid'.format(trait, trait), sep="\t", index=False, header=False)
    ## next step: motifbreakR
    sys.exit()
## process motifbreakR output
mb_out = pd.read_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.anno'.format(trait, trait), sep="\t")
# retain only strong effect TF and remove errornous TF name sush as "T"
mb_out = mb_out[(mb_out["effect"] == "strong") & (mb_out["tf"].str.len() > 1)]
mb_out = mb_out.drop_duplicates()
tf_counts_df = mb_out['tf'].value_counts().to_frame().reset_index().rename(columns={"index":"TF","tf":"count"})

mb_out_control = pd.read_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.control.anno'.format(trait, trait), sep="\t")
mb_out_control = mb_out_control[(mb_out_control["effect"] == "strong") & (mb_out["tf"].str.len() > 1)]
tf_counts_control_df = mb_out_control['tf'].value_counts().to_frame().reset_index().rename(columns={"index":"TF","tf":"count"})

tf_count_compare = pd.merge(tf_counts_df, tf_counts_control_df, on="TF", how="outer").rename(columns={"count_x":"sig","count_y":"control"})
tf_count_compare["diff"] = (tf_count_compare["sig"] - tf_count_compare["control"])/(tf_count_compare["sig"] + tf_count_compare["control"])
# significantly enriched motifs
tf_compare_sig = tf_count_compare[((tf_count_compare["sig"] > tf_count_compare["sig"].median()) & \
                                   (tf_count_compare["diff"] > tf_count_compare["diff"].median())) | \
                                    ((tf_count_compare["sig"] > tf_count_compare["sig"].median()) & (tf_count_compare["control"] == "nan"))]
tf_compare_sig[['sig','control']] = tf_compare_sig[['sig','control']].astype(int)


# enrichr input
tf_compare_sig[["TF"]].to_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.enrichr.bed'.format(trait, trait), sep="\t", index=False, header=False)
tf_compare_sig = tf_compare_sig.rename(columns={"sig":"{}".format(trait)})
tf_compare_sig[["TF","{}".format(trait)]].to_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.for_merging.bed'.format(trait, trait), sep="\t", index=False, header=True)

shared_motif = pd.read_csv(r'C:\Users\libin\R_projects\motifbreakR\shared_motif', sep="\t")
tf_compare_sig_disease_specific = tf_compare_sig[~tf_compare_sig["TF"].isin(shared_motif["TF"])]
tf_compare_sig_disease_specific[["TF"]].to_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.enrichr.specific.bed'.format(trait, trait), sep="\t", index=False, header=False)

# plot enrichment heatmap
plt.figure(figsize=(7,20))
sns.heatmap(tf_compare_sig.set_index("TF")[["{}".format(trait)]],cmap='viridis')
matplotlib.rc('xtick', labelsize=15) 
# plt.savefig(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.heatmap_new.pdf'.format(trait, trait), transparent=True)

plt.figure(figsize=(7,20))
sns.heatmap(tf_compare_sig_disease_specific.set_index("TF")[["{}".format(trait)]],cmap='viridis')
matplotlib.rc('xtick', labelsize=15) 
# plt.savefig(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.heatmap_new_specific.pdf'.format(trait, trait), transparent=True)

tf_compare_sig_disease_specific_part = tf_compare_sig_disease_specific[tf_compare_sig_disease_specific["{}".format(trait)] > 4]
plt.figure(figsize=(12,15))
sns.heatmap(tf_compare_sig_disease_specific_part.set_index("TF")[["{}".format(trait)]],cmap='inferno')
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 
plt.yticks(rotation=0)
plt.savefig(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.heatmap_new_specific_part.pdf'.format(trait, trait), transparent=True)


# SNPs that disrupt enriched motifs
mb_out = mb_out.rename(columns={"tf":'TF'})
# collapse on SNP
mb_out_merged = (mb_out \
                 .groupby(['rsid'], as_index=False) \
                 .agg(",".join))
print("SNPs overlapping TF: ", mb_out_merged.shape[0])
tf_snps = mb_out_merged[["rsid","TF"]]
# collapse on TF
mb_out_merged_TF = (mb_out \
                 .groupby(['TF'], as_index=False) \
                 .agg(",".join))


### !!!!!!!!!!!!!!!!! ------------------------------- motif analyses end here ---------------------------------


## ---------------------------- processing intersecting results -----------------------------
# peak intersection
peak_annotations = []
cs_out_sig = pd.merge(cs_out, sig_snps, on=["rsid"], how="inner")
cs_out_sig = cs_out_sig[cs_out_sig["prob"] > 0.1]

# plot annotation summary
value_counts_df = cs_out_sig['annotation'].value_counts().to_frame().reset_index().rename(columns={"index":"annotation","annotation":"count"})
plt.figure(figsize=(7,5))
sns.barplot(x="annotation", y="count", data=value_counts_df, palette="Pastel1")
plt.savefig(r'C:\Users\libin\R_projects\variant_ann\{}\{}_snp.chipseek.output.anno.pdf'.format(trait, trait), transparent=True)
value_counts_df.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_annotation_summary.csv'.format(trait, trait), sep=",", index=False, header=True)


for cell_type in cell_types:
    atac_peak_intersect = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_{}_peaks.intersect'.format(trait, trait, cell_type), sep="\t", names=\
                                  ["snp_chr","snp_start","snp_end","rsid","snpID",
                                   "{}_peak_chr".format(cell_type),"{}_peak_start".format(cell_type),"{}_peak_end".format(cell_type),"{}_peak_pvalue".format(cell_type),"{}_peakID".format(cell_type)])

# merge on snp
    atac_peak_intersect_merge_snp = (atac_peak_intersect \
                                     .groupby(['snpID','rsid','snp_chr'], as_index=False) \
                                     .agg(",".join))
    peak_annotations.append(atac_peak_intersect_merge_snp)
# reformat cs_out for merging
    
    cs_out_sig = cs_out_sig.rename(columns={"seqnames":"snp_chr", "ID":"snpID"}) 
    anno_peak = reduce(lambda left,right: pd.merge(left, right, on=["snp_chr","snpID","rsid"], how="outer").fillna("0"), peak_annotations)
    anno_cs_peak = pd.merge(cs_out_sig, anno_peak,on=["snp_chr","snpID","rsid"], how="outer").fillna("0")
# interaction intersect
interaction_annotations = []
for cell_type in cell_types:
    interactions_intersect_lh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_{}_lh.interaction'.format(trait, trait, cell_type), sep="\t", names=\
                                      ["snp_chr","snp_start","snp_end","rsid","snpID",
                                       "{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"score","{}_lh_interactionID".format(cell_type)])
    interactions_intersect_rh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_{}_rh.interaction'.format(trait, trait, cell_type), sep="\t", names=\
                                      ["snp_chr","snp_start","snp_end","rsid","snpID",
                                       "{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"score","{}_rh_interactionID".format(cell_type)])

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
# drop variants that were not eligible for fine-mapping
anno_cs_peak_interaction = anno_cs_peak_interaction[anno_cs_peak_interaction["prob"] != "0"]
# merge TF motif annotations
anno_cs_peak_tf_interaction = pd.merge(anno_cs_peak_interaction,tf_snps, on="rsid", how="left").fillna("0")

# final cleaning & format before downstream filtering analysis
anno_cs_peak_tf_interaction = anno_cs_peak_tf_interaction[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'distanceToTSS',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'astrocyte_peakID',
 'motor_peakID',
 'cortical_peakID',
 'hippocampal_peakID',
 'astrocyte_lh_interactionID',
 'astrocyte_rh_interactionID',
 'motor_lh_interactionID',
 'motor_rh_interactionID',
 'cortical_lh_interactionID',
 'cortical_rh_interactionID',
 'hippocampal_lh_interactionID',
 'hippocampal_rh_interactionID']]

anno_cs_peak_tf_interaction[['start','end','distanceToTSS','maf','beta','se','z','prob']] = \
anno_cs_peak_tf_interaction[['start','end','distanceToTSS','maf','beta','se','z','prob']].astype(float)
# write file
anno_cs_peak_tf_interaction.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_anno_cs_peak_tf_interaction'.format(trait, trait, trait), sep="\t", index=False, header=True)
# non-coding variants
anno_cs_peak_tf_interaction_nc = anno_cs_peak_tf_interaction[anno_cs_peak_tf_interaction["annotation"] != "Exon"]

# take the ones fall into chromatin accessible region
anno_cs_peak_tf_interaction_nc_open_all = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["astrocyte_peakID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["motor_peakID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["cortical_peakID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["hippocampal_peakID"] != "0"))]
        
prob_compare = pd.DataFrame()
prob_compare["{}_Accessible".format(trait)] =  anno_cs_peak_tf_interaction_nc_open_all["prob"].to_numpy()
prob_control = anno_cs_peak_tf_interaction_nc[["prob"]].sample(n=anno_cs_peak_tf_interaction_nc_open_all.shape[0], random_state=1)
prob_compare["{}_control".format(trait)] = prob_control["prob"].to_numpy()
prob_compare.to_csv(r'C:\Users\libin\Desktop\thesis\{}_prob_compare'.format(trait), sep="\t", index=False, header=True)

plt.figure()
sns.violinplot(x="variable", y="value", data=prob_compare.melt(), palette="Set2", cut=0)
stats.ttest_ind(prob_compare["{}_Accessible".format(trait)],prob_compare["{}_control".format(trait)])

# TODO： whether variants in open chromatin regions are more enriched in terms of prob
# t-test !!!!!!!!!!!!!!!

# open chromatin & overlap at least one interacting bins
anno_cs_peak_tf_interaction_nc_open_interaction_all = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["astrocyte_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["astrocyte_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["astrocyte_rh_interactionID"] != "0"))) \
        | \
        ((anno_cs_peak_tf_interaction_nc["motor_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["motor_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["motor_rh_interactionID"] != "0"))) \
        | \
        ((anno_cs_peak_tf_interaction_nc["cortical_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["cortical_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["cortical_rh_interactionID"] != "0"))) \
        | \
        ((anno_cs_peak_tf_interaction_nc["hippocampal_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["hippocampal_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["hippocampal_rh_interactionID"] != "0")))]
      
 
### ----------------------- cell-type_specific analyses ------------------------------     
# --- excitory neurons ---
cell_type = "cortical"
# process promoter annotation first
promoter_interaction_lh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_lh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_lh_pro_chr".format(cell_type),"{}_lh_pro_start".format(cell_type),"{}_lh_pro_end".format(cell_type),
                                                     "{}_lh_gene_id".format(cell_type),"{}_lh_gene_name".format(cell_type),"{}_lh_gene_type".format(cell_type),
                                                     "{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)])
promoter_interaction_rh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_rh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_rh_pro_chr".format(cell_type),"{}_rh_pro_start".format(cell_type),"{}_rh_pro_end".format(cell_type),
                                                     "{}_rh_gene_id".format(cell_type),"{}_rh_gene_name".format(cell_type),"{}_rh_gene_type".format(cell_type),
                                                     "{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)])
# merge on interaction ID
promoter_interaction_lh = (promoter_interaction_lh \
.groupby(["{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))                           
promoter_interaction_rh = (promoter_interaction_rh \
.groupby(["{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))

anno_cs_peak_tf_interaction_nc_open_interaction_cortical = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["cortical_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["cortical_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["cortical_rh_interactionID"] != "0")))] \
        .drop(['astrocyte_peakID', 'motor_peakID','hippocampal_peakID',
               'astrocyte_lh_interactionID','astrocyte_rh_interactionID',
               'motor_lh_interactionID','motor_rh_interactionID',
               'hippocampal_lh_interactionID','hippocampal_rh_interactionID'],axis=1)
# expand on lh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_cortical = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_cortical.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_cortical.columns.drop("cortical_lh_interactionID", 1).tolist())
    ["cortical_lh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"cortical_lh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_cortical.columns]
    )
# expand on rh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_cortical = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_cortical.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_cortical.columns.drop("cortical_rh_interactionID", 1).tolist())
    ["cortical_rh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"cortical_rh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_cortical.columns]
    )
# integrate promoter annotation
anno_cs_peak_tf_interaction_nc_open_interaction_cortical = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_cortical, promoter_interaction_lh,
                                                                    on = "cortical_lh_interactionID" , how="left").fillna("0")
anno_cs_peak_tf_interaction_nc_open_interaction_cortical = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_cortical, promoter_interaction_rh,
                                                                    on = "cortical_rh_interactionID", how="left").fillna("0")
# SNPs overlapping with lh bins
# !!!!!!!!!!!!!! "cortical_lh_start" = 0 means no promoter overlap this interacting bin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
anno_cs_peak_tf_interaction_nc_open_interaction_cortical_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_cortical[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'cortical_peakID',
 'cortical_lh_interactionID',
 'cortical_lh_start',
 'cortical_lh_end']].rename(columns={"cortical_lh_interactionID" : "cortical_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_cortical_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_cortical_lh[anno_cs_peak_tf_interaction_nc_open_interaction_cortical_lh["cortical_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_cortical_lh[['start','end','cortical_lh_start','cortical_lh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_cortical_lh[['start','end','cortical_lh_start','cortical_lh_end']].astype(int)

cortical_promoter_interaction_rh = promoter_interaction_rh.rename(columns={"cortical_rh_interactionID" : "cortical_interactionID"})[["cortical_interactionID",'cortical_rh_chr','cortical_rh_start','cortical_rh_end','cortical_score', 'cortical_rh_gene_id','cortical_rh_gene_name','cortical_rh_gene_type']]

cortical_lh_snp_rh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_cortical_lh, cortical_promoter_interaction_rh, how="left", on=['cortical_interactionID']).fillna("0").drop_duplicates()
cortical_lh_snp_rh_interaction[['cortical_rh_start','cortical_rh_end']] = cortical_lh_snp_rh_interaction[['cortical_rh_start','cortical_rh_end']].astype(int)
cortical_lh_snp_rh_interaction_promoter = cortical_lh_snp_rh_interaction[cortical_lh_snp_rh_interaction["cortical_rh_start"] != 0] 

anno_cs_peak_tf_interaction_nc_open_interaction_cortical_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_cortical[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'cortical_peakID',
 'cortical_rh_interactionID',
 'cortical_rh_start',
 'cortical_rh_end']].rename(columns={"cortical_rh_interactionID" : "cortical_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_cortical_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_cortical_rh[anno_cs_peak_tf_interaction_nc_open_interaction_cortical_rh["cortical_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_cortical_rh[['start','end','cortical_rh_start','cortical_rh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_cortical_rh[['start','end','cortical_rh_start','cortical_rh_end']].astype(int)

cortical_promoter_interaction_lh = promoter_interaction_lh.rename(columns={"cortical_lh_interactionID" : "cortical_interactionID"})[["cortical_interactionID",'cortical_lh_chr','cortical_lh_start','cortical_lh_end','cortical_score', 'cortical_lh_gene_id','cortical_lh_gene_name','cortical_lh_gene_type']]

cortical_rh_snp_lh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_cortical_rh, cortical_promoter_interaction_lh, how="left", on=['cortical_interactionID']).fillna("0").drop_duplicates()
cortical_rh_snp_lh_interaction[['cortical_lh_start','cortical_lh_end']] = cortical_rh_snp_lh_interaction[['cortical_lh_start','cortical_lh_end']].astype(int)
# SNPs overlapping left hand bin of interactions whose right hand bins overlap with enhancers :)
cortical_rh_snp_lh_interaction_promoter = cortical_rh_snp_lh_interaction[cortical_rh_snp_lh_interaction["cortical_lh_start"] != 0] 
# write file 
cortical_rh_snp_lh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_cortical_rh_snp_lh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)
cortical_lh_snp_rh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_cortical_lh_snp_rh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)

cortical_snp_interaction_promoter_all = pd.concat([cortical_rh_snp_lh_interaction_promoter,cortical_lh_snp_rh_interaction_promoter], axis=0, ignore_index=True)
# most credible set of SNPs 
cortical_credible_variants_final = cortical_snp_interaction_promoter_all[['rsid']]
# GO anaysis round two
cortical_target_genes_final_list = cortical_snp_interaction_promoter_all["cortical_lh_gene_name"].tolist() + cortical_snp_interaction_promoter_all["cortical_rh_gene_name"].tolist()
# remove nan
cortical_target_genes_final_list = [str(i).split(",") for i in cortical_target_genes_final_list if isinstance(i, str)]
# flattern list & remove duplicates
cortical_target_genes_final_list = list(set(list(itertools.chain.from_iterable(cortical_target_genes_final_list))))

anno_cs_peak_tf_interaction_nc_open_interaction_cortical_only = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["cortical_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["cortical_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["cortical_rh_interactionID"] != "0")) & \
        (anno_cs_peak_tf_interaction_nc["astrocyte_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["motor_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["hippocampal_peakID"] == "0"))]

    
# --- astrocyte ---
cell_type = "astrocyte"
promoter_interaction_lh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_lh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_lh_pro_chr".format(cell_type),"{}_lh_pro_start".format(cell_type),"{}_lh_pro_end".format(cell_type),
                                                     "{}_lh_gene_id".format(cell_type),"{}_lh_gene_name".format(cell_type),"{}_lh_gene_type".format(cell_type),
                                                     "{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)])
promoter_interaction_rh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_rh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_rh_pro_chr".format(cell_type),"{}_rh_pro_start".format(cell_type),"{}_rh_pro_end".format(cell_type),
                                                     "{}_rh_gene_id".format(cell_type),"{}_rh_gene_name".format(cell_type),"{}_rh_gene_type".format(cell_type),
                                                     "{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)])
# merge on interaction ID
promoter_interaction_lh = (promoter_interaction_lh \
.groupby(["{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))                           
promoter_interaction_rh = (promoter_interaction_rh \
.groupby(["{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))

anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["astrocyte_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["astrocyte_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["astrocyte_rh_interactionID"] != "0")))] \
        .drop(['cortical_peakID', 'motor_peakID','hippocampal_peakID',
               'cortical_lh_interactionID','cortical_rh_interactionID',
               'motor_lh_interactionID','motor_rh_interactionID',
               'hippocampal_lh_interactionID','hippocampal_rh_interactionID'],axis=1)
# expand on lh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte.columns.drop("astrocyte_lh_interactionID", 1).tolist())
    ["astrocyte_lh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"astrocyte_lh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte.columns]
    )
# expand on rh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte.columns.drop("astrocyte_rh_interactionID", 1).tolist())
    ["astrocyte_rh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"astrocyte_rh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte.columns]
    )
# integrate promoter annotation
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte, promoter_interaction_lh,
                                                                    on = "astrocyte_lh_interactionID" , how="left").fillna("0")
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte, promoter_interaction_rh,
                                                                    on = "astrocyte_rh_interactionID", how="left").fillna("0")
# SNPs overlapping with lh bins
# !!!!!!!!!!!!!! "astrocyte_lh_start" = 0 means no promoter overlap this interacting bin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'astrocyte_peakID',
 'astrocyte_lh_interactionID',
 'astrocyte_lh_start',
 'astrocyte_lh_end']].rename(columns={"astrocyte_lh_interactionID" : "astrocyte_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_lh[anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_lh["astrocyte_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_lh[['start','end','astrocyte_lh_start','astrocyte_lh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_lh[['start','end','astrocyte_lh_start','astrocyte_lh_end']].astype(int)

astrocyte_promoter_interaction_rh = promoter_interaction_rh.rename(columns={"astrocyte_rh_interactionID" : "astrocyte_interactionID"})[["astrocyte_interactionID",'astrocyte_rh_chr','astrocyte_rh_start','astrocyte_rh_end','astrocyte_score', 'astrocyte_rh_gene_id','astrocyte_rh_gene_name','astrocyte_rh_gene_type']]

astrocyte_lh_snp_rh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_lh, astrocyte_promoter_interaction_rh, how="left", on=['astrocyte_interactionID']).fillna("0").drop_duplicates()
astrocyte_lh_snp_rh_interaction[['astrocyte_rh_start','astrocyte_rh_end']] = astrocyte_lh_snp_rh_interaction[['astrocyte_rh_start','astrocyte_rh_end']].astype(int)
astrocyte_lh_snp_rh_interaction_promoter = astrocyte_lh_snp_rh_interaction[astrocyte_lh_snp_rh_interaction["astrocyte_rh_start"] != 0] 

anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'astrocyte_peakID',
 'astrocyte_rh_interactionID',
 'astrocyte_rh_start',
 'astrocyte_rh_end']].rename(columns={"astrocyte_rh_interactionID" : "astrocyte_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_rh[anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_rh["astrocyte_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_rh[['start','end','astrocyte_rh_start','astrocyte_rh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_rh[['start','end','astrocyte_rh_start','astrocyte_rh_end']].astype(int)

astrocyte_promoter_interaction_lh = promoter_interaction_lh.rename(columns={"astrocyte_lh_interactionID" : "astrocyte_interactionID"})[["astrocyte_interactionID",'astrocyte_lh_chr','astrocyte_lh_start','astrocyte_lh_end','astrocyte_score', 'astrocyte_lh_gene_id','astrocyte_lh_gene_name','astrocyte_lh_gene_type']]

astrocyte_rh_snp_lh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_rh, astrocyte_promoter_interaction_lh, how="left", on=['astrocyte_interactionID']).fillna("0").drop_duplicates()
astrocyte_rh_snp_lh_interaction[['astrocyte_lh_start','astrocyte_lh_end']] = astrocyte_rh_snp_lh_interaction[['astrocyte_lh_start','astrocyte_lh_end']].astype(int)
# SNPs overlapping left hand bin of interactions whose right hand bins overlap with enhancers :)
astrocyte_rh_snp_lh_interaction_promoter = astrocyte_rh_snp_lh_interaction[astrocyte_rh_snp_lh_interaction["astrocyte_lh_start"] != 0] 
# write file 
astrocyte_rh_snp_lh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_astrocyte_rh_snp_lh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)
astrocyte_lh_snp_rh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_astrocyte_lh_snp_rh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)

astrocyte_snp_interaction_promoter_all = pd.concat([astrocyte_rh_snp_lh_interaction_promoter,astrocyte_lh_snp_rh_interaction_promoter], axis=0, ignore_index=True)
# most credible set of SNPs 
astrocyte_credible_variants_final = astrocyte_snp_interaction_promoter_all[['rsid']]
# GO anaysis round two
astrocyte_target_genes_final_list = astrocyte_snp_interaction_promoter_all["astrocyte_lh_gene_name"].tolist() + astrocyte_snp_interaction_promoter_all["astrocyte_rh_gene_name"].tolist()
# remove nan
astrocyte_target_genes_final_list = [str(i).split(",") for i in astrocyte_target_genes_final_list if isinstance(i, str)]
# flattern list & remove duplicates
astrocyte_target_genes_final_list = list(set(list(itertools.chain.from_iterable(astrocyte_target_genes_final_list))))

anno_cs_peak_tf_interaction_nc_open_interaction_astrocyte_only = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["astrocyte_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["astrocyte_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["astrocyte_rh_interactionID"] != "0")) & \
        (anno_cs_peak_tf_interaction_nc["cortical_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["motor_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["hippocampal_peakID"] == "0"))]


# --- hippocampal ---
cell_type = "hippocampal"
# process promoter annotation first
promoter_interaction_lh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_lh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_lh_pro_chr".format(cell_type),"{}_lh_pro_start".format(cell_type),"{}_lh_pro_end".format(cell_type),
                                                     "{}_lh_gene_id".format(cell_type),"{}_lh_gene_name".format(cell_type),"{}_lh_gene_type".format(cell_type),
                                                     "{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)])
promoter_interaction_rh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_rh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_rh_pro_chr".format(cell_type),"{}_rh_pro_start".format(cell_type),"{}_rh_pro_end".format(cell_type),
                                                     "{}_rh_gene_id".format(cell_type),"{}_rh_gene_name".format(cell_type),"{}_rh_gene_type".format(cell_type),
                                                     "{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)])
# merge on interaction ID
promoter_interaction_lh = (promoter_interaction_lh \
.groupby(["{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))                           
promoter_interaction_rh = (promoter_interaction_rh \
.groupby(["{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))

anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["hippocampal_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["hippocampal_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["hippocampal_rh_interactionID"] != "0")))] \
        .drop(['cortical_peakID', 'motor_peakID','astrocyte_peakID',
               'cortical_lh_interactionID','cortical_rh_interactionID',
               'motor_lh_interactionID','motor_rh_interactionID',
               'astrocyte_lh_interactionID','astrocyte_rh_interactionID'],axis=1)

# expand on lh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal.columns.drop("hippocampal_lh_interactionID", 1).tolist())
    ["hippocampal_lh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"hippocampal_lh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal.columns]
    )
# expand on rh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal.columns.drop("hippocampal_rh_interactionID", 1).tolist())
    ["hippocampal_rh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"hippocampal_rh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal.columns]
    )
# integrate promoter annotation
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal, promoter_interaction_lh,
                                                                    on = "hippocampal_lh_interactionID" , how="left").fillna("0")
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal, promoter_interaction_rh,
                                                                    on = "hippocampal_rh_interactionID", how="left").fillna("0")
# SNPs overlapping with lh bins
# !!!!!!!!!!!!!! "hippocampal_lh_start" = 0 means no promoter overlap this interacting bin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'hippocampal_peakID',
 'hippocampal_lh_interactionID',
 'hippocampal_lh_start',
 'hippocampal_lh_end']].rename(columns={"hippocampal_lh_interactionID" : "hippocampal_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_lh[anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_lh["hippocampal_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_lh[['start','end','hippocampal_lh_start','hippocampal_lh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_lh[['start','end','hippocampal_lh_start','hippocampal_lh_end']].astype(int)

hippocampal_promoter_interaction_rh = promoter_interaction_rh.rename(columns={"hippocampal_rh_interactionID" : "hippocampal_interactionID"})[["hippocampal_interactionID",'hippocampal_rh_chr','hippocampal_rh_start','hippocampal_rh_end','hippocampal_score', 'hippocampal_rh_gene_id','hippocampal_rh_gene_name','hippocampal_rh_gene_type']]

hippocampal_lh_snp_rh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_lh, hippocampal_promoter_interaction_rh, how="left", on=['hippocampal_interactionID']).fillna("0").drop_duplicates()
hippocampal_lh_snp_rh_interaction[['hippocampal_rh_start','hippocampal_rh_end']] = hippocampal_lh_snp_rh_interaction[['hippocampal_rh_start','hippocampal_rh_end']].astype(int)
hippocampal_lh_snp_rh_interaction_promoter = hippocampal_lh_snp_rh_interaction[hippocampal_lh_snp_rh_interaction["hippocampal_rh_start"] != 0] 

anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'hippocampal_peakID',
 'hippocampal_rh_interactionID',
 'hippocampal_rh_start',
 'hippocampal_rh_end']].rename(columns={"hippocampal_rh_interactionID" : "hippocampal_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_rh[anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_rh["hippocampal_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_rh[['start','end','hippocampal_rh_start','hippocampal_rh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_rh[['start','end','hippocampal_rh_start','hippocampal_rh_end']].astype(int)

hippocampal_promoter_interaction_lh = promoter_interaction_lh.rename(columns={"hippocampal_lh_interactionID" : "hippocampal_interactionID"})[["hippocampal_interactionID",'hippocampal_lh_chr','hippocampal_lh_start','hippocampal_lh_end','hippocampal_score', 'hippocampal_lh_gene_id','hippocampal_lh_gene_name','hippocampal_lh_gene_type']]

hippocampal_rh_snp_lh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_rh, hippocampal_promoter_interaction_lh, how="left", on=['hippocampal_interactionID']).fillna("0").drop_duplicates()
hippocampal_rh_snp_lh_interaction[['hippocampal_lh_start','hippocampal_lh_end']] = hippocampal_rh_snp_lh_interaction[['hippocampal_lh_start','hippocampal_lh_end']].astype(int)
# SNPs overlapping left hand bin of interactions whose right hand bins overlap with enhancers :)
hippocampal_rh_snp_lh_interaction_promoter = hippocampal_rh_snp_lh_interaction[hippocampal_rh_snp_lh_interaction["hippocampal_lh_start"] != 0] 
# write file 
hippocampal_rh_snp_lh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_hippocampal_rh_snp_lh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)
hippocampal_lh_snp_rh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_hippocampal_lh_snp_rh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)

hippocampal_snp_interaction_promoter_all = pd.concat([hippocampal_rh_snp_lh_interaction_promoter,hippocampal_lh_snp_rh_interaction_promoter], axis=0, ignore_index=True)
# most credible set of SNPs 
hippocampal_credible_variants_final = hippocampal_snp_interaction_promoter_all[['rsid']]
# GO anaysis round two
hippocampal_target_genes_final_list = hippocampal_snp_interaction_promoter_all["hippocampal_lh_gene_name"].tolist() + hippocampal_snp_interaction_promoter_all["hippocampal_rh_gene_name"].tolist()
# remove nan
hippocampal_target_genes_final_list = [str(i).split(",") for i in hippocampal_target_genes_final_list if isinstance(i, str)]
# flattern list & remove duplicates
hippocampal_target_genes_final_list = list(set(list(itertools.chain.from_iterable(hippocampal_target_genes_final_list))))

nno_cs_peak_tf_interaction_nc_open_interaction_hippocampal_only = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["hippocampal_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["hippocampal_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["hippocampal_rh_interactionID"] != "0")) & \
        (anno_cs_peak_tf_interaction_nc["cortical_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["motor_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["astrocyte_peakID"] == "0"))]
# write file


cell_type = "motor"
# process promoter annotation first
promoter_interaction_lh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_lh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_lh_pro_chr".format(cell_type),"{}_lh_pro_start".format(cell_type),"{}_lh_pro_end".format(cell_type),
                                                     "{}_lh_gene_id".format(cell_type),"{}_lh_gene_name".format(cell_type),"{}_lh_gene_type".format(cell_type),
                                                     "{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)])
promoter_interaction_rh = pd.read_csv(r'C:\Users\libin\R_projects\variant_ann\{}_rh.promoter'.format(cell_type), sep="\t",
                                            names = ["{}_rh_pro_chr".format(cell_type),"{}_rh_pro_start".format(cell_type),"{}_rh_pro_end".format(cell_type),
                                                     "{}_rh_gene_id".format(cell_type),"{}_rh_gene_name".format(cell_type),"{}_rh_gene_type".format(cell_type),
                                                     "{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)])
# merge on interaction ID
promoter_interaction_lh = (promoter_interaction_lh \
.groupby(["{}_lh_chr".format(cell_type),"{}_lh_start".format(cell_type),"{}_lh_end".format(cell_type),"{}_score".format(cell_type),"{}_lh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))                           
promoter_interaction_rh = (promoter_interaction_rh \
.groupby(["{}_rh_chr".format(cell_type),"{}_rh_start".format(cell_type),"{}_rh_end".format(cell_type),"{}_score".format(cell_type),"{}_rh_interactionID".format(cell_type)], as_index=False) \
.agg(",".join))

anno_cs_peak_tf_interaction_nc_open_interaction_motor = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["motor_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["motor_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["motor_rh_interactionID"] != "0")))] \
        .drop(['cortical_peakID', 'hippocampal_peakID','astrocyte_peakID',
               'cortical_lh_interactionID','cortical_rh_interactionID',
               'hippocampal_lh_interactionID','hippocampal_lh_interactionID',
               'astrocyte_lh_interactionID','astrocyte_rh_interactionID'],axis=1)

# expand on lh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_motor = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_motor.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_motor.columns.drop("motor_lh_interactionID", 1).tolist())
    ["motor_lh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"motor_lh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_motor.columns]
    )

# expand on rh interaction bin
anno_cs_peak_tf_interaction_nc_open_interaction_motor = \
    (anno_cs_peak_tf_interaction_nc_open_interaction_motor.set_index(
    anno_cs_peak_tf_interaction_nc_open_interaction_motor.columns.drop("motor_rh_interactionID", 1).tolist())
    ["motor_rh_interactionID"].str.split(",", expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:"motor_rh_interactionID"}) 
    .loc[:, anno_cs_peak_tf_interaction_nc_open_interaction_motor.columns]
    )
# integrate promoter annotation
anno_cs_peak_tf_interaction_nc_open_interaction_motor = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_motor, promoter_interaction_lh,
                                                                    on = "motor_lh_interactionID" , how="left").fillna("0")
anno_cs_peak_tf_interaction_nc_open_interaction_motor = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_motor, promoter_interaction_rh,
                                                                    on = "motor_rh_interactionID", how="left").fillna("0")
# SNPs overlapping with lh bins
# !!!!!!!!!!!!!! "motor_lh_start" = 0 means no promoter overlap this interacting bin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
anno_cs_peak_tf_interaction_nc_open_interaction_motor_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_motor[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'motor_peakID',
 'motor_lh_interactionID',
 'motor_lh_start',
 'motor_lh_end']].rename(columns={"motor_lh_interactionID" : "motor_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_motor_lh = \
anno_cs_peak_tf_interaction_nc_open_interaction_motor_lh[anno_cs_peak_tf_interaction_nc_open_interaction_motor_lh["motor_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_motor_lh[['start','end','motor_lh_start','motor_lh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_motor_lh[['start','end','motor_lh_start','motor_lh_end']].astype(int)

motor_promoter_interaction_rh = promoter_interaction_rh.rename(columns={"motor_rh_interactionID" : "motor_interactionID"})[["motor_interactionID",'motor_rh_chr','motor_rh_start','motor_rh_end','motor_score', 'motor_rh_gene_id','motor_rh_gene_name','motor_rh_gene_type']]

motor_lh_snp_rh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_motor_lh, motor_promoter_interaction_rh, how="left", on=['motor_interactionID']).fillna("0").drop_duplicates()
motor_lh_snp_rh_interaction[['motor_rh_start','motor_rh_end']] = motor_lh_snp_rh_interaction[['motor_rh_start','motor_rh_end']].astype(int)
motor_lh_snp_rh_interaction_promoter = motor_lh_snp_rh_interaction[motor_lh_snp_rh_interaction["motor_rh_start"] != 0] 

anno_cs_peak_tf_interaction_nc_open_interaction_motor_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_motor[
['snp_chr',
 'start',
 'end',
 'rsid',
 'snpID',
 'annotation',
 'geneId',
 'promoter',
 'five_UTR',
 'three_UTR',
 'intron',
 'exon',
 'maf',
 'beta',
 'se',
 'z',
 'prob',
 'TF',
 'motor_peakID',
 'motor_rh_interactionID',
 'motor_rh_start',
 'motor_rh_end']].rename(columns={"motor_rh_interactionID" : "motor_interactionID"})

anno_cs_peak_tf_interaction_nc_open_interaction_motor_rh = \
anno_cs_peak_tf_interaction_nc_open_interaction_motor_rh[anno_cs_peak_tf_interaction_nc_open_interaction_motor_rh["motor_interactionID"] != "0"]
anno_cs_peak_tf_interaction_nc_open_interaction_motor_rh[['start','end','motor_rh_start','motor_rh_end']] = anno_cs_peak_tf_interaction_nc_open_interaction_motor_rh[['start','end','motor_rh_start','motor_rh_end']].astype(int)

motor_promoter_interaction_lh = promoter_interaction_lh.rename(columns={"motor_lh_interactionID" : "motor_interactionID"})[["motor_interactionID",'motor_lh_chr','motor_lh_start','motor_lh_end','motor_score', 'motor_lh_gene_id','motor_lh_gene_name','motor_lh_gene_type']]

motor_rh_snp_lh_interaction = pd.merge(anno_cs_peak_tf_interaction_nc_open_interaction_motor_rh, motor_promoter_interaction_lh, how="left", on=['motor_interactionID']).fillna("0").drop_duplicates()
motor_rh_snp_lh_interaction[['motor_lh_start','motor_lh_end']] = motor_rh_snp_lh_interaction[['motor_lh_start','motor_lh_end']].astype(int)
# SNPs overlapping left hand bin of interactions whose right hand bins overlap with enhancers :)
motor_rh_snp_lh_interaction_promoter = motor_rh_snp_lh_interaction[motor_rh_snp_lh_interaction["motor_lh_start"] != 0] 
# write file 
motor_rh_snp_lh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_motor_rh_snp_lh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)
motor_lh_snp_rh_interaction_promoter.to_csv(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\{}_motor_lh_snp_rh_interaction_promoter'.format(trait, trait, trait), sep="\t", index=False, header=True)

motor_snp_interaction_promoter_all = pd.concat([motor_rh_snp_lh_interaction_promoter,motor_lh_snp_rh_interaction_promoter], axis=0, ignore_index=True)
# most credible set of SNPs 
motor_credible_variants_final = motor_snp_interaction_promoter_all[['rsid']]
# GO anaysis round two
motor_target_genes_final_list = motor_snp_interaction_promoter_all["motor_lh_gene_name"].tolist() + motor_snp_interaction_promoter_all["motor_rh_gene_name"].tolist()
# remove nan
motor_target_genes_final_list = [str(i).split(",") for i in motor_target_genes_final_list if isinstance(i, str)]
# flattern list & remove duplicates
motor_target_genes_final_list = list(set(list(itertools.chain.from_iterable(motor_target_genes_final_list))))

anno_cs_peak_tf_interaction_nc_open_interaction_motor_only = anno_cs_peak_tf_interaction_nc[
        ((anno_cs_peak_tf_interaction_nc["motor_peakID"] != "0") & \
        ((anno_cs_peak_tf_interaction_nc["motor_lh_interactionID"] != "0") | \
        (anno_cs_peak_tf_interaction_nc["motor_rh_interactionID"] != "0")) & \
        (anno_cs_peak_tf_interaction_nc["cortical_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["hippocampal_peakID"] == "0") & \
        (anno_cs_peak_tf_interaction_nc["astrocyte_peakID"] == "0"))]     
# TODO: write file

# combine four cell types
# SNPs that 1) overlap open chromatin region, 2) overlap interacting bins,
# 3) interacting bins overlap with promoters
credible_variant_set_combined = pd.concat([motor_credible_variants_final, cortical_credible_variants_final,
                                astrocyte_credible_variants_final, hippocampal_credible_variants_final]).drop_duplicates(keep="first").reset_index().drop(["index"], axis=1)
putative_target_genes_combined = list(set(motor_target_genes_final_list + astrocyte_target_genes_final_list + \
                        hippocampal_target_genes_final_list + cortical_target_genes_final_list))
putative_target_genes_combined_df = pd.DataFrame()
putative_target_genes_combined_df["genes"] = putative_target_genes_combined
putative_target_genes_combined_df.to_csv(r'C:\Users\libin\R_projects\variant_ann\{}\{}_target_genes.set'.format(trait, trait), sep="\t", index=False, header=False)

print(trait, "fine-mapped", anno_cs_peak_tf_interaction.shape[0])
print(trait, "non_coding", anno_cs_peak_tf_interaction_nc.shape[0])
print(trait, "chromatin_accessible", anno_cs_peak_tf_interaction_nc_open_all.shape[0])
print(trait, "final credible set: ", credible_variant_set_combined.shape[0])

# TODO : whether credible sets are more enriched in terms of prob









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

# tf_counts_df_sig = tf_counts_df[tf_counts_df["count"] > tf_counts_df["count"].mean()]
# tf_counts_control_df_sig = tf_counts_control_df[tf_counts_control_df["count"] > tf_counts_control_df["count"].mean()]
# tf_count_compare['diff'] = tf_count_compare['diff'].apply("int64")

# tf_count_compare_sig = pd.merge(tf_counts_df_sig, tf_counts_control_df_sig, on="TF", how="outer").rename(columns={"count_x":"sig","count_y":"control"})
# tf_count_compare_sig["diff"] = (tf_count_compare_sig["sig"] - tf_count_compare_sig["control"])/(tf_count_compare_sig["sig"] + tf_count_compare_sig["control"])

# tf_compare_sig['sig'] = tf_compare_sig['sig'].apply("int64")
# tf_compare_sig['control'] = tf_compare_sig['control'].apply("int64")

# anno_cs_peak_tf_interaction_nc_open_interaction = anno_cs_peak_tf_interaction_nc_open[
#        (anno_cs_peak_tf_interaction_nc_open["astrocyte_lh_interactionID"] != "0") | \
#        (anno_cs_peak_tf_interaction_nc_open["astrocyte_rh_interactionID"] != "0") | \
#        (anno_cs_peak_tf_interaction_nc_open["motor_lh_interactionID"] != "0") | \
#        (anno_cs_peak_tf_interaction_nc_open["motor_rh_interactionID"] != "0") | \
#        (anno_cs_peak_tf_interaction_nc_open["cortical_lh_interactionID"] != "0") | \
#        (anno_cs_peak_tf_interaction_nc_open["cortical_rh_interactionID"] != "0") | \
#        (anno_cs_peak_tf_interaction_nc_open["hippocampal_lh_interactionID"] != "0") | \
#        (anno_cs_peak_tf_interaction_nc_open["hippocampal_rh_interactionID"] != "0")]