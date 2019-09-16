# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:07:16 2019

@author: libin
"""

from pathlib import Path
import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from functools import reduce
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

traits=["SCZ", "ASD","AD", "ADHD", "BP"]


for_plot = pd.DataFrame()
estimated_LD = pd.read_csv(r'C:\Users\libin\Desktop\thesis\nygcresearch-ldetect-data-ac125e47bf7f\nygcresearch-ldetect-data-ac125e47bf7f\EUR\fourier_ls-all.bed', sep="\t")
estimated_LD["size"] = (estimated_LD[' stop'] - estimated_LD[' start ']).apply("int64")
sns.violinplot(estimated_LD['size'])

for_plot["estimated"] = estimated_LD["size"]

# plotting basic statistics
for trait in traits:
    size = pd.read_csv(r'C:\Users\libin\Desktop\thesis\{}_count_size'.format(trait), sep="\t", names=["{}".format(trait)]).dropna()
    size["{}".format(trait)] = size["{}".format(trait)].apply("int64")
    # sns.violinplot(size)
    for_plot["{}".format(trait)] = size["{}".format(trait)]
    # for_plot["{}".format(trait)] = for_plot["{}".format(trait)].apply("int64")
   
plt.figure(figsize=(15,10))
sns.violinplot(x="variable", y="value", data=for_plot.melt(), palette="Set2", cut=0)
matplotlib.rc('xtick', labelsize=40) 
    

# plot shared motif
motif_df_list = []
for trait in traits:
    disease_motif = pd.read_csv(r'C:\Users\libin\R_projects\motifbreakR\{}\{}_motif.for_merging.bed'.format(trait,trait), sep="\t")
    motif_df_list.append(disease_motif)
motif_df_join = anno_peak = reduce(lambda left,right: pd.merge(left, right, on=["TF"], how="inner"), motif_df_list)    

plt.figure(figsize=(10,10))
sns.heatmap(motif_df_join.set_index("TF"),cmap='viridis')
matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 
plt.yticks(rotation=0)
plt.savefig(r'C:\Users\libin\Desktop\thesis\motif.heatmap_shared.pdf', transparent=True)
motif_df_join[["TF"]].to_csv(r'C:\Users\libin\R_projects\motifbreakR\shared_motif', sep="\t", index=False, header=True)
    
# plot prob compare
prob_compare_df_all = pd.DataFrame()
for trait in traits:
    prob_compare = pd.read_csv(r'C:\Users\libin\Desktop\thesis\{}_prob_compare'.format(trait), sep="\t")
    prob_compare_df_all = prob_compare_df_all.append(prob_compare, ignore_index=True)
    
plt.figure(figsize=(21,5))
sns.violinplot(x="variable", y="value", data=prob_compare_df_all.melt().dropna(), 
               palette=["#ef475d", "#ef475d",  "#1a94bc", "#1a94bc","#f9d770", "#f9d770", "#29b7cb", "#29b7cb", "#fca106", "#fca106", ], cut=0).set(ylim=(0.1,1.01))
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

# for trait in traits:
#    size_prob = pd.read_csv(r'C:\Users\libin\Desktop\thesis\{}_prob_stat'.format(trait), sep="\t", names=["SNP_number", "prob"])
#    sns.scatterplot(x="SNP_number", y="prob", data=size_prob, color="blue")











### ---------------- deprecated code below ---------------------
# fig.1c: # of regions & posterior bins
# can't really make this plot because SNPs overlap among regions
# snp_path = Path(r'C:\Users\libin\UCSF\FINEMAP\{}_set\{}_snp\snp'.format(trait, trait))
# snp_file_path_list = [str(fp) for fp in snp_path.glob('**/*_block_*.snp')]
# snp_file_name_list = ["".join(re.findall(r'(.+).snp',fp.name)) for fp in snp_path.glob('**/*_block_*.snp')]
# list to store all the snp dataframes for integrating later
# all_snps = []
# for file in snp_file_path_list:
#    filename = snp_file_name_list[snp_file_path_list.index(file)]
#    snp_df = pd.read_csv(file, delim_whitespace=True)
#    snp_df = snp_df[["rsid", 'chromosome', 'position', 'maf', 'beta', 'se', 'z', 'prob']]   
#    filename_lst = [filename for _ in range(0, snp_df.shape[0])]
#    snp_df["loci"] = filename_lst
#    all_snps.append(snp_df)
# merge all df in the list