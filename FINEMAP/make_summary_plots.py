# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:07:16 2019

@author: libin
"""

from pathlib import Path
import re
import pandas as pd

# fig.1c: # of regions & posterior bins
# can't really make this plot because SNPs overlap among regions
snp_path = Path(r'C:\Users\libin\UCSF\FINEMAP\AD_set\AD_snp\snp')
snp_file_path_list = [str(fp) for fp in snp_path.glob('**/*_block_*.snp')]
snp_file_name_list = ["".join(re.findall(r'(.+).snp',fp.name)) for fp in snp_path.glob('**/*_block_*.snp')]
# list to store all the snp dataframes for integrating later
all_snps = []
for file in snp_file_path_list:
    filename = snp_file_name_list[snp_file_path_list.index(file)]
    snp_df = pd.read_csv(file, delim_whitespace=True)
    snp_df = snp_df[["rsid", 'chromosome', 'position', 'maf', 'beta', 'se', 'z', 'prob']]   
    filename_lst = [filename for _ in range(0, snp_df.shape[0])]
    snp_df["loci"] = filename_lst
    all_snps.append(snp_df)
# merge all df in the list
