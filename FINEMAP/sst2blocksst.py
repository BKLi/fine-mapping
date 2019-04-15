"""
Step one of pipeline.
Takes in chromosome-wide summary statistics containing significant SNPs, groups SNPs into blocks by distance,
then includes LD buddies of each SNP in each group.
Outputs block-wide summary statistics as input for maf2finemap.py

"""

import sys
# import json
import pandas as pd
import numpy as np
from statistics import mean
from collections import OrderedDict


def sst2blocksst(filt_sst, chr_sst, header, distance, ldFile, prefix, cutoff):
    # input SST should be lifted-over and split by chromosome first
    # forms of SST varies; change import step accordingly
    # prefix:usually chromosome identifier
    # cutoff: LD expanding cutoff
    # header : string, header of sst file
    cutoff = int(cutoff)
    header_list = header.split(" ")
    ld = pd.read_table(ldFile, delim_whitespace=True)
    sst = pd.read_table(filt_sst, header=None, sep="\t",
                        names=header_list)
    chr_sst = pd.read_table(chr_sst, header=None, sep="\t",
                            names=header_list)

    if sst.shape[0] == 0:
        print("empty file")
        sys.exit(0)

    group_ids = (sst["BP"] > (sst["BP"].shift() + distance)).cumsum()
    grouped_sst = sst.groupby(group_ids)

    locus_sizes = []
    for k, sub_sst in grouped_sst:

        # print(sub_sst.shape[0])

        outfile = "{}_block_{}.sst".format(prefix, str(k + 1))
        # file = open(outfile, "w+")
        position = sub_sst["BP"].tolist()

        sub_id = []  # list to store position of LD-buddies

        for id in position:
            # id may be snp1 or snp2 in LD matrix file
            sub_id.append(id)
            sub_ld_1 = ld[(ld["SNP1"] == id) & (ld["R2"] > cutoff)]
            sub_id.extend(sub_ld_1["SNP2"].tolist())
            sub_ld_2 = ld[(ld["SNP2"] == id) & (ld["R2"] > cutoff)]
            sub_id.extend(sub_ld_2["SNP1"].tolist())
        # remove duplicates from list
        uniq_sub_id = list(OrderedDict.fromkeys(sub_id))
        uniq_sorted_sub_id = sorted(uniq_sub_id)

        # intersect each LD block with sst
        uid_df = pd.DataFrame({"BP": uniq_sorted_sub_id})
        merged = pd.merge(uid_df, chr_sst, how='inner', on=["BP"])
        # print(merged.head())
        print("number of SNPs: ", merged.shape[0])
        # file.write("{}\n".format(uid))
        merged.to_csv(outfile, sep="\t", index=False)

        # get length of each locus
        locus_sizes.append(merged["BP"].max() - merged["BP"].min())
    print(locus_sizes)
    print("average: ", mean(locus_sizes))


# SSTtoLDblock("C:\\Users\\libin\\Desktop\\tmp\\MDD_chr1_filtered.sst", "CHR BP SNP A1 A2 FRQ_A FRQ_U INFO OR SE P",
# 1000000, "C:\\Users\\libin\\Desktop\\tmp\\chr1_europe_0.2_1000000.txt", "chr1")
sst2blocksst(filt_sst=sys.argv[1], chr_sst=sys.argv[2], header=sys.argv[3], distance=int(sys.argv[4]),
             ldFile=sys.argv[5], prefix=sys.argv[6])
