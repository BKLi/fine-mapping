import pandas as pd
import itertools
import operator
import numpy as np
from statistics import mean
from collections import OrderedDict
import sys

'''
take in after-filtering chromosome wide sst files, group SNPs into blocks by distance.
output block-wide sst files as input for MAFtoFINEMAP.py

'''
# def SSTtoLDblock(sstFile, header, distance, ldFile, prefix):
# input SST should be lifted-over and split by chromosome first
# forms of SST varies; change import step accordingly


def SSTtoLDblock(sstFile, header, distance, ldFile, prefix):

    # header = "CHR BP SNP A1 A2 FRQ_A FRQ_U INFO OR SE P"
    # distance = 1000000
    # prefix = "MDD"
    # sstFile = "C:\\Users\\libin\\Desktop\\tmp\\MDD_chr22_filtered.sst"
    # ldFile = "C:\\Users\libin\Desktop\\tmp\\chr22_all_pop.ld"
    header_list = header.split(" ")
    ld = pd.read_table(ldFile, delim_whitespace=True)
    # print(ld.head())
    sst = pd.read_table(sstFile, header=None, sep="\t",
                        names=header_list)

    # in case some chromosomes do not contain index SNPs
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
        ID = sub_sst["BP"].tolist()
        sub_id = []  # list to store ID of LD-buddies

        for id in ID:
            sub_id.append(id)
            sub_ld_1 = ld[(ld["BP_A"] == id) & (ld["R2"] > 0.6)]
            sub_id.extend(sub_ld_1["BP_B"].tolist())
            sub_ld_2 = ld[(ld["BP_B"] == id) & (ld["R2"] > 0.6)]
            sub_id.extend(sub_ld_2["BP_A"].tolist())

        uniq_sub_id = list(OrderedDict.fromkeys(sub_id))
        uniq_sorted_sub_id = sorted(uniq_sub_id)

        # intersect each LD block with sst
        uid_df = pd.DataFrame({"BP": uniq_sorted_sub_id})
        merged = pd.merge(uid_df, sst, how='inner', on=["BP"])
        # print(merged.head())
        # print(merged.shape[0])
        # file.write("{}\n".format(uid))
        merged.to_csv(outfile, sep="\t", index=False)

        # get length of each locus
        locus_sizes.append(sub_sst["BP"].max() - sub_sst["BP"].min())

    print(locus_sizes)
    print("average: ", mean(locus_sizes))


# SSTtoLDblock("C:\\Users\\libin\\Desktop\\tmp\\MDD_chr1_filtered.sst", "CHR BP SNP A1 A2 FRQ_A FRQ_U INFO OR SE P",
# 1000000, "C:\\Users\\libin\\Desktop\\tmp\\chr1_europe_0.2_1000000.txt", "chr1")
SSTtoLDblock(sstFile=sys.argv[1], header=sys.argv[2], distance=int(sys.argv[3]), ldFile=sys.argv[4], prefix=sys.argv[5])
