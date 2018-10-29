import pandas as pd
import itertools
import operator
import numpy as np
from statistics import mean
from collections import OrderedDict
import sys


def SSTtoLDblock(sstFile, distance, ldFile, prefix):
    # input SST should be lifted-over and split by chromosome first
    # forms of SST varies; change import step accordingly

    ld = pd.read_table(ldFile, delim_whitespace=True)
    sst = pd.read_table(sstFile, header=None, sep="\t",
                        names=["CHR", "BP", "SNP", "A1", "A2", "FRQ_A", "FRQ_U", "INFO", "OR", "SE", "P", "report"])

    group_ids = (sst["BP"] > (sst["BP"].shift() + distance)).cumsum()
    grouped_sst = sst.groupby(group_ids)

    locus_sizes = []
    for k, sub_sst in grouped_sst:
        outfile = "{}_block_{}".format(prefix, str(k + 1))
        file = open(outfile, "w+")
        ID = sub_sst["BP"].tolist()
        sub_id = []  # list to store ID of LD-buddies

        for id in ID:
            sub_id.append(id)
            sub_ld_1 = ld[(ld["SNP1"] == id) & (ld["R2"] > 0.6)]
            sub_id.extend(sub_ld_1["SNP2"].tolist())
            sub_ld_2 = ld[(ld["SNP2"] == id) & (ld["R2"] > 0.6)]
            sub_id.extend(sub_ld_2["SNP1"].tolist())

        uniq_sub_id = list(OrderedDict.fromkeys(sub_id))
        uniq_sorted_sub_id = sorted(uniq_sub_id)

        for uid in uniq_sorted_sub_id:
            file.write("{}\n".format(uid))

        # get length of each locus
        locus_sizes.append(sub_sst["BP"].max() - sub_sst["BP"].min())

    print(locus_sizes)
    print("average: ", mean(locus_sizes))


# SSTtoLDblock("C:\\Users\libin\Desktop\MDD_chr1_filtered.sst", 1000000,
             # "C:\\Users\libin\Desktop\chr1_europe_0.2_1000000.txt", "chr1")
SSTtoLDblock(sstFile=sys.argv[1], distance=int(sys.argv[2]), ldFile=sys.argv[3], prefix=sys.argv[4])
