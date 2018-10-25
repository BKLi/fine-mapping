import pandas as pd
import itertools
import operator
import numpy as np
from statistics import mean
import sys


def splitSSTbyDistance(sstFile, distance, ldFile):
    # input SST should be lifted-over and split by chromosome first
    # forms of SST varies; change import step accordingly

    ld = pd.read_table(ldFile, delim_whitespace=True)
    sst = pd.read_table(sstFile, header=None, sep="\t",
                        names=["CHR", "BP", "SNP", "A1", "A2", "FRQ_A", "FRQ_U", "INFO", "OR", "SE", "P", "report"])

    group_ids = (sst["BP"] > (sst["BP"].shift() + distance)).cumsum()
    grouped_sst = sst.groupby(group_ids)

    locus_sizes = []
    for k, sub_sst in grouped_sst:
        Pmin = sub_sst["P"].min()
        IDmin = int(sub_sst[sub_sst["P"] == Pmin]["BP"])
        sub_ld_1 = ld[(ld["SNP1"] == IDmin) & (ld["R2"] > 0.6)]
        sub_id = sub_ld_1["SNP2"].tolist()
        sub_ld_2 = ld[(ld["SNP2"] == IDmin) & (ld["R2"] > 0.6)]
        sub_id.extend(sub_ld_2["SNP1"].tolist())
        print(sub_id)
        print(len(sub_id))
        # get length of each locus
        locus_sizes.append(sub_sst["BP"].max() - sub_sst["BP"].min())

    print(locus_sizes)
    print("average: ", mean(locus_sizes))


splitSSTbyDistance("C:\\Users\libin\Desktop\MDD_chr21.sst", 1000000,
                   "C:\\Users\libin\Desktop\chr21_europe_0.2_1000000.txt")
# splitSSTbyDistance(sstFile=sys.argv[1], distance=int(sys.argv[2]), chromosome=sys.argv[3])
