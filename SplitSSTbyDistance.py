import pandas as pd
import itertools
import operator
import numpy as np
from statistics import mean
import sys


def splitSSTbyDistance(sstFile, distance, chromosome):
    # input SST should be lifted-over and split by chromosome first
    # forms of SST varies; change import step accordingly

    # import sst as dataframe
    sst = pd.read_table(sstFile, header=None, sep="\t",
                        names=["CHR", "BP", "SNP", "A1", "A2", "FRQ_A", "FRQ_U", "INFO", "OR", "SE", "P", "report"])

    # position = sst["BP"].tolist()
    # print(position[:10])
    # dif = list(itertools.starmap(operator.sub, zip(position[1:], position)))
    # splitPoint = [dif.index(i) for i in dif if(i > distance)]
    # print(splitPoint[:10])
    # print(len(splitPoint))

    group_ids = (sst["BP"] > (sst["BP"].shift() + distance)).cumsum()
    grouped_sst = sst.groupby(group_ids)

    locus_sizes = []
    for k, sub_sst in grouped_sst:
        print(sub_sst.head())
        outfile = "{}_sub_set_{}".format(chromosome, str(k+1))
        # print(outfile)
        file = open(outfile, "w+")
        sub_sst.to_csv(file, sep="\t", index=False)
        # get length of each locus
        locus_sizes.append(sub_sst["BP"].max() - sub_sst["BP"].min())

    print(locus_sizes)
    print("average: ", mean(locus_sizes))


splitSSTbyDistance("C:\\Users\libin\Desktop\MDD_chr21.sst", 1000000, 1)
# splitSSTbyDistance(sstFile=sys.argv[1], distance=int(sys.argv[2]), chromosome=sys.argv[3])
