"""
Step3 of the pipeline. Write ld matrix files (one for each block)
"""

import numpy as np
import pandas as pd
import sys, os
import glob
import re
import time
import math

start = time.clock()


def ldToFINEMAP(ld_file_folder, z_file_folder, output_folder):

    # ld_file_folder = "C:\\Users\libin\Desktop\\tmp\\"
    # z_file_folder = "C:\\Users\libin\Desktop\\tmp\\"
    # output_folder = "C:\\Users\libin\Desktop\\tmp\\"

    list_of_zfile = glob.glob(z_file_folder + "*.z")
    print("number of blocks: ", len(list_of_zfile))

    for z in list_of_zfile:
        filename = "".join(os.path.split(z)[-1])[:-2]
        chromosome = re.findall(r"(chr\d+)_.+", filename)[0]
        print(filename, chromosome)
        # read in zfiles
        fzfile = pd.read_table(z, delim_whitespace=True)
        sidList = fzfile["position"].tolist()
        # print(sidList)
        print('read in: ', z)

        # initialize FINEMAP input matrix
        dim = fzfile.shape[0]
        if dim == 1:
            print("only one SNP found")
            sys.exit(0)

        initial_matrix = np.eye(dim, dtype=int)  # return matrix with 1 on diagnal and 0 elsewhere
        matrix_listed = initial_matrix.tolist()
        print("dim: ", dim)
        print("Done: Matrix initialization")

        # read in LD matraix as pandas dataframe
        ldFile = ld_file_folder + "{}_europe_0.2_1000000.txt".format(chromosome)
        # print(ldFile)

        ld_df = pd.read_table(ldFile, delim_whitespace=True)
        ld_df = pd.DataFrame(ld_df)
        # multi-indexing
        ldindex = ld_df.set_index(["SNP1", "SNP2"])["R2"]
        # convert to dict for future use
        ld_dict = ldindex.to_dict()
        print("Done: multi-index")

        FINEMAPin = output_folder + "{}.ld".format(filename)
        with open(FINEMAPin, "w+") as outfile:

            # list of SNPs in fz, in terms of BP position(not ID)
            for id1 in sidList:
                # check if SNP exists in first column of LD file
                if int(id1) in set(ld_df["SNP1"].tolist()):
                    idindex = sidList.index(id1)
                    for j in range(1, dim - idindex):
                        # iterate all the SNPs below SNP1 in sst file
                        id2 = sidList[idindex + j]
                        # check whether pair included in LD file
                        if (int(id1), int(id2)) in ld_dict:
                            ld_value = ldindex[int(id1)][int(id2)]
                            sqrt_ld_value = math.sqrt(ld_value)

                            # fill in matrix
                            matrix_listed[idindex][idindex+j] = sqrt_ld_value
                            matrix_listed[idindex+j][idindex] = sqrt_ld_value
                            print("Done: ", id1, id2)

            print("Writing File...")
            for line in matrix_listed:
                outfile.write(" ".join(str(i) for i in line))
                outfile.write("\n")

    end = time.clock()
    print(end - start)


ldToFINEMAP(ld_file_folder=sys.argv[1], z_file_folder=sys.argv[2], output_folder=sys.argv[3])

