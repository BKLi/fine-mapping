import numpy as np
import pandas as pd
import sys
import time

start = time.clock()


def ldToFINEMAP(ldFile, fzFile, FINEMAPIn):

    # read in LD matraix as pandas dataframe
    ld_df = pd.read_table(ldFile, delim_whitespace=True)
    ld_df = pd.DataFrame(ld_df)
    # multi-indexing
    ldindex = ld_df.set_index(["SNP1", "SNP2"])["R2"]
    # convert to dict for future use
    ld_dict = ldindex.to_dict()
    print("Done: multi-index")

    with open(fzFile) as fzfile:
        # initialize FINEMAP input matrix
        fzfile = fzfile.readlines()[1:]
        dim = len(fzfile)
        initial_matrix = np.eye(dim, dtype=int)  # return matrix with 1 on diagnal and 0 elsewhere
        matrix_listed = initial_matrix.tolist()
        print("dim: ", dim)
        print("Done: Matrix initialization")

        with open(FINEMAPIn, "w+") as outfile:
            sidList = []
            # list of SNPs in fz, in terms of BP position(not ID)
            for fz in fzfile:
                s = fz.strip().split()
                sid = s[2]
                sidList.append(sid)
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


ldToFINEMAP(ldFile=sys.argv[1], fzFile=sys.argv[2], FINEMAPIn=sys.argv[3])

