import numpy as np
import pandas as pd
import sys
import time

start = time.clock()


def ldToCaviar(ldFile, sstFile, caviarIn):

    # read in LD matraix as pandas dataframe
    ld_df = pd.read_table(ldFile, delim_whitespace=True)
    ld_df = pd.DataFrame(ld_df)
    # multi-indexing
    ldindex = ld_df.set_index(["SNP1", "SNP2"])["R2"]
    # convert to dict for future use
    ld_dict = ldindex.to_dict()
    print("Done: multi-index")

    with open(sstFile) as sstfile:
        # initialize Caviar input matrix
        sstfile = sstfile.readlines()
        dim = len(sstfile)
        initial_matrix = np.eye(dim, dtype=int)  # return matrix with 1 on diagnal and 0 elsewhere
        matrix_listed = initial_matrix.tolist()
        print("Done: Matrix initialization")

        with open(caviarIn, "w+") as outfile:
            sidList = []
            # list of SNPs in sst
            for sst in sstfile:
                s = sst.strip().split()
                sid = s[1]
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

                            # fill in matrix
                            matrix_listed[idindex][idindex+j] = ld_value
                            matrix_listed[idindex+j][idindex] = ld_value
                            print("Done: ", id1, id2)

            print("Writing File...")
            for line in matrix_listed:
                oline = '\t'.join(str(i) for i in line)
                # output = f"{oline}\n"
                # outfile.write("\t".join(str(i) for i in line))
                # outfile.write("\n")
                # outfile.write(output)
                print(f"{oline}\n", file=outfile)

    end = time.clock()
    print(end - start)


'''if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage : ldToCaviar")
        sys.exit()'''

ldToCaviar(ldFile="C:\\Users\libin\Desktop\\testLD", sstFile="C:\\Users\libin\Desktop\\testSST", caviarIn="C:\\Users\libin\Desktop\\ldOut")
    #ldToCaviar(ldFile=sys.argv[1], sstFile=sys.argv[2], caviarIn=sys.argv[3])

