import sys


def rmMissSnp(ldFile, inSst, outSst):
    with open(ldFile) as ld:
        ldDict = {}
        for pair in ld:
            spair = pair.strip().split()
            snp1 = spair[0]
            snp2 = spair[1]
            ldDict[snp1] = ""
            ldDict[snp2] = ""
        print(len(ldDict))
        print(list(ldDict.keys())[0], list(ldDict.keys())[1])
        with open(inSst) as sst:
            with open(outSst, "w+") as outfile:
                i = 0
                j = 0
                for s in sst:
                    ss = s.strip().split()
                    i += 1
                    snp = ss[1]
                    if snp in ldDict:
                        j += 1
                        outfile.write("\t".join(ss))
                        outfile.write("\n")
    print(i)
    print(j)
    print(j/i)


rmMissSnp(ldFile=sys.argv[1], inSst=sys.argv[2], outSst=sys.argv[3])