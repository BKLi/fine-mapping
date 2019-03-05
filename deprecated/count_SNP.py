# This script retain all the SNPs present in 1000g
# Apply this to LD matrix generated with 1000g


def count_snp(ld, snpALL, missingSNP):
    count = 0
    with open(ld) as infile:
        with open(snpALL) as snpFile:
            with open(missingSNP, "w+") as outfile:
                snpList = []
                infile_read = infile.readlines()[1:]
                for line in infile_read:
                    sline = line.strip().split()
                    snpA = sline[2]
                    snpList.append(snpA)
                snpSet = set(snpList)
                # store all the reported SNPs in LD matrix
                uniqSNP = {}
                for snp in snpSet:
                    uniqSNP[snp] = ""
                for i in snpFile:
                    si = i.strip()
                    if si not in uniqSNP:
                        count += 1
                        outfile.write(si)
                        outfile.write("\n")

            print(count)
            print(len(snpList))
            print(len(snpSet))


count_snp("C:\\Users\libin\Desktop\MDD_full_modified.ld", "C:\\Users\libin\Desktop\MDD_full_SNPlist",
          "C:\\Users\libin\Desktop\missingSNPs")
