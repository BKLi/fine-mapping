import pandas as pd
import math
import sys


def inputFINEMAP(sstFile, MAFfile, FINEMAPinput):

    sst = pd.read_table(sstFile, header=None, sep="\t",
                        names=["CHR", "BP", "SNP", "A1", "A2", "FRQ_A", "FRQ_U", "INFO", "OR", "SE", "P", "report"])

    MAF = pd.read_table(MAFfile, sep="\s+")
    # eliminate SNPs missing in MAFfile
    merged = pd.merge(sst, MAF, how='inner', on=["SNP"])
    # compute effect size and add one column
    merged["beta"] = merged["OR"].map(lambda x: math.log(x))
    # reform output file
    col_to_keep = ["SNP", "CHR_x", "BP", "A1_x", "A2_x", "MAF", "beta", "SE"]
    reformed = merged[col_to_keep]
    reformed = reformed.rename(columns={"CHR_x": "chromosome", "A1_x": "allele1", "A2_x": "allele2",
                                        "SNP": "rsid", "MAF": "maf", "SE": "se", "BP": "position"})

    reformed.to_csv(FINEMAPinput, sep=" ", index=False)


inputFINEMAP(sstFile=sys.argv[1], MAFfile=sys.argv[2], FINEMAPinput=sys.argv[3])
