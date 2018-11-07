import pandas as pd
import math
import sys


def inputFINEMAP(sstFile, header,  MAFfile, FINEMAPinput):

    header_list = header.split(" ")

    sst = pd.read_table(sstFile, sep="\t")
    print(sst)

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
    if reformed.shape[0] != 1:
        reformed.to_csv(FINEMAPinput, sep=" ", index=False)


inputFINEMAP("C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.sst", header="CHR BP SNP A1 A2 FRQ_A FRQ_U INFO OR SE P",
             MAFfile="C:\\Users\\libin\\Desktop\\tmp\\MDD_chr1.frq", FINEMAPinput="C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.zfile")
# inputFINEMAP(sstFile=sys.argv[1], header=sys.argv[2], MAFfile=sys.argv[3], FINEMAPinput=sys.argv[4])
