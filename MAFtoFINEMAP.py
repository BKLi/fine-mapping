import pandas as pd
import math
import sys, os
import glob
import re


def inputFINEMAP(sstFolder, MAFfileFolder, maf_prefix, FINEMAPinputFolder):
    # FINEMAPinputFolder = "C:\\Users\\libin\\Desktop\\tmp"
    # sstFolder = "C:\\Users\\libin\\Desktop\\tmp"
    # MAFfileFolder = "C:\\Users\\libin\\Desktop\\tmp"
    # maf_prefix = 'MDD'

    if not MAFfileFolder.endswith("\\"):
        MAFfileFolder = MAFfileFolder + "\\"

    if not sstFolder.endswith("\\"):
        sstFolder = sstFolder + "\\"

    if not FINEMAPinputFolder.endswith("\\"):
        FINEMAPinputFolder = FINEMAPinputFolder + "\\"

    list_of_sst = glob.glob(sstFolder + "*block*.sst")
    print(list_of_sst)

    single_snps = []

    for sst_file in list_of_sst:
        sst = pd.read_table(sst_file, sep="\t")
        sst_name = "".join(os.path.split(sst_file)[-1])
        chromosome = re.findall(r"(chr\d+)_.+", sst_name)
        maf_file = "{}{}_{}.frq".format(MAFfileFolder, maf_prefix, chromosome[0])
        # print(maf_file)

        MAF = pd.read_table(maf_file, sep="\s+")
        # eliminate SNPs missing in MAFfile
        merged = pd.merge(sst, MAF, how='inner', on=["SNP"])
        # compute effect size and add one column
        merged["beta"] = merged["OR"].map(lambda x: math.log(x))
        # reform output file
        col_to_keep = ["SNP", "CHR_x", "BP", "A1_x", "A2_x", "MAF", "beta", "SE"]
        reformed = merged[col_to_keep]
        reformed = reformed.rename(columns={"CHR_x": "chromosome", "A1_x": "allele1", "A2_x": "allele2",
                                            "SNP": "rsid", "MAF": "maf", "SE": "se", "BP": "position"})

        finemap_input = FINEMAPinputFolder + sst_name[:-4] + ".zfile"
        # print(finemap_input)
        if reformed.shape[0] != 1:
            reformed.to_csv(finemap_input, sep=" ", index=False)
        else:
            single_snps.append(reformed.iloc[0].values.tolist())

    with open(FINEMAPinputFolder + "single_snps", "w+") as outfile:
        for i in single_snps:
            si = [str(s) for s in i]
            outfile.write("{}\n".format("\t".join(si)))


# inputFINEMAP("C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.sst", header="CHR BP SNP A1 A2 FRQ_A FRQ_U INFO OR SE P",
# MAFfile="C:\\Users\\libin\\Desktop\\tmp\\MDD_chr1.frq",
# FINEMAPinput="C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.zfile")
inputFINEMAP(sstFolder=sys.argv[1], maf_prefix=sys.argv[2], MAFfileFolder=sys.argv[3], FINEMAPinputFolder=sys.argv[4])
