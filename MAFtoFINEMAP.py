"""
Step 2 of pipeline.
Takes in block-wide summary statistics fom step one, incorporating pre-calculated MAF (.frq file from plink).
Outputs z file required as input for FINEMAP.
maf_prefix is usually the name of disease/phenotype

"""

import pandas as pd
import math
import sys, os
import glob
import re


def inputfinemap(sstFolder, mafFolder, maf_prefix, FINEMAPinputFolder):

    # FINEMAPinputFolder = "C:\\Users\\libin\\Desktop\\tmp\\"
    # sstFolder = "C:\\Users\\libin\\Desktop\\tmp\\"
    # MAFfileFolder = "C:\\Users\\libin\\Desktop\\tmp\\"
    # maf_prefix = 'MDD'

    # the following doesn't work in linux
    # if not MAFfileFolder.endswith("\\"):
        # MAFfileFolder = MAFfileFolder + "\\"

    # if not sstFolder.endswith("\\"):
        # sstFolder = sstFolder + "\\"

    # if not FINEMAPinputFolder.endswith("\\"):
        # FINEMAPinputFolder = FINEMAPinputFolder + "\\"

    list_of_sst = glob.glob(sstFolder + "*block*.sst")

    # print(list_of_sst)
    single_snps = []
    double_snps = []

    for sst_file in list_of_sst:
        sst = pd.read_table(sst_file, sep="\t")
        sst_name = "".join(os.path.split(sst_file)[-1])
        chromosome = re.findall(r"(chr\d+)_.+", sst_name)
        maf_file = "{}{}_{}.frq".format(mafFolder, maf_prefix, chromosome[0])
        print(maf_file)

        MAF = pd.read_table(maf_file, sep="\s+")
        # eliminate SNPs missing in MAFfile
        merged = pd.merge(sst, MAF, how='inner', on=["SNP"])
        # compute effect size and add one column
        if "BETA" in list(merged):
            pass
        else:
            merged["BETA"] = merged["OR"].map(lambda x: math.log(x))
        print(merged.head())
        # reform output file
        if "MAF" in list(sst):
            col_to_keep = ["SNP", "CHR_x", "BP", "A1_x", "A2_x", "MAF_y", "BETA", "SE"]
            reformed = merged[col_to_keep]
            reformed = reformed.rename(columns={"CHR_x": "chromosome", "A1_x": "allele1", "A2_x": "allele2",
                                                "SNP": "rsid", "MAF_y": "maf", "SE": "se", "BP": "position",
                                                "BETA": "beta"})
        else:
            col_to_keep = ["SNP", "CHR_x", "BP", "A1_x", "A2_x", "MAF", "BETA", "SE"]
            reformed = merged[col_to_keep]
            reformed = reformed.rename(columns={"CHR_x": "chromosome", "A1_x": "allele1", "A2_x": "allele2",
                                                "SNP": "rsid", "MAF": "maf", "SE": "se", "BP": "position", "BETA": "beta"})

        # to unify format of credible sets files
        if "MAF" in list(sst):
            reformed_set_col = ["BP", "CHR_x", "SNP", "A1_x", "A2_x", "MAF_y", "BETA", "SE", "P"]
        else:
            reformed_set_col = ["BP", "CHR_x", "SNP", "A1_x", "A2_x", "MAF", "BETA", "SE", "P"]
        reformed_set = merged[reformed_set_col]

        finemap_input = FINEMAPinputFolder + sst_name[:-4] + ".z"
        # print(finemap_input)

        if reformed.empty:
            # when SNPs in block do not exist in 1000G freq files
            single_snps.append(sst.iloc[0].values.tolist())
        elif reformed.shape[0] == 1:
            # when block only contains single SNP
            single_snps.append(reformed_set.iloc[0].values.tolist())
        elif reformed.shape[0] == 2:
            # FINEMAP cannot finemap blocks with only 2 SNPs
            double_snps.append(reformed_set.iloc[0].values.tolist())
            double_snps.append(reformed_set.iloc[1].values.tolist())
        elif reformed.shape[0] > 2:
            # write zfile
            reformed.to_csv(finemap_input, sep=" ", index=False)

    with open(FINEMAPinputFolder + "single_snps.set", "w+") as outfile:
        for i in single_snps:
            si = [str(s) for s in i]
            outfile.write("{}\n".format("\t".join(si)))

    with open(FINEMAPinputFolder + "double_snps.set", "w+") as outfile:
        for j in double_snps:
            sj = [str(s) for s in j]
            outfile.write("{}\n".format("\t".join(sj)))


# inputFINEMAP("C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.sst", header="CHR BP SNP A1 A2 FRQ_A FRQ_U INFO OR SE P",
# MAFfile="C:\\Users\\libin\\Desktop\\tmp\\MDD_chr1.frq",
# FINEMAPinput="C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.zfile")
inputfinemap(sstFolder=sys.argv[1], mafFolder=sys.argv[2], maf_prefix=sys.argv[3], FINEMAPinputFolder=sys.argv[4])
