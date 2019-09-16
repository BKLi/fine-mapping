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


def inputfinemap(sstFolder, mafFolder, maf_prefix, FINEMAPinputFolder, chromosome):

    chromosome = "chr" + chromosome
    print(chromosome)

    list_of_sst = glob.glob(sstFolder + "{}_block*.sst".format(chromosome))

    # print(list_of_sst)
    missed_snps = []
    single_snps = []
    double_snps = []

    # count how many block sst files in total
    sst_count = 0
    # read in sst and find correspondent maf file
    for sst_file in list_of_sst:
        sst_count += 1
        # store name of the sst file
        sst_name = "".join(os.path.split(sst_file)[-1])
        print("Merging {}...".format(sst_name))
        sst = pd.read_csv(sst_file, sep="\t")
        sst_size = sst.shape[0]
        print("{} SNPs found in summary statistics".format(sst_size))
        # chromosome = re.findall(r"(chr\d+)_.+", sst_name)
        maf_file = "{}{}_{}.frq".format(mafFolder, maf_prefix, chromosome)
        maf = pd.read_csv(maf_file, sep="\s+")
        maf_size = maf.shape[0]
        print("{} SNPs found in {}".format(maf_size, maf_file))
        # eliminate SNPs missing in MAFfile
        merged = pd.merge(sst, maf, how='inner', on=["SNP"])
        # compute effect size (if not present in original sst) and add one column
        if "BETA" in list(merged):
            pass
        else:
            print("No BETA found in sst, calculating...")
            merged["BETA"] = merged["OR"].map(lambda x: math.log(x))
        merged_size = merged.shape[0]
        print("{} SNPs overlap between SNPs in block {} and 1000G".format(merged_size, sst_name))

        # reform output file
        # ignore pre-existing maf if available, use plink results instead.
        if "MAF" in list(sst):
            col_to_keep = ["SNP", "CHR_x", "BP", "A1_x", "A2_x", "MAF_y", "BETA", "SE"]
            reformed = merged[col_to_keep]
            reformed = reformed.rename(columns={"CHR_x": "chromosome",
                                                "A1_x": "allele1",
                                                "A2_x": "allele2",
                                                "SNP": "rsid",
                                                "MAF_y": "maf",
                                                "SE": "se",
                                                "BP": "position",
                                                "BETA": "beta"})
        else:
            col_to_keep = ["SNP", "CHR_x", "BP", "A1_x", "A2_x", "MAF", "BETA", "SE"]
            reformed = merged[col_to_keep]
            reformed = reformed.rename(columns={"CHR_x": "chromosome",
                                                "A1_x": "allele1",
                                                "A2_x": "allele2",
                                                "SNP": "rsid",
                                                "MAF":"maf",
                                                "SE": "se",
                                                "BP": "position",
                                                "BETA": "beta"})

        # to unify format of credible sets files
        if "MAF" in list(sst):
            reformed_set_col = ["BP", "CHR_x", "SNP", "A1_x", "A2_x", "MAF_y", "BETA", "SE", "P"]
        else:
            reformed_set_col = ["BP", "CHR_x", "SNP", "A1_x", "A2_x", "MAF", "BETA", "SE", "P"]
        reformed_set = merged[reformed_set_col]

        # create z files
        finemap_input = FINEMAPinputFolder + sst_name[:-4] + ".z"
        # print(finemap_input)
        reformed = reformed[reformed["maf"] != 0]
        if reformed.empty:
            # when SNPs in block do not exist in 1000G freq files
            # not tested yet
            for i in range(sst.shape[0]):
                missed_snps.append(sst.iloc[i].values.tolist())
        elif reformed.shape[0] == 1:
            # when block only contains one SNP after merging
            single_snps.append(reformed_set.iloc[0].values.tolist())
        elif reformed.shape[0] == 2:
            # FINEMAP cannot map blocks with only 2 SNPs
            reformed_set["block"] = sst_name
            double_snps.append(reformed_set.iloc[0].values.tolist())
            double_snps.append(reformed_set.iloc[1].values.tolist())
        elif reformed.shape[0] > 2:
            # write zfile
            reformed.to_csv(finemap_input, sep=" ", index=False)

    with open(FINEMAPinputFolder + "missed_snps_{}.set".format(chromosome), "w+") as outfile:
        for i in missed_snps:
            si = [str(s) for s in i]
            outfile.write("{}\n".format("\t".join(si)))

    with open(FINEMAPinputFolder + "single_snps_{}.set".format(chromosome), "w+") as outfile:
        for q in single_snps:
            sq = [str(s) for s in q]
            outfile.write("{}\n".format("\t".join(sq)))

    with open(FINEMAPinputFolder + "double_snps_{}.set".format(chromosome), "w+") as outfile:
        for j in double_snps:
            sj = [str(s) for s in j]
            outfile.write("{}\n".format("\t".join(sj)))


# inputFINEMAP("C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.sst", header="CHR BP SNP A1 A2 FRQ_A FRQ_U INFO OR SE P",
# MAFfile="C:\\Users\\libin\\Desktop\\tmp\\MDD_chr1.frq",
# FINEMAPinput="C:\\Users\\libin\\Desktop\\tmp\\chr1_block_2.zfile")
inputfinemap(sstFolder=sys.argv[1], mafFolder=sys.argv[2], maf_prefix=sys.argv[3], FINEMAPinputFolder=sys.argv[4], chromosome=sys.argv[5])

# --------deprecated codes below--------
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
    # print(sstFolder + "{}_block*.sst".format(chromosome))