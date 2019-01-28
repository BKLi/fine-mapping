import pandas as pd
import glob
import os
import sys

'''
Extract credible set from config file (FINEMAP output) and intersect with sst file to grep full info
'''

pd.set_option('display.max_columns', 500)


def extract_snp(config_folder, sst_folder):
# config_folder = "C:\\Users\libin\Scripts\\fineMapping\\fine-mapping\\test\\"
# sst_folder = "C:\\Users\libin\Scripts\\fineMapping\\fine-mapping\\test\\"
    list_of_config = glob.glob(config_folder + "*.config")
    for config_file in list_of_config:
        filename = "".join(os.path.split(config_file)[-1][:-7])
        sst_file = "{}{}.sst".format(sst_folder, filename)
        print(filename)
        # config_file = "C:\\Users\libin\Desktop\\chr1_block_2.config"
        config = pd.read_table(config_file, delim_whitespace=True)
        sst = pd.read_table(sst_file, delim_whitespace=True)
        snp_list = " ".join(config.head(1)["config"].tolist()).split(",")
        prob = config.head(1)["prob"].tolist()
        prob_list = prob*len(snp_list)
        pre_merge = pd.DataFrame()
        block_name = []
        [block_name.append(filename) for i in range(len(snp_list))]
        pre_merge["BLOCK"] = block_name
        pre_merge["SNP"] = snp_list
        pre_merge["PROB"] = prob_list
        merged = pd.merge(pre_merge, sst, how="inner", on=['SNP'])
        # print(merged)
        # col_ordered = ["CHR", "BP", "SNP", "PROB", "P", "BLOCK", "A1", "A2", "Z", "MAF", "BETA", "SE"]
        # merged = merged[col_ordered]
        # print(merged)
        credible_set = "{}{}-{}.set".format(config_folder, filename, str(prob))
        merged.to_csv(credible_set, sep="\t", index=False)
    # os.system("grep -f {} {} > {}".format(credible_set, sst_file, set_file))


extract_snp(config_folder=sys.argv[1], sst_folder=sys.argv[2])
