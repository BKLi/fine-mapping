"""
Final step of pipeline. Extracts credible set from config file (FINEMAP output), intersecting with summary statistics
to grep complete information. Output one file for each block.
"""

import pandas as pd
import glob
import os
import sys


# pd.set_option('display.max_columns', 500)
def extract_snp(config_folder, sst_folder):
    # config_folder = "C:\\Users\libin\Scripts\\fineMapping\\fine-mapping\\test\\"
    # sst_folder = "C:\\Users\libin\Scripts\\fineMapping\\fine-mapping\\test\\"
    list_of_config = glob.glob(config_folder + "*.config")
    stat_file = open("prob_stat", "w+")
    stat_list = []
    for config_file in list_of_config:
        filename = "".join(os.path.split(config_file)[-1][:-7])
        sst_file = "{}{}.sst".format(sst_folder, filename)
        print(filename)
        # config_file = "C:\\Users\libin\Desktop\\chr1_block_2.config"
        config = pd.read_table(config_file, delim_whitespace=True)
        print("config read")
        sst = pd.read_table(sst_file, delim_whitespace=True)
        snp_list = " ".join(config.head(1)["config"].tolist()).split(",")
        # record the number of SNPs in top configuration
        snp_count = len(snp_list)
        print(snp_count)
        # prob_value = config.head(1)["prob"]
        # print(prob_value)
        prob = config.head(1)["prob"].tolist()
        prob_list = prob*len(snp_list)
        pre_merge = pd.DataFrame()
        block_name = []
        [block_name.append(filename) for _ in range(len(snp_list))]
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

        prob_str = "".join(str(i) for i in prob)
        stat = "{}\t{}".format(str(snp_count), str(prob_str))
        stat_list.append(stat)
    for s in stat_list:
        stat_file.write("{}\n".format(s))

    # os.system("grep -f {} {} > {}".format(credible_set, sst_file, set_file))


extract_snp(config_folder=sys.argv[1], sst_folder=sys.argv[2])
