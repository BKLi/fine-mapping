import pandas as pd


def expand_ld_block(subgroup, ldFile):
    sst = pd.read_table(subgroup, sep="\t")
    ld = pd.read_table(ldFile, delim_whitespace=True)
    # print(sst)
    Pmin = sst["P"].min()
    IDmin = int(sst[sst["P"] == Pmin]["BP"])
    sub_ld_1 = ld[(ld["SNP1"] == IDmin) & (ld["R2"] > 0.6)]
    sub_id = sub_ld_1["SNP2"].tolist()
    sub_ld_2 = ld[(ld["SNP2"] == IDmin) & (ld["R2"] > 0.6)]
    sub_id.extend(sub_ld_2["SNP1"].tolist())
    print(sub_id)


expand_ld_block("C:\\Users\libin\Scripts\\fineMapping\\fine-mapping\\test\\1_sub_set_1",
                "C:\\Users\libin\Desktop\chr21_europe_0.2_1000000.txt")


