import sys
import glob
import re
import pandas as pd

'''
write masterfiles required by FINEMAP, meanwhile generating scripts for running FINEMAP.
'''


def generateMasterfile(zfile_folder, ld_folder, output_folder, n_sample, mem):

    # zfile_folder = "C:\\Users\\libin\\Desktop\\tmp\\"
    # ld_folder = "C:\\Users\\libin\\Desktop\\tmp\\"
    # output_folder = "C:\\Users\\libin\\Desktop\\tmp\\"
    list_of_zfile = glob.glob(zfile_folder + "*.z")
    # n_sample = 40000

    # print(list_of_zfile)

    multi_com = []
    # double_com = []

    for zf in list_of_zfile:
        zf_df = pd.read_table(zf, delim_whitespace=True)
        n_zf_row = zf_df.shape[0]
        # if zf_df.shape[0] > 2:
        basename = re.findall(r"([A-Za-z0-9_]+)\..+", zf)[0]
        # print(basename)
        filename = "{}.masterfile".format(basename)
        # print(filename)
        if n_zf_row > 2 and n_zf_row < 200:
            n_causal = n_zf_row
            mcom = "finemap --sss --in-files {} --log --n-causal-snps {}".format(filename, n_causal)
            multi_com.append(mcom)
        elif n_zf_row >= 200:
            n_causal = round(n_zf_row*0.2)
            mcom = "finemap --sss --in-files {} --log --n-causal-snps {}".format(filename, n_causal)
            multi_com.append(mcom)
        # elif n_zf_row == 2:
            # snps = ",".join(zf_df["rsid"].tolist())
            # dcom = "finemap --config --in-files {} --rsids {}".format(filename, snps)
            # double_com.append(dcom)
        with open(output_folder + filename, "w+") as outfile:
            zfile = zf
            ldFile = "{}{}.ld".format(ld_folder, basename)
            snp = "{}.snp".format(basename)
            config = "{}.config".format(basename)
            log = "{}.log".format(basename)
            # cred = "{}.cred".format(basename)

            outfile.write("z;ld;snp;config;log;n_samples\n")
            outfile.write("{};{};{};{};{};{}".format(zfile, ldFile, snp, config, log, n_sample))

        with open(output_folder + "FINEMAP_multi_SNP.sh", "w+") as multiout:
            multiout.write('''#!/bin/bash\n#$ -l h_rt=12:0:0\n#$ -l mem_free={}G\n#$ -S /bin/bash\n#$ -cwd\n#$ -j y\n#$ -r y\n\n'''.format(mem))
            multiout.write("export PATH=/scrapp2/bingkun/FINEMAP/finemap_v1.3_x86_64:$PATH\n")
            for i in multi_com:
                multiout.write("{}\n".format(i))
        # with open(output_folder + "FINEMAP_doub_SNP.snp", "w+") as doubout:
            # for j in double_com:
                # doubout.write("{}\n".format(j))


generateMasterfile(zfile_folder=sys.argv[1], ld_folder=sys.argv[2], output_folder=sys.argv[3], n_sample=sys.argv[4], mem=sys.argv[5])
