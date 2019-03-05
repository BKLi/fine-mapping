"""
Write input files for FINEMAP. Output masterfiles(one for each block) and wrapper sh file.
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--prefix", type=str, help="Prefix of blockfile names", required=True)
parser.add_argument("-chr", "--chrom", type=str, help="Chromosome number", required=True)
parser.add_argument("-z", "--zDir", type=str, help="Path to zfile directory", required=True)
parser.add_argument("-ld", "--ldDir", type=str, help="Path to ld matrix directory", required=True)
parser.add_argument("-n", "--nSamples", type=str, help="GWAS sample size", required=True)
args = parser.parse_args()

def writeinput():
    prefix = args.prefix
    chrom = args.chrom
    zDir = args.zDir
    ldDir = args.ldDir
    nSamples = args.nSamples

    filename = "{}_chr{}".format(prefix, chrom)
    with open(filename + ".masterfile", "w+") as outfile:

        zfile = "{}{}.z".format(zDir, filename)
        ldFile = "{}{}.ld".format(ldDir, filename)
        snp = "{}.snp".format(filename)
        config = "{}.config".format(filename)
        log = "{}.log".format(filename)

        outfile.write("z;ld;snp;config;log;n_samples\n")
        outfile.write("{};{};{};{};{};{}".format(zfile, ldFile, snp, config, log, nSamples))


writeinput()
