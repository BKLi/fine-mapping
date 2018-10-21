import sys


def generateMasterfile(prefix, chr, zDir, ldDir, nSamples):

    filename = "{}_chr{}".format(prefix, chr)
    with open(filename + ".masterfile", "w+") as outfile:

        zfile = "{}{}.z".format(zDir, filename)
        ldFile = "{}{}.ld".format(ldDir, filename)
        snp = "{}.snp".format(filename)
        config = "{}.config".format(filename)
        log = "{}.log".format(filename)

        outfile.write("z;ld;snp;config;log;n_samples\n")
        outfile.write("{};{};{};{};{};{}".format(zfile, ldFile, snp, config, log, nSamples))


generateMasterfile(prefix=sys.argv[1], chr=sys.argv[2], zDir=sys.argv[3], ldDir=sys.argv[4], nSamples=sys.argv[5])
