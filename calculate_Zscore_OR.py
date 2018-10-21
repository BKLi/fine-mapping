import sys
import math


def calcZscore(summary, zscore):
    with open(summary) as infile:
        with open(zscore, "w+") as outfile:
            lines = infile.readlines()
            for line in lines[1:]:
                sline = line.strip().split()
                OR = float(sline[8])
                SE = float(sline[9])
                Zscore = (math.log(OR))/SE
                sline.append(str(Zscore))
                outline = sline[1]+"\t"+sline[-1]+"\n"
                outfile.write(outline)


calcZscore(summary=sys.argv[1], zscore=sys.argv[2])



