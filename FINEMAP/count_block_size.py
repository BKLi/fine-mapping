import sys
import glob
import pandas as pd


def countsize(z_file_folder):
    list_of_zfile = glob.glob(z_file_folder + "*.z")
    count = []
    for z in list_of_zfile:
        print(z)
        z_df = pd.read_csv(z, delim_whitespace=True)
        pos_min = z_df["position"].min()
        pos_max = z_df["position"].max()
        pos_size = pos_max - pos_min
        count.append(pos_size)
    with open("count_size", "w+") as outfile:
        for i in count:
            outfile.write("{}\n".format(i))


countsize(z_file_folder=sys.argv[1])

