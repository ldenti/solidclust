import sys
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate

def parse_time(fn):
    ram, time = 0, 0
    for line in open(fn):
        line = line.strip("\t\n ")
        print(line)
        if line.startswith("Maximum"):
            ram = round(int(line.split(" ")[5])/1024/1024, 1)
        elif line.startswith("Elapsed"):
            time = line.split(" ")[7]
    return ram, time

def main():
    wd = sys.argv[1]
    data = []
    Ws = ["0.1", "0.25", "0.33", "0.5", "0.9"]
    for txt in glob.glob(os.path.join(wd, "logs", "*", "*.time")):
        print(txt)
        dataset, fname = txt.split("/")[-2:]
        print(dataset, fname)
        mode, t, w = None, None, None
        split = fname.split(".")
        if len(split) == 7:
            _, mode, t, td, w, wd, _ = split
            t = float(t[1:] + "." + td)
            w = float(w[1:] + "." + wd)
        else:
            _, mode, t, td, w, _ = split
            t = float(t[1:] + "." + td)
            w = float(w[1:])
        ram, time = parse_time(txt)
        data.append([dataset, mode, t, w, time, ram])

    df = pd.DataFrame(
        data,
        columns=[
            "dataset",
            "mode",
            "t",
            "w",
            "Time",
            "RAM",
        ],
    )
    df.sort_values(["t", "w"], ascending=[True, True], inplace=True)
    #ont_df.to_csv('ont.csv', index=False)
    #pb_df.to_csv('pacbio.csv', index=False)

    print(tabulate(df, showindex=False))
    print(tabulate(df, tablefmt="latex", showindex=False))

    # g = sns.FacetGrid(df, col="t", row="dataset")
    # g.map(sns.lineplot, "w", "V")
    # plt.savefig("x.png")


if __name__ == "__main__":
    main()
