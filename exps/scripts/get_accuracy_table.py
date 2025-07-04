import sys
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate


def main():
    wd = sys.argv[1]
    data = []
    # Ws = ["0.1", "0.25", "0.33", "0.5", "0.9"]
    for txt in glob.glob(os.path.join(wd, "*", "results.*.txt")):
        print(txt)
        dataset, fname = txt.split("/")[-2:]
        # print(dataset, fname)
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
        for i, line in enumerate(open(txt)):
            if i < 5:
                continue
            V, c, h, ARI, _, _, non_singl, singl, *rest = line.strip("\n").split(",")
            V = round(float(V), 3)
            c = round(float(c), 3)
            h = round(float(h), 3)
            ARI = round(float(ARI), 3)
            non_singl = int(non_singl)
            singl = int(singl)
            data.append([dataset, mode, t, w, V, c, h, ARI, non_singl, singl])
            # if mode != "weighted":
            #     for w in Ws:
            #         data.append([dataset, mode, t, w, V, c, h, ARI, non_singl, singl])
    df = pd.DataFrame(
        data,
        columns=[
            "dataset",
            "mode",
            "t",
            "w",
            "V",
            "c",
            "h",
            "ARI",
            "nonsingl",
            "singl",
        ],
    )
    df["Clusters"] = df["singl"] + df["nonsingl"]
    df.sort_values(["t", "w"], ascending=[True, True], inplace=True)

    # pb_df = df[df["dataset"] == "PB"].drop(labels=["dataset"], axis=1)
    # ont_df = df[df["dataset"] == "ONT"].drop(labels=["dataset"], axis=1)
    # ont_df.to_csv("ont.csv", index=False)
    # pb_df.to_csv("pacbio.csv", index=False)
    # print(tabulate(pb_df, tablefmt="latex", showindex=False))
    # print(tabulate(ont_df, tablefmt="latex", showindex=False))

    print(tabulate(df, tablefmt="latex", showindex=False))

    # g = sns.FacetGrid(df, col="t", row="dataset")
    # g.map(sns.lineplot, "w", "V")
    # plt.savefig("x.png")


if __name__ == "__main__":
    main()
