import sys
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate
from matplotlib.lines import Line2D

sns.set_style("whitegrid")


def main():
    wd = sys.argv[1]
    data = []
    Ws = ["0.1", "0.25", "0.33", "0.5", "0.9"]
    for txt in glob.glob(os.path.join(wd, "PB", "results.*.txt")):
        # print(txt)
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
        if t == 0.75:
            continue
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

    print(tabulate(df, tablefmt="latex", showindex=False))

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True, sharey=True)
    for col, t in enumerate(df["t"].unique()):
        print(t)
        subdf_w = df[(df["t"] == t) & (df["mode"] == "weighted")]
        subdf = df[(df["t"] == t) & (df["mode"] == "original")]

        axes[col].axhline(
            y=subdf["V"].iloc[0],
            linestyle="--",
            color="salmon",
            linewidth=2,
        )
        axes[col].axhline(
            y=subdf["ARI"].iloc[0],
            linestyle="--",
            color="forestgreen",
            linewidth=2,
        )

        sns.pointplot(
            data=subdf_w,
            x="w",
            y="V",
            ax=axes[col],
            color="salmon",
            linewidth=2,
        )
        sns.pointplot(
            data=subdf_w,
            x="w",
            y="ARI",
            ax=axes[col],
            color="forestgreen",
            linewidth=2,
        )

        axes[col].set_title(r"$\tau$" + "=" + str(t))
        axes[col].set_xlabel(r"$\sigma$")
        axes[col].set_ylim(-0.1, 1.1)
        axes[col].set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        if col == 0:
            axes[col].set_ylabel("Value")

    custom_lines = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="salmon",
            linewidth=2,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="forestgreen",
            linewidth=2,
        ),
        Line2D(
            [0],
            [0],
            color="salmon",
            linestyle="--",
            linewidth=2,
        ),
        Line2D(
            [0],
            [0],
            color="forestgreen",
            linestyle="--",
            linewidth=2,
        ),
    ]

    leg = plt.legend(
        custom_lines,
        ["V (SC)", "ARI (SC)", "V (IOC3)", "ARI (IOC3)"],
        ncol=2,
        # bbox_to_anchor=(0.15, -0.35),
        loc="lower left",
        borderaxespad=0,
    )
    fig.add_artist(leg)

    plt.tight_layout()
    plt.show()
    # plt.savefig("x.png")


if __name__ == "__main__":
    main()
