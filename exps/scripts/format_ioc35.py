import sys
from Bio import SeqIO


def main():
    fq_fn = sys.argv[1]
    txt_fn = sys.argv[2]

    ids = []
    for read in SeqIO.parse(fq_fn, "fastq"):
        ids.append(read.id)
        if len(ids) % 50000 == 0:
            print(len(ids), file=sys.stderr)

    for line in open(txt_fn):
        if line.startswith("read ID"):
            continue
        ridx, cidx = line.strip("\n").split(",")
        print(cidx, ids[ridx], sep="\t")


def main2():
    fq_fn = sys.argv[1]
    txt_fn = sys.argv[2]

    clusters = {}
    for line in open(txt_fn):
        if line.startswith("read ID"):
            continue
        ridx, cidx = line.strip("\n").split(",")
        ridx, cidx = int(ridx), int(cidx)
        clusters[ridx] = clusters[ridx] + [cidx] if ridx in clusters else [cidx]
    print("Parsed clusters", file=sys.stderr)
    for i, read in enumerate(SeqIO.parse(fq_fn, "fastq")):
        if i % 50000 == 0:
            print(i, file=sys.stderr)
        if i in clusters:
            for cluster in clusters[i]:
                print(cluster, read.id, sep="\t")


if __name__ == "__main__":
    main2()
