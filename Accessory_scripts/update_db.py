#!/usr/bin/env python3

import argparse
import os
import sys

def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "-i",
        required = True,
        help="input dir"
        )
    p.add_argument(
        "-o",
        required = True,
        help="output dir"
        )

    return p.parse_args()

def read_alleles(path: str):
    locus = path.split("/")[-1].split(".")[0]
    if locus.lower() in {"neua", "neuah"}:
        locus = "neuA_neuAH"
    alleles = []
    with open(path) as fin:
        for line in fin.readlines()[1:]:
            seq, num = line.strip().split(",")
            alleles.append(f">{locus}_{num}\n{seq}\n")
    return alleles

def read_sbt(sbt_file: str, outfile: str):
    content = ["st\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA_neuAH\n"]
    with open(sbt_file) as fin:
        for line in fin.readlines()[1:]:
            content.append(line.replace(",", "\t"))
    
    with open(outfile, "w") as fout:
        fout.write("".join(content))
        

def main():
    args = parse_args()
    if not os.path.isdir(args.o) and not os.path.exists(args.o):
        os.makedirs(args.o)
    all_alleles = []
    for file in os.listdir(args.i):
        path = args.i + file
        if file == "sbt.csv":
            pass
            read_sbt(path, f"{args.o.rstrip('/')}/lpneumophila.txt")
        else:
            all_alleles += read_alleles(path)
    
    with open(f"{args.o.rstrip('/')}/all_loci.fasta", "w") as fout:
        fout.write("".join(all_alleles))
    

if __name__ == "__main__":
    main()
