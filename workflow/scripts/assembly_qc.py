#!/usr/bin/env python3
import sys
from collections import Counter

def read_fasta_lengths_and_gc(path: str):
    lengths = []
    gc = 0
    at = 0
    n = 0
    cur = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur > 0:
                    lengths.append(cur)
                    cur = 0
                continue
            s = line.upper()
            cur += len(s)
            gc += s.count("G") + s.count("C")
            at += s.count("A") + s.count("T")
            n  += s.count("N")
    if cur > 0:
        lengths.append(cur)

    bp = sum(lengths)
    gc_pct = (gc / (gc + at) * 100.0) if (gc + at) > 0 else 0.0
    return lengths, bp, gc_pct

def n50(lengths):
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    half = total / 2
    c = 0
    for L in lengths:
        c += L
        if c >= half:
            return L
    return 0

def infer_accession(path: str):
    # Espera .../GCF_xxx/GCF_xxx.fna
    base = path.split("/")[-1]
    if base.endswith(".fna"):
        return base[:-4]
    return base

def main():
    fna_paths = sys.argv[1:]
    if not fna_paths:
        sys.stderr.write("Usage: assembly_qc.py <genome1.fna> <genome2.fna> ...\n")
        sys.exit(1)

    out = sys.stdout
    out.write("\t".join(["accession","path","contigs","bp","N50","gc_percent"]) + "\n")

    for p in fna_paths:
        acc = infer_accession(p)
        lengths, bp, gc_pct = read_fasta_lengths_and_gc(p)
        out.write("\t".join([
            acc,
            p,
            str(len(lengths)),
            str(bp),
            str(n50(lengths)),
            f"{gc_pct:.2f}"
        ]) + "\n")

if __name__ == "__main__":
    main()
