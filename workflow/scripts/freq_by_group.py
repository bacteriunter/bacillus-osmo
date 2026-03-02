#!/usr/bin/env python3
import pandas as pd
import os

MAT = "results/matrices/osmo_presence.tsv"
GEN = "configs/genomes.tsv"
OUT_FREQ  = "results/tables/osmo_freq_by_group.tsv"
OUT_COUNT = "results/tables/osmo_counts_by_group.tsv"

def main():
    mat = pd.read_csv(MAT, sep="\t")
    gen = pd.read_csv(GEN, sep="\t")

    gen = gen[["Accession","Group"]].copy()
    gen["Group"] = gen["Group"].astype(str).str.strip()

    # merge
    df = mat.merge(gen, left_on="accession", right_on="Accession", how="inner")
    if df.empty:
        raise SystemExit("[ERROR] No hubo match entre accessions de matriz y configs/genomes.tsv")

    kos = [c for c in mat.columns if c != "accession"]

    # conteos
    counts = df.groupby("Group")[kos].sum(numeric_only=True).reset_index()

    # denominadores (n genomas por grupo)
    denom = df.groupby("Group")["accession"].nunique().to_dict()

    # frecuencias
    freq = counts.copy()
    for g in freq["Group"]:
        n = denom.get(g, 0)
        if n == 0:
            freq.loc[freq["Group"]==g, kos] = 0.0
        else:
            freq.loc[freq["Group"]==g, kos] = (freq.loc[freq["Group"]==g, kos] / n) * 100.0

    os.makedirs("results/tables", exist_ok=True)
    counts.to_csv(OUT_COUNT, sep="\t", index=False)
    freq.to_csv(OUT_FREQ, sep="\t", index=False)

    print(f"[OK] counts: {OUT_COUNT}")
    print(f"[OK] freq:   {OUT_FREQ}")
    print("[INFO] N por grupo:", denom)

if __name__ == "__main__":
    main()
