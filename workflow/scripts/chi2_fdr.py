#!/usr/bin/env python3
import os
import pandas as pd

MAT = "results/matrices/osmo_presence.tsv"
GEN = "configs/genomes.tsv"
OUT = "results/tables/osmo_chi2_results.tsv"

def bh_fdr(pvals):
    # Benjamini-Hochberg
    m = len(pvals)
    order = sorted(range(m), key=lambda i: (pvals[i] if pvals[i]==pvals[i] else 1.0))
    adj = [None]*m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        p = pvals[i]
        if p != p:  # NaN
            adj[i] = float("nan")
            continue
        q = min(prev, p * m / (m - rank + 1))
        adj[i] = q
        prev = q
    return adj

def main():
    try:
        from scipy.stats import chi2_contingency
    except Exception as e:
        raise SystemExit("[ERROR] Falta scipy. Instala con: conda install -c conda-forge scipy\n" + str(e))

    mat = pd.read_csv(MAT, sep="\t")
    gen = pd.read_csv(GEN, sep="\t")[["Accession","Group"]]
    df = mat.merge(gen, left_on="accession", right_on="Accession", how="inner")
    if df.empty:
        raise SystemExit("[ERROR] No hubo match entre matriz y genomes.tsv")

    kos = [c for c in mat.columns if c != "accession"]
    groups = ["A","B","C"]

    rows = []
    for ko in kos:
        # tabla 3x2: filas=grupos, cols=[presente, ausente]
        table = []
        for g in groups:
            sub = df[df["Group"] == g][ko]
            present = int(sub.sum())
            total = int(sub.shape[0])
            absent = total - present
            table.append([present, absent])

        # si no hay variación, chi2 no aplica bien
        total_present = sum(r[0] for r in table)
        total_absent  = sum(r[1] for r in table)
        if total_present == 0 or total_absent == 0:
            chi2 = float("nan")
            p = float("nan")
        else:
            chi2, p, dof, exp = chi2_contingency(table)

        rows.append({
            "KO": ko,
            "chi2": chi2,
            "p_value": p
        })

    res = pd.DataFrame(rows)
    res["p_adj"] = bh_fdr(res["p_value"].tolist())
    res["significant"] = res["p_adj"].apply(lambda x: False if pd.isna(x) else x < 0.05)

    os.makedirs("results/tables", exist_ok=True)
    res.to_csv(OUT, sep="\t", index=False)

    print(f"[OK] {OUT} (KOs: {len(res)})")

if __name__ == "__main__":
    main()
