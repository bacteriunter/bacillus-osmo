#!/usr/bin/env python3
import glob
import os
import re
from collections import defaultdict

# Este script asume exactamente estas rutas (las mismas del Snakefile/config):
TBL_GLOB = "results/hmmer/*/*.tbl"
OSMO_KOS = "data/hmm/osmo_kos.txt"
KO_LIST  = "data/kofam/ko_list"
OUT_MAT  = "results/matrices/osmo_presence.tsv"
OUT_LONG = "results/tables/osmo_hits_long.tsv"

def read_kos(path):
    kos=[]
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            kos.append(line)
    return kos

def read_kofam_thresholds(ko_list_path):
    """
    ko_list: columnas separadas por espacios/tab.
    Formato típico: KO  threshold  ...  (descripción)
    Tomamos:
      KO = col1
      threshold = col2 (float)
    """
    thr={}
    with open(ko_list_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            parts=re.split(r"\s+", line)
            if len(parts) < 2:
                continue
            ko=parts[0]
            try:
                thr[ko]=float(parts[1])
            except:
                continue
    return thr

def infer_accession_from_tbl(path):
    # results/hmmer/GCF_xxx/GCF_xxx.tbl
    base=os.path.basename(path)
    if base.endswith(".tbl"):
        return base[:-4]
    return base

def parse_tbl(tbl_path):
    """
    HMMER --tblout:
    target name  target acc  query name  query acc  E-value score bias ...
    En tus archivos, el KO está en la columna 3 (query name), ej: K00133
    y el bit-score en la columna 6.
    """
    hits=[]
    with open(tbl_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts=line.strip().split()
            if len(parts) < 6:
                continue
            target = parts[0]
            ko     = parts[2]
            evalue = parts[4]
            score  = parts[5]
            try:
                score_f=float(score)
            except:
                continue
            hits.append((target, ko, score_f, evalue))
    return hits

def main():
    kos = read_kos(OSMO_KOS)
    if not kos:
        raise SystemExit("[ERROR] data/hmm/osmo_kos.txt está vacío.")

    thr = read_kofam_thresholds(KO_LIST)
    # si algún KO no tiene umbral, lo dejamos como None y NO filtra por score
    # (pero lo reportamos igual en la tabla larga)
    tbls = sorted(glob.glob(TBL_GLOB))
    if not tbls:
        raise SystemExit("[ERROR] No se encontraron .tbl en results/hmmer/*/*.tbl")

    os.makedirs(os.path.dirname(OUT_MAT), exist_ok=True)
    os.makedirs(os.path.dirname(OUT_LONG), exist_ok=True)

    # presence[acc][ko]=0/1
    presence = { }
    long_rows = []

    for tbl in tbls:
        acc = infer_accession_from_tbl(tbl)
        presence.setdefault(acc, {k:0 for k in kos})

        for target, ko, score, evalue in parse_tbl(tbl):
            if ko not in presence[acc]:
                # ignorar KOs fuera del panel
                continue

            t = thr.get(ko, None)
            passed = (score >= t) if t is not None else True

            # tabla larga
            long_rows.append((acc, ko, target, score, evalue, t if t is not None else "", "True" if passed else "False"))

            # presencia si pasa threshold
            if passed:
                presence[acc][ko] = 1

    # escribir matriz
    with open(OUT_MAT, "w", encoding="utf-8") as out:
        out.write("\t".join(["accession"] + kos) + "\n")
        for acc in sorted(presence.keys()):
            out.write("\t".join([acc] + [str(presence[acc][k]) for k in kos]) + "\n")

    # escribir tabla larga
    with open(OUT_LONG, "w", encoding="utf-8") as out:
        out.write("\t".join(["accession","KO","target","score","evalue","threshold","pass"]) + "\n")
        for r in long_rows:
            out.write("\t".join(map(str,r)) + "\n")

    print(f"[OK] Matriz: {OUT_MAT} (genomas: {len(presence)}, KOs: {len(kos)})")
    print(f"[OK] Hits:   {OUT_LONG} (filas: {len(long_rows)})")

if __name__ == "__main__":
    main()
