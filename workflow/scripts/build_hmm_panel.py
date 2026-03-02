#!/usr/bin/env python3
import argparse
import os
import sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--osmo_kos", required=True, help="Ruta a osmo_kos.txt (1 KO por línea)")
    ap.add_argument("--profiles", required=True, help="Directorio con perfiles KOfam (*.hmm)")
    ap.add_argument("--out", required=True, help="Salida: HMM panel concatenado")
    args = ap.parse_args()

    # leer KOs
    kos = []
    with open(args.osmo_kos, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            kos.append(line)

    if not kos:
        sys.stderr.write("[ERROR] osmo_kos.txt está vacío.\n")
        sys.exit(1)

    missing = []
    paths = []
    for ko in kos:
        p = os.path.join(args.profiles, f"{ko}.hmm")
        if not os.path.exists(p):
            missing.append(p)
        else:
            paths.append(p)

    if missing:
        sys.stderr.write("[ERROR] Faltan perfiles HMM:\n" + "\n".join(missing) + "\n")
        sys.exit(2)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # concatenar
    with open(args.out, "w", encoding="utf-8") as out:
        for p in paths:
            with open(p, "r", encoding="utf-8", errors="ignore") as h:
                out.write(h.read())
                if not out.tell() or not str(out).endswith("\n"):
                    out.write("\n")

    print(f"[OK] Panel HMM creado: {args.out} (perfiles: {len(paths)})")

if __name__ == "__main__":
    main()
