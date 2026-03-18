#!/usr/bin/env python3
import os
import pandas as pd

def main():
    chi_path = "results/tables/osmo_chi2_results.tsv"
    freq_path = "results/tables/osmo_freq_by_group.tsv"
    out_path = "results/tables/table_S3_significant_KOs_qvalues.tsv"

    # Leer tablas
    chi = pd.read_csv(chi_path, sep="\t")
    freq = pd.read_csv(freq_path, sep="\t")

    # Limpiar nombres de columnas
    chi.columns = [c.strip() for c in chi.columns]
    freq.columns = [c.strip() for c in freq.columns]

    # Asegurar tipos numéricos
    chi["p_value"] = pd.to_numeric(chi["p_value"], errors="coerce")
    chi["p_adj"] = pd.to_numeric(chi["p_adj"], errors="coerce")

    # Filtrar KOs significativos
    sig = chi[chi["significant"] == True].copy()

    if sig.empty:
        raise SystemExit("[ERROR] No significant KOs found (q<0.05). Cannot generate Table S3.")

    # Pasar frecuencias a formato largo
    freq_long = freq.melt(id_vars=["Group"], var_name="KO", value_name="frequency_percent")

    # Pivotear para tener una columna por grupo
    freq_wide = freq_long.pivot(index="KO", columns="Group", values="frequency_percent").reset_index()

    # Renombrar columnas
    rename_map = {}
    if "A" in freq_wide.columns:
        rename_map["A"] = "Group_A_frequency_percent"
    if "B" in freq_wide.columns:
        rename_map["B"] = "Group_B_frequency_percent"
    if "C" in freq_wide.columns:
        rename_map["C"] = "Group_C_frequency_percent"
    freq_wide = freq_wide.rename(columns=rename_map)

    # Unir con resultados estadísticos
    out = sig.merge(freq_wide, on="KO", how="left")

    # Ordenar por q-value ascendente
    out = out.sort_values("p_adj", ascending=True)

    # Reordenar columnas si existen
    desired_cols = [
        "KO",
        "chi2",
        "p_value",
        "p_adj",
        "significant",
        "Group_A_frequency_percent",
        "Group_B_frequency_percent",
        "Group_C_frequency_percent",
    ]
    out = out[[c for c in desired_cols if c in out.columns]]

    # Crear carpeta y guardar
    os.makedirs("results/tables", exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    print(f"[OK] Table S3 saved to {out_path}")
    print(f"[INFO] Significant KOs: {len(out)}")

if __name__ == "__main__":
    main()
