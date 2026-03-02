#!/usr/bin/env python3
import os
import pandas as pd

def main():
    chi_path  = "results/tables/osmo_chi2_results.tsv"
    freq_path = "results/tables/osmo_freq_by_group.tsv"
    out_path  = "results/figures/figure2_enrichment_significant_KOs_final.png"

    chi  = pd.read_csv(chi_path, sep="\t")
    freq = pd.read_csv(freq_path, sep="\t")

    # Filtrar significativos (q<0.05)
    sig = chi[chi["significant"] == True].copy()
    if sig.empty:
        raise SystemExit("[ERROR] No hay KOs significativos (q<0.05). No se puede generar figura 2.")

    # Preparar datos en formato largo
    sig_kos = sig["KO"].tolist()
    f = freq[freq["Group"].isin(["A","B","C"])].copy()
    f_long = f.melt(id_vars=["Group"], value_vars=sig_kos, var_name="KO", value_name="freq")
    # unimos q-values
    f_long = f_long.merge(sig[["KO","p_adj"]], on="KO", how="left")

    # Orden por q-value ascendente
    ko_order = sig.sort_values("p_adj")["KO"].tolist()

    # Plot
    import matplotlib.pyplot as plt
    import seaborn as sns

    os.makedirs("results/figures", exist_ok=True)

    plt.figure(figsize=(10, max(4, 0.35*len(ko_order))))
    sns.set_style("whitegrid")

    ax = sns.barplot(
        data=f_long,
        y="KO",
        x="freq",
        hue="Group",
        order=ko_order,
        palette={"A":"#4C78A8","B":"#F58518","C":"#54A24B"}
    )

    ax.set_xlabel("Frecuencia de presencia (%)")
    ax.set_ylabel("KO (significativos, FDR q<0.05)")
    ax.legend(title="Grupo", loc="lower right")

    # No imprimimos q-values en la figura (como ya decidiste).
    # Los q-values van en la tabla S3.
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"[OK] Figura 2 guardada en {out_path}")
    print("[INFO] KOs significativos:", ", ".join(ko_order))

if __name__ == "__main__":
    main()
