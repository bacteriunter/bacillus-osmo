# bacillus-osmo (Snakemake)

Reproducible workflow for osmoadaptation gene profiling in Bacillus-related genomes.

## Install

```bash
conda env create -f envs/environment.yaml
conda activate bacillus-osmo
```

## Run

```bash
snakemake --cores 4
```

## Main Outputs

- results/figures/Figure1_Osmoadaptation_Heatmap.png
- results/figures/Figure2_Enrichment_Osmoadaptation.png
- results/tables/Table_S1_GeneCounts_ByGroup.tsv
- results/tables/Table_S2_GeneFrequency_ByGroup.tsv
- results/tables/Table_S3_Significant_KOs_FDR.tsv
