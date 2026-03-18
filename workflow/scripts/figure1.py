import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Inputs
presence_file = "results/matrices/osmo_presence.tsv"
genomes_file = "configs/genomes.tsv"
output_file = "results/figures/figure1_heatmap_osmo_technical_v3.png"

# Load data
presence = pd.read_csv(presence_file, sep="\t")
genomes = pd.read_csv(genomes_file, sep="\t")

# Normalize column names
presence.columns = [c.strip() for c in presence.columns]
genomes.columns = [c.strip() for c in genomes.columns]

# Set index
presence = presence.set_index("accession")
genomes = genomes.set_index("Accession")

# Keep only shared genomes
shared = [acc for acc in genomes.index if acc in presence.index]
presence = presence.loc[shared]
genomes = genomes.loc[shared]

# Order genomes by lineage/group
group_order = ["A", "B", "C"]
genomes["Group"] = pd.Categorical(genomes["Group"], categories=group_order, ordered=True)
genomes = genomes.sort_values(["Group", "Note"])
presence = presence.loc[genomes.index]

# KO order by functional module
ko_modules = {
    "Proline": ["K00133", "K00931", "K01799", "K00147", "K11910"],
    "Ectoine": ["K06718", "K06719", "K06720", "K06721"],
    "Betaine": ["K00129", "K00130", "K03464", "K00965"],
    "Trehalose": ["K00697", "K05375", "K01231", "K01230"],
    "Transport": ["K05845", "K05846", "K05847", "K02000", "K02001", "K02002", "K03312", "K03313", "K03314"],
    "Ion homeostasis": ["K03315", "K03316", "K03317", "K03318", "K01546", "K01547", "K01548", "K03455"],
}

ordered_kos = []
module_boundaries = []
module_centers = []

start = 0
for module, kos in ko_modules.items():
    kos_present = [ko for ko in kos if ko in presence.columns]
    ordered_kos.extend(kos_present)
    end = start + len(kos_present)
    module_boundaries.append(end)
    module_centers.append((start + end - 1) / 2)
    start = end

presence = presence[ordered_kos]

# Convert to numeric safely
presence = presence.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

# Plot
fig, ax = plt.subplots(figsize=(16, 9))
cmap = ListedColormap(["#f2f2f2", "#123b7a"])  # 0 = absent, 1 = present

ax.imshow(
    presence.values,
    aspect="auto",
    interpolation="none",
    cmap=cmap,
    vmin=0,
    vmax=1
)

# X axis
ax.set_xticks(range(len(ordered_kos)))
ax.set_xticklabels(ordered_kos, rotation=90, fontsize=9)

# Y axis
ax.set_yticks([])
ax.set_ylabel("Genomes grouped by lineage", fontsize=11)
ax.set_xlabel("KEGG Ortholog grouped by functional module", fontsize=11)

# Vertical separators between modules
for b in module_boundaries[:-1]:
    ax.axvline(b - 0.5, color="black", linewidth=1.2)

# Horizontal separators between groups
group_sizes = genomes["Group"].value_counts().reindex(group_order).fillna(0).astype(int)
cum = 0
for g in group_order[:-1]:
    cum += group_sizes[g]
    ax.axhline(cum - 0.5, color="black", linewidth=1.2)

# Module labels on top
for center, module in zip(module_centers, ko_modules.keys()):
    ax.text(
        center,
        -0.9,
        module,
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold"
    )

# Clean frame
for spine in ax.spines.values():
    spine.set_visible(False)

# Layout and save
plt.tight_layout(rect=[0, 0, 1, 0.97])
os.makedirs(os.path.dirname(output_file), exist_ok=True)
plt.savefig(output_file, dpi=600, bbox_inches="tight")
plt.close()
