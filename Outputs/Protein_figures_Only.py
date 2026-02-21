# -*- coding: utf-8 -*-
"""
Publication-quality 3x3 protein functional class figure
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# =================================================
# STYLE
# =================================================
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# =================================================
# INPUT
# =================================================
BASE_DIR = r"C:\PlotRevSample\ProteinsForPlots"
CYTOKINE_FILE = os.path.join(BASE_DIR, "Cytokines.csv")
FUNCTION_FILE = os.path.join(BASE_DIR, "UPdatedtry.csv")

cyto = pd.read_csv(CYTOKINE_FILE)
func = pd.read_csv(FUNCTION_FILE)

cyto = cyto.rename(columns={"Observed  direction": "Direction"})
func = func.rename(columns={"Panel_Name": "Element"})
cyto["Direction"] = cyto["Direction"].astype(str).str.strip().str.capitalize()

# =================================================
# MERGE
# =================================================
need = ["Element","Major_Functional_Class"]
func_cols = [c for c in need if c in func.columns]
func2 = func[func_cols]

merged_raw = cyto.merge(func2, on="Element", how="left")
merged_raw = merged_raw.dropna(subset=["Major_Functional_Class"])

# =================================================
# DEDUP
# =================================================
merged_dedup = (
    merged_raw
    .groupby(["Pat","Element","Direction","Major_Functional_Class"], as_index=False)
    .size()
)

# =================================================
# COUNT
# =================================================
summary = (
    merged_dedup
    .groupby(["Pat","Major_Functional_Class","Direction"])
    .size()
    .reset_index(name="Count")
)

summary["SignedCount"] = np.where(
    summary["Direction"]=="Decrease",
    -summary["Count"],
    summary["Count"]
)

# =================================================
# SCALE
# =================================================
max_abs = summary["SignedCount"].abs().max()
global_max = 1 if pd.isna(max_abs) else int(np.ceil(max_abs))

LABEL_STEP = 30
X_LIMIT = max(LABEL_STEP,int(np.ceil(global_max/LABEL_STEP)*LABEL_STEP))
label_ticks = np.arange(-X_LIMIT,X_LIMIT+LABEL_STEP,LABEL_STEP)

COLOR_INCREASE="#0072B2"
COLOR_DECREASE="#D55E00"

ALL_FUNCTIONS=[
"Immune cell activation",
"Immune cell amplification",
"Immune cell positioning",
"Tissue destruction",
"Tissue construction / repair",
"Resolution",
"Bone remodeling"
]

# =================================================
# 3x3 FIGURE
# =================================================
patients = sorted(summary["Pat"].unique(), reverse=True)

fig, axes = plt.subplots(3, 3, figsize=(14, 13))
axes = axes.flatten()

for i, patient in enumerate(patients[:9]):

    ax = axes[i]
    dfp = summary[summary["Pat"] == patient]

    inc = (
        dfp[dfp["Direction"] == "Increase"]
        .set_index("Major_Functional_Class")["Count"]
        .reindex(ALL_FUNCTIONS, fill_value=0)
    )

    dec = (
        dfp[dfp["Direction"] == "Decrease"]
        .set_index("Major_Functional_Class")["Count"]
        .reindex(ALL_FUNCTIONS, fill_value=0) * -1
    )

    plot_df = pd.DataFrame({"Decrease": dec, "Increase": inc}).loc[ALL_FUNCTIONS]

    plot_df.plot(
        kind="barh",
        ax=ax,
        width=0.55,
        color=[COLOR_DECREASE, COLOR_INCREASE],
        edgecolor="none",
        legend=False
    )

    # Zero line
    ax.axvline(0,color="black",linewidth=0.8)

    # Limits
    ax.set_xlim(-X_LIMIT, X_LIMIT)
    ax.set_xticks(label_ticks)

    # Only bottom row shows numbers
    if i < 6:
        ax.set_xticklabels([])
    else:
        ax.set_xticklabels([str(abs(int(x))) for x in label_ticks])

    # Only left column shows functional names
    if i % 3 != 0:
        ax.set_yticklabels([])

    # remove labels
    ax.set_xlabel("")
    ax.set_ylabel("")

    # minimal spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # subtle title
    ax.set_title(str(patient), pad=4)

# hide empty panels
for j in range(len(patients),9):
    axes[j].axis("off")

# =================================================
# GLOBAL LABEL
# =================================================
fig.supxlabel("Number of proteins", y=0.035, fontsize=12)

plt.subplots_adjust(
    left=0.20,
    right=0.98,
    top=0.96,
    bottom=0.08,
    wspace=0.15,
    hspace=0.22
)

plt.savefig(os.path.join(BASE_DIR,"Combined_3x3_Publication.png"),dpi=600)
plt.close()

print("Publication figure exported")