import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re

def safe_filename(name):
    """Remove characters not allowed in Windows filenames."""
    name = str(name)
    name = re.sub(r'[\\/*?:"<>|]', "_", name)
    name = name.replace(" ", "_")
    return name

# =================================================
# STYLE SETTINGS (Journal quality)
# =================================================
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "axes.linewidth": 1.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# =================================================
# INPUT & OUTPUT DIRECTORY
# =================================================
BASE_DIR = r"C:\PlotRevSample\ProteinsForPlots"

CYTOKINE_FILE = os.path.join(BASE_DIR, "PatData.csv")
FUNCTION_FILE = os.path.join(BASE_DIR, "ProteinGoAnnotation.csv")

# =================================================
# LOAD DATA
# =================================================
cyto = pd.read_csv(CYTOKINE_FILE)
func = pd.read_csv(FUNCTION_FILE)

cyto = cyto.rename(columns={"Observed  direction": "Direction"})
func = func.rename(columns={"Panel_Name": "Element"})

merged = cyto.merge(func, on="Element", how="left")

# =================================================
# FIXED ORDER
# =================================================
ALL_FUNCTIONS = [
 "Immune cell activation",
 "Immune cell amplification",
 "Immune cell positioning",
 "Tissue destruction",
 "Tissue construction / repair",
 "Resolution",
 "Bone remodeling"
]

summary = (
    merged
    .groupby(["Pat", "Major_Functional_Class", "Direction"])
    .size()
    .reset_index(name="Count")
)

# =================================================
# GLOBAL SCALE
# =================================================
summary["SignedCount"] = np.where(
    summary["Direction"] == "Decrease",
    -summary["Count"],
    summary["Count"]
)

global_max = int(np.ceil(summary["SignedCount"].abs().max()))

LABEL_STEP = 50
GRID_STEP = 5

X_LIMIT = int(np.ceil(global_max / LABEL_STEP) * LABEL_STEP)

label_ticks = np.arange(-X_LIMIT, X_LIMIT + LABEL_STEP, LABEL_STEP)
grid_ticks = np.arange(-X_LIMIT, X_LIMIT + GRID_STEP, GRID_STEP)

# =================================================
# COLORS (color-blind safe)
# =================================================
COLOR_INCREASE = "#0072B2"
COLOR_DECREASE = "#D55E00"

# =================================================
# PLOT PER PATIENT
# =================================================
for patient in summary["Pat"].unique():

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

    plot_df = pd.DataFrame({
        "Decrease": dec,
        "Increase": inc
    }).loc[ALL_FUNCTIONS]

    # =================================================
    # FIGURE
    # =================================================
    fig, ax = plt.subplots(figsize=(7.2, 3.6))

    plot_df.plot(
        kind="barh",
        ax=ax,
        width=0.45,
        color=[COLOR_DECREASE, COLOR_INCREASE],
        edgecolor="none",
        legend=False   # <<< prevents legend creation
    )

    # ----- LOCK SCALE -----
    ax.set_autoscale_on(False)
    ax.set_xlim(-X_LIMIT, X_LIMIT)
    ax.margins(y=0.02)

    # ----- ticks -----
    ax.set_xticks(label_ticks)
    ax.set_xticklabels([str(abs(int(x))) for x in label_ticks])

    # ----- center line -----
    ax.axvline(0, color="black", linewidth=1.2)

    # ----- neutral zone -----
    ax.axvspan(-GRID_STEP/2, GRID_STEP/2, color="black", alpha=0.00, zorder=0)

    # ----- labels -----
    ax.set_xlabel("Number of Proteins")
    ax.set_ylabel("")
    ax.set_title(patient, weight="bold")

    # ----- clean spines -----
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.tick_params(axis="y", pad=10)

    plt.tight_layout()

    # =================================================
    # SAVE
    # =================================================
    safe_patient = safe_filename(patient)
    png_out = os.path.join(BASE_DIR, f"{safe_patient}_Figure.png")
    pdf_out = os.path.join(BASE_DIR, f"{safe_patient}_Figure.pdf")

    plt.savefig(png_out, dpi=600, bbox_inches="tight")
    plt.savefig(pdf_out, bbox_inches="tight")
    plt.close(fig)

print("Publication-quality figures exported (NO LEGENDS).")

