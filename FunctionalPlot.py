# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 13:17:58 2026

@author: newfaculty
"""

# -*- coding: utf-8 -*-
"""
Correct counting + elements per bar
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re

def safe_filename(name):
    name = str(name)
    name = re.sub(r'[\\/*?:"<>|]', "_", name)
    name = name.replace(" ", "_")
    return name

# =================================================
# STYLE SETTINGS
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
# MERGE (raw trace)
# =================================================
need = ["Element","Major_Functional_Class","TERM","Symbol"]
func_cols = [c for c in need if c in func.columns]
func2 = func[func_cols]

merged_raw = cyto.merge(func2, on="Element", how="left")
merged_raw = merged_raw.dropna(subset=["Major_Functional_Class"])

raw_out = os.path.join(BASE_DIR, "01_Raw_Merged_Proteins_With_TERMs.csv")
merged_raw.to_csv(raw_out, index=False)
print("Saved:", raw_out)

# =================================================
# DEDUP LOGIC (THE FIX)
# Each element counted ONCE per functional class
# =================================================
group_keys = ["Pat","Element","Direction","Major_Functional_Class"]

agg_dict = {}
if "TERM" in merged_raw.columns:
    agg_dict["TERM"] = lambda x: "; ".join(sorted(set(x.dropna().astype(str))))
if "Symbol" in merged_raw.columns:
    agg_dict["Symbol"] = lambda x: "; ".join(sorted(set(x.dropna().astype(str))))

merged_dedup = merged_raw.groupby(group_keys, as_index=False).agg(agg_dict)

dedup_out = os.path.join(BASE_DIR, "02_Deduped_Proteins_UniquePerFunction.csv")
merged_dedup.to_csv(dedup_out, index=False)
print("Saved:", dedup_out)

# =================================================
# COUNTING (based ONLY on dedup table)
# =================================================
summary = (
    merged_dedup
    .groupby(["Pat","Major_Functional_Class","Direction"])
    .size()
    .reset_index(name="Count")
)

summary["SignedCount"] = np.where(summary["Direction"]=="Decrease",-summary["Count"],summary["Count"])

# =================================================
# TABLE 3 — FINAL RESULT WITH ELEMENT MEMBERS
# =================================================
bar_elements = (
    merged_dedup
    .groupby(["Pat","Major_Functional_Class","Direction"])["Element"]
    .apply(lambda x: "; ".join(sorted(set(x))))
    .reset_index(name="Elements")
)

table3 = summary.merge(bar_elements,on=["Pat","Major_Functional_Class","Direction"],how="left")

table3_out = os.path.join(BASE_DIR,"03_BarCounts_ByFunction_ByDirection.csv")
table3.sort_values(["Pat","Major_Functional_Class","Direction"]).to_csv(table3_out,index=False)
print("Saved:",table3_out)

# =================================================
# SCALE
# =================================================
max_abs = summary["SignedCount"].abs().max()
global_max = 1 if pd.isna(max_abs) or max_abs==0 else int(np.ceil(max_abs))

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
# PER PATIENT EXPORT + PLOT
# =================================================
per_patient_dir=os.path.join(BASE_DIR,"04_PerPatient_BarCounts")
os.makedirs(per_patient_dir,exist_ok=True)

for patient in summary["Pat"].unique():

    safe_patient=safe_filename(patient)
    dfp=summary[summary["Pat"]==patient]

    inc=dfp[dfp["Direction"]=="Increase"].set_index("Major_Functional_Class")["Count"].reindex(ALL_FUNCTIONS,fill_value=0)
    dec=dfp[dfp["Direction"]=="Decrease"].set_index("Major_Functional_Class")["Count"].reindex(ALL_FUNCTIONS,fill_value=0)*-1

    plot_df=pd.DataFrame({"Decrease":dec,"Increase":inc}).loc[ALL_FUNCTIONS]

    long_df=plot_df.reset_index().melt(id_vars="Major_Functional_Class",var_name="Direction",value_name="Count")
    long_df=long_df[long_df["Count"]!=0]
    long_df["Count"]=long_df["Count"].abs()
    long_df.to_csv(os.path.join(per_patient_dir,f"{safe_patient}_BarCounts.csv"),index=False)

    # plot
    fig,ax=plt.subplots(figsize=(7.2,3.6))
    plot_df.plot(kind="barh",ax=ax,width=0.45,color=[COLOR_DECREASE,COLOR_INCREASE],edgecolor="none")

    ax.set_xlim(-X_LIMIT,X_LIMIT)
    ax.set_xticks(label_ticks)
    ax.set_xticklabels([str(abs(int(x))) for x in label_ticks])
    ax.axvline(0,color="black",linewidth=1.2)

    ax.set_xlabel("Number of Proteins")
    ax.set_ylabel("")
    ax.set_title(patient,weight="bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.legend(title="Direction",frameon=False,loc="lower right")

    plt.tight_layout()
    plt.savefig(os.path.join(BASE_DIR,f"{safe_patient}_Figure.png"),dpi=600,bbox_inches="tight")
    plt.savefig(os.path.join(BASE_DIR,f"{safe_patient}_Figure.pdf"),bbox_inches="tight")
    plt.close(fig)

print("All figures and corrected tables exported successfully.")
