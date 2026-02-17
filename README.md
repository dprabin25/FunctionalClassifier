# FunctionalClassifier

GO-based protein classification and longitudinal immune response visualization

This pipeline converts multiplex immune mediator measurements (cytokines, chemokines, growth factors, and MMPs) into interpretable biological mechanisms and visualizes patient-specific immune response patterns.

---

## Overview

The workflow has two independent stages:

| Stage  | Script              | Purpose                                                                                           |
| ------ | ------------------- | ------------------------------------------------------------------------------------------------- |
| Step 1 | `Classifier.R`      | Generates biological functional annotation reference using Gene Ontology (GO: Biological Process) |
| Step 2 | `FunctionalPlot.py` | Visualizes longitudinal patient immune response                                                   |

The R step prepares a reusable biological reference table.
The Python step performs the actual analysis on patient data.

---

## Biological Concept

Proteins are grouped into biological functional processes using GO Biological Process (GO-BP) terms and then visualized based on directional change.

```
Increase  → Right side
Decrease  → Left side
```

This produces a symmetric biological response profile for each patient.

---

## Step 1 — Functional Annotation (R)

No input file is required.

The R script generates the annotation reference automatically.

### Output directory

```
Working dir/
 ├── UniqueFunctional.csv
 └── ProteinGoAnnotation.csv
```

| File                    | Description                                  |
| ----------------------- | -------------------------------------------- |
| UniqueFunctional.csv    | Unique panel → functional class mapping      |
| ProteinGoAnnotation.csv | Full GO annotation with all biological terms |

Functional classes produced:

| Class                        | Interpretation                               |
| ---------------------------- | -------------------------------------------- |
| Immune cell activation       | Cytokine signaling & inflammatory activation |
| Immune cell positioning      | Chemotaxis & recruitment                     |
| Immune cell amplification    | Expansion & survival                         |
| Resolution                   | Anti-inflammatory regulation                 |
| Tissue destruction           | Proteolysis & cell death                     |
| Tissue construction / repair | Healing & matrix organization                |
| Bone remodeling              | Osteoclast / osteoblast pathways             |

The Python script counts proteins based on these functional groups.

---

### R Requirements

Tested on: **R ≥ 4.2**

CRAN packages:

```
dplyr
```

Bioconductor packages:

```
AnnotationDbi
org.Hs.eg.db
GO.db
```

### Installation

Run inside R:

```r
install.packages("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "GO.db"
))
```

### Usage

Run:

```
Rscript Classifier.R
```

---

## Step 2 — Patient Data Mapping and Visualization (Python)

This script converts longitudinal biomarker changes into functional immune response profiles using the annotation table from Step 1.

The Python script **does not perform GO analysis**.
It only maps and counts biological responses.

---

### Required Input Files

#### 1) Patient Data (`PatData.csv`)

| Pat      | Element | Observed direction |
| -------- | ------- | ------------------ |
| Patient1 | IL-6    | Increase           |
| Patient1 | MMP-8   | Decrease           |
| Patient2 | MCP-1   | Increase           |

Rules:

* `Element` must exactly match `Panel_Name` in annotation file
* Each row = one observation
* Any two directional labels can be used (Increase/Decrease recommended)

Examples of allowed labels:

```
Increase / Decrease
Up / Down
Elevated / Reduced
Treatment / Baseline
Positive / Negative
```

---

#### 2) Functional Annotation (`ProteinGoAnnotation.csv`)

Generated from `Classifier.R`

Required columns:

| Panel_Name | SYMBOL | TERM | Major_Functional_Class |

This acts as the biological dictionary linking proteins to functional mechanisms.

---

### Python Requirements

Tested on: **Python ≥ 3.9**

Core libraries:

```
pandas
numpy
matplotlib
```

Install:

```bash
pip install pandas numpy matplotlib
```

---

### Usage

```bash
python FunctionalPlot.py
```

---

## Output

One figure is generated per patient in the same directory.

| File                   | Description                       |
| ---------------------- | --------------------------------- |
| PatientName_Figure.png | High-resolution image (600 dpi)   |
| PatientName_Figure.pdf | Vector publication-quality figure |

---

## Interpretation

Right side → Dominant biological processes
Left side → Suppressed biological processes

Each bar represents:

```
Number of proteins belonging to a functional class
```

---

## Citation

Miura S., Dawadi P., Tobin R.M., Frias-Lopez J., Kantarci A., Teles F.
**Uncovering Periodontal Ecosystem Complexity with Sample Trees (2026).**
Under review.

---

## License

BSD 3-Clause License © 2026
Authors and University of Mississippi

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the BSD-3 conditions are met.
