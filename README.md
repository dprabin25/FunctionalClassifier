# FunctionalClassifier

This repository provides a small, focused data pipeline for turning per‑patient protein expression calls into a set of well-documented CSV tables. These tables quantify, for each patient, how many proteins increase or decrease within broad functional classes. Horizontal diverging bar plots are optional visualizations built on top of these tables; the **CSVs are the primary scientific outputs**.

---

## What This Repository Provides

- A standardized workflow for going from raw patient–protein calls to:
  - A merged, annotated protein table per patient.
  - A deduplicated table ensuring each protein is counted once per functional class.
  - A “Table 3”–style summary of counts by patient, functional class, and direction.
  - Per‑patient, plot‑ready bar‑count tables.
- Optional per‑patient bar plots that visualize these counts, but do not replace the underlying tables.

The design assumes that downstream users may want to replot or reanalyze the data, so every plot can be reconstructed from the exported CSVs.

---

## Conceptual Workflow

At a high level, the pipeline does the following:

1. **Start from two input tables**
   - A patient–protein table indicating which proteins increased or decreased for each patient.
   - A lookup table mapping each protein to a functional class (and, optionally, ontology terms or symbols).

2. **Merge and annotate**
   - Join the patient–protein calls with the functional lookup so that each observation carries a functional class and optional ontology information.

3. **Deduplicate within functional classes**
   - Ensure that a given protein is counted at most once per patient, direction (increase/decrease), and functional class.
   - Preserve any available ontology or symbol information as semicolon-separated lists.

4. **Summarize to functional-class counts**
   - For each patient, functional class, and direction, count how many unique proteins fall into that cell.
   - Record both the count and a signed version (negative for decreases) to support diverging bar plots.
   - Keep track of which proteins contribute to each bar as a separate “elements” list.

5. **Produce per-patient, plot-ready tables**
   - For each patient, create a compact table listing:
     - Functional class.
     - Direction (increase/decrease).
     - The number of proteins in that category.
   - These per‑patient tables are suitable for plotting in any environment (R, Python, etc.).

6. **Generate optional figures**
   - Create horizontal diverging bar plots that simply visualize the counts in the summary tables.
   - Save one figure per patient as PNG and PDF.
   - Treat these as convenience outputs; the authoritative results are the CSVs.

---

## Inputs

The pipeline expects two CSV files in a user-specified base directory:

1. **Patient–protein file (`PatData.csv`)**

   Contains per‑patient protein direction calls.

   - Columns:
     - `Pat`: patient identifier.
     - `Element`: protein or panel element name.
     - `Observed  direction`: direction of change (e.g., “Increase” or “Decrease”).

2. **Protein–function lookup (`ProteinGOAnnotation.csv`)**

   Maps proteins to functional classes and optional ontology metadata.

   - Required columns:
     - `Panel_Name`: protein or panel element name.
     - `Major_Functional_Class`: broad functional category.

   - Optional columns:
     - `TERM`: ontology term(s) or descriptions.
     - `Symbol`: gene/protein symbol.

The only strict requirement is that the protein identifiers in both files match so that the join is meaningful.

---

## Core CSV Outputs (Backbone Tables)

These outputs are the heart of the repository. The figure files are derived products.

### 1. `01_Raw_Merged_Proteins_With_TERMs.csv`

A fully merged table linking patient–protein calls to functional classes and optional ontology information.

- Typical columns:
  - `Pat`
  - `Element`
  - `Direction` (standardized from the input column)
  - `Major_Functional_Class`
  - `TERM` (if present)
  - `Symbol` (if present)

Use this file when you want to see exactly how each input row was annotated.

---

### 2. `02_Deduped_Proteins_UniquePerFunction.csv`

A deduplicated table ensuring that each protein is counted once per patient, direction, and functional class.

- Typical columns:
  - `Pat`
  - `Element`
  - `Direction`
  - `Major_Functional_Class`
  - `TERM` (aggregated, if present)
  - `Symbol` (aggregated, if present)

This is the appropriate table to use when you want to avoid within-class double counting of proteins.

---

### 3. `03_BarCounts_ByFunction_ByDirection.csv` (Table 3)

A summary table that provides the basis for all bar plots.

- Columns:
  - `Pat`: patient ID.
  - `Major_Functional_Class`: functional category.
  - `Direction`: “Increase” or “Decrease”.
  - `Count`: number of unique proteins in this category.
  - `SignedCount`: same number but negative for decreases (used only for plotting).
  - `Elements`: semicolon-separated list of the contributing proteins.

This is the main quantitative result: one row corresponds to one bar in the conceptual figure.

---

### 4. `04_PerPatient_BarCounts/{Patient}_BarCounts.csv`

A tidy, per‑patient file that is easy to feed into plotting tools.

- Columns:
  - `Major_Functional_Class`
  - `Direction`
  - `Count` (absolute)

This file is intentionally small and “plot-ready,” making it simple to reproduce or customize figures outside this repository.

---

## Figures (Derived Outputs)

For each patient, the pipeline optionally creates:

- `{Patient}_Figure.png`
- `{Patient}_Figure.pdf`

These are horizontal diverging bar plots showing, by functional class:

- Bars to the right (increased proteins).
- Bars to the left (decreased proteins).
- X‑axis labeled with absolute counts.

They are **not** required to use or interpret the data. Everything needed for analysis is contained in the CSVs listed above.


---

## Installation and Execution (Brief)

The technical details are deliberately minimal here.

1. Install Python (3.8+) and the standard scientific packages (`pandas`, `numpy`, `matplotlib`).
2. Place your two input files (Patient Data and Annotation Table) in your chosen base directory.
3. Run the main script (`FunctionalPlot.py`) from that directory.
4. Use the generated CSVs as your primary outputs; use the figures only as visual summaries.

---

# Reference
[1] Miura S, Dawadi P, Tobin RM, Frias-Lopez J, Kantarci A, Teles F. **Uncovering Periodontal Ecosystem Complexity with Sample Trees.**. 2026. Under review

# Copyright 2025, Authors and University of Mississippi
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
