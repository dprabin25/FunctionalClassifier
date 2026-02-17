# FunctionalClassifier

Deterministic functional annotation of multiplex immune mediator panels using Gene Ontology (GO: Biological Process) terms with curated biological fallback rules.

This tool converts heterogeneous biomarker measurements (cytokines and chemokines) into interpretable biological mechanisms. 

# Functional Classes Produced
| Class                        | Interpretation                               |
| ---------------------------- | -------------------------------------------- |
| Immune cell activation       | Cytokine signaling & inflammatory activation |
| Immune cell positioning      | Chemotaxis & recruitment                     |
| Immune cell amplification    | Expansion & survival                         |
| Resolution                   | Regulatory / anti-inflammatory responses     |
| Tissue destruction           | Proteolysis & cell death                     |
| Tissue construction / repair | Healing & matrix organization                |
| Bone remodeling              | Osteoclast / osteoblast pathways             |


# Method Overview
Panel name → HGNC gene symbol mapping

GO Biological Process retrieval (org.Hs.eg.db)

Pattern-based biological classification

Curated correction for incomplete annotation (e.g., IFNL2)

Export long and unique mapping tables

No statistical modeling or machine learning is used.
Classification is deterministic and reproducible.

# Requirements
Tested on: R ≥ 4.2

Cran packages: dplyr

Bioconductor packages: AnnotationDbi, org.Hs.eg.db, GO.db

# Input files
No input files required

# Output files
1. CSV file 1: Full GO annotation (long format)
2. CSV file 2: Unique panel → functional class mapping

# Usage
Install the required packages and then run the required Rscript (Classifier.R).

 `install.packages("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "GO.db"
))`



# Limitations
1. Not enrichment analysis
2. Not pathway overrepresentation testing
3. Not predictive modeling
4. Designed for interpretation, not inference

# Citation
[1] Sayaka Miura, Prabin Dawadi, Ryan M Tobin, Jorge Frias-Lopez, Alpdogan Kantarci, Flavia Teles. Uncovering Periodontal Ecosystem Complexity with Sample Trees (2026). Under review. 


# Copyright 2026, Authors and University of Mississippi 
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
