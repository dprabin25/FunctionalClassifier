# ============================
# Functional Grouping for Cytokines / Chemokines / Growth factors / MMPs
# GO BP–based + curated fallback for IFNL2
# CSV OUTPUT ONLY
# Outputs saved to: C:/PlotRevSample
# ============================

# ----------------------------
# OUTPUT DIRECTORY
# ----------------------------
OUTDIR <- "C:/PlotRevSample"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# ----------------------------
# LIBRARIES (LOAD ONLY)
# ----------------------------
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(dplyr)
})

# ----------------------------
# STEP 1 — PANEL NAME → GENE SYMBOL MAP
# ----------------------------
panel_map <- data.frame(
  Panel_Name = c(
    "HSP70","MMP-8","Eotaxin-2","BCA-1","MCP-4","I-309","IL-16","6CKine","LIF","SCF","TSLP",
    "IL-21","IL-23","TRAIL","SDF-1a.b","ENA-78","MIP-1d","IL-28A","EGF","Eotaxin","TGF-A",
    "G-CSF","Flt-3L","GM-CSF","Fractalkine","IFN-A","IFN-y","GRO","IL-10","MCP-3",
    "IL-12p40","MDC","IL-12p70","PDGF-AA","IL-13","IL-15","sCD40L","IL-17A","IL-1RA",
    "IL-1A","IL-1B","IL-2","IL-6","IL-7","IL-8","IP-10","MCP-1","MIP-1A","MIP-1B",
    "RANTES","TNF-A","VEGF","MMP-1","MMP-2","MMP-9","MMP-10","RANKL","GDF-15","DKK1",
    "Periostin","TRAP5","OPG","YKL40"
  ),
  SYMBOL = c(
    "HSPA1A","MMP8","CCL24","CXCL13","CCL13","CCL1","IL16","CCL21","LIF","KITLG","TSLP",
    "IL21","IL23A","TNFSF10","CXCL12","CXCL5","CCL15","IFNL2","EGF","CCL11","TGFA",
    "CSF3","FLT3LG","CSF2","CX3CL1","IFNA1","IFNG","CXCL1","IL10","CCL7",
    "IL12B","CCL22","IL12A","PDGFA","IL13","IL15","CD40LG","IL17A","IL1RN",
    "IL1A","IL1B","IL2","IL6","IL7","CXCL8","CXCL10","CCL2","CCL3","CCL4",
    "CCL5","TNF","VEGFA","MMP1","MMP2","MMP9","MMP10","TNFSF11","GDF15","DKK1",
    "POSTN","ACP5","TNFRSF11B","CHI3L1"
  ),
  stringsAsFactors = FALSE
)

genes <- unique(panel_map$SYMBOL)

# ----------------------------
# STEP 2 — GET GO BIOLOGICAL PROCESS TERMS
# ----------------------------
go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes,
  columns = c("GO", "ONTOLOGY"),
  keytype = "SYMBOL"
)

bp <- go %>%
  filter(ONTOLOGY == "BP") %>%
  mutate(TERM = Term(GOTERM[GO])) %>%
  filter(!is.na(TERM)) %>%
  inner_join(panel_map, by = "SYMBOL")

# ----------------------------
# STEP 3 — FUNCTIONAL CLASS ASSIGNMENT
# ----------------------------
assign_major_class <- function(term) {
  t <- tolower(term)
  
  if (grepl(
    "activation|cytokine production|immune response activation|
     positive regulation of immune response|innate immune response|jak-stat",
    t
  ))
    return("Immune cell activation")
  
  if (grepl("chemotaxis|migration|recruitment|homing|positioning", t))
    return("Immune cell positioning")
  
  if (grepl("proliferation|expansion|clonal|amplification|survival|hematopoiesis", t))
    return("Immune cell amplification")
  
  if (grepl("negative regulation|tolerance|resolution|anti-inflammatory|immune suppression", t))
    return("Resolution")
  
  if (grepl("matrix degradation|bone resorption|apoptotic|cell death|tissue destruction|cartilage", t))
    return("Tissue destruction")
  
  if (grepl(
    "epithelial|wound healing|tissue repair|extracellular matrix organization|
     fibrosis|angiogenesis|regeneration",
    t
  ))
    return("Tissue construction / repair")
  
  if (grepl("bone remodeling|osteoclast|osteoblast|mineralization", t))
    return("Bone remodeling")
  
  return(NA)
}

bp$Major_Functional_Class <- sapply(bp$TERM, assign_major_class)

# ----------------------------
# STEP 3.5 — CURATED BIOLOGICAL FALLBACK (CRITICAL FIX)
# IFNL2 has innate immune BP via IBA evidence (not in org.Hs.eg.db)
# ----------------------------
bp <- bp %>%
  mutate(
    Major_Functional_Class = ifelse(
      SYMBOL == "IFNL2" & is.na(Major_Functional_Class),
      "Immune cell activation",
      Major_Functional_Class
    )
  )

# ----------------------------
# STEP 4 — LONG FORMAT OUTPUT
# ----------------------------
classified_long <- bp %>%
  filter(!is.na(Major_Functional_Class)) %>%
  distinct(Panel_Name, SYMBOL, TERM, Major_Functional_Class)

# ----------------------------
# STEP 5 — SAVE LONG CSV
# ----------------------------
write.csv(
  classified_long,
  file.path(OUTDIR, "Updatedtry.csv"),
  row.names = FALSE
)

# ----------------------------
# STEP 6 — UNIQUE PANEL × FUNCTION CLASS
# ----------------------------
classified_unique <- classified_long %>%
  distinct(Panel_Name, Major_Functional_Class)

# ----------------------------
# STEP 7 — SAVE UNIQUE CSV
# ----------------------------
write.csv(
  classified_unique,
  file.path(OUTDIR, "abc_Unique.csv"),
  row.names = FALSE
)

# ----------------------------
# DONE
# ----------------------------
cat("\nFILES SAVED TO:", OUTDIR, "\n")
cat(" - Updatedtry.csv (all GO BP terms)\n")
cat(" - abc_Unique.csv (unique functional classes)\n")
