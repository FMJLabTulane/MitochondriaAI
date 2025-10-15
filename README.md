# MitochondriaAI: Sex-Based Analysis of Mitochondrial Function in Pancreatic Beta Cells

## Author
**Fahd Qadir (Dragonmasterx87)**  
FMJ Lab, Tulane University School of Medicine

**Date**: November 16, 2022  
**R Version**: 4.2.1 (2019-12-12) 'Funny-Looking Kid'

---

## Overview

This repository contains the comprehensive analytical pipeline for investigating sex-based differences in mitochondrial function and quality in pancreatic beta cells from individuals with and without type 2 diabetes (T2D). The analysis integrates single-cell RNA sequencing data with machine learning approaches to identify mitochondrial signatures associated with disease progression.

### Key Features

- **Single-cell RNA-seq analysis** of pancreatic beta cells
- **Sex-stratified differential expression** analysis
- **Mitochondrial quality index (MQI)** calculation and prediction using Kolmogorov-Arnold Networks (KAN)
- **PINK1-mediated mitophagy** pathway analysis
- **Algorithme.AI disease scoring** integration and validation
- **DHT treatment effects** on beta cell mitochondrial function
- **External dataset integration** (GSE217775 PINK1 KO mice)
- Comprehensive visualization suite (UMAP, heatmaps, Sankey diagrams, volcano plots)

---

## System Requirements

### R Environment
- **R version**: 4.2.1 or higher
- **RStudio**: Recommended
- **Rtools**: Compatible version with your R installation
- **Python**: conda environment via reticulate

### Operating System
- Windows (primary development)
- Linux/Mac (compatible with path adjustments)

---

## Installation

### 1. CRAN Packages

```r
install.packages(c(
  'ggplot2', 'cowplot', 'Matrix', 'ggridges', 'ggrepel', 'dplyr', 
  'plotly', 'clustree', 'patchwork', 'future', 'devtools', 'rlang', 
  'pROC', 'harmony', 'SoupX', 'tidyverse', 'viridis', 'circlize', 
  'scCustomize', 'archive', 'R.utils', 'qs', 'fastmap', 'torch',
  'readr', 'ggpubr', 'pheatmap', 'gridExtra', 'broom', 'scales', 'ggbreak'
))
```

### 2. Bioconductor Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.18")

BiocManager::install(c(
  'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
  'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
  'SummarizedExperiment', 'batchelor', 'HDF5Array',
  'terra', 'ggrastr', 'EnhancedVolcano', 'DoubletFinder', 
  'glmGamPoi', 'GOSemSim', 'org.Hs.eg.db', 'AnnotationHub',
  'GenomeInfoDb', 'MeSHDbi', 'clusterProfiler', 'dittoSeq', 
  'escape', 'ComplexHeatmap', 'DropletUtils', 'Nebulosa', 
  'hdf5r', 'scDblFinder', 'JASPAR2020', 'TFBSTools', 
  'motifmatchr', 'chromVAR', 'EnrichmentBrowser',
  'BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86',
  'GenomicRanges', 'Gviz', 'rtracklayer', 'DOSE'
))
```

### 3. GitHub Packages (Version-Specific)

```r
# Seurat ecosystem (v4.3.0.1 for Pando compatibility)
remotes::install_version(package = 'Seurat', version = '4.3.0.1')
devtools::install_version("SeuratObject", version = "4.1.4")
devtools::install_version("Signac", version = "1.11.0")

# Development versions
devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)

# Additional tools
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("gaospecial/ggVennDiagram")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
devtools::install_github('quadbiolab/Pando')
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')
```

### 4. Verify Installation

```r
# Check package versions
packageVersion("clusterProfiler")
packageVersion("Seurat")
packageVersion("Signac")
packageVersion("monocle3")
```

---

## Data Requirements

### Input Files

| File | Description | Format |
|------|-------------|--------|
| `processed_rna.qs` | Pre-processed Seurat object with beta cells | QS |
| `result.csv` | Algorithme.AI disease confidence scores | CSV |
| `Donor_Summary_186.csv` | Clinical metadata (HbA1c, demographics) | CSV |
| `AUROC.csv` | Model performance metrics | CSV |
| `interpretability_metrics.csv` | Model interpretability comparison | CSV |
| `pancreas.integrated.rds` | DHT treatment dataset (optional) | RDS |
| GSE217775 files | PINK1 KO mouse data (optional) | TXT |

### Directory Structure

```
Project/
├── DATA/
│   ├── QC_supplemental/
│   │   └── interpretability_metrics.csv
│   ├── DGE/
│   │   └── wilcox/
│   │       ├── by_beta_sex/
│   │       ├── cross_cluster/
│   │       ├── beta_vs_beta/
│   │       ├── T2D_vs_beta1/
│   │       └── pseudoperterbseq/
│   ├── ORA/
│   │   ├── by_beta_sex/
│   │   │   ├── UP/
│   │   │   └── DOWN/
│   │   ├── beta_vs_beta/
│   │   ├── T2D_vs_beta1/
│   │   └── ptseqbeta1_vs_t2d/
│   ├── seurat_objects/
│   │   ├── beta_cells.qs
│   │   └── combined.qs
│   ├── KANN/
│   │   ├── qc_images/
│   │   ├── fit_all_state.pt
│   │   └── input_scaler.qs
│   ├── PINK1_discovery/
│   │   └── result_with_all_features.csv
│   ├── Dataset_GSE217775/
│   │   ├── GSE217775_RAW/
│   │   ├── dge/
│   │   └── ora/
│   ├── AUROC.csv
│   ├── result.csv
│   └── Donor_Summary_186.csv
├── images/
└── REnvironments/
    └── Analytical_Pipeline.R
```

---

## Analysis Workflow

### 1. Data Loading and Preprocessing

```r
# Load Seurat object
processed_rna <- qread("path/to/processed_rna.qs")

# Load Algorithme.AI scores
df <- read_csv("path/to/result.csv")

# Subset and merge
subset_rna <- subset(processed_rna, cells = df$barcode)
# Add metadata
for (col in colnames(df)) {
  subset_rna[[col]] <- df[[col]]
}
```

### 2. Beta Cell Subtyping

**Four beta cell subtypes identified (β1-β4):**

- **β1**: Healthy/functional beta cells (baseline)
- **β2**: Intermediate dysfunction
- **β3**: Stressed beta cells
- **β4**: Most dysfunctional

```r
# Recluster beta cells
beta_cells <- NormalizeData(beta_cells)
beta_cells <- FindVariableFeatures(beta_cells)
beta_cells <- ScaleData(beta_cells)
beta_cells <- RunPCA(beta_cells, verbose = FALSE)
beta_cells <- RunUMAP(beta_cells, dims = 1:8)
beta_cells <- FindNeighbors(beta_cells, dims = 1:8)
beta_cells <- FindClusters(beta_cells, resolution = 0.2)

# Collapse clusters into subtypes
beta_cells$collapsed_cluster <- case_when(
  beta_cells$seurat_clusters %in% c(0, 6) ~ "A",  # β1
  beta_cells$seurat_clusters %in% c(1, 5) ~ "B",  # β3
  beta_cells$seurat_clusters %in% c(2, 4) ~ "C",  # β4
  beta_cells$seurat_clusters == 3 ~ "D"           # β2
)
```

### 3. Algorithme.AI Validation

**Performance Metrics:**
- ROC curve analysis
- Confusion matrices (stratified by sex and cell type)
- Correlation with HbA1c
- AUROC comparison with baseline models (Random Forest, Gradient Boosting)

```r
# ROC analysis
library(pROC)
roc_obj <- roc(subset_rna$is_T2D, subset_rna$algorithme_confidence)
plot(roc_obj, print.auc = TRUE)

# Optimal threshold
opt_thresh <- coords(roc_obj, "best", ret = "threshold")
```

### 4. Differential Gene Expression (DGE)

**Three comparison strategies:**

#### A. Within-Status (ND or T2D)
Compare β2-4 vs β1 within healthy or diseased state

```r
# Example: β2 vs β1 in ND, stratified by sex
for (sex in c("M", "F")) {
  results <- FindMarkers(
    object = beta_cells,
    ident.1 = paste0("β2_", sex, "_ND"),
    ident.2 = paste0("β1_", sex, "_ND"),
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.137504  # ~10% change
  )
}
```

#### B. Cross-Status
Compare T2D β-clusters vs ND β1 (reference)

```r
# β2 T2D vs β1 ND
for (sex in c("M", "F")) {
  results <- FindMarkers(
    object = beta_cells,
    ident.1 = paste0("β2_T2D_", sex),
    ident.2 = paste0("β1_ND_", sex),
    test.use = "wilcox"
  )
}
```

#### C. Disease Effect (T2D vs ND)
Within each β-subtype

```r
# β2 T2D vs β2 ND
results <- FindMarkers(
  object = beta_cells,
  ident.1 = "β2_T2D",
  ident.2 = "β2_ND"
)
```

**Test Parameters:**
- Test: Wilcoxon rank-sum
- `min.pct`: 0.1 (10% minimum expression)
- `logfc.threshold`: 0.137504 (natural log, ~10% fold-change)
- Multiple testing: Bonferroni correction

### 5. Over-Representation Analysis (ORA)

```r
library(gprofiler2)

# Separate UP and DOWN genes
up_genes <- filter(deg_results, p_val_adj < 0.05, avg_log2FC > 0)$gene
down_genes <- filter(deg_results, p_val_adj < 0.05, avg_log2FC < 0)$gene

# Run GO enrichment
go_results <- gost(
  up_genes, 
  organism = "hsapiens",
  correction_method = "fdr",
  sources = "GO"
)
```

**Automated ORA Pipeline:**
- Processes all DGE results
- Splits by UP/DOWN regulation
- FDR < 0.05 threshold
- Removes redundant columns
- Saves to organized directory structure

### 6. Module Score Analysis

**Pathway gene sets:**

```r
gene_sets <- list(
  UPR = c("HSPA5", "XBP1", "ATF4", "ATF6", "DDIT3", ...),
  Mitophagy = c("PINK1", "PRKN", "BNIP3", "SQSTM1", ...),
  OXPHOS = c("NDUFA9", "COX4I1", "ATP5A1", ...),
  Stress_Death = c("BAX", "BCL2", "CASP3", "TP53", ...),
  Inflammation = c("NFKB1", "TNF", "IL1B", "CXCL10", ...)
)

# Calculate scores
for (name in names(gene_sets)) {
  beta_cells <- AddModuleScore(
    beta_cells, 
    features = list(gene_sets[[name]]),
    name = name
  )
}
```

**Correlation with disease score:**
```r
# Per-donor aggregation
df_long <- beta_cells@meta.data %>%
  group_by(Library, Sex, diabetes_status) %>%
  summarise(mean_ModuleScore = mean(Module1, na.rm = TRUE))

# Correlation test
corr_tbl <- df_long %>%
  group_by(Sex) %>%
  summarise(
    r = cor(mean_DiseaseScore, mean_ModuleScore),
    p = cor.test(mean_DiseaseScore, mean_ModuleScore)$p.value
  )
```

### 7. PINK1-Centered Analysis

#### A. Mitophagy Gene Expression

**40+ mitophagy genes analyzed:**
- Core machinery: PINK1, PRKN, SQSTM1, OPTN
- Receptors: BNIP3, BNIP3L, FUNDC1, PHB2
- Regulators: ULK1, TBK1, MFN2, FIS1
- Autophagy: MAP1LC3A/B/C, ATG5, ATG7, BECN1

```r
mitophagy_genes <- c(
  "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", ...
)

# DotPlot visualization
DotPlot(
  beta_cells,
  features = mitophagy_genes,
  group.by = "beta_sex_diab",
  split.by = "Sex"
)
```

#### B. Pseudo-Perturbation Analysis

**Strategy**: Compare PINK1-high vs PINK1-low cells within ND β1

```r
# Define PINK1 status
pink1_counts <- FetchData(beta1_nd, "PINK1", slot = "counts")[,1]
beta1_nd$PINK1_status <- "Other"
beta1_nd$PINK1_status[pink1_counts >= 2] <- "PINK1_pos"
beta1_nd$PINK1_status[pink1_counts == 0] <- "PINK1_neg"

# Differential expression
de_results <- FindMarkers(
  beta1_nd,
  ident.1 = "PINK1_pos",
  ident.2 = "PINK1_neg"
)
```

#### C. Interaction Models

**Identify genes that modify PINK1 effects on T2D:**

```r
# Logistic regression with interaction
fit_interaction <- function(gene) {
  formula <- as.formula(paste0("is_T2D ~ PINK1 * ", gene))
  model <- glm(formula, data = mito_algo, family = binomial)
  
  # Extract interaction term
  tidy(model) %>% filter(term == paste0("PINK1:", gene))
}

# Run for all mitophagy genes
results <- lapply(mitophagy_genes, fit_interaction) %>%
  bind_rows() %>%
  arrange(p.value)
```

**Visualization:**
- Volcano plots with FDR correction
- Sex-stratified analysis
- Label significant interactions

### 8. DHT Treatment Analysis

**Dataset**: Male pancreatic beta cells treated with DHT (androgen)

```r
# Load DHT dataset
dht_dataset <- readRDS("pancreas.integrated.rds")
dht_male_beta <- subset(dht_dataset, sex == "Male" & 
                        celltype %in% c("Beta INS-hi", "Beta INS-low"))

# Map to reference
anchors <- FindTransferAnchors(
  reference = beta_cells,
  query = dht_male_beta,
  normalization.method = "SCT"
)

dht_male_beta <- MapQuery(
  anchorset = anchors,
  query = dht_male_beta,
  reference = beta_cells,
  refdata = list(beta_cluster = "beta_cluster")
)
```

**Comparisons:**
- EtOH (vehicle control) vs DHT[10nM]
- Mitophagy module scores
- Pathway ratio analysis (OXPHOS/Glycolysis/TCA)

### 9. GSE217775 Analysis (PINK1 KO Mice)

**Dataset**: Kidney tissue from PINK1 knockout mice
- Wildtype vs PINK1 KO
- 4 months vs 24 months
- 2 biological replicates per group

```r
# Read raw count files
counts_matrix <- [combined from individual sample files]

# Differential expression (limma-voom)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~0 + group)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)

# Contrasts
contrast.matrix <- makeContrasts(
  KO_4M_vs_WT_4M = KO_4M - WT_4M,
  KO_24M_vs_WT_24M = KO_24M - WT_24M,
  levels = design
)
```

**Key analyses:**
- Mitochondrial gene expression heatmaps
- Pathway module scores
- Metabolic ratios (OXPHOS/TCA/Glycolysis)
- Age-dependent effects

### 10. Kolmogorov-Arnold Network (KAN) for MQI Prediction

**Mitochondrial Quality Index (MQI) Formula:**

```
MQI = Σ(good_signals) - Σ(stress_penalties)

Good signals:
  + Mitophagy activity (weight: 1.0)
  + OXPHOS capacity (weight: 1.0)
  + UPRmt (weight: 0.7)

Stress penalties:
  - TCA overactivity (weight: 0.2)
  - MT-DNA fraction (weight: 0.5)
  - OXPHOS-Glycolysis imbalance (weight: 0.5)
```

**KAN Architecture:**

```r
# Input features: ~50-140 mitochondrial regulators
regulators <- c(
  mitophagy_core, dynamics, biogenesis, 
  uprmt_proteostasis, antioxidant, signaling_tfs
)

# Model structure
KANNet <- nn_module(
  initialize = function(in_features, hidden = 48, n_bins = 10) {
    self$layer1 <- KANLayer(in_features, hidden, n_bins)
    self$act1 <- nn_relu()
    self$layer2 <- KANLayer(hidden, 1, n_bins)
  }
)

# Training
fit <- train_kan(
  X, y, 
  epochs = 300, 
  lr = 2e-3,
  hidden = 48,
  n_bins = 10
)
```

**Spline-based activation:**
- B-spline basis functions with learnable coefficients
- 10 knots per input dimension
- Captures non-linear relationships

**Model artifacts:**
- `fit_all_state.pt`: PyTorch state dict
- `input_scaler.qs`: Feature normalization (μ, σ)
- `imp_all.qs`: Feature importance scores

---

## Visualization

### 1. UMAP Projections

```r
# Disease score by sex and diabetes status
p <- FeaturePlot(
  beta_cells, 
  features = "algorithme_confidence",
  split.by = "Sex_Disease"
) + scale_color_viridis_c()
```

### 2. Sankey Diagrams

**Interactive flow diagrams showing pathway enrichment:**

```r
make_sankey_plot <- function(
  my_comparison_colors,
  my_pathway_colors,
  my_ora_dir,
  my_pathways
) {
  # [Implementation in script]
  # Creates alluvial flow from comparisons → pathways
  # Color by comparison, size by -log10(FDR)
}
```

**Features:**
- Comparison → Pathway connections
- Flow width proportional to significance
- Color-coded by comparison groups
- Legend with dot size scale

### 3. DotPlots

**Gene expression heatmaps:**

```r
# Sex-aware DotPlot
DotPlot(
  beta_cells,
  features = gene_list,
  group.by = "cluster_sex_diabetes",
  dot.scale = 6,
  scale = TRUE
) + 
  facet_grid(. ~ gene_category + Sex) +
  scale_fill_gradient2(low = "white", high = "red")
```

**Features:**
- Dot size = % cells expressing
- Color = scaled average expression
- Faceted by pathway and sex

### 4. Volcano Plots

**PINK1 interaction effects:**

```r
volcano_plot(
  results,
  title = "PINK1 × Gene Interactions",
  xlim = c(-10, 5),
  ylim = c(0, 20)
)
```

**Features:**
- X-axis: Interaction coefficient (log-odds)
- Y-axis: -log10(p-value)
- Colors: FDR < 0.05 (red), p < 0.05 (yellow), NS (grey)
- Labels: All significant genes with `ggrepel`

### 5. Heatmaps

**ComplexHeatmap visualization:**

```r
# Mitophagy genes grouped by function
Heatmap(
  expr_mat,
  col = colorRamp2(c(-1.5, 0, 1.5), 
                   c("dodgerblue4", "white", "firebrick4")),
  cluster_rows = TRUE,
  split = gene_groups,
  left_annotation = row_annot
)
```

**Features:**
- Z-scored expression
- Hierarchical clustering within groups
- Functional category annotations
- Sex-aware column grouping

### 6. Correlation Plots

**Module scores vs Disease score or HbA1c:**

```r
# Per-donor correlation
ggplot(df, aes(x = mean_HbA1c, y = mean_ModuleScore)) +
  geom_point(aes(color = Sex, shape = Disease), size = 3) +
  geom_smooth(method = "lm", aes(color = Sex)) +
  stat_cor(aes(color = Sex))
```

---

## Key Outputs

### Figure Files

| Figure | Description | Format |
|--------|-------------|--------|
| UMAP projections | Disease scores by sex/diabetes | PNG/PDF |
| Sankey diagrams | Pathway enrichment flows | PNG/PDF |
| Volcano plots | PINK1 interactions | PNG/PDF |
| Heatmaps | Gene expression patterns | PNG/PDF |
| Correlation plots | Module scores vs clinical | PNG/PDF |
| Violin plots | Distribution comparisons | PNG/PDF |
| DotPlots | Multi-gene expression | PNG/PDF |

### Data Tables

| File | Description |
|------|-------------|
| `*_results.csv` | DGE results per comparison |
| `gost_*_UP.csv` | GO enrichment (upregulated) |
| `gost_*_DOWN.csv` | GO enrichment (downregulated) |
| `AUROC.csv` | Model performance metrics |
| `beta_cells_metadata.csv` | Annotated cell metadata |

### Model Artifacts

| File | Description |
|------|-------------|
| `fit_all_state.pt` | KAN model weights (PyTorch) |
| `input_scaler.qs` | Feature normalization parameters |
| `imp_all.qs` | Feature importance scores |
| `beta_cells.qs` | Annotated Seurat object |
| `combined.qs` | Merged beta cells + DHT dataset |

---

## Color Schemes

### Beta Cell Subtypes
```r
beta_cols <- c(
  "β1" = "#B22222",  # firebrick
  "β2" = "#DA70D6",  # orchid
  "β3" = "#FF8C00",  # darkorange
  "β4" = "#1E90FF"   # dodgerblue
)
```

### Sex
```r
sex_cols <- c(
  "M" = "#1976D2",  # blue
  "F" = "#C2185B"   # pink
)
```

### Disease Status
```r
disease_cols <- c(
  "ND"  = "grey65",
  "T2D" = "firebrick"
)
```

### Treatment
```r
treatment_cols <- c(
  "Untreated" = "grey65",
  "DHT[10nM]" = "firebrick4",
  "EtOH"      = "dodgerblue4"
)
```

---

## Reproducibility

### Set Seed

```r
set.seed(1234)  # Used throughout analysis

# For torch models
torch::torch_manual_seed(42)
```

### Session Info

```r
sessionInfo()
# R version 4.2.1 (2019-12-12)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64
```

### Package Versions

```r
packageVersion("Seurat")        # 4.3.0.1
packageVersion("Signac")        # 1.11.0
packageVersion("clusterProfiler") # [version]
packageVersion("ComplexHeatmap") # [version]
```

## 🧠 Algorithme.ai Models

Algorithme.ai introduces a new paradigm in machine learning — **lookalike-based reasoning**.  
The following Python modules demonstrate two model classes built around this concept.

---

### 🐍 Part I — `tulane.py` → `SnakeClassifier`
![SnakeClassifier](Algorithm.Ai/Python%20Models/Part1/SnakeClassifier.jpeg)

The **SnakeClassifier** model establishes the foundation of the lookalike rationale.  
It learns by identifying patterns and extrapolating confidence from similar data points.

```python
from tulane import (
    make_population,
    extrapolate,
    my_function,
    datapoint_to_csv,
    my_function_audit
)
from random import choice

# Generate a list of dictionaries from a CSV file
population = make_population("train.csv")

# Extrapolate confidence percentages from 'train.csv' → 'result.csv'
extrapolate("train.csv")

# Select a random element from the population
item = choice(population)

# Print the CSV of lookalikes for the chosen datapoint
print(datapoint_to_csv(item))

# Retrieve elements in favor and against the datapoint
favorable, against = my_function_audit(item)

# Compute the confidence percentage for the specific datapoint
confidence = my_function(item)
print("Confidence:", confidence)
```

---

### 🦢 Part II — `tulane.py` → `BlackSwanClassifier`
![BlackSwanClassifier](Algorithm.Ai/Python%20Models/Part2/BlackSwanClassifier.jpeg)
The **BlackSwanClassifier** extends the SnakeClassifier by incorporating **top feature selection** and **feature importance weighting**, enhancing predictive precision.

```python
from tulane import (
    make_population,
    extrapolate,
    my_function,
    datapoint_to_csv,
    my_function_audit,
    extrapolate_with_top_features,
    extrapolate_with_all_features
)
from random import choice

# Generate a list of dictionaries from a CSV file
population = make_population("train.csv")

# Standard extrapolation
extrapolate("train.csv")

# Select a random element from the population
item = choice(population)

# Print the CSV of lookalikes for the chosen datapoint
print(datapoint_to_csv(item))

# Retrieve elements in favor and against the datapoint
favorable, against = my_function_audit(item)

# Compute the confidence percentage for the specific datapoint
confidence = my_function(item)
print("Confidence:", confidence)

# Perform extrapolation using top-ranked features
extrapolate_with_top_features("train.csv")

# Perform extrapolation using all available features
extrapolate_with_all_features("train.csv")
```

---

### 🧩 Summary

| Model | Description | Key Functions |
|--------|--------------|----------------|
| **SnakeClassifier** | Base model introducing lookalike-based reasoning | `make_population`, `extrapolate`, `my_function` |
| **BlackSwanClassifier** | Enhanced model with top feature selection and broader extrapolation logic | `extrapolate_with_top_features`, `extrapolate_with_all_features` |

---

## Troubleshooting

### Common Issues

#### 1. Memory Errors
```r
# Increase future global size
options(future.globals.maxSize = 8000 * 1024^2)
```

#### 2. Missing Genes
```r
# Check gene presence before module scores
present_genes <- intersect(gene_list, rownames(seurat_obj))
if (length(present_genes) == 0) {
  warning("No genes found in dataset")
}
```

#### 3. Python Configuration
```r
# Configure reticulate
library(reticulate)
use_condaenv("r-reticulate")
py_config()
```

#### 4. UTF-8 Encoding
```r
# For β symbols in Windows
Sys.setlocale("LC_ALL", "English")
```

---

## Citation

If you use this pipeline in your research, please cite:

```
[Paper citation to be added upon publication]

Qadir, F. et al. (2025). Sex-based analysis of mitochondrial function in 
pancreatic beta cells reveals PINK1-dependent mechanisms in type 2 diabetes.
[Journal Name]. [DOI]
```

---

## Contact

**Fahd Qadir (Dragonmasterx87)**  
FMJ Lab  
Tulane University School of Medicine  
Email: [Add email address]

---

## Acknowledgments

- **Algorithme.AI team** for disease scoring algorithm development
- **Human Pancreas Analysis Program (HPAP)** for donor samples and data access
- **GEO contributors** for GSE217775 PINK1 KO dataset
- **Barko et al.** for DHT treatment dataset
- **Tulane HPC** for computational resources

---

## License

[Specify license - e.g., MIT, GPL-3, Apache 2.0]

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2022-11-16 | Initial pipeline development |
| 1.1 | [Date] | Added KAN implementation |
| 1.2 | [Date] | Integrated GSE217775 analysis |
| 1.3 | [Date] | Added DHT treatment mapping |

---

**Last Updated**: October 13, 2025  
**Pipeline Version**: 1.3  
**Status**: Production
