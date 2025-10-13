# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

# LOAD LIBRARIES ####
# Restart Rstudio or R
# make sure a Rtools version is installed that supports your current versin of R
install.packages('ggplot2')
install.packages('cowplot')
install.packages('Matrix')
install.packages('ggridges')
install.packages('ggrepel')
install.packages('dplyr')
#install.packages('Seurat')
install.packages('plotly')
install.packages('clustree')
install.packages('patchwork')
install.packages('future')
install.packages("devtools")
install.packages("rlang")
install.packages("pROC")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'), force = TRUE)

BiocManager::install(c("EnhancedVolcano", "DoubletFinder", "glmGamPoi",
                       "GOSemSim", "org.Hs.eg.db", "AnnotationHub",
                       "GenomeInfoDb", "MeSHDbi", "clusterProfiler",
                       "dittoSeq", "escape", "ComplexHeatmap", "DropletUtils", 
                       "Nebulosa", "hdf5r", "scDblFinder", "JASPAR2020",
                       "TFBSTools", "motifmatchr", "GreenleafLab/chromVAR",
                       "EnrichmentBrowser"),
                     dependencies = T, force = TRUE)

BiocManager::install("clusterProfiler")

# Pando compatible
remotes::install_version(package = 'Seurat', version = package_version('4.3.0.1'))
devtools::install_version("SeuratObject", version = "4.1.4", repos = "http://cran.us.r-project.org")
devtools::install_version("Signac", version = "1.11.0", repos = "http://cran.us.r-project.org")

BiocManager::install("scDblFinder", dependencies = T, force = TRUE)
BiocManager::install("HDO.db")

# install Seurat from Github (automatically updates sctransform)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
install.packages("devtools")
devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("gaospecial/ggVennDiagram")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("harmony")
install.packages('SoupX')
install.packages('tidyverse')
install.packages("viridis")
install.packages("circlize")
install.packages("scCustomize")
install.packages("archive")
install.packages("R.utils")
install.packages("qs")
install.packages("ggseqlogo", force = TRUE)
install.packages('fastmap', force = TRUE)
devtools::install_github('quadbiolab/Pando')
usethis::create_github_token()

#devtools::install_github("stuart-lab/signac", ref = "develop", force = TRUE)
BiocManager::install(c("GenomicRanges", force = TRUE))
BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')

# Run the following code once you have Seurat installed
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(leiden)
  library(stringr)
  library(hdf5r)
  library(SoupX)
  library(Rcpp)
  library(cowplot)
  library(Matrix)
  library(ggridges)
  library(dplyr)
  library(tidyverse)
  library(data.table)
  library(reticulate)
  library(Seurat)
  library(monocle3)
  library(harmony)
  library(Signac)
  library(scDblFinder)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(plotly)
  library(clustree)
  library(patchwork)
  library(future)
  library(DoubletFinder)
  library(EnhancedVolcano)
  library(glmGamPoi)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(AnnotationHub)
  library(MeSHDbi)
  library(clusterProfiler)
  library(DOSE)
  library(dittoSeq)
  library(escape)
  library(EnrichmentBrowser)
  library(viridisLite)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  #library(scCustomize)
  library(Nebulosa)
  library(DropletUtils)
  library(ggvenn)
  library(ggVennDiagram)
  library(devtools)
  library(R.utils)
  library(qs)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(motifmatchr)
  library(chromVAR)
  #library(SeuratWrappers)
  library(cicero)
  library(BiocParallel)
  #library(ggseqlogo)
  library(DESeq2)
  library(Pando)
  # Set global environment parameter par-proc
  # options(future.globals.maxSize = 8000 * 1024^2)
  library(ragg)
  library(gprofiler2)
  set.seed(1234)
})
set.seed(1234)

# Python env
if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/mqadir/AppData/Local/r-miniconda/envs/r-reticulate",Sys.getenv()["PATH"],sep=";"))
py_config()

# WD
setwd(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\WD)")
(WD <- getwd())


# Check package versions
packageVersion("clusterProfiler")
packageVersion("dittoSeq")
packageVersion("escape")
packageVersion("seurat")
packageVersion("signac")
packageVersion("EnrichmentBrowser")
packageVersion("monocle3")
packageVersion("cicero")
packageVersion("scDblFinder")

#Fig1 Supplemental
# --- Packages ---
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)   # for arranging plots in one row
})

# --- Input / Output paths ---
infile <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/QC_supplemental/interpretability_metrics.csv"
outfile <- sub("\\.csv$", "_5panel.pdf", infile)

# --- Load (robust to comma- or tab-delimited) ---
df <- tryCatch(read_csv(infile, show_col_types = FALSE),
               error = function(e) readr::read_tsv(infile, show_col_types = FALSE))

# Order models as shown in your table
model_levels <- c("Snake", "BlackSwan", "RandomForest", "GradientBoosting")
df <- df %>% mutate(model = factor(model, levels = model_levels))

# --- Helper: make a clean bar plot for a metric ---
make_bar <- function(data, col, title, zero_line = FALSE) {
  p <- ggplot(data, aes(x = model, y = .data[[col]])) +
    geom_col(width = 0.7) +
    geom_text(aes(label = ifelse(abs(.data[[col]]) < 10,
                                 sprintf("%.3f", .data[[col]]),
                                 sprintf("%.1f", .data[[col]]))),
              vjust = ifelse(data[[col]] >= 0, -0.35, 1.2), size = 3) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  if (zero_line) p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p
}

# --- Build five plots ---
p1 <- make_bar(df, "sparsity_frac", "Sparsity Fraction")
p2 <- make_bar(df, "stability_kendall", "Stability (Kendall \u03C4)")
p3 <- make_bar(df, "monotonic_consistency", "Monotonic Consistency")
p4 <- make_bar(df, "rule_complexity", "Rule Complexity")
p5 <- make_bar(df, "surrogate_R2", "Surrogate R\u00B2", zero_line = TRUE)

# --- Arrange in one row and save ---
panel <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 5)
ggsave(outfile, panel, width = 20, height = 4, units = "in")

# Also print to the viewer if running interactively
print(panel)


# Load dataset
processed_rna <- qread(r"(F:\2.SexbasedStudyCurrent\QS files\processed_rna.qs)")

# Load required libraries
library(Seurat)
library(dplyr)
library(readr)

# Step 1: Read the CSV file (adjust the path to your actual file location)
df <- read_csv(R"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\result.csv)")  # Update this path

# Step 2: Extract the barcodes from the CSV
barcodes_to_keep <- df$barcode

# Step 3: Subset the Seurat object using those barcodes
subset_rna <- subset(processed_rna, cells = barcodes_to_keep)

# Optional: Check how many cells were retained
length(Cells(subset_rna))

# Step 1: Set barcodes as rownames (if not already)
df <- df %>% column_to_rownames("barcode")

# Step 2: Reorder rows to match Seurat cell names
df <- df[Cells(subset_rna), ]

# Step 3: Make sure it's a proper data.frame (not tibble)
df <- as.data.frame(df)

# Step 4: Add metadata properly
for (col in colnames(df)) {
  subset_rna[[col]] <- df[[col]]
}

FeaturePlot(subset_rna, features = "algorithme_confidence", reduction = "umap") +
  scale_color_viridis_c(option = "A") +
  ggtitle("UMAP: Algorithmme.ai Disease Score")

VlnPlot(subset_rna, features = "algorithme_confidence", pt.size = 0, group.by = "celltype_sex", split.by = "diabetes_status") +
  ggtitle("Disease Score by Cell Type")

subset_rna$score_group <- case_when(
  subset_rna$algorithme_confidence >= 0.75 ~ "High",
  subset_rna$algorithme_confidence <= 0.25 ~ "Low",
  TRUE ~ "Mid"
)
table(subset_rna$score_group)

processed_rna <- NULL

# Step 1: Make sure metadata is clean
table(subset_rna$Sex)
table(subset_rna$diabetes_status)

# Step 2: Subset each group
male_nd     <- subset(subset_rna, subset = Sex == "M" & diabetes_status == "ND")
male_t2d    <- subset(subset_rna, subset = Sex == "M" & diabetes_status == "T2D")
female_nd   <- subset(subset_rna, subset = Sex == "F" & diabetes_status == "ND")
female_t2d  <- subset(subset_rna, subset = Sex == "F" & diabetes_status == "T2D")

# Step 3: Generate FeaturePlots for each
p1 <- FeaturePlot(male_nd, features = "algorithme_confidence", reduction = "umap") + ggtitle("Male - ND")
p2 <- FeaturePlot(male_t2d, features = "algorithme_confidence", reduction = "umap") + ggtitle("Male - T2D")
p3 <- FeaturePlot(female_nd, features = "algorithme_confidence", reduction = "umap") + ggtitle("Female - ND")
p4 <- FeaturePlot(female_t2d, features = "algorithme_confidence", reduction = "umap") + ggtitle("Female - T2D")

# Step 4: Combine and display
(p1 | p2) / (p3 | p4)

# Combine plots
combined_plot <- (p1 | p2) / (p3 | p4)

# Save as high-resolution PNG (or TIFF if needed)
agg_png("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/images/umap_disease_sex.png", width = 3000, height = 3000, res = 300)
print(combined_plot)
dev.off()

library(ggplot2)

# Violin + boxplot overlay
ggplot(subset_rna@meta.data, aes(x = factor(is_T2D), y = algorithme_confidence)) +
  geom_violin(fill = "lightgray") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  labs(x = "Actual T2D Status", y = "Algorithmme.ai Disease Score") +
  theme_minimal()

library(pROC)

# Create ROC object
roc_obj <- roc(subset_rna$is_T2D, subset_rna$algorithme_confidence)

# Plot
plot(roc_obj, print.auc = TRUE, col = "#2C3E50")

# Optimal threshold
opt_thresh <- coords(roc_obj, "best", ret = "threshold")
# Extract numeric value from the threshold result
thresh_value <- as.numeric(opt_thresh)
opt_thresh

# Use the extracted numeric threshold for comparison
subset_rna$predicted_T2D <- ifelse(subset_rna$algorithme_confidence >= thresh_value, 1, 0)

# Step 2: Confusion matrix
table(Predicted = subset_rna$predicted_T2D, Actual = subset_rna$is_T2D)

# Step 3: Accuracy
mean(subset_rna$predicted_T2D == subset_rna$is_T2D)

library(ggplot2)
library(dplyr)

# Step 1: Create group-wise confusion matrix (Predicted vs Actual)
conf_df <- subset_rna@meta.data %>%
  mutate(
    Predicted = factor(predicted_T2D, levels = c(0, 1), labels = c("ND", "T2D")),
    Actual = factor(is_T2D, levels = c(0, 1), labels = c("ND", "T2D"))
  ) %>%
  group_by(Sex, celltype_qadir, Actual, Predicted) %>%
  summarise(Count = n(), .groups = "drop")

# Step 2: Calculate percentages *within each actual class per group*
conf_df <- conf_df %>%
  group_by(Sex, celltype_qadir, Actual) %>%
  mutate(Percent = round(Count / sum(Count) * 100, 1)) %>%
  ungroup()

# Step 3: Plot faceted confusion matrix
ggplot(conf_df, aes(x = Actual, y = Predicted, fill = Percent)) +
  geom_tile(color = "black") +
  geom_text(aes(label = paste0(Percent, "%")), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#2C3E50") +
  labs(
    title = "Stratified Confusion Matrix by Sex and Cell Type",
    x = "Actual T2D",
    y = "Predicted T2D",
    fill = "Percent"
  ) +
  facet_grid(Sex ~ celltype_qadir, switch = "both") +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.position = "right"
  )

library(dplyr)
library(ggplot2)
library(ggpubr)

# Step 1: Calculate actual vs predicted % T2D per celltype and sex
agg_df <- subset_rna@meta.data %>%
  group_by(Sex, celltype_qadir) %>%
  summarise(
    actual_pct_T2D = mean(is_T2D == 1) * 100,
    predicted_pct_T2D = mean(predicted_T2D == 1) * 100,
    .groups = "drop"
  )

# Step 2: Correlation plot
# Define your color palette
celltype_colors <- c(
  "beta" = "dodgerblue3",
  "beta+alpha" = "turquoise2",
  "alpha" = "lightseagreen",
  "cycling_endo" = "darkseagreen2",
  "epsilon" = "khaki2",
  "gamma" = "springgreen4",
  "delta" = "chartreuse3",
  "beta+delta" = "burlywood3",
  "ductal" = "darkorange2",
  "acinar" = "salmon3",
  "activated_stellate" = "orange",
  "quiescent_stellate" = "salmon",
  "endothelial" = "red",
  "macrophages" = "magenta3",
  "lymphocyte" = "orchid1",
  "mast" = "red4",
  "schwann" = "grey30"
)

# Libraries
library(ggplot2)
library(ggpubr)
library(dplyr)

# Step 1: Compute slope and intercept per sex
coefs <- agg_df %>%
  group_by(Sex) %>%
  summarise(
    intercept = coef(lm(predicted_pct_T2D ~ actual_pct_T2D))[1],
    slope = coef(lm(predicted_pct_T2D ~ actual_pct_T2D))[2],
    .groups = "drop"
  )

# Step 2: Clip regression lines to visible region: x = [0, 100]
segments <- coefs %>%
  rowwise() %>%
  mutate(
    x_start = max(0, (0 - intercept) / slope),
    x_end = min(100, (100 - intercept) / slope),
    y_start = slope * x_start + intercept,
    y_end = slope * x_end + intercept
  )

# Step 3: Plot
ggplot(agg_df, aes(x = actual_pct_T2D, y = predicted_pct_T2D, color = celltype_qadir)) +
  geom_point(size = 3) +
  geom_text(aes(label = celltype_qadir), hjust = 0, vjust = 1.5, size = 5) +
  geom_segment(
    data = segments,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    color = "black",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  facet_wrap(~Sex) +
  scale_color_manual(values = celltype_colors) +
  stat_cor(
    method = "pearson",
    aes(group = Sex),
    color = "black",
    label.x = 5, label.y = 95,
    size = 5,
    label.sep = "\n"
  ) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Actual % T2D Cells (per Cell Type & Sex)",
    y = "Predicted % T2D Cells (Algorithm)",
    color = "Cell Type",
    title = "Correlation of Predicted vs Actual T2D Assignment"
  ) +
  theme_minimal(base_size = 12)

# AUROC
library(tidyverse)
library(ggpubr)

# Load and reshape
auroc <- read_csv("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/AUROC.csv")

df_long <- auroc %>%
  pivot_longer(cols = c("Algorithm", "Random Forest", "Gradient Boosting"),
               names_to = "Model", values_to = "AUROC") %>%
  mutate(Model = factor(Model, levels = c("Algorithm", "Gradient Boosting", "Random Forest")))

# Get best model per group
winner_df <- df_long %>%
  group_by(Sex, ancestry, celltype_qadir) %>%
  slice_max(order_by = AUROC, n = 1, with_ties = FALSE) %>%
  ungroup()

# Compute summary stats
summary_df <- winner_df %>%
  group_by(Sex, Model) %>%
  summarise(mean_auroc = mean(AUROC),
            sem = sd(AUROC) / sqrt(n()), .groups = "drop")

# Compute p-values per Sex
comparisons <- list(
  c("Algorithm", "Gradient Boosting"),
  c("Algorithm", "Random Forest"),
  c("Gradient Boosting", "Random Forest")
)

p_df <- compare_means(AUROC ~ Model, group.by = "Sex", data = winner_df,
                      method = "t.test", comparisons = comparisons)

# Add custom y-positions per comparison
p_df <- p_df %>%
  mutate(y.position = case_when(
    group1 == "Algorithm" & group2 == "Gradient Boosting" ~ 0.87,
    group1 == "Algorithm" & group2 == "Random Forest" ~ 0.91,
    TRUE ~ 0.95
  ))

# Plot
model_colors <- c("Algorithm" = "purple", "Gradient Boosting" = "gold", "Random Forest" = "green3")

ggplot() +
  geom_col(data = summary_df,
           aes(x = Model, y = mean_auroc, fill = Model),
           width = 0.6, alpha = 0.85) +
  geom_errorbar(data = summary_df,
                aes(x = Model, ymin = mean_auroc - sem, ymax = mean_auroc + sem),
                width = 0.2, linewidth = 0.8) +
  geom_jitter(data = winner_df,
              aes(x = Model, y = AUROC),
              color = "black", size = 1.5, width = 0.15) +
  geom_text(data = summary_df,
            aes(x = Model, y = mean_auroc + sem + 0.01, 
                label = sprintf("%.3f", mean_auroc)),
            size = 4.5, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = p_df,
            aes(x = group1, label = p.signif, y = y.position),
            vjust = -0.5, size = 5, inherit.aes = FALSE) +
  scale_fill_manual(values = model_colors) +
  facet_wrap(~Sex) +
  labs(
    title = "Mean AUROC ± SEM per Model by Sex (Only for Winning Subgroups)",
    subtitle = "p-values from pairwise t-tests; each dot is a winning group",
    y = "Mean AUROC", x = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

library(tidyverse)

# Load data
auroc <- read_csv("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/AUROC.csv")

# Step 1: Reshape
df_long <- auroc %>%
  pivot_longer(cols = c("Algorithm", "Random Forest", "Gradient Boosting"),
               names_to = "Model", values_to = "AUROC") %>%
  mutate(Model = factor(Model, levels = c("Algorithm", "Gradient Boosting", "Random Forest")))

# Step 2: Get winning model per subgroup
winner_df <- df_long %>%
  group_by(Sex, ancestry, celltype_qadir) %>%
  slice_max(order_by = AUROC, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 3: Compute percent of test data within each sex
winner_df <- winner_df %>%
  group_by(Sex) %>%
  mutate(scaled_share = `Share of test data` / sum(`Share of test data`) * 100) %>%
  ungroup()

# Step 4: Summary stats
summary_df <- winner_df %>%
  group_by(Sex, Model) %>%
  summarise(
    Percent = sum(scaled_share),
    SEM = sd(scaled_share) / sqrt(n()),
    .groups = "drop"
  )

# Step 5: Plot
model_colors <- c("Algorithm" = "purple", "Gradient Boosting" = "gold", "Random Forest" = "green3")

ggplot(summary_df, aes(x = Model, y = Percent, fill = Model)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = Percent - SEM, ymax = Percent + SEM),
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            vjust = -0.8, size = 5) +
  scale_fill_manual(values = model_colors) +
  facet_wrap(~Sex) +
  labs(
    title = "Weighted-Winner Summary by Sex",
    subtitle = "Each bar reflects only groups where the model performed best",
    y = "Weighted Share of Test Set (%)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Lollip
# Load packages
library(tidyverse)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)

# Load CSV from model input
df <- read_csv("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/result.csv")

# Define metadata columns to remove
meta_cols <- c(
  "barcode", "bmi", "age", "diabetes_status", "Chemistry", "ancestry", 
  "is_T2D", "Sex", "celltype_qadir",
  "algorithme_confidence", "random_forest_confidence", "gradient_boosting_confidence"
)

# Subset to feature columns only
feature_df <- df %>% select(-all_of(meta_cols))

# Get top 20 variable genes by variance
top_genes <- feature_df %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "variance") %>%
  arrange(desc(variance)) %>%
  slice_head(n = 10) %>%
  pull(gene)

# Load Seurat object
processed_rna <- qread("F:/2.SexbasedStudyCurrent/QS files/processed_rna.qs")

# Make sure relevant metadata exists
Idents(processed_rna) <- "celltype_qadir"

# Split by Sex and facet by diabetes_status using patchwork
# Females
DotPlot(
  subset(processed_rna, subset = Sex == "F"), 
  features = top_genes, 
  group.by = "celltype_qadir", 
  split.by = "diabetes_status", 
  cols = c("grey90", "darkorchid4"),
  dot.scale = 3  # 👈 reduce dot size (default is 6)
) +
  labs(title = "Females") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Males
DotPlot(
  subset(processed_rna, subset = Sex == "M"), 
  features = top_genes, 
  group.by = "celltype_qadir", 
  split.by = "diabetes_status", 
  cols = c("grey90", "steelblue4"),
  dot.scale = 3  # 👈 same here
) +
  labs(title = "Males") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

# PLOTTING corelation
# Columns: Donor_ID, HbA1c, Disease_Score, Sex
# Load HbA1c from CSV — extract correct donor ID format
clinical_df <- read_csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\Donor_Summary_186.csv)") %>%
  mutate(Library = donor_ID) %>%
  select(Library, hba1c)

# Step 2: Add manually missing HbA1c values (with full donor names)
hba1c_missing <- tribble(
  ~Library,        ~hba1c,
  "HP2022801",     5.5,
  "HP2024001",     5.4,
  "SAMN15877725",  5.6,
  "HP2024001",     5.4,
  "HP2031401",     5.4,
  "HP2105501",     5.6,
  "HP2106201",     5.3,
  "HP2107001",     5.1,
  "HP2107901",     5.2,
  "HP2108601",     5.1,
  "HP2108901",     5.9,
  "HP2110001",     5.5,
  "HP2123201",     5.3,
  "HP2132801",     5.5,
  "HP2202101",     5.5
)

# Combine all known HbA1c values
hba1c_full <- bind_rows(clinical_df, hba1c_missing)

# Get all donor scores (no filtering, all Library IDs intact)
donor_scores <- subset_rna@meta.data %>%
  group_by(Library, Sex, diabetes_status) %>%
  summarise(
    Disease_Score = mean(algorithme_confidence, na.rm = TRUE),
    .groups = "drop"
  )

# Merge based on full Library ID
merged_df <- donor_scores %>%
  left_join(hba1c_full, by = "Library") %>%
  filter(!is.na(hba1c))  # Keep rows with valid HbA1c

# Plot
ggplot(merged_df, aes(x = hba1c, y = Disease_Score)) +
  geom_point(aes(color = diabetes_status), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  stat_cor(method = "pearson", size = 5, na.rm = TRUE, parse = FALSE) +
  scale_color_manual(values = c("ND" = "black", "T2D" = "firebrick")) +
  labs(
    title = "Panel G: Disease Score Correlates with HbA1c (All Donors)",
    x = "HbA1c (%)",
    y = "Mean Disease Score (Per Donor)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# ALL
donor_scores_all <- subset_rna@meta.data %>%
  group_by(Library) %>%
  summarise(
    Disease_Score = mean(algorithme_confidence, na.rm = TRUE),
    Sex = dplyr::first(Sex),
    diabetes_status = dplyr::first(diabetes_status),
    .groups = "drop"
  )
# Auto-map using last 3 digits
donor_scores_all <- donor_scores_all %>%
  mutate(donor_ID = paste0("HPAP-", substr(Library, 6, 8)))
merged_all <- donor_scores_all %>%
  left_join(clinical_df, by = "donor_ID")
merged_all_hba1c <- merged_all %>%
  filter(!is.na(hba1c))
ggplot(merged_all_hba1c, aes(x = hba1c, y = Disease_Score)) +
  geom_point(aes(color = diabetes_status), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  stat_cor(method = "pearson", size = 5, parse = FALSE, na.rm = TRUE) +
  scale_color_manual(values = c("ND" = "grey60", "T2D" = "firebrick")) +
  labs(
    title = "Disease Score Correlates with HbA1c (All Donors)",
    x = "HbA1c (%)",
    y = "Mean Disease Score (Per Donor)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# beta cells
processed_rna <- NULL
beta_cells <- subset(subset_rna, subset = celltype_qadir == "beta")
summary(beta_cells$algorithme_confidence)
hist(beta_cells$algorithme_confidence, breaks = 50,
     col = "dodgerblue", main = "β-cell Disease Score Distribution", xlab = "Disease Score")
beta_cells$score_group <- cut(beta_cells$algorithme_confidence,
                              breaks = quantile(beta_cells$algorithme_confidence, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE),
                              labels = c("Low", "Mid", "High"),
                              include.lowest = TRUE)
DimPlot(beta_cells, group.by = "score_group", pt.size = 0.2) +
  scale_color_manual(values = c("Low" = "grey80", "Mid" = "orange", "High" = "firebrick")) +
  labs(title = "β-cell Stratification by Disease Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#Recluster
beta_cells <- NormalizeData(beta_cells)
beta_cells <- FindVariableFeatures(beta_cells)
beta_cells <- ScaleData(beta_cells)
beta_cells <- RunPCA(beta_cells, verbose = FALSE)
beta_cells <- RunUMAP(beta_cells, dims = 1:8)
beta_cells <- FindNeighbors(beta_cells, dims = 1:8)
beta_cells <- FindClusters(beta_cells, resolution = 0.2)  # try 0.4–0.8

#Plotting
DimPlot(beta_cells, group.by = "seurat_clusters", label = TRUE, pt.size = 0.2)
FeaturePlot(beta_cells, features = "algorithme_confidence", 
            cols = c("grey90", "darkred")) +
  labs(title = "Disease Score in β-cell Clusters")
FeaturePlot(beta_cells, features = "algorithme_confidence", 
            cols = c("grey90", "darkred")) +
  labs(title = "Disease Score in β-cell Clusters")
DimPlot(beta_cells, group.by = "Sex", pt.size = 0.2) +
  scale_color_manual(values = c("M" = "deepskyblue3", "F" = "orchid")) +
  labs(title = "β-cell UMAP Colored by Sex") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
DimPlot(beta_cells, group.by = "diabetes_status", pt.size = 0.2) +
  scale_color_manual(values = c("ND" = "grey50", "T2D" = "firebrick")) +
  labs(title = "β-cells Colored by Diabetes Status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
DotPlot(beta_cells, features = c("PDX1", "MAFA", "INS"), group.by = "seurat_clusters") +
  scale_color_gradient(low = "grey90", high = "darkorchid4") +
  labs(title = "DotPlot of Key β-cell Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))


beta_cells$collapsed_cluster <- case_when(
  beta_cells$seurat_clusters %in% c(0, 6) ~ "A",
  beta_cells$seurat_clusters %in% c(1, 5) ~ "B",
  beta_cells$seurat_clusters %in% c(2, 4) ~ "C",
  beta_cells$seurat_clusters == 3 ~ "D",
  TRUE ~ "Other"
)
beta_cells$collapsed_cluster <- factor(beta_cells$collapsed_cluster, levels = c("A", "B", "C", "D"))

library(ggplot2)

# Extract UMAP coordinates + cluster group
umap_df <- as.data.frame(Embeddings(beta_cells, "umap"))
umap_df$cluster <- beta_cells$collapsed_cluster

# Plot
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.01, alpha = 0.8) +  # EPS-safe, clean vector points
  scale_color_manual(values = c("A" = "firebrick", "B" = "darkorange", 
                                "C" = "dodgerblue", "D" = "orchid")) +
  labs(title = "β-cell States (Collapsed Clusters)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # white plot panel
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.8),  # black box
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/images/Figure2_Collapsed_BetaClusters.png",
       width = 6, height = 5, dpi = 600, bg = "white")


# Step 1: Subset beta cells by sex & diabetes status
male_nd     <- subset(beta_cells, subset = Sex == "Male"   & diabetes_status == "ND")
male_t2d    <- subset(beta_cells, subset = Sex == "Male"   & diabetes_status == "T2D")
female_nd   <- subset(beta_cells, subset = Sex == "Female" & diabetes_status == "ND")
female_t2d  <- subset(beta_cells, subset = Sex == "Female" & diabetes_status == "T2D")

# Step 2: FeaturePlots for disease score
p1 <- FeaturePlot(male_nd, features = "algorithme_confidence", reduction = "umap") +
  ggtitle("Male - ND") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

p2 <- FeaturePlot(male_t2d, features = "algorithme_confidence", reduction = "umap") +
  ggtitle("Male - T2D") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

p3 <- FeaturePlot(female_nd, features = "algorithme_confidence", reduction = "umap") +
  ggtitle("Female - ND") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

p4 <- FeaturePlot(female_t2d, features = "algorithme_confidence", reduction = "umap") +
  ggtitle("Female - T2D") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

# Step 3: Combine with patchwork
library(patchwork)
combined_plot <- (p1 | p2) / (p3 | p4)

# Step 4: Save high-res image
agg_png("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/images/umap_disease_sex.png",
        width = 3000, height = 3000, res = 300)
print(combined_plot)
dev.off()

library(ggplot2)
library(dplyr)

# Label clusters β1–β4
beta_cells$beta_cluster <- recode(beta_cells$collapsed_cluster,
                                  "A" = "β1", "D" = "β2", "B" = "β3", "C" = "β4")
beta_cells$beta_cluster <- factor(beta_cells$beta_cluster, levels = c("β1", "β2", "β3", "β4"))

# Combine sex and diabetes into group
beta_cells$group <- paste0(beta_cells$Sex, "_", beta_cells$diabetes_status)

# Build plotting dataframe
score_df <- beta_cells@meta.data %>%
  select(beta_cluster, group, algorithme_confidence)

# Plot: Split by β cluster
ggplot(score_df, aes(x = 1, y = algorithme_confidence, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, width = 0.9, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", alpha = 0.9) +
  facet_grid(group ~ beta_cluster, switch = "y") +  # 4 rows (F_ND, F_T2D, M_ND, M_T2D)
  scale_fill_manual(values = c(
    "M_ND"   = "#5BAFBF",
    "M_T2D"  = "#228B22",
    "F_ND"   = "#E69F00",
    "F_T2D"  = "#B22222"
  )) +
  labs(
    x = "",
    y = "Disease Score",
    title = "Panel B: β-cell Disease Score by Cluster, Sex, and T2D Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

#
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(tibble)

# Step 1: Aggregate beta cell disease scores by donor and cluster
donor_scores_beta <- beta_cells@meta.data %>%
  group_by(Library, beta_cluster) %>%
  summarise(
    Disease_Score = mean(algorithme_confidence, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Load clinical HbA1c data using exact Seurat donor names (no gsub truncation!)
clinical_df <- read_csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\Donor_Summary_186.csv)") %>%
  mutate(Library = donor_ID) %>%
  select(Library, hba1c)

# Step 3: Add missing hba1c values (with full names)
hba1c_missing <- tribble(
  ~Library,        ~hba1c,
  "HP2022801",     5.5,
  "HP2024001",     5.4,
  "SAMN15877725",  5.6,
  "HP2031401",     5.4,
  "HP2105501",     5.6,
  "HP2106201",     5.3,
  "HP2107001",     5.1,
  "HP2107901",     5.2,
  "HP2108601",     5.1,
  "HP2108901",     5.9,
  "HP2110001",     5.5,
  "HP2123201",     5.3,
  "HP2132801",     5.5,
  "HP2202101",     5.5,
  "HP2121601",     5.8
)

# Step 4: Combine full hba1c data
hba1c_full <- bind_rows(clinical_df, hba1c_missing)

# Step 5: Merge disease scores with hba1c
merged_df <- left_join(donor_scores_beta, hba1c_full, by = "Library") %>%
  filter(!is.na(hba1c))  # remove NAs

# Step 6: Plot regression per beta cluster
ggplot(merged_df, aes(x = hba1c, y = Disease_Score, color = beta_cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.9) +
  stat_cor(
    method = "pearson",
    size = 4,
    aes(group = beta_cluster),
    label.x.npc = "left",
    label.y.npc = "top",
    parse = FALSE  # 🔒 properly disables expression parsing
  ) +
  scale_color_manual(values = c(
    "β1" = "firebrick",   # darkred
    "β2" = "orchid",      # orange
    "β3" = "darkorange",  # blue
    "β4" = "dodgerblue"   # purple
  )) +
  labs(
    title = "Panel C: β-cell Disease Score vs HbA1c by Cluster",
    x = "HbA1c (%)",
    y = "Mean Disease Score (Per Donor)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

# Beta cell genes
# ============================================
# DotPlot by Gene Category and Sex — Sex-aware aggregation
# ============================================

suppressPackageStartupMessages({
  library(qs)       # for qread
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- 0) Load object & checks ---
beta_cells <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)")
stopifnot(all(c("beta_cluster","diabetes_status","Sex") %in% colnames(beta_cells@meta.data)))
DefaultAssay(beta_cells) <- "RNA"

# --- 1) Define gene categories (no overlaps) ---
gene_categories_raw <- list(
  "Insulin Secretion"   = c("INS", "PCSK1", "MAFA", "VAMP2", "SYN1"),
  "Glucose Sensing"     = c("GCK", "SLC2A2", "GLP1R", "KCNJ11", "ABCC8", "KCNJ3"),
  "Beta Cell Identity"  = c("PDX1", "ISL1", "ARX", "PAX4", "HNF1A", "HNF4A", "RFX6"),
  "Beta Cell Stress"    = c("DDIT3", "ATF4", "GNAS"),
  "Vesicle Exocytosis"  = c("STX1A", "CHGA", "TTR"),
  "Endocrine Signaling" = c("SST", "GAD1", "GAD2", "VEGFA")
)

# --- 2) Keep only genes present in object; warn on drops; drop empty categories ---
filter_present <- function(gvec) intersect(unique(gvec), rownames(beta_cells))
gene_categories <- lapply(gene_categories_raw, filter_present)

dropped <- setdiff(unlist(gene_categories_raw), unlist(gene_categories))
if (length(dropped)) message("Warning: dropped (not found): ", paste(dropped, collapse = ", "))

gene_categories <- gene_categories[vapply(gene_categories, length, 1L) > 0]
stopifnot(length(gene_categories) > 0)

genes_to_plot <- unlist(gene_categories, use.names = FALSE)

# --- 3) Normalize Sex and build Sex-aware grouping ---
beta_cells$Sex <- as.character(beta_cells$Sex)
beta_cells$Sex <- ifelse(beta_cells$Sex %in% c("Male","M"), "M",
                         ifelse(beta_cells$Sex %in% c("Female","F"), "F", beta_cells$Sex))

beta_cells$group_id_base <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, sep = "_")
beta_cells$group_id_sex  <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, beta_cells$Sex, sep = "_")

# --- 4) Build DotPlot with Sex-aware aggregation ---
p_raw <- DotPlot(
  beta_cells,
  features  = genes_to_plot,
  group.by  = "group_id_sex",
  dot.scale = 6,
  scale     = TRUE,
  col.min   = -1,
  col.max   =  1
)
dotdata <- p_raw$data

# --- 5) Restore group parts & annotate categories ---
dotdata <- dotdata %>%
  tidyr::separate(id, into = c("beta_cluster","diabetes_status","Sex"), sep = "_", remove = FALSE) %>%
  mutate(
    group_base = paste(beta_cluster, diabetes_status, sep = "_"),
    Sex = factor(Sex, levels = c("M","F"), labels = c("Male","Female"))
  )

gene_category_map <- stack(gene_categories)
colnames(gene_category_map) <- c("gene","gene_category")
dotdata <- left_join(dotdata, gene_category_map, by = c("features.plot" = "gene")) %>%
  mutate(gene_category = factor(gene_category, levels = names(gene_categories)))

# --- 6) Set β-group ordering (β4..β1; T2D then ND — match your original) ---
dotdata$group_base <- factor(
  dotdata$group_base,
  levels = c("β4_T2D", "β3_T2D", "β2_T2D",
             "β4_ND",  "β3_ND",  "β2_ND",  "β1_ND")
)

# Preserve the gene order (as listed within categories)
dotdata$features.plot <- factor(dotdata$features.plot, levels = genes_to_plot)

# --- 7) Plot (facet by category × Sex) ---
fig_gc_sexaware <- ggplot(dotdata, aes(x = features.plot, y = group_base, size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, color = "black") +
  # With scale=TRUE above, the natural midpoint is ~0
  scale_fill_gradient2(low = "white", high = "red", midpoint = -1) +
  facet_grid(. ~ gene_category + Sex, scales = "free_x", space = "free_x") +
  labs(
    title = "DotPlot by Gene Category and Sex",
    x = "Gene",
    y = "β-cell Subtype × Diabetes",
    fill = "Avg Expr (scaled)",
    size = "% Expressed"
  ) +
  theme_light() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 10, face = "bold", colour = "black"),
    axis.text.y   = element_text(size = 10, face = "bold", colour = "black"),
    strip.text.x  = element_text(size = 12, face = "bold"),
    plot.title    = element_text(size = 14, face = "bold"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(fig_gc_sexaware)

#DGE AND ORA
# DGE AND ORA
# Load object
beta_cells <- qread("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/seurat_objects/beta_cells.qs")

# Composite identity
beta_cells$beta_sex_diab <- paste(
  beta_cells$beta_cluster,
  beta_cells$Sex,
  beta_cells$diabetes_status,
  sep = "_"
)
Idents(beta_cells) <- "beta_sex_diab"

# Unique values
beta_clusters <- unique(beta_cells$beta_cluster)
sexes <- unique(beta_cells$Sex)
cell_types <- unique(beta_cells$beta_sex_diab)

# Output directory
output_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/cross_cluster"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- (A) Within-status comparisons ----
for (sex in sexes) {
  for (status in c("ND", "T2D")) {
    for (target in setdiff(beta_clusters, "β1")) {
      ident_1 <- paste(target, sex, status, sep = "_")
      ident_2 <- paste("β1", sex, status, sep = "_")
      comp_name <- paste0(target, "_vs_β1_", sex, "_", status)
      if (ident_1 %in% cell_types & ident_2 %in% cell_types) {
        message("Testing: ", ident_1, " vs ", ident_2)
        results <- tryCatch({
          FindMarkers(
            object = beta_cells,
            ident.1 = ident_1,
            ident.2 = ident_2,
            test.use = "wilcox",
            min.pct = 0.1,
            logfc.threshold = 0.137504, # 10% change
            pseudocount.use = 1,
            assay = "RNA",
            group.by = "beta_sex_diab",
            only.pos = FALSE,
            return.thresh = 0.1
          )
        }, error = function(e) {
          message("Error in comparison: ", comp_name, " : ", e$message)
          NULL
        })
        if (!is.null(results)) {
          write.csv(results, file = file.path(output_dir, paste0(comp_name, "_results.csv")), row.names = TRUE)
        }
      } else {
        message("Skipping (no data): ", ident_1, " or ", ident_2)
      }
    }
  }
}

# ---- (B) Cross-status: T2D clusters vs β1 ND ----
for (sex in sexes) {
  nd_ref <- paste("β1", sex, "ND", sep = "_")
  for (target in setdiff(beta_clusters, "β1")) {
    t2d_cluster <- paste(target, sex, "T2D", sep = "_")
    comp_name <- paste0(target, "_T2D_vs_β1_", sex, "_ND")
    if (t2d_cluster %in% cell_types & nd_ref %in% cell_types) {
      message("Testing: ", t2d_cluster, " vs ", nd_ref)
      results <- tryCatch({
        FindMarkers(
          object = beta_cells,
          ident.1 = t2d_cluster,
          ident.2 = nd_ref,
          test.use = "wilcox",
          min.pct = 0.1,
          logfc.threshold = 0.137504, # 10% change
          pseudocount.use = 1,
          assay = "RNA",
          group.by = "beta_sex_diab",
          only.pos = FALSE,
          return.thresh = 0.1
        )
      }, error = function(e) {
        message("Error in comparison: ", comp_name, " : ", e$message)
        NULL
      })
      if (!is.null(results)) {
        write.csv(results, file = file.path(output_dir, paste0(comp_name, "_results.csv")), row.names = TRUE)
      }
    } else {
      message("Skipping (no data): ", t2d_cluster, " or ", nd_ref)
    }
  }
}


#DGE
# 1. Load Object & Set Composite Identity
# Rememebr we will use beta_cells as the suerat object
beta_cells <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)")

# Create composite group: beta_cluster + Sex + diabetes_status (e.g., β1_M_ND)
beta_cells$beta_sex_diab <- paste(
  beta_cells$beta_cluster,
  beta_cells$Sex,
  beta_cells$diabetes_status,
  sep = "_"
)

Idents(beta_cells) <- "beta_sex_diab"
cell_types <- unique(beta_cells$beta_sex_diab)

# 2. Build ND vs T2D Comparisons for Each Subtype × Sex
comparisons <- list()
beta_clusters <- unique(beta_cells$beta_cluster)
sexes <- unique(beta_cells$Sex)

for (clust in beta_clusters) {
  for (sex in sexes) {
    nd_id  <- paste(clust, sex, "ND", sep = "_")
    t2d_id <- paste(clust, sex, "T2D", sep = "_")
    if (nd_id %in% cell_types & t2d_id %in% cell_types) {
      comp_name <- paste(clust, sex, "T2D_vs_ND", sep = "_")
      comparisons[[comp_name]] <- c(t2d_id, nd_id)
    }
  }
}

# 3. Run FindMarkers for Each Comparison
output_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/by_beta_sex"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (comp_name in names(comparisons)) {
  ident.1 <- comparisons[[comp_name]][1] # T2D
  ident.2 <- comparisons[[comp_name]][2] # ND
  print(paste("Testing:", ident.1, "vs", ident.2))
  results <- tryCatch({
    FindMarkers(
      object = beta_cells,
      ident.1 = ident.1,
      ident.2 = ident.2,
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.137504, # 10% change
      pseudocount.use = 1,
      assay = "RNA",
      group.by = "beta_sex_diab",
      only.pos = FALSE,
      return.thresh = 0.1
    )
  }, error = function(e) {
    print(paste("Error in comparison:", comp_name, ":", e$message))
    return(NULL)
  })
  if (!is.null(results)) {
    results$comparison <- comp_name
    filename <- file.path(output_dir, paste0(comp_name, "_results.csv"))
    write.csv(results, file = filename, row.names = TRUE)
  }
}


# List the DGE folders to analyze
dge_folders <- list(
  by_beta_sex = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/by_beta_sex",
  beta_vs_beta = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/cross_cluster/beta_vs_beta",
  T2D_vs_beta1 = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/cross_cluster/T2D_vs_beta1"
)

# Set base ORA output directory
ora_base <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ORA"

for (type in names(dge_folders)) {
  input_dir <- dge_folders[[type]]
  output_up <- file.path(ora_base, type, "UP")
  output_down <- file.path(ora_base, type, "DOWN")
  
  # Create output directories if not exist
  dir.create(output_up, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_down, showWarnings = FALSE, recursive = TRUE)
  
  # List CSV files
  dge_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for (file in dge_files) {
    # Read the file
    dat <- read.csv(file, row.names = 1)
    # Standardize avg_log2FC column if needed
    if (!"avg_log2FC" %in% names(dat)) {
      logfc_col <- grep("log2FC", names(dat), value = TRUE)
      if (length(logfc_col) == 1) names(dat)[names(dat) == logfc_col] <- "avg_log2FC"
    }
    if (!"p_val_adj" %in% names(dat)) {
      adjp_col <- grep("adj", names(dat), value = TRUE)
      if (length(adjp_col) == 1) names(dat)[names(dat) == adjp_col] <- "p_val_adj"
    }
    # Get gene names
    dat$gene <- rownames(dat)
    # Genes up/down
    up <- filter(dat, p_val_adj < 0.05 & avg_log2FC > 0)$gene
    down <- filter(dat, p_val_adj < 0.05 & avg_log2FC < 0)$gene
    
    # GO for upregulated
    if (length(up) > 0) {
      go_up <- gost(up, organism = "hsapiens", significant = TRUE, user_threshold = 0.05,
                    correction_method = "fdr", domain_scope = "annotated", sources = "GO", evcodes = TRUE)
      go_up_res <- as.data.frame(go_up$result)
      # Optional: rename FDR col
      if ("p_value" %in% names(go_up_res)) names(go_up_res)[names(go_up_res) == "p_value"] <- "hypergeometric FDR"
      # Remove columns (if present)
      drop_cols <- c("parents", "source_order", "effective_domain_size", "query", "precision", "recall", "evidence_codes")
      go_up_res <- select(go_up_res, -any_of(drop_cols))
      # Write
      outfile <- file.path(output_up, paste0(basename(file)))
      write.csv(go_up_res, outfile, row.names = FALSE)
    }
    # GO for downregulated
    if (length(down) > 0) {
      go_down <- gost(down, organism = "hsapiens", significant = TRUE, user_threshold = 0.05,
                      correction_method = "fdr", domain_scope = "annotated", sources = "GO", evcodes = TRUE)
      go_down_res <- as.data.frame(go_down$result)
      if ("p_value" %in% names(go_down_res)) names(go_down_res)[names(go_down_res) == "p_value"] <- "hypergeometric FDR"
      go_down_res <- select(go_down_res, -any_of(drop_cols))
      outfile <- file.path(output_down, paste0(basename(file)))
      write.csv(go_down_res, outfile, row.names = FALSE)
    }
  }
}

cat("GO/ORA completed for all DGE folders!\n")

# SANKEY
# ---- FUNCTION ----
make_sankey_plot <- function(
    ora_dir,
    comp_names,
    pathways = NULL,         # Character vector of pathways (priority)
    top_n_pathways = NULL,   # Integer: take top N pathways (by sum(-log10(FDR)))
    top_n_genes = 5,
    comparison_colors = NULL,
    pathway_colors = NULL,
    plot_width = 1000,
    plot_height = 800
) {
  # Load required libraries
  require(readr)
  require(dplyr)
  require(tidyr)
  require(networkD3)
  require(htmlwidgets)
  require(stringr)
  
  # --- Build full dataframe ---
  all_rows <- list()
  for (fname in comp_names) {
    fpath <- file.path(ora_dir, fname)
    if (!file.exists(fpath)) {
      warning(sprintf("File not found: %s", fpath))
      next
    }
    df <- read_csv(fpath, show_col_types = FALSE)
    if (!"intersection" %in% colnames(df)) next
    df <- df %>%
      filter(!is.na(`hypergeometric FDR`), `hypergeometric FDR` < 0.05, !is.na(intersection), intersection != "")
    if (nrow(df) == 0) next
    for (i in seq_len(nrow(df))) {
      this <- df[i,]
      genes <- str_split(this$intersection, ",")[[1]] %>% str_trim()
      genes <- head(genes, top_n_genes)
      all_rows[[length(all_rows)+1]] <- data.frame(
        Comparison = gsub("_results.csv", "", fname),
        Pathway = this$term_name,
        FDR = this$`hypergeometric FDR`,
        Genes = I(list(genes)),
        stringsAsFactors = FALSE
      )
    }
  }
  sankey_df <- bind_rows(all_rows)
  if (nrow(sankey_df) == 0) stop("No data for selected comparisons/files!")
  
  # --- Filter for pathways ---
  if (!is.null(pathways)) {
    sankey_df <- sankey_df %>% filter(Pathway %in% pathways)
    if (nrow(sankey_df) == 0) stop("No data for selected pathway(s)!")
  } else if (!is.null(top_n_pathways)) {
    path_scores <- sankey_df %>%
      group_by(Pathway) %>%
      summarize(total_score = sum(-log10(FDR)), .groups = "drop") %>%
      arrange(desc(total_score)) %>%
      slice_head(n = top_n_pathways)
    sankey_df <- sankey_df %>% filter(Pathway %in% path_scores$Pathway)
  }
  if (nrow(sankey_df) == 0) stop("No data for pathway selection!")
  
  # --- Build nodes ---
  comparison_names <- unique(sankey_df$Comparison)
  pathway_names <- unique(sankey_df$Pathway)
  gene_names <- unique(unlist(sankey_df$Genes))
  nodes <- data.frame(
    name = c(comparison_names, pathway_names, gene_names),
    group = NA,
    stringsAsFactors = FALSE
  )
  
  # --- Assign groups and colors ---
  # Colors: fallback if not provided
  if (is.null(comparison_colors)) {
    pal <- RColorBrewer::brewer.pal(8, "Dark2")
    comparison_colors <- setNames(rep(pal, length.out = length(comparison_names)), comparison_names)
  }
  if (is.null(pathway_colors)) {
    pal2 <- RColorBrewer::brewer.pal(8, "Set1")
    pathway_colors <- setNames(rep(pal2, length.out = length(pathway_names)), pathway_names)
  }
  
  nodes$group <- nodes$name
  for (gene in gene_names) {
    which_rows <- sapply(sankey_df$Genes, function(glist) gene %in% glist)
    pathway <- sankey_df$Pathway[which_rows][1]
    nodes$group[nodes$name == gene] <- pathway
  }
  all_groups <- unique(nodes$group)
  color_vec <- c(comparison_colors, pathway_colors)
  fallback_colors <- rep("#BDBDBD", sum(!(all_groups %in% names(color_vec))))
  color_vec <- c(color_vec[all_groups[all_groups %in% names(color_vec)]], fallback_colors)
  names(color_vec) <- all_groups
  colJS <- paste0(
    "d3.scaleOrdinal().domain([", 
    paste(shQuote(names(color_vec)), collapse=","), 
    "]).range([", 
    paste(shQuote(color_vec), collapse=","), 
    "])"
  )
  
  # --- Build links ---
  links1 <- sankey_df %>%
    mutate(
      source = match(Comparison, nodes$name) - 1,
      target = match(Pathway, nodes$name) - 1,
      value = -log10(FDR) + 1,
      LinkGroup = Comparison
    ) %>%
    select(source, target, value, LinkGroup)
  links2 <- sankey_df %>%
    select(Pathway, Genes) %>%
    unnest(Genes) %>%
    mutate(
      source = match(Pathway, nodes$name) - 1,
      target = match(Genes, nodes$name) - 1,
      value = 1,
      LinkGroup = Pathway
    ) %>%
    select(source, target, value, LinkGroup)
  links <- bind_rows(links1, links2)
  
  # --- Make Sankey plot object ---
  sankey <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    NodeGroup = "group",
    LinkGroup = "LinkGroup",
    fontSize = 16,
    nodeWidth = 30,
    width = plot_width,
    height = plot_height,
    sinksRight = FALSE,
    colourScale = JS(colJS)
  )
  return(sankey)
}

# PAss paths
my_comparison_colors <- c(
  "β2_M_T2D_vs_ND" = "#123456",
  "β2_F_T2D_vs_ND" = "#C2185B",
  "β3_M_T2D_vs_ND" = "#388E3C",
  "β3_F_T2D_vs_ND" = "#FBC02D",
  "β4_M_T2D_vs_ND" = "#7B1FA2",
  "β4_F_T2D_vs_ND" = "#E64A19"
)

my_pathway_colors <- c(
  "cytosolic ribosome" = "#e74c3c",
  "ribosomal subunit" = "#8D6E63",
  "cytoplasmic translation" = "#43A047",
  "cytoplasm" = "#1976D2",
  "response to peptide hormone" = "#F57C00"
)

# Example vector of pathway names
my_pathways <- c(
  "cytosolic ribosome",
  "ribosomal subunit",
  "cytoplasmic translation",
  "cytoplasm",
  "response to peptide hormone"
)

my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/by_beta_sex/DOWN"
my_comp_names <- c(
  "β2_M_T2D_vs_ND_results.csv",
  "β2_F_T2D_vs_ND_results.csv",
  "β3_M_T2D_vs_ND_results.csv",
  "β3_F_T2D_vs_ND_results.csv",
  "β4_M_T2D_vs_ND_results.csv",
  "β4_F_T2D_vs_ND_results.csv"
)

my_sankey <- make_sankey_plot(
  ora_dir = my_ora_dir,
  comp_names = my_comp_names,
  pathways = my_pathways,     # <--- Only this pathway will show
  comparison_colors = my_comparison_colors,
  pathway_colors = my_pathway_colors,
  top_n_genes = 5,
  plot_width = 1200,
  plot_height = 800
)
my_sankey

# DOTPLOTS
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)

# ---- INPUT ----
#GROUP ============== β_vs_β1 UP
base_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/beta_vs_beta/UP"

up_files <- c(
  "β2_vs_β1_F_ND_results.csv",
  "β2_vs_β1_M_ND_results.csv",
  "β3_vs_β1_F_ND_results.csv",
  "β3_vs_β1_M_ND_results.csv",
  "β4_vs_β1_F_ND_results.csv",
  "β4_vs_β1_M_ND_results.csv"
)
comp_names <- c(
  "β2_F_vs_β1_F_ND",
  "β2_M_vs_β1_M_ND",
  "β3_F_vs_β1_F_ND",
  "β3_M_vs_β1_M_ND",
  "β4_F_vs_β1_F_ND",
  "β4_M_vs_β1_M_ND"
)

# Read & combine
up_results <- map2_dfr(
  up_files, comp_names,
  ~ read.csv(file.path(base_dir, .x)) %>% mutate(comparison = .y)
)

# Select pathways
selected_pathways <- c("translation", "DNA-templated transcription", "positive regulation of RNA splicing",
                       "programmed cell death", "apoptotic process", "autophagy", 
                       "cell cycle", "response to hypoxia",
                       "lipid homeostasis", "cholesterol biosynthetic process",
                       "focal adhesion", "growth",
                       "glycogen metabolic process", "pyruvate transport",
                       "integrated stress response signaling", "cellular response to unfolded protein", "endoplasmic reticulum unfolded protein response"
                       )

# Category assignment (edit as needed)
up_results <- up_results %>%
  mutate(category = case_when(
    term_name %in% c("translation", "DNA-templated transcription", "positive regulation of RNA splicing") ~ "RNA/Protein Production",
    term_name %in% c("programmed cell death", "apoptotic process", "autophagy") ~ "Cell Death",
    term_name %in% c("cell cycle", "response to hypoxia") ~ "Stemness",
    term_name %in% c("lipid homeostasis", "cholesterol biosynthetic process") ~ "Lipid Homeostasis",
    term_name %in% c("focal adhesion", "growth") ~ "ECM",
    term_name %in% c("glycogen metabolic process", "pyruvate transport") ~ "Energy Metabolism",
    term_name %in% c("integrated stress response signaling", "cellular response to unfolded protein", "endoplasmic reticulum unfolded protein response") ~ "Misfolded Protein Response",
    TRUE ~ "Other"
  ))

filtered_results <- up_results %>% filter(term_name %in% selected_pathways)

# Faceting and color style (DotPlot-style)
filtered_results$comparison <- factor(filtered_results$comparison, levels = comp_names)
filtered_results$category <- factor(filtered_results$category, levels = unique(filtered_results$category))
filtered_results$term_name <- factor(filtered_results$term_name, levels = rev(selected_pathways))

# "DotPlot"-style color
my_blue <- "dodgerblue3"
my_red  <- "firebrick2"

p <- ggplot(filtered_results, aes(x = comparison, y = term_name, size = intersection_size, fill = hypergeometric.FDR)) +
  geom_point(shape = 21, color = "black", stroke = 0.4) +
  scale_fill_gradient2(
    low = my_red, high = my_blue,
    midpoint = 0.05, name = "Hypergeometric FDR"
  ) +
  scale_size_continuous(range = c(3.5, 8), name = "Intersection Size") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +  # << vertical split!
  labs(
    title = "UPregulated Pathways: β2–4 vs β1 (ND)",
    x = "Comparison",
    y = "Pathway Name"
  ) +
  theme_light(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", colour = "black", size = 13),
    axis.text.y = element_text(size = 13, face = "bold", colour = "black"),
    strip.text.y.left = element_text(angle = 0, size = 15, face = "bold"), # strip on left
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.position = "right",
    panel.spacing.y = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p

# ---- INPUT ----
#GROUP ============== β_vs_β1 DOWN
base_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/beta_vs_beta/DOWN"

up_files <- c(
  "β2_vs_β1_F_ND_results.csv",
  "β2_vs_β1_M_ND_results.csv",
  "β3_vs_β1_F_ND_results.csv",
  "β3_vs_β1_M_ND_results.csv",
  "β4_vs_β1_F_ND_results.csv",
  "β4_vs_β1_M_ND_results.csv"
)
comp_names <- c(
  "β2_F_vs_β1_F_ND",
  "β2_M_vs_β1_M_ND",
  "β3_F_vs_β1_F_ND",
  "β3_M_vs_β1_M_ND",
  "β4_F_vs_β1_F_ND",
  "β4_M_vs_β1_M_ND"
)

# Read & combine
up_results <- map2_dfr(
  up_files, comp_names,
  ~ read.csv(file.path(base_dir, .x)) %>% mutate(comparison = .y)
)

# Select pathways
selected_pathways <- c("electron transport chain", "oxidative phosphorylation", "ATP biosynthetic process", "mitochondrion organization",
                       "response to unfolded protein", "'de novo' protein folding", "ERAD pathway",
                       "insulin secretion", "insulin processing", "response to glucose", "response to glucagon",
                       "canonical glycolysis", "pyruvate metabolic process", "tricarboxylic acid cycle", "glucose homeostasis", "calcium ion transport",
                       "fatty acid beta-oxidation", "glycolipid metabolic process"
)

# Category assignment (edit as needed)
up_results <- up_results %>%
  mutate(category = case_when(
    term_name %in% c("electron transport chain", "oxidative phosphorylation", "ATP biosynthetic process", "mitochondrion organization") ~ "ETC",
    term_name %in% c("response to unfolded protein", "'de novo' protein folding", "ERAD pathway") ~ "Protein Folding",
    term_name %in% c("insulin secretion", "insulin processing", "response to glucose", "response to glucagon") ~ "Insulin Secretion",
    term_name %in% c("canonical glycolysis", "pyruvate metabolic process", "tricarboxylic acid cycle", "glucose homeostasis", "calcium ion transport") ~ "Metabolism",
    term_name %in% c("fatty acid beta-oxidation", "glycolipid metabolic process") ~ "Fat Oxidation",
    TRUE ~ "Other"
  ))

filtered_results <- up_results %>% filter(term_name %in% selected_pathways)

# Faceting and color style (DotPlot-style)
filtered_results$comparison <- factor(filtered_results$comparison, levels = comp_names)
filtered_results$category <- factor(filtered_results$category, levels = unique(filtered_results$category))
filtered_results$term_name <- factor(filtered_results$term_name, levels = rev(selected_pathways))

# "DotPlot"-style color
my_blue <- "dodgerblue3"
my_red  <- "firebrick2"

p <- ggplot(filtered_results, aes(x = comparison, y = term_name, size = intersection_size, fill = hypergeometric.FDR)) +
  geom_point(shape = 21, color = "black", stroke = 0.4) +
  scale_fill_gradient2(
    low = my_red, high = my_blue,
    midpoint = 0.05, name = "Hypergeometric FDR"
  ) +
  scale_size_continuous(range = c(3.5, 8), name = "Intersection Size") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +  # << vertical split!
  labs(
    title = "UDOWNregulated Pathways: β2–4 vs β1 (ND)",
    x = "Comparison",
    y = "Pathway Name"
  ) +
  theme_light(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", colour = "black", size = 13),
    axis.text.y = element_text(size = 13, face = "bold", colour = "black"),
    strip.text.y.left = element_text(angle = 0, size = 15, face = "bold"), # strip on left
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.position = "right",
    panel.spacing.y = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p


# # ---- INPUT ----
# #GROUP ============== βT2D_vs_βND DOWN
# library(networkD3)
# library(dplyr)
# library(tidyr)
# library(readr)
# library(stringr)
# library(scales)
# 
# make_sankey_plot <- function(
#     ora_dir,
#     comp_names,
#     pathways = NULL,
#     top_n_pathways = NULL,
#     show_genes = FALSE,
#     comparison_colors = NULL,
#     plot_width = 1000,
#     plot_height = 800
# ) {
#   all_rows <- list()
#   for (fname in comp_names) {
#     fpath <- file.path(ora_dir, fname)
#     if (!file.exists(fpath)) next
#     df <- read_csv(fpath, show_col_types = FALSE)
#     if (!"intersection" %in% colnames(df)) next
#     df <- df %>% filter(!is.na(`hypergeometric FDR`), `hypergeometric FDR` < 0.05, !is.na(intersection), intersection != "")
#     if (nrow(df) == 0) next
#     for (i in seq_len(nrow(df))) {
#       this <- df[i,]
#       genes <- str_split(this$intersection, ",")[[1]] %>% str_trim()
#       all_rows[[length(all_rows) + 1]] <- data.frame(
#         Comparison = gsub("_results.csv", "", fname),
#         Pathway = this$term_name,
#         FDR = this$`hypergeometric FDR`,
#         Genes = I(list(genes)),
#         stringsAsFactors = FALSE
#       )
#     }
#   }
#   sankey_df <- bind_rows(all_rows)
#   if (nrow(sankey_df) == 0) stop("No data for selected comparisons/files!")
#   
#   # Filter for pathways
#   if (!is.null(pathways)) {
#     sankey_df <- sankey_df %>% filter(Pathway %in% pathways)
#     if (nrow(sankey_df) == 0) stop("No data for selected pathway(s)!")
#   } else if (!is.null(top_n_pathways)) {
#     path_scores <- sankey_df %>%
#       group_by(Pathway) %>%
#       summarize(total_score = sum(-log10(FDR)), .groups = "drop") %>%
#       arrange(desc(total_score)) %>%
#       slice_head(n = top_n_pathways)
#     sankey_df <- sankey_df %>% filter(Pathway %in% path_scores$Pathway)
#   }
#   if (nrow(sankey_df) == 0) stop("No data for pathway selection!")
#   
#   comparison_names <- unique(sankey_df$Comparison)
#   pathway_names <- unique(sankey_df$Pathway)
#   
#   # Assign comparison node colors
#   if (is.null(comparison_colors)) {
#     pal <- RColorBrewer::brewer.pal(8, "Dark2")
#     comparison_colors <- setNames(rep(pal, length.out = length(comparison_names)), comparison_names)
#   }
#   comp_colors <- unname(comparison_colors[comparison_names])
#   
#   # Assign pathway node colors by FDR
#   pathway_fdr <- sankey_df %>% group_by(Pathway) %>% summarize(FDR = min(FDR), .groups = "drop")
#   col_func <- col_numeric(c("dodgerblue3", "white", "firebrick2"),
#                           domain = range(-log10(pathway_fdr$FDR), na.rm = TRUE))
#   path_colors <- col_func(-log10(pathway_fdr$FDR))
#   
#   # Key: use hex color *as group*, so node group = node color!
#   nodes <- data.frame(
#     name = c(comparison_names, pathway_names),
#     group = c(comp_colors, path_colors),  # <--- Group = Color!
#     color = c(comp_colors, path_colors),
#     stringsAsFactors = FALSE
#   )
#   
#   # Assemble links as usual
#   links <- sankey_df %>%
#     mutate(
#       source = match(Comparison, nodes$name) - 1,
#       target = match(Pathway, nodes$name) - 1,
#       value = -log10(FDR) + 1,
#       LinkGroup = Comparison
#     ) %>%
#     select(source, target, value, LinkGroup)
#   
#   # JS: color scale domain and range both hex codes!
#   hexes <- unique(nodes$group)
#   colJS <- paste0(
#     "d3.scaleOrdinal().domain([",
#     paste0("'", hexes, "'", collapse = ","),
#     "]).range([",
#     paste0("'", hexes, "'", collapse = ","),
#     "])"
#   )
#   
#   # PLOT
#   sankeyNetwork(
#     Links = links,
#     Nodes = nodes,
#     Source = "source",
#     Target = "target",
#     Value = "value",
#     NodeID = "name",
#     NodeGroup = "group",   # Group is now color!
#     LinkGroup = "LinkGroup",
#     fontSize = 16,
#     nodeWidth = 30,
#     width = plot_width,
#     height = plot_height,
#     sinksRight = FALSE,
#     colourScale = JS(colJS)
#   )
# }
# 
# # ==== USAGE ====
# my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/by_beta_sex/DOWN"
# my_comp_names <- c(
#   "β2_M_T2D_vs_ND_results.csv",
#   "β2_F_T2D_vs_ND_results.csv",
#   "β3_M_T2D_vs_ND_results.csv",
#   "β3_F_T2D_vs_ND_results.csv",
#   "β4_M_T2D_vs_ND_results.csv",
#   "β4_F_T2D_vs_ND_results.csv"
# )
# my_comparison_colors <- c(
#   "β2_M_T2D_vs_ND" = "#123456",
#   "β2_F_T2D_vs_ND" = "#C2185B",
#   "β3_M_T2D_vs_ND" = "#388E3C",
#   "β3_F_T2D_vs_ND" = "#FBC02D",
#   "β4_M_T2D_vs_ND" = "#7B1FA2",
#   "β4_F_T2D_vs_ND" = "#E64A19"
# )
# my_pathways <- c(
#   "cytosolic ribosome",
#   "ribosomal subunit",
#   "cytoplasmic translation",
#   "cytoplasm",
#   "response to peptide hormone"
# )
# 
# # ==== RUN ====
# my_sankey <- make_sankey_plot(
#   ora_dir = my_ora_dir,
#   comp_names = my_comp_names,
#   pathways = my_pathways,
#   show_genes = FALSE,
#   comparison_colors = my_comparison_colors,
#   plot_width = 1200,
#   plot_height = 800
# )
# my_sankey



# MAKE REPRODUCIBLE SANKEY PLOTS
# Install if needed:
# remotes::install_github("davidsjoberg/ggsankey")
# install.packages(c("readr", "dplyr", "ggalluvial", "ggplot2", "cowplot"))
make_sankey_plot <- function(my_comparison_colors, my_pathway_colors, my_ora_dir, my_pathways) {
  library(readr)
  library(dplyr)
  library(ggalluvial)
  library(ggplot2)
  library(cowplot)
  library(scales)
  
  nukefix <- function(x) {
    x <- gsub("β", "b", x)
    x <- trimws(x)
    x <- iconv(x, to = "ASCII//TRANSLIT")
    x <- tolower(x)
    gsub("\\s+", " ", x)
  }
  
  # --- Get file names ---
  my_comp_names <- list.files(my_ora_dir, pattern = "_results\\.csv$", full.names = FALSE)
  
  all_links <- list()
  for (fname in my_comp_names) {
    comp_raw <- gsub("_results.csv", "", fname)
    comp <- nukefix(comp_raw)
    fpath <- file.path(my_ora_dir, fname)
    if (!file.exists(fpath)) next
    df <- read_csv(fpath, show_col_types = FALSE)
    if (!"term_name" %in% colnames(df)) next
    df <- df %>% filter(!is.na(`hypergeometric FDR`), `hypergeometric FDR` < 0.05, !is.na(term_name))
    df$Pathway <- nukefix(df$term_name)
    for (p in my_pathways) {
      pw <- nukefix(p)
      match_row <- df[df$Pathway == pw,]
      if (nrow(match_row) == 0) next
      all_links[[length(all_links)+1]] <- data.frame(
        Comparison = comp,
        Pathway = pw,
        Value = -log10(match_row$`hypergeometric FDR`[1]) + 1,
        stringsAsFactors = FALSE
      )
    }
  }
  
  sankey_df <- dplyr::bind_rows(all_links)
  sankey_df$Comparison <- factor(sankey_df$Comparison, levels = names(my_comparison_colors))
  sankey_df$Pathway <- factor(sankey_df$Pathway, levels = nukefix(my_pathways))
  
  # --- Main Sankey plot ---
  p_main <- ggplot(sankey_df, aes(axis1 = Comparison, axis2 = Pathway, y = Value)) +
    geom_alluvium(aes(fill = Comparison), width = 1/12, alpha = 0.8) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_discrete(limits = c("Comparison", "Pathway"), expand = c(.1, .1)) +
    scale_fill_manual(
      values = c(my_comparison_colors, my_pathway_colors),
      breaks = names(my_comparison_colors)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold")
    )
  
  # --- Dot legend ---
  fdr_values <- sankey_df$Value
  legend_vals <- pretty(fdr_values, n = 5)
  legend_vals <- sort(unique(legend_vals[legend_vals > 0]), decreasing = TRUE)
  
  rescale_size <- function(x, to = c(3, 10)) {
    rng <- range(x, na.rm = TRUE)
    scales::rescale(x, to = to, from = rng)
  }
  
  dot_sizes <- rescale_size(legend_vals, to = c(3, 10))
  
  legend_df <- data.frame(
    x = 1,
    y = seq(50, 38, length.out = length(dot_sizes)),
    y_label = seq(50, 38, length.out = length(dot_sizes)),
    size = dot_sizes,
    label = paste0("-log10(FDR)+1 = ", legend_vals)
  )
  
  p_legend <- ggplot(legend_df) +
    geom_point(aes(x = x, y = y, size = size), shape = 21, fill = "steelblue", color = "black", stroke = 0.25) +
    geom_text(aes(x = x + 0.4, y = y_label, label = label), hjust = 0, vjust = 0.5, size = 4) +
    theme_void() +
    coord_cartesian(clip = "off") +
    scale_size_identity() +
    scale_x_continuous(limits = c(0.9, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(plot.margin = margin(5, 20, 5, 5))
  
  # --- Combine ---
  cowplot::plot_grid(
    p_main,
    p_legend,
    rel_widths = c(4.2, 1.2),
    nrow = 1,
    axis = "none",
    align = "none"
  )
}

## βT2D vs ND DOWN
## Load prereqs
# --- Colors ---
my_comparison_colors <- c(
  "b2_m_t2d_vs_nd" = "#274472",
  "b3_m_t2d_vs_nd" = "#406E8E",
  "b4_m_t2d_vs_nd" = "#6998B7",
  "b2_f_t2d_vs_nd" = "#C97D60",
  "b3_f_t2d_vs_nd" = "#B85C5C",
  "b4_f_t2d_vs_nd" = "#8C3333"
)

my_pathway_colors <- c(
  # Protein synthesis & ribosome (matte red)
  "translation"                        = "#C0392B",  # Matte muted red
  "rrna processing"                    = "#C0392B",
  
  # Oxphos & mitochondrial (matte orange)
  "electron transport chain"           = "#D35400",  # Matte orange
  "oxidative phosphorylation"          = "#D35400",
  "atp synthesis coupled electron transport" = "#D35400",
  "proton transmembrane transport"     = "#D35400",
  
  # Antioxidant/ROS (matte green)
  "cellular oxidant detoxification"    = "#27AE60",  # Muted green (same as before, but fits matte)
  "antioxidant activity"               = "#27AE60",
  
  # Insulin/hormone signaling (matte blue)
  "response to insulin"                = "#2980B9",  # Muted blue (less saturated than "#1976D2")
  "insulin receptor signaling pathway" = "#2980B9",
  "hormone secretion"                  = "#2980B9",
  "insulin processing"                 = "#2980B9",
  "insulin secretion"                  = "#2980B9",
  
  # Lipid metabolism (matte purple)
  "lipid metabolic process"            = "#7D3C98",  # Matte muted purple
  "response to fatty acid"             = "#7D3C98",
  
  # Glucose/glycemic (matte brown)
  "glucose homeostasis"                = "#8D6748",  # Matte brown
  "response to glucagon"               = "#8D6748",
  "glucose metabolic process"          = "#8D6748"
)

nukefix <- function(x) {
  x <- gsub("β", "b", x)
  x <- trimws(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("\\s+", " ", x)
}

# --- Input ---
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/by_beta_sex/DOWN"
my_comp_names <- c(
  "β2_M_T2D_vs_ND_results.csv",
  "β3_M_T2D_vs_ND_results.csv",
  "β4_M_T2D_vs_ND_results.csv",
  "β2_F_T2D_vs_ND_results.csv",
  "β3_F_T2D_vs_ND_results.csv",
  "β4_F_T2D_vs_ND_results.csv"
)

my_pathways <- c(
  "translation",
  "rRNA processing",
  "electron transport chain",
  "oxidative phosphorylation",
  "ATP synthesis coupled electron transport",
  "proton transmembrane transport",
  "cellular oxidant detoxification",
  "antioxidant activity",
  "response to insulin",
  "insulin receptor signaling pathway",
  "hormone secretion",
  "insulin processing",
  "insulin secretion",
  "lipid metabolic process",
  "response to fatty acid",
  "glucose homeostasis",
  "response to glucagon",
  "glucose metabolic process"
)

# RUN plot
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/by_beta_sex/DOWN",
  my_pathways = my_pathways
)

## βT2D vs ND UP
## Load prereqs
# --- Colors ---
my_comparison_colors <- c(
  #"b2_m_t2d_vs_nd" = "#274472",
  "b3_m_t2d_vs_nd" = "#406E8E",
  "b4_m_t2d_vs_nd" = "#6998B7",
  "b2_f_t2d_vs_nd" = "#C97D60",
  "b3_f_t2d_vs_nd" = "#B85C5C",
  "b4_f_t2d_vs_nd" = "#8C3333"
)

my_pathway_colors <- c(
  # Protein synthesis & ribosome (matte red)
  "translation"                        = "#C0392B",  # Matte muted red
  "protein metabolic process"                    = "#C0392B",
  
  # Mitophagy (matte brown)
  "macroautophagy"                = "#8D6748",
  "autophagy of mitochondrion"   = "#8D6748",
  "autophagy"                    = "#8D6748",
  
  # Oxphos & mitochondrial (matte orange)
  "atp biosynthetic process"              = "#D35400",
  "oxidative phosphorylation"             = "#D35400",
  "proton transmembrane transport"        = "#D35400",
  "aerobic electron transport chain"      = "#D35400",
  
  # UPR (matte green)
  "erad pathway"                         = "#27AE60",
  "perk-mediated unfolded protein response" = "#27AE60",
  "protein folding"                      = "#27AE60",
  "response to unfolded protein"         = "#27AE60",
  "response to endoplasmic reticulum stress" = "#27AE60",
  "response to topologically incorrect protein" = "#27AE60",
  "er overload response"                = "#27AE60",
  "positive regulation of proteolysis"  = "#27AE60",
  
  # Stress and Cell Death (matte blue)
  "intrinsic apoptotic signaling pathway" = "#2980B9",
  "response to oxidative stress"          = "#2980B9",
  "cell death"                            = "#2980B9",
  "apoptotic process"                     = "#2980B9",
  "cellular response to stress"           = "#2980B9",
  "response to stress"                    = "#2980B9",
  
  # Inflammation (matte purple)
  "p38MAPK cascade"                                 = "#7D3C98",
  "cellular response to interleukin-7"              = "#7D3C98",
  "response to virus"                               = "#7D3C98",
  "cellular response to interleukin-4"              = "#7D3C98",
  "interleukin-11-mediated signaling pathway"       = "#7D3C98",
  "positive regulation of NF-kappaB transcription factor activity" = "#7D3C98",
  "inflammatory response"                           = "#7D3C98",
  "tumor necrosis factor production"                = "#7D3C98",
  
  # Dedifferentiation / Stem-like Reprogramming (matte teal)
  "epithelial to mesenchymal transition"            = "#148F77",  # Matte teal
  "positive regulation of cell fate commitment"     = "#148F77",
  "neurogenesis"                                    = "#148F77",
  "chromatin remodeling"                            = "#148F77"
)


nukefix <- function(x) {
  x <- gsub("β", "b", x)
  x <- trimws(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("\\s+", " ", x)
}

# --- Input ---
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/by_beta_sex/UP"
my_comp_names <- c(
  #"β2_M_T2D_vs_ND_results.csv",
  "β3_M_T2D_vs_ND_results.csv",
  "β4_M_T2D_vs_ND_results.csv",
  "β2_F_T2D_vs_ND_results.csv",
  "β3_F_T2D_vs_ND_results.csv",
  "β4_F_T2D_vs_ND_results.csv"
)

my_pathways <- c(
  # Protein synthesis & ribosome
  "translation",
  "protein metabolic process",
  
  # Mitophagy
  "macroautophagy",
  "autophagy of mitochondrion",
  "autophagy",
  
  # Oxphos & mitochondrial
  "atp biosynthetic process",
  "oxidative phosphorylation",
  "proton transmembrane transport",
  "aerobic electron transport chain",
  
  # UPR
  "erad pathway",
  "perk-mediated unfolded protein response",
  "protein folding",
  "response to unfolded protein",
  "response to endoplasmic reticulum stress",
  "response to topologically incorrect protein",
  "er overload response",
  "positive regulation of proteolysis",
  
  # Stress and Cell Death
  "intrinsic apoptotic signaling pathway",
  "response to oxidative stress",
  "cell death",
  "apoptotic process",
  "cellular response to stress",
  "response to stress",
  
  # Inflammation
  "p38MAPK cascade",
  "cellular response to interleukin-7",
  "response to virus",
  "cellular response to interleukin-4",
  "interleukin-11-mediated signaling pathway",
  "positive regulation of NF-kappaB transcription factor activity",
  "inflammatory response",
  "tumor necrosis factor production",
  
  # Dedifferentiation / Stem-like
  "epithelial to mesenchymal transition",
  "positive regulation of cell fate commitment",
  "neurogenesis",
  "chromatin remodeling"
)

# RUN plot
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/by_beta_sex/UP",
  my_pathways = my_pathways
)

# Module scores and Gene expression
# Beta cell genes
# ============================================
# Module scores & Gene expression — All Pathways, Sex-aware DotPlot
# ============================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- 0) Load object ---
beta_cells <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)")
DefaultAssay(beta_cells) <- "RNA"

# --- 1) Define ALL pathways ---
all_genes_raw <- list(
  UPR                 = c("HSPA5", "XBP1", "ATF6", "ATF4", "DDIT3"),
  #Mitophagy           = c("PINK1", "PRKN", "BNIP3", "SQSTM1", "MAP1LC3B"),
  #OXPHOS              = c("NDUFS1", "SDHB", "UQCRC2", "COX4I1", "ATP5F1B")
  Antioxidant         = c("SOD2", "GPX4", "PRDX3", "NQO1", "HMOX1"),
  beta_cell_identity   = c("KCNJ3","MAFA","GLP1R","SYN1","VAMP2"),
  Glucose_Homeostasis = c("GCK","SLC2A2","PCSK1","KCNJ11","ABCC8"),
  #Stress_Death        = c("BAX","BCL2","CASP3","TP53","BBC3"),
  Protein_Synthesis   = c("EIF2S1","EIF4E","RPS6","RPLP0","RPL10A"),
  Inflammation        = c("NFKB1","RELA","TNF","IL1B","CXCL10")
  #Stemness            = c("SOX9","EPCAM","PROM1","KLF4","LGR5"),
  #Lipid_Metabolism    = c("SREBF1","ACACA","FASN","CPT1A","PPARG")
)

# --- 2) Filter to genes present in object ---
filter_present <- function(gvec) intersect(unique(gvec), rownames(beta_cells))
gene_categories <- lapply(all_genes_raw, filter_present)
dropped <- setdiff(unlist(all_genes_raw), unlist(gene_categories))
if (length(dropped)) message("Warning: dropped (not found): ", paste(dropped, collapse = ", "))

genes_to_plot <- unlist(gene_categories, use.names = FALSE)

# --- 3) Sex-cleaning & grouping ---
beta_cells$Sex <- as.character(beta_cells$Sex)
beta_cells$Sex <- ifelse(beta_cells$Sex %in% c("Male","M"), "M",
                         ifelse(beta_cells$Sex %in% c("Female","F"), "F", beta_cells$Sex))

beta_cells$group_id_base <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, sep = "_")
beta_cells$group_id_sex  <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, beta_cells$Sex, sep = "_")

# --- 4) Build DotPlot ---
p_raw <- DotPlot(
  beta_cells,
  features  = genes_to_plot,
  group.by  = "group_id_sex",
  dot.scale = 6,
  scale     = TRUE,
  col.min   = -1,
  col.max   =  1
)
dotdata <- p_raw$data

# --- 5) Restore group parts + annotate categories ---
dotdata <- dotdata %>%
  tidyr::separate(id, into = c("beta_cluster","diabetes_status","Sex"), sep = "_", remove = FALSE) %>%
  mutate(
    group_base = paste(beta_cluster, diabetes_status, sep = "_"),
    Sex = factor(Sex, levels = c("M","F"), labels = c("Male","Female"))
  )

gene_category_map <- stack(gene_categories)
colnames(gene_category_map) <- c("gene","gene_category")
dotdata <- left_join(dotdata, gene_category_map, by = c("features.plot" = "gene")) %>%
  mutate(gene_category = factor(gene_category, levels = names(gene_categories)))

# --- 6) Keep β2–β4 (like your previous script) ---
keep_groups <- c("β2_ND","β2_T2D","β3_ND","β3_T2D","β4_ND","β4_T2D")
dotdata <- dotdata %>% filter(group_base %in% keep_groups)

dotdata$group_base <- factor(dotdata$group_base,
                             levels = c("β2_ND","β2_T2D","β3_ND","β3_T2D","β4_ND","β4_T2D")
)
dotdata$features.plot <- factor(dotdata$features.plot, levels = genes_to_plot)

# --- 7) Plot ---
fig_modules_all <- ggplot(dotdata, aes(x = features.plot, y = group_base, size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradient2(low = "white", high = "red", midpoint = -1) +
  facet_grid(. ~ gene_category + Sex, scales = "free_x", space = "free_x") +
  labs(
    title = "All Pathway Gene Sets (β2–β4): Sex-aware DotPlot",
    x = "Gene", y = "β Subtype × Diabetes",
    fill = "Avg Expr (scaled)", size = "% Expressed"
  ) +
  theme_light() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9, face = "bold", colour = "black"),
    axis.text.y   = element_text(size = 10, face = "bold", colour = "black"),
    strip.text.x  = element_text(size = 11, face = "bold"),
    plot.title    = element_text(size = 14, face = "bold"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(fig_modules_all)


# Load required package
library(Seurat)

# --- 1. Define gene lists per pathway ---
gene_sets <- list(
  # Oxidative phosphorylation / ETC (Complex I–V; concise, beta-cell friendly core)
  OXPHOS = c(
    # Complex I (NDUFs)
    "NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFA9","NDUFA10",
    # Complex II
    "SDHA","SDHB","SDHC","SDHD",
    # Complex III
    "UQCRC1","UQCRC2","UQCRQ","CYC1",
    # Cytochrome c
    "CYCS",
    # Complex IV
    "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX7A2","COX7B","COX8A",
    # Complex V (ATP synthase)
    "ATP5F1A","ATP5F1B","ATP5F1C","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5PF","ATP5PB"
  ),
  
  # Antioxidant / ROS detox
  Antioxidant = c(
    "SOD2","SOD1","CAT","GPX1","GPX4","GSR","TXN","TXNRD1","PRDX1","PRDX2","PRDX3","PRDX5",
    "NQO1","HMOX1","GCLC","GCLM","KEAP1","NFE2L2"
  ),
  
  # Insulin signaling (receptor → PI3K → AKT → mTOR / FOXO; β-cell relevant)
  Insulin_Signaling = c(
    "INS","MAFA",
    "PCSK1","SYN1",   # regulatory + catalytic most stable in β-cells
    "VAMP2","ABCC8",
    "GCK","GLP1R",
    "KCNJ11", "KCNJ3",
    "SLC2A2", "CHGA",
    "STX1A", "TTR"
  ),
  
  # Glucose homeostasis / stimulus–secretion coupling (β-cell core)
  Glucose_Homeostasis = c(
    "INS","INS-IGF2","PDX1","MAFA","NKX6-1","ISL1","GCK","PCSK1","PCSK2","SLC2A2",
    "KCNJ11","ABCC8","FXYD2","SLC30A8","GLUD1","DLK1","RFX6","CHGA","CHGB"
  ),
  
  UPR = c(
    "HSPA5", "XBP1", "ATF4", "ATF6", "DDIT3", "EIF2AK3", "HERPUD1", "DNAJB9", "HYOU1", "PDIA6", 
    "HSP90B1", "SELENOS", "CLGN", "UGGT2", "CANX", "ERLEC1", "CCDC47", "TMEM259", "HSP90AB1", 
    "DNAJC3", "HSP90AA1", "TMBIM6", "CREBRF", "HSPA4", "MANF", "UFL1", "PPP1R15A", "PDIA3", 
    "PDIA4", "PDIA5", "EDEM1", "SYVN1", "DERL1", "ASNS", "TRIB3"
  ),
  
  Stress_Death = c(
    "BAX", "BCL2", "CASP3", "CASP7", "CASP9", "BBC3", "GADD45A", "FOS", "JUN", "TP53", 
    "BAD", "BID", "MCL1", "PMAIP1", "CYCS", "DAP", "DUSP1", "TXNIP"
  ),
  
  Inflammation = c(
    "NFKB1", "NFKBIA", "TNF", "IL6", "IL1B", "CXCL8", "CCL2", "ICAM1", "PTGS2", "IRF1",
    "STAT1", "STAT3", "IRAK1", "MYD88", "RELA", "TLR4", "TLR2", "IFIT1", "IFIT3", "CXCL10"
  ),
  
  Mitophagy = c(
    "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "SQSTM1", "MAP1LC3B", "OPTN", "TOMM20", 
    "VDAC1", "MFN2", "DNM1L", "RHOT1", "NIPSNAP1", "NIPSNAP2", "TBC1D15", "ATG5", "ATG7", "ULK1"
  ),
  
  # Protein synthesis / translation (initiation & elongation; representative RPs)
  Protein_Synthesis = c(
    # Initiation & control
    "EIF2S1","EIF2B1","EIF2B2","EIF2B3","EIF2B4","EIF2B5",
    "EIF4E","EIF4A1","EIF4A2","EIF4G1","EIF4G2","EIF4EBP1","EIF4EBP2",
    "EIF3A","EIF3B","EIF3C","EIF3D","EIF3E","EIF3F","EIF3G","EIF3H","EIF3I","EIF3J","EIF3K","EIF3L","EIF3M",
    "PABPC1","EEF2","EEF1A1","EEF1A2","EEF1B2","EEF1D","EEF1G",
    # A representative, compact panel of ribosomal proteins (avoid listing all ~80)
    "RPS3","RPS6","RPS8","RPS10","RPS12","RPS14","RPS18","RPS19","RPS20","RPS24","RPS27A",
    "RPL3","RPL7","RPL10","RPL10A","RPL11","RPL13","RPL13A","RPL18","RPLP0","RPLP1","RPLP2"
  ),
  
  # Stemness / ductal-progenitor & Notch/YAP signaling
  Stemness = c(
    "SOX9","EPCAM","KRT19","KRT8","KRT18","MUC1","PROM1","CD44","ITGA6","ITGB1",
    "ALDH1A1","KLF4","LGR5","MKI67",
    "HNF1B","MMP7","KRT7","KRT23",
    # Notch / progenitor
    "NOTCH1","NOTCH2","NOTCH3","JAG1","JAG2","DLL1","DLL4","HES1","HEY1",
    # Hippo/YAP axis
    "YAP1","WWTR1","TEAD1","TEAD4","CTGF"
  ),
  
  # Lipid metabolism (uptake, synthesis, storage, oxidation, cholesterol)
  Lipid_Metabolism = c(
    # Uptake & trafficking
    "CD36","FABP4","FABP5","LDLR","SCARB1","LPL","SLC27A1","SLC27A4",
    # Lipogenesis / desaturation
    "SREBF1","MLXIPL","ACACA","ACACB","FASN","SCD","ELOVL6","FADS1","FADS2","GPAM","DGAT1","DGAT2","AGPAT2",
    # Storage / lipolysis
    "PLIN2","PLIN5","LIPE","PNPLA2","MGLL","ATGL","CIDEC",
    # β-oxidation (mitochondrial & peroxisomal)
    "CPT1A","CPT1B","CPT2","ACADM","ACADL","ACADVL","HADHA","HADHB","ECHS1","ACOX1",
    # Regulators
    "PPARG","PPARA","PPARD","PPARGC1A",
    # Cholesterol synthesis & efflux
    "HMGCR","HMGCS1","MVK","PMVK","MVD","FDPS","FDFT1","SQLE","LSS","DHCR7","DHCR24","ABCA1","ABCG1",
    # Lipoproteins
    "APOE","APOA1"
  )
)


# --- 2. Calculate Module Scores ---
# This adds new metadata columns: UPR1, Stress_Death1, etc.
for (name in names(gene_sets)) {
  beta_cells <- AddModuleScore(
    object = beta_cells,
    features = list(gene_sets[[name]]),
    name = name
  )
}

# --- 3. Confirm new columns in metadata ---
colnames(beta_cells@meta.data)

# List of module score features
features <- c(
  #"UPR1",
  #"Stress_Death1",
  "Inflammation1",
  #"Mitophagy1",
  #"OXPHOS1",
  #"Antioxidant1",
  #"Insulin_Signaling1",
  #"Glucose_Homeostasis1",
  "Protein_Synthesis1",
  #"Stemness1",
  "Lipid_Metabolism1"
)

# Prepare data
plot_df <- beta_cells@meta.data %>%
  select(all_of(features), diabetes_status, Sex) %>%
  pivot_longer(
    cols = all_of(features),
    names_to = "Pathway",
    values_to = "Score"
  )

# Violin + boxplot + stat test (Wilcoxon)
library(ggpubr)  # for stat_compare_means
ggplot(plot_df, aes(x = diabetes_status, y = Score, fill = diabetes_status)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.4) +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  facet_grid(Pathway ~ Sex, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4) +
  scale_fill_manual(values = c("ND" = "#1F77B4", "T2D" = "#D62728")) +
  labs(x = "", y = "Module Score", fill = "Diabetes Status") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# Reshape metadata to long format
plot_df <- beta_cells@meta.data %>%
  select(all_of(features), diabetes_status, Sex) %>%
  pivot_longer(
    cols = all_of(features),
    names_to = "Pathway",
    values_to = "Score"
  )

# Run Wilcoxon tests: compare ND vs T2D for each Pathway within each Sex
pval_df <- plot_df %>%
  group_by(Pathway, Sex) %>%
  summarise(
    p = wilcox.test(Score ~ diabetes_status)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "bonferroni"),
    log10FDR = -log10(p.adj),
    p.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# View result
print(pval_df)

####################
# Figure 4
####################

## βT2D vs β1 ND DOWN
## Load prereqs
make_sankey_plot <- function(my_comparison_colors, my_pathway_colors, my_ora_dir, my_pathways) {
  library(readr)
  library(dplyr)
  library(ggalluvial)
  library(ggplot2)
  library(cowplot)
  library(scales)
  
  nukefix <- function(x) {
    x <- gsub("β", "b", x)
    x <- trimws(x)
    x <- iconv(x, to = "ASCII//TRANSLIT")
    x <- tolower(x)
    gsub("\\s+", " ", x)
  }
  
  # --- Get file names ---
  my_comp_names <- list.files(my_ora_dir, pattern = "_results\\.csv$", full.names = FALSE)
  
  all_links <- list()
  for (fname in my_comp_names) {
    comp_raw <- gsub("_results.csv", "", fname)
    comp <- nukefix(comp_raw)
    fpath <- file.path(my_ora_dir, fname)
    if (!file.exists(fpath)) next
    df <- read_csv(fpath, show_col_types = FALSE)
    if (!"term_name" %in% colnames(df)) next
    df <- df %>% filter(!is.na(`hypergeometric FDR`), `hypergeometric FDR` < 0.05, !is.na(term_name))
    df$Pathway <- nukefix(df$term_name)
    for (p in my_pathways) {
      pw <- nukefix(p)
      match_row <- df[df$Pathway == pw,]
      if (nrow(match_row) == 0) next
      all_links[[length(all_links)+1]] <- data.frame(
        Comparison = comp,
        Pathway = pw,
        Value = -log10(match_row$`hypergeometric FDR`[1]) + 1,
        stringsAsFactors = FALSE
      )
    }
  }
  
  sankey_df <- dplyr::bind_rows(all_links)
  sankey_df$Comparison <- factor(sankey_df$Comparison, levels = names(my_comparison_colors))
  sankey_df$Pathway <- factor(sankey_df$Pathway, levels = nukefix(my_pathways))
  
  # --- Main Sankey plot ---
  p_main <- ggplot(sankey_df, aes(axis1 = Comparison, axis2 = Pathway, y = Value)) +
    geom_alluvium(aes(fill = Comparison), width = 1/12, alpha = 0.8) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_discrete(limits = c("Comparison", "Pathway"), expand = c(.1, .1)) +
    scale_fill_manual(
      values = c(my_comparison_colors, my_pathway_colors),
      breaks = names(my_comparison_colors)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold")
    )
  
  # --- Dot legend ---
  fdr_values <- sankey_df$Value
  legend_vals <- pretty(fdr_values, n = 5)
  legend_vals <- sort(unique(legend_vals[legend_vals > 0]), decreasing = TRUE)
  
  rescale_size <- function(x, to = c(3, 10)) {
    rng <- range(x, na.rm = TRUE)
    scales::rescale(x, to = to, from = rng)
  }
  
  dot_sizes <- rescale_size(legend_vals, to = c(3, 10))
  
  legend_df <- data.frame(
    x = 1,
    y = seq(50, 38, length.out = length(dot_sizes)),
    y_label = seq(50, 38, length.out = length(dot_sizes)),
    size = dot_sizes,
    label = paste0("-log10(FDR)+1 = ", legend_vals)
  )
  
  p_legend <- ggplot(legend_df) +
    geom_point(aes(x = x, y = y, size = size), shape = 21, fill = "steelblue", color = "black", stroke = 0.25) +
    geom_text(aes(x = x + 0.4, y = y_label, label = label), hjust = 0, vjust = 0.5, size = 4) +
    theme_void() +
    coord_cartesian(clip = "off") +
    scale_size_identity() +
    scale_x_continuous(limits = c(0.9, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(plot.margin = margin(5, 20, 5, 5))
  
  # --- Combine ---
  cowplot::plot_grid(
    p_main,
    p_legend,
    rel_widths = c(4.2, 1.2),
    nrow = 1,
    axis = "none",
    align = "none"
  )
}

# --- Colors ---
my_comparison_colors <- c(
  "b2_t2d_vs_b1_m_nd" = "#274472",
  "b3_t2d_vs_b1_m_nd" = "#406E8E",
  "b4_t2d_vs_b1_m_nd" = "#6998B7",
  "b2_t2d_vs_b1_f_nd" = "#C97D60",
  "b3_t2d_vs_b1_f_nd" = "#B85C5C",
  "b4_t2d_vs_b1_f_nd" = "#8C3333"
)

my_pathway_colors <- c(
  # Protein synthesis & ribosome (matte red)
  "translation"                        = "#C0392B",  # Matte muted red
  "rrna processing"                    = "#C0392B",
  
  # Oxphos & mitochondrial (matte orange)
  "electron transport chain"           = "#D35400",  # Matte orange
  "oxidative phosphorylation"          = "#D35400",
  "atp synthesis coupled electron transport" = "#D35400",
  "proton transmembrane transport"     = "#D35400",
  
  # Antioxidant/ROS (matte green)
  "cellular oxidant detoxification"    = "#27AE60",  # Muted green (same as before, but fits matte)
  "antioxidant activity"               = "#27AE60",
  
  # Insulin/hormone signaling (matte blue)
  "response to insulin"                = "#2980B9",  # Muted blue (less saturated than "#1976D2")
  "insulin receptor signaling pathway" = "#2980B9",
  "hormone secretion"                  = "#2980B9",
  "insulin processing"                 = "#2980B9",
  "insulin secretion"                  = "#2980B9",
  
  # Lipid metabolism (matte purple)
  "lipid metabolic process"            = "#7D3C98",  # Matte muted purple
  "response to fatty acid"             = "#7D3C98",
  
  # Glucose/glycemic (matte brown)
  "glucose homeostasis"                = "#8D6748",  # Matte brown
  "response to glucagon"               = "#8D6748",
  "glucose metabolic process"          = "#8D6748"
)

nukefix <- function(x) {
  x <- gsub("β", "b", x)
  x <- trimws(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("\\s+", " ", x)
}

# --- Input ---
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/T2D_vs_beta1/DOWN"
my_comp_names <- c(
  "β2_T2D_vs_β1_M_ND_results.csv",
  "β3_T2D_vs_β1_M_ND_results.csv",
  "β4_T2D_vs_β1_M_ND_results.csv",
  "β2_T2D_vs_β1_F_ND_results.csv",
  "β3_T2D_vs_β1_F_ND_results.csv",
  "β4_T2D_vs_β1_F_ND_results.csv"
)

my_pathways <- c(
  "translation",
  "rRNA processing",
  "electron transport chain",
  "oxidative phosphorylation",
  "ATP synthesis coupled electron transport",
  "proton transmembrane transport",
  "cellular oxidant detoxification",
  "antioxidant activity",
  "response to insulin",
  "insulin receptor signaling pathway",
  "hormone secretion",
  "insulin processing",
  "insulin secretion",
  "lipid metabolic process",
  "response to fatty acid",
  "glucose homeostasis",
  "response to glucagon",
  "glucose metabolic process"
)

# RUN plot
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/T2D_vs_beta1/DOWN",
  my_pathways = my_pathways
)


## βT2D vs β1 ND UP
# --- Colors ---
my_comparison_colors <- c(
  "b2_t2d_vs_b1_m_nd" = "#274472",
  "b3_t2d_vs_b1_m_nd" = "#406E8E",
  "b4_t2d_vs_b1_m_nd" = "#6998B7",
  "b2_t2d_vs_b1_f_nd" = "#C97D60",
  "b3_t2d_vs_b1_f_nd" = "#B85C5C",
  "b4_t2d_vs_b1_f_nd" = "#8C3333"
)

my_pathway_colors <- c(
  # Protein synthesis & ribosome (matte red)
  "translation"                        = "#C0392B",  # Matte muted red
  "protein metabolic process"                    = "#C0392B",
  
  # Mitophagy (matte brown)
  "macroautophagy"                = "#8D6748",
  "autophagy of mitochondrion"   = "#8D6748",
  "autophagy"                    = "#8D6748",
  
  # Oxphos & mitochondrial (matte orange)
  "atp biosynthetic process"              = "#D35400",
  "oxidative phosphorylation"             = "#D35400",
  "proton transmembrane transport"        = "#D35400",
  "aerobic electron transport chain"      = "#D35400",
  
  # UPR (matte green)
  "erad pathway"                         = "#27AE60",
  "perk-mediated unfolded protein response" = "#27AE60",
  "protein folding"                      = "#27AE60",
  "response to unfolded protein"         = "#27AE60",
  "response to endoplasmic reticulum stress" = "#27AE60",
  "response to topologically incorrect protein" = "#27AE60",
  "er overload response"                = "#27AE60",
  "positive regulation of proteolysis"  = "#27AE60",
  
  # Stress and Cell Death (matte blue)
  "intrinsic apoptotic signaling pathway" = "#2980B9",
  "response to oxidative stress"          = "#2980B9",
  "cell death"                            = "#2980B9",
  "apoptotic process"                     = "#2980B9",
  "cellular response to stress"           = "#2980B9",
  "response to stress"                    = "#2980B9",
  
  # Inflammation (matte purple)
  "p38MAPK cascade"                                 = "#7D3C98",
  "cellular response to interleukin-7"              = "#7D3C98",
  "response to virus"                               = "#7D3C98",
  "cellular response to interleukin-4"              = "#7D3C98",
  "interleukin-11-mediated signaling pathway"       = "#7D3C98",
  "positive regulation of NF-kappaB transcription factor activity" = "#7D3C98",
  "inflammatory response"                           = "#7D3C98",
  "tumor necrosis factor production"                = "#7D3C98",
  
  # Dedifferentiation / Stem-like Reprogramming (matte teal)
  "epithelial to mesenchymal transition"            = "#148F77",  # Matte teal
  "positive regulation of cell fate commitment"     = "#148F77",
  "neurogenesis"                                    = "#148F77",
  "chromatin remodeling"                            = "#148F77"
)

nukefix <- function(x) {
  x <- gsub("β", "b", x)
  x <- trimws(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("\\s+", " ", x)
}

# --- Input ---
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/T2D_vs_beta1/UP"
my_comp_names <- c(
  "β2_T2D_vs_β1_M_ND_results.csv",
  "β3_T2D_vs_β1_M_ND_results.csv",
  "β4_T2D_vs_β1_M_ND_results.csv",
  "β2_T2D_vs_β1_F_ND_results.csv",
  "β3_T2D_vs_β1_F_ND_results.csv",
  "β4_T2D_vs_β1_F_ND_results.csv"
)

my_pathways <- c(
  # Protein synthesis & ribosome
  "translation",
  "protein metabolic process",
  
  # Mitophagy
  "macroautophagy",
  "autophagy of mitochondrion",
  "autophagy",
  
  # Oxphos & mitochondrial
  "atp biosynthetic process",
  "oxidative phosphorylation",
  "proton transmembrane transport",
  "aerobic electron transport chain",
  
  # UPR
  "erad pathway",
  "perk-mediated unfolded protein response",
  "protein folding",
  "response to unfolded protein",
  "response to endoplasmic reticulum stress",
  "response to topologically incorrect protein",
  "er overload response",
  "positive regulation of proteolysis",
  
  # Stress and Cell Death
  "intrinsic apoptotic signaling pathway",
  "response to oxidative stress",
  "cell death",
  "apoptotic process",
  "cellular response to stress",
  "response to stress",
  
  # Inflammation
  "p38MAPK cascade",
  "cellular response to interleukin-7",
  "response to virus",
  "cellular response to interleukin-4",
  "interleukin-11-mediated signaling pathway",
  "positive regulation of NF-kappaB transcription factor activity",
  "inflammatory response",
  "tumor necrosis factor production",
  
  # Dedifferentiation / Stem-like
  "epithelial to mesenchymal transition",
  "positive regulation of cell fate commitment",
  "neurogenesis",
  "chromatin remodeling"
)

# RUN plot
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/T2D_vs_beta1/UP",
  my_pathways = my_pathways
)

# Module Scores
# Load required package
library(Seurat)
beta_cells <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)")

# --- 1. Define gene lists per pathway ---
gene_sets <- list(
  UPR = c(
    "HSPA5", "XBP1", "ATF4", "ATF6", "DDIT3", "EIF2AK3", "HERPUD1", "DNAJB9", "HYOU1", "PDIA6", 
    "HSP90B1", "SELENOS", "CLGN", "UGGT2", "CANX", "ERLEC1", "CCDC47", "TMEM259", "HSP90AB1", 
    "DNAJC3", "HSP90AA1", "TMBIM6", "CREBRF", "HSPA4", "MANF", "UFL1", "PPP1R15A", "PDIA3", 
    "PDIA4", "PDIA5", "EDEM1", "SYVN1", "DERL1", "ASNS", "TRIB3"
  ),
  
  Stress_Death = c(
    "BAX", "BCL2", "CASP3", "CASP7", "CASP9", "BBC3", "GADD45A", "FOS", "JUN", "TP53", 
    "BAD", "BID", "MCL1", "PMAIP1", "CYCS", "DAP", "DUSP1", "TXNIP"
  ),
  
  Inflammation = c(
    "NFKB1", "NFKBIA", "TNF", "IL6", "IL1B", "CXCL8", "CCL2", "ICAM1", "PTGS2", "IRF1",
    "STAT1", "STAT3", "IRAK1", "MYD88", "RELA", "TLR4", "TLR2", "IFIT1", "IFIT3", "CXCL10"
  ),
  
  Mitophagy = c(
    "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "SQSTM1", "MAP1LC3B", "OPTN", "TOMM20", 
    "VDAC1", "MFN2", "DNM1L", "RHOT1", "NIPSNAP1", "NIPSNAP2", "TBC1D15", "ATG5", "ATG7", "ULK1"
  ),
  
  OXPHOS = c( # Complex I–V (structural core)
    "NDUFA1","NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
    "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11",
    "NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFV3",
    "SDHA","SDHB","SDHC","SDHD",
    "UQCRC1","UQCRC2","UQCRB","UQCRQ","UQCRH","UQCRFS1","CYC1",
    "CYCS",
    "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7B","COX8A",
    "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5PB","ATP5PD","ATP5PF","ATP5ME","ATP5MC1","ATP5MC2","ATP5MC3","ATP5MD"
  ),
  
  Antioxidant = c(
    "SOD1","SOD2","CAT","GPX1","GPX4","GSR","TXN","TXNRD1","PRDX1","PRDX2","PRDX3","PRDX5",
    "NQO1","HMOX1","GCLC","GCLM","KEAP1","NFE2L2","G6PD","PGD","TKT","TALDO1","IDH1","IDH2","ME1","ME2","SLC7A11"
  ),
  
  Protein_Synthesis = c(
    "EIF2S1","EIF2B1","EIF2B2","EIF2B3","EIF2B4","EIF2B5",
    "EIF4E","EIF4A1","EIF4G1","EIF4EBP1","EIF3A","EIF3B","EIF3C","EIF3D","EIF3E","EIF3F","EIF3G","EIF3H","EIF3I","EIF3J","EIF3K","EIF3L","EIF3M",
    "PABPC1","EEF2","EEF1A1","EEF1B2",
    "RPS3","RPS6","RPS10","RPS18","RPS19","RPS24","RPS27A",
    "RPL3","RPL7","RPL10","RPL10A","RPL11","RPL13","RPL13A","RPL18","RPLP0","RPLP1","RPLP2"
  ),
  
  Lipid_Metabolism = c(
    "CD36","FABP4","FABP5","LDLR","SCARB1","LPL","SLC27A1","SLC27A4",
    "SREBF1","ACACA","FASN","SCD","ELOVL6","FADS1","FADS2","GPAM","DGAT1","DGAT2","AGPAT2",
    "PLIN2","PLIN5","PNPLA2","LIPE","MGLL",
    "CPT1A","CPT2","ACADM","ACADL","ACADVL","HADHA","HADHB","ECHS1","ACOX1",
    "PPARG","PPARA","PPARD","PPARGC1A",
    "HMGCR","HMGCS1","MVK","PMVK","MVD","FDPS","FDFT1","SQLE","LSS","DHCR7","DHCR24","ABCA1","ABCG1",
    "APOE","APOA1"
  ),
  
  Glucose_Homeostasis = c(
    "PDX1","NKX6-1","NEUROD1","ISL1","MAFA","RFX6",
    "INS","IAPP","PCSK1","PCSK2","CPE","CHGA","CHGB","SLC30A8",
    "SLC2A2","GCK",
    "KCNJ11","ABCC8","KCNQ1","KCNK16",
    "CACNA1C","CACNA1D","CACNB2","CACNA2D1","RYR2","ATP2A2","ATP2B2","SLC8A1",
    "SNAP25","STX1A","STXBP1","VAMP2","SYT7","RIMS2","UNC13A","DOC2B","PCLO","RAB3A","SYTL4"
  )
)


# --- 2. Calculate Module Scores ---
# This adds new metadata columns: UPR1, Stress_Death1, etc.
for (name in names(gene_sets)) {
  beta_cells <- AddModuleScore(
    object = beta_cells,
    features = list(gene_sets[[name]]),
    name = name
  )
}

# --- 3. Confirm new columns in metadata ---
colnames(beta_cells@meta.data)

# List of module score features
features <- c("UPR1", "Stress_Death1", "Inflammation1", "Mitophagy1", "OXPHOS1", "Antioxidant1", "Protein_Synthesis1", "Lipid_Metabolism1", "Glucose_Homeostasis1")

beta_cells@meta.data <- beta_cells@meta.data %>%
  mutate(
    plot_group = case_when(
      beta_cluster == "β1" & diabetes_status == "ND" ~ "ND_β1",
      diabetes_status == "T2D" ~ paste0("T2D_", beta_cluster),
      TRUE ~ NA_character_
    )
  )

# Optional: Set factor levels for nice plotting order
beta_cells@meta.data$plot_group <- factor(
  beta_cells@meta.data$plot_group,
  levels = c("ND_β1", "T2D_β2", "T2D_β3", "T2D_β4")
)

plot_df <- beta_cells@meta.data %>%
  filter(!is.na(plot_group)) %>%
  select(plot_group, Sex, UPR1, Stress_Death1, Inflammation1, Mitophagy1, OXPHOS1, Antioxidant1, Protein_Synthesis1, Lipid_Metabolism1, Glucose_Homeostasis1) %>%
  pivot_longer(cols = starts_with(c("UPR", "Stress_Death", "Inflammation", "Mitophagy", "OXPHOS", "Antioxidant", "Protein_Synthesis", "Lipid_Metabolism", "Glucose_Homeostasis")),
               names_to = "Module", values_to = "Score")

# Clean module names if needed
plot_df$Module <- gsub("1$", "", plot_df$Module)

ggplot(plot_df, aes(x = plot_group, y = Score, fill = Sex)) +
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, alpha = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black",
               position = position_dodge(width = 0.9), lwd = 0.3) +
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("M" = "#0073C2", "F" = "#E31963")) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top") +
  labs(x = NULL, y = "Module Score", title = "Key Pathway Module Scores by Sex")

library(dplyr)
library(tidyr)
library(purrr)

# 1. Gather required columns into long format
module_scores <- c("UPR1", "Stress_Death1", "Inflammation1", "Mitophagy1", "OXPHOS1", "Antioxidant1", "Protein_Synthesis1", "Lipid_Metabolism1", "Glucose_Homeostasis1")

plot_df <- beta_cells@meta.data %>%
  select(plot_group, Sex, all_of(module_scores)) %>%
  pivot_longer(
    cols = all_of(module_scores),
    names_to = "Pathway",
    values_to = "Score"
  ) %>%
  mutate(Pathway = gsub("1$", "", Pathway))


# 2. Define the T2D groups to compare to ND_β1
t2d_groups <- c("T2D_β1", "T2D_β2", "T2D_β3", "T2D_β4")

# 3. Wilcoxon test for each Pathway, Sex, T2D group vs ND_β1 (same sex)
pval_df <- map_dfr(
  t2d_groups,
  function(grp) {
    plot_df %>%
      filter(plot_group %in% c("ND_β1", grp)) %>%
      group_by(Pathway, Sex) %>%
      summarise(
        group = grp,
        p = tryCatch(
          wilcox.test(Score ~ plot_group)$p.value,
          error = function(e) NA_real_
        ),
        mean_T2D = mean(Score[plot_group == grp], na.rm = TRUE),
        mean_ND = mean(Score[plot_group == "ND_β1"], na.rm = TRUE),
        median_T2D = median(Score[plot_group == grp], na.rm = TRUE),
        median_ND = median(Score[plot_group == "ND_β1"], na.rm = TRUE),
        diff_mean = mean_T2D - mean_ND,
        diff_median = median_T2D - median_ND,
        .groups = "drop"
      )
  }
)


# 4. Adjust p-values (FDR) for all comparisons
# 1. Compute and annotate log10FDR (Inf safe)
pval_df <- pval_df %>%
  group_by(Pathway, Sex) %>%
  mutate(
    p.adj = p.adjust(p, method = "bonferroni"),
    # Replace zeros with smallest possible value to avoid Inf
    p.adj.nonzero = ifelse(p.adj == 0, .Machine$double.xmin, p.adj),
    log10FDR = -log10(p.adj.nonzero),
    p.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()

# 2. Replace Inf with highest finite value
max_finite <- max(pval_df$log10FDR[is.finite(pval_df$log10FDR)], na.rm = TRUE)

pval_df <- pval_df %>%
  mutate(
    log10FDR = ifelse(is.infinite(log10FDR), max_finite, log10FDR)
  )

# Optional: For plot annotation, you can add a label column:
pval_df <- pval_df %>%
  mutate(
    log10FDR_label = ifelse(!is.finite(-log10(p.adj)), paste0(">", round(max_finite, 1)), as.character(round(log10FDR, 2)))
  )


# 5. View result
print(pval_df, n=32)

library(ggplot2)

# Clean up group column for X axis (remove T2D_ prefix if you want, or not)
pval_df$Cluster <- gsub("T2D_", "", pval_df$group)

# Set consistent color mapping
sex_colors <- c("F" = "#C2185B", "M" = "#1976D2")  # Use your preferred palette

# Option 1: Facet by Pathway (each module is a facet)
ggplot(
  pval_df %>% filter(!is.na(log10FDR)), 
  aes(x = Cluster, y = log10FDR, fill = Sex)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  facet_wrap(~ Pathway, nrow = 1) +
  scale_fill_manual(values = sex_colors) +
  labs(
    x = "Beta Cell Cluster (T2D vs ND_β1)", 
    y = expression(-log[10]~"(FDR)"),
    title = expression("Pathway Significance (" * -log[10]~"FDR) by Cluster and Sex")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    panel.grid.major.x = element_blank()
  ) +
  geom_text(
    aes(label = p.signif), 
    position = position_dodge(width = 0.7), 
    vjust = -0.2, 
    size = 5, 
    color = "black"
  )

# UMAP
beta_cells@meta.data

library(Seurat)
library(patchwork)

# Make a combined group column in metadata
beta_cells$Sex_Disease <- paste0(beta_cells$Sex, "_", beta_cells$diabetes_status)

# Optionally, make it a factor for plotting order:
beta_cells$Sex_Disease <- factor(
  beta_cells$Sex_Disease,
  levels = c("F_ND", "F_T2D", "M_ND", "M_T2D")
)

# List of your module score columns
module_scores <- c("UPR1", "Stress_Death1", "Inflammation1", "Mitophagy1", "OXPHOS1", "Antioxidant1", "Protein_Synthesis1", "Lipid_Metabolism1", "Glucose_Homeostasis1")

# Loop through each module and plot, split by Sex_Disease
plots <- lapply(module_scores, function(feat) {
  FeaturePlot(
    object = beta_cells,
    features = feat,
    reduction = "umap",
    split.by = "Sex_Disease",
    pt.size = 0.5,
    min.cutoff = "q10",
    max.cutoff = "q90"
  ) + labs(title = gsub("1$", "", feat))
})

wrap_plots(plots, ncol = 2)


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Module columns and pretty labels
module_scores <- c("UPR1", "Stress_Death1", "Inflammation1", "Mitophagy1", "OXPHOS1", "Antioxidant1", "Protein_Synthesis1", "Lipid_Metabolism1", "Glucose_Homeostasis1")
module_labels <- c(
  UPR1                 = "UPR",
  Stress_Death1        = "Stress / Death",
  Inflammation1        = "Inflammation",
  Mitophagy1           = "Mitophagy",
  OXPHOS1              = "OXPHOS (ETC)",
  Antioxidant1         = "Antioxidant",
  Protein_Synthesis1   = "Protein synthesis",
  Lipid_Metabolism1    = "Lipid metabolism",
  Glucose_Homeostasis1 = "Glucose homeostasis"
)
disease_score_col <- "algorithme_confidence"

# 1. Prepare long dataframe
df_long <- beta_cells@meta.data %>%
  filter(!is.na(Sex), !is.na(.data[[disease_score_col]])) %>%
  pivot_longer(
    cols = all_of(module_scores),
    names_to = "Module",
    values_to = "ModuleScore"
  ) %>%
  filter(!is.na(ModuleScore)) %>%
  mutate(Module = module_labels[Module])

# 2. Calculate correlation and (raw) p, then FDR for each Module × Sex
corr_tbl <- df_long %>%
  group_by(Module, Sex) %>%
  summarise(
    r = cor(.data[[disease_score_col]], ModuleScore, method = "pearson"),
    p = cor.test(.data[[disease_score_col]], ModuleScore, method = "pearson")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(p.adj = p.adjust(p, method = "fdr"),
         p.label = ifelse(p.adj < 2.2e-16, "<2.2e-16", signif(p.adj, 2)),
         cor.label = paste0("R = ", sprintf("%.2f", r), "\nFDR = ", p.label)
  )

# 3. Join back to df_long for plotting
df_long <- df_long %>%
  left_join(corr_tbl %>% select(Module, Sex, cor.label), by = c("Module", "Sex"))

# 4. Plot
p <- ggplot(df_long, aes_string(x = disease_score_col, y = "ModuleScore")) +
  geom_point(aes(color = Sex), alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", lwd = 0.5) +
  facet_grid(Module ~ Sex, scales = "free_y") +
  # Only one annotation per facet
  geom_text(
    data = distinct(df_long, Module, Sex, cor.label),
    aes(x = -Inf, y = Inf, label = cor.label),
    hjust = -0.03, vjust = 1.1, size = 4, fontface = "bold"
  ) +
  labs(
    x = "Algorithm.AI Disease Score",
    y = "Module Score",
    title = "Module Score vs Disease Score by Sex\n(R = Pearson, FDR-adjusted p)"
  ) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = c("F" = "#C51B7D", "M" = "#2B83BA")) +
  theme(strip.text = element_text(face = "bold", size = 11))

p


library(dplyr)
library(tidyr)
library(ggplot2)

# Set column names
donor_col <- "Library"
module_scores <- c("UPR1", "Stress_Death1", "Inflammation1", "Mitophagy1", "OXPHOS1", "Antioxidant1", "Protein_Synthesis1", "Lipid_Metabolism1", "Glucose_Homeostasis1")
module_labels <- c(
  UPR1                 = "UPR",
  Stress_Death1        = "Stress / Death",
  Inflammation1        = "Inflammation",
  Mitophagy1           = "Mitophagy",
  OXPHOS1              = "OXPHOS (ETC)",
  Antioxidant1         = "Antioxidant",
  Protein_Synthesis1   = "Protein synthesis",
  Lipid_Metabolism1    = "Lipid metabolism",
  Glucose_Homeostasis1 = "Glucose homeostasis"
)
disease_score_col <- "algorithme_confidence"
sex_col <- "Sex"
disease_col <- "diabetes_status"

# 1. Collapse to donor × sex × disease × module (mean values)
df_long <- beta_cells@meta.data %>%
  filter(
    !is.na(.data[[sex_col]]), 
    !is.na(.data[[disease_score_col]]), 
    !is.na(.data[[donor_col]]), 
    !is.na(.data[[disease_col]])
  ) %>%
  pivot_longer(
    cols = all_of(module_scores),
    names_to = "Module",
    values_to = "ModuleScore"
  ) %>%
  filter(!is.na(ModuleScore)) %>%
  mutate(Module = module_labels[Module]) %>%
  group_by(Module, Sex = .data[[sex_col]], Disease = .data[[disease_col]], Donor = .data[[donor_col]]) %>%
  summarise(
    mean_ModuleScore = mean(ModuleScore, na.rm = TRUE),
    mean_DiseaseScore = mean(.data[[disease_score_col]], na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# 2. Correlation and FDR per Module × Sex
corr_tbl <- df_long %>%
  group_by(Module, Sex) %>%
  summarise(
    r = cor(mean_DiseaseScore, mean_ModuleScore, method = "pearson"),
    p = cor.test(mean_DiseaseScore, mean_ModuleScore, method = "pearson")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    #pval
    p.label = ifelse(p < 2.2e-16, "<2.2e-16", signif(p, 2)),
    cor.label = paste0("R=", sprintf("%.2f", r), "\np=", p.label)
    #FDR
    # p.adj = p.adjust(p, method = "fdr"),
    # p.label = ifelse(p.adj < 2.2e-16, "<2.2e-16", signif(p.adj, 2)),
    # cor.label = paste0("R=", sprintf("%.2f", r), "\nFDR=", p.label)
  )

# 3. Calculate y.label per Module (top of range)
max_y_tbl <- df_long %>%
  group_by(Module) %>%
  summarise(y_max = max(mean_ModuleScore, na.rm = TRUE))

corr_tbl <- corr_tbl %>%
  left_join(max_y_tbl, by = "Module") %>%
  mutate(
    y.label = y_max + 0.05 * abs(y_max),
    x.label = -Inf
  )

# 4. Plot: color by Sex, shape by Disease
ggplot(df_long, aes(x = mean_DiseaseScore, y = mean_ModuleScore, color = Sex, shape = Disease)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Sex, color = Sex), linetype = "solid") +
  geom_text(
    data = corr_tbl,
    aes(x = x.label, y = y.label, label = cor.label, color = Sex),
    hjust = -0.1, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE
  ) +
  facet_wrap(~ Module, scales = "free_y", ncol = 3) +
  labs(
    x = "Mean Algorithm.AI Disease Score (per donor)",
    y = "Mean Module Score (per donor)",
    title = "Correlation of Module Scores with Disease Score by Donor\nColored by Sex, Shape by Disease",
    color = "Sex",
    shape = "Disease"
  ) +
  theme_minimal(base_size = 14)

# Select and arrange key columns
corr_tbl %>%
  select(Module, Sex, r, p, n) %>%
  mutate(
    R = sprintf("%.2f", r),
    pval = formatC(p, format = "e", digits = 2)
  ) %>%
  select(Module, Sex, R, pval, n) %>%
  arrange(Module, Sex) %>%
  print(n = Inf)


# Module scores for Figure 3
# --- Load required libraries ---
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# --- 1. Define gene sets for beta cell pathways ---
gene_sets <- list(
  Insulin_Secretion = c("INS", "PCSK1", "PCSK2", "CPE", "CHGA", "CHGB", "VAMP2", "SYT7", "SLC30A8", "G6PC2", "PAM", "SCG2", "SCG5", "ABCC8", "KCNJ11", "IDE", "UCN3", "PDX1", "MAFA"),
  Glucose_Sensing = c("GCK", "SLC2A2", "ABCC8", "KCNJ11", "GPD2", "PDHA1", "PC", "PKM", "GLUD1", "SDHA", "NDUFA5"),
  Beta_Cell_Development = c("PDX1", "NKX6-1", "MAFA", "ISL1", "NEUROD1", "MNX1", "HNF1A", "HNF4A", "RFX6", "FOXA2", "GLIS3", "SOX9", "PAX6", "ARX", "ONECUT1", "GATA6", "PAX4", "SIX2"),
  Beta_Cell_Survival = c("BCL2", "MCL1", "BAX", "CASP3", "CASP7", "TP53", "DUSP1", "TXNIP", "SOD2", "PRDX1", "HSPA5", "ATF4", "ATF6", "DDIT3", "HERPUD1", "PPP1R15A", "XBP1", "NR4A1", "MANF"),
  Vesicle_Exocytosis = c("VAMP2", "STX1A", "SYT7", "RAB3A", "RAB27A", "UNC13A", "CPLX1", "SNAP25", "DOC2B", "DNM1", "STXBP1", "UNC13B", "NSF"),
  Endocrine_Signaling = c("INS", "GCG", "SST", "GLP1R", "INSR", "GIPR", "FFAR1", "CXCR4", "ADRA2A", "VIPR1", "VIPR2", "SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", "P2RY1", "CALCR", "ADCY5", "PRKACB")
)

# --- 2. Add module scores to the Seurat object ---
for (name in names(gene_sets)) {
  beta_cells <- AddModuleScore(
    object = beta_cells,
    features = list(gene_sets[[name]]),
    name = name
  )
}

# --- 3. Create cluster label for plotting (both ND and T2D for stats) ---
beta_cells@meta.data <- beta_cells@meta.data %>%
  mutate(
    plot_group = case_when(
      diabetes_status == "ND" ~ paste0("ND_", beta_cluster),
      diabetes_status == "T2D" ~ paste0("T2D_", beta_cluster),
      TRUE ~ NA_character_
    )
  )

# --- 4. Set consistent order for cluster labels (edit if needed) ---
cluster_levels <- c("ND_β1", "ND_β2", "ND_β3", "ND_β4", "T2D_β1", "T2D_β2", "T2D_β3", "T2D_β4")
beta_cells@meta.data$plot_group <- factor(
  beta_cells@meta.data$plot_group,
  levels = cluster_levels
)

# --- 5. Gather module scores into long format ---
module_scores <- paste0(names(gene_sets), "1")
plot_df <- beta_cells@meta.data %>%
  filter(!is.na(plot_group)) %>%
  select(plot_group, Sex, all_of(module_scores)) %>%
  pivot_longer(
    cols = all_of(module_scores),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(Module = gsub("1$", "", Module))

# --- 6. Violin + boxplot: Only ND clusters (for visual QC etc) ---
ggplot(plot_df %>% filter(grepl("^ND_", plot_group)),
       aes(x = plot_group, y = Score, fill = Sex)) +
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, alpha = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black",
               position = position_dodge(width = 0.9), lwd = 0.3) +
  facet_wrap(~ Module, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("M" = "#0073C2", "F" = "#E31963")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top") +
  labs(x = NULL, y = "Module Score", title = "Beta Cell Functional Pathway Scores by Sex (ND Only)")

# --- 7. Wilcoxon test: Each ND cluster vs ND_β1 ---
nd_groups <- c("ND_β2", "ND_β3", "ND_β4")
pval_df <- map_dfr(
  nd_groups,
  function(grp) {
    plot_df %>%
      filter(plot_group %in% c("ND_β1", grp)) %>%
      group_by(Module, Sex) %>%
      summarise(
        group = grp,
        p = tryCatch(
          wilcox.test(Score ~ plot_group)$p.value,
          error = function(e) NA_real_
        ),
        mean_other = mean(Score[plot_group == grp], na.rm = TRUE),
        mean_ND = mean(Score[plot_group == "ND_β1"], na.rm = TRUE),
        diff_mean = mean_other - mean_ND,
        .groups = "drop"
      )
  }
)

# --- 8. FDR correction, -log10(FDR), annotation for plotting ---
pval_df <- pval_df %>%
  group_by(Module, Sex) %>%
  mutate(
    p.adj = p.adjust(p, method = "bonferroni"),
    p.adj.nonzero = ifelse(p.adj == 0, .Machine$double.xmin, p.adj),
    log10FDR = -log10(p.adj.nonzero),
    p.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()

max_finite <- max(pval_df$log10FDR[is.finite(pval_df$log10FDR)], na.rm = TRUE)
pval_df <- pval_df %>%
  mutate(log10FDR = ifelse(is.infinite(log10FDR), max_finite, log10FDR))

pval_df$Cluster <- gsub("T2D_", "", pval_df$group)

# --- 9. Barplot: FDR for each pathway, cluster, sex ---
sex_colors <- c("F" = "#C2185B", "M" = "#1976D2")
ggplot(
  pval_df %>% filter(!is.na(log10FDR)),
  aes(x = Cluster, y = log10FDR, fill = Sex)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  facet_wrap(~ Module, nrow = 2) +
  scale_fill_manual(values = sex_colors) +
  labs(
    x = "Beta Cell Cluster (T2D vs ND_β1)",
    y = expression(-log[10]~"(FDR)"),
    title = expression("Functional Pathway Significance (" * -log[10]~"FDR) by Cluster and Sex")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    panel.grid.major.x = element_blank()
  ) +
  geom_text(
    aes(label = p.signif),
    position = position_dodge(width = 0.7),
    vjust = -0.2,
    size = 5,
    color = "black"
  )

# 2. Replace Inf with highest finite value
max_finite <- max(pval_df$log10FDR[is.finite(pval_df$log10FDR)], na.rm = TRUE)

pval_df <- pval_df %>%
  mutate(
    log10FDR = ifelse(is.infinite(log10FDR), max_finite, log10FDR)
  )

# Optional: For plot annotation, you can add a label column:
pval_df <- pval_df %>%
  mutate(
    log10FDR_label = ifelse(!is.finite(-log10(p.adj)), paste0(">", round(max_finite, 1)), as.character(round(log10FDR, 2)))
  )


# 5. View result
print(pval_df, n=32)

# Cor plots
# --- Load libraries ---
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. Define your custom gene sets ---
gene_sets <- list(
  Insulin_Secretion = c("INS", "PCSK1", "PCSK2", "CPE", "CHGA", "CHGB", "VAMP2", "SYT7", "SLC30A8", "G6PC2", "PAM", "SCG2", "SCG5", "ABCC8", "KCNJ11", "IDE", "UCN3", "PDX1", "MAFA"),
  Glucose_Sensing = c("GCK", "SLC2A2", "ABCC8", "KCNJ11", "GPD2", "PDHA1", "PC", "PKM", "GLUD1", "SDHA", "NDUFA5"),
  Beta_Cell_Development = c("PDX1", "NKX6-1", "MAFA", "ISL1", "NEUROD1", "MNX1", "HNF1A", "HNF4A", "RFX6", "FOXA2", "GLIS3", "SOX9", "PAX6", "ARX", "ONECUT1", "GATA6", "PAX4", "SIX2"),
  Beta_Cell_Survival = c("BCL2", "MCL1", "BAX", "CASP3", "CASP7", "TP53", "DUSP1", "TXNIP", "SOD2", "PRDX1", "HSPA5", "ATF4", "ATF6", "DDIT3", "HERPUD1", "PPP1R15A", "XBP1", "NR4A1", "MANF"),
  Vesicle_Exocytosis = c("VAMP2", "STX1A", "SYT7", "RAB3A", "RAB27A", "UNC13A", "CPLX1", "SNAP25", "DOC2B", "DNM1", "STXBP1", "UNC13B", "NSF"),
  Endocrine_Signaling = c("INS", "GCG", "SST", "GLP1R", "INSR", "GIPR", "FFAR1", "CXCR4", "ADRA2A", "VIPR1", "VIPR2", "SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", "P2RY1", "CALCR", "ADCY5", "PRKACB")
)

# --- 2. Add module scores for each pathway ---
for (name in names(gene_sets)) {
  beta_cells <- AddModuleScore(
    object = beta_cells,
    features = list(gene_sets[[name]]),
    name = name
  )
}

# --- 3. Prepare module score and meta data ---
donor_col <- "Library"                # or your donor/sample column
sex_col <- "Sex"
disease_col <- "diabetes_status"
disease_score_col <- "algorithme_confidence"

# List module score columns (created by AddModuleScore)
module_scores <- paste0(names(gene_sets), "1")
module_labels <- setNames(names(gene_sets), module_scores) # Pretty labels

# --- 4. Collapse to mean per donor × sex × disease × module ---
df_long <- beta_cells@meta.data %>%
  filter(
    !is.na(.data[[sex_col]]), 
    !is.na(.data[[disease_score_col]]), 
    !is.na(.data[[donor_col]]), 
    !is.na(.data[[disease_col]])
  ) %>%
  pivot_longer(
    cols = all_of(module_scores),
    names_to = "Module",
    values_to = "ModuleScore"
  ) %>%
  filter(!is.na(ModuleScore)) %>%
  mutate(Module = module_labels[Module]) %>%
  group_by(Module, Sex = .data[[sex_col]], Disease = .data[[disease_col]], Donor = .data[[donor_col]]) %>%
  summarise(
    mean_ModuleScore = mean(ModuleScore, na.rm = TRUE),
    mean_DiseaseScore = mean(.data[[disease_score_col]], na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# --- 5. Correlation and p per module × sex ---
corr_tbl <- df_long %>%
  group_by(Module, Sex) %>%
  summarise(
    r = cor(mean_DiseaseScore, mean_ModuleScore, method = "pearson"),
    p = cor.test(mean_DiseaseScore, mean_ModuleScore, method = "pearson")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    p.label = ifelse(p < 2.2e-16, "<2.2e-16", signif(p, 2)),
    cor.label = paste0("R=", sprintf("%.2f", r), "\np=", p.label)
  )

# --- 6. Calculate label position for annotation ---
max_y_tbl <- df_long %>%
  group_by(Module) %>%
  summarise(y_max = max(mean_ModuleScore, na.rm = TRUE))

corr_tbl <- corr_tbl %>%
  left_join(max_y_tbl, by = "Module") %>%
  mutate(
    y.label = y_max + 0.05 * abs(y_max),
    x.label = -Inf
  )

# --- 7. Plot: color by Sex, shape by Disease (ND/T2D), facet by module ---
p <- ggplot(df_long, aes(x = mean_DiseaseScore, y = mean_ModuleScore, color = Sex, shape = Disease)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Sex, color = Sex), linetype = "solid") +
  geom_text(
    data = corr_tbl,
    aes(x = x.label, y = y.label, label = cor.label, color = Sex),
    hjust = -0.1, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE
  ) +
  facet_wrap(~ Module, scales = "free_y", nrow = 1) +   # <-- key change here!
  scale_color_manual(values = c("M" = "#1976D2", "F" = "#C2185B")) +
  labs(
    x = "Mean Algorithm.AI Disease Score (per donor)",
    y = "Mean Module Score (per donor)",
    title = "Correlation of Functional Module Scores with Disease Score by Donor\nColor: Sex, Shape: Disease",
    color = "Sex",
    shape = "Disease"
  ) +
  theme_minimal(base_size = 14)

print(p)


# --- 8. Print summary of R and p-values (legend/table ready) ---
corr_tbl %>%
  select(Module, Sex, r, p, n) %>%
  mutate(
    R = sprintf("%.2f", r),
    pval = formatC(p, format = "e", digits = 2)
  ) %>%
  select(Module, Sex, R, pval, n) %>%
  arrange(Module, Sex) %>%
  print(n = Inf)

# HBa1c plotting
# --- 1. Load libraries ---
library(dplyr)
library(readr)
library(tibble)

# --- 2. Load clinical HbA1c data ---
clinical_df <- read_csv("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Donor_Summary_186.csv") %>%
  mutate(Library = as.character(donor_ID)) %>%
  select(Library, hba1c)

# --- 3. Manually add missing donors ---
hba1c_missing <- tribble(
  ~Library,        ~hba1c,
  "HP2022801",     5.5,
  "HP2024001",     5.4,
  "SAMN15877725",  5.6,
  "HP2031401",     5.4,
  "HP2105501",     5.6,
  "HP2106201",     5.3,
  "HP2107001",     5.1,
  "HP2107901",     5.2,
  "HP2108601",     5.1,
  "HP2108901",     5.9,
  "HP2110001",     5.5,
  "HP2123201",     5.3,
  "HP2132801",     5.5,
  "HP2202101",     5.5,
  "HP2121601",     5.8
)

# --- 4. Combine full HbA1c data (ensure Library is character for both) ---
hba1c_full <- bind_rows(
  clinical_df %>% mutate(Library = as.character(Library)),
  hba1c_missing %>% mutate(Library = as.character(Library))
)

# --- 5. Add HbA1c to Seurat meta.data ---
meta <- beta_cells@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(Library = as.character(Library)) %>%
  left_join(hba1c_full, by = "Library") %>%
  tibble::column_to_rownames("cell_id")
beta_cells@meta.data <- meta

# --- 6. QC: How many cells have non-missing HbA1c? ---
cat("Non-missing HbA1c per cell:\n")
print(table(!is.na(beta_cells@meta.data$hba1c)))

# Donors in meta.data with missing HbA1c
beta_cells@meta.data %>%
  as.data.frame() %>%
  filter(is.na(hba1c)) %>%
  pull(Library) %>%
  unique()

cat("\nExample of donors with HbA1c:\n")
print(beta_cells@meta.data %>% filter(!is.na(hba1c)) %>% count(Library, hba1c) %>% head(10))

cat("\nExample of donors with NA HbA1c:\n")
print(beta_cells@meta.data %>% filter(is.na(hba1c)) %>% count(Library))

# --- 7. Check all unique Libraries that are NA after join ---
cat("\nDonors (Library) with missing HbA1c annotation:\n")
print(unique(beta_cells@meta.data$Library[is.na(beta_cells@meta.data$hba1c)]))


library(dplyr)
library(tidyr)
library(ggplot2)

# Define your columns
donor_col <- "Library"
hba1c_col <- "hba1c"          # <--- change if your column is named differently
sex_col <- "Sex"
# ---- Define module score column names and pretty labels ----

# This matches the output column names from AddModuleScore (always "name1" for each gene set)
module_scores <- c(
  "Insulin_Secretion1",
  "Glucose_Sensing1",
  "Beta_Cell_Development1",
  "Beta_Cell_Survival1",
  "Vesicle_Exocytosis1",
  "Endocrine_Signaling1"
)

# Human-readable labels for plotting, in the same order
module_labels <- c(
  Insulin_Secretion1      = "Insulin Secretion & Processing",
  Glucose_Sensing1        = "Glucose Sensing & Metabolism",
  Beta_Cell_Development1  = "Beta Cell Development & Differentiation",
  Beta_Cell_Survival1     = "Beta Cell Survival",
  Vesicle_Exocytosis1     = "Vesicle Trafficking & Exocytosis",
  Endocrine_Signaling1    = "Endocrine Signaling & Communication"
)

# Get mean HbA1c per donor (if needed)
hba1c_per_donor <- beta_cells@meta.data %>%
  filter(!is.na(.data[[donor_col]]), !is.na(.data[[hba1c_col]])) %>%
  group_by(Donor = .data[[donor_col]]) %>%
  summarise(mean_HbA1c = mean(.data[[hba1c_col]], na.rm = TRUE), .groups = "drop")

# Collapse module scores per donor × sex × module
df_long <- beta_cells@meta.data %>%
  filter(
    !is.na(Sex), 
    !is.na(hba1c), 
    !is.na(Library), 
    !is.na(diabetes_status)
  ) %>%
  pivot_longer(
    cols = paste0(names(gene_sets), "1"),
    names_to = "Module",
    values_to = "ModuleScore"
  ) %>%
  filter(!is.na(ModuleScore)) %>%
  mutate(Module = gsub("1$", "", Module)) %>%
  group_by(Module, Sex, Disease = diabetes_status, Donor = Library) %>%
  summarise(
    mean_ModuleScore = mean(ModuleScore, na.rm = TRUE),
    mean_HbA1c = mean(hba1c, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# --- Correlation and p-value ---
corr_tbl <- df_long %>%
  group_by(Module, Sex) %>%
  summarise(
    r = cor(mean_HbA1c, mean_ModuleScore, method = "pearson"),
    p = cor.test(mean_HbA1c, mean_ModuleScore, method = "pearson")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    p.label = ifelse(p < 2.2e-16, "<2.2e-16", signif(p, 2)),
    cor.label = paste0("R=", sprintf("%.2f", r), "\np=", p.label)
  )

# Calculate y.label for annotation
max_y_tbl <- df_long %>%
  group_by(Module) %>%
  summarise(y_max = max(mean_ModuleScore, na.rm = TRUE))

corr_tbl <- corr_tbl %>%
  left_join(max_y_tbl, by = "Module") %>%
  mutate(
    y.label = y_max + 0.05 * abs(y_max),
    x.label = -Inf
  )

# --- Plot ---
ggplot(df_long, aes(x = mean_HbA1c, y = mean_ModuleScore, color = Sex, shape = Disease)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Sex, color = Sex), linetype = "solid") +
  geom_text(
    data = corr_tbl,
    aes(x = x.label, y = y.label, label = cor.label, color = Sex),
    hjust = -0.1, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE
  ) +
  facet_wrap(~ Module, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("M" = "#1976D2", "F" = "#C2185B")) +
  scale_shape_manual(values = c("ND" = 16, "T2D" = 17)) +   # 16: circle, 17: triangle
  labs(
    x = "Mean HbA1c (per donor)",
    y = "Mean Module Score (per donor)",
    title = "Correlation of Beta Cell Module Scores with HbA1c by Donor",
    color = "Sex",
    shape = "Disease"
  ) +
  theme_minimal(base_size = 14)

# ================================
# Figure 5 — Mitophagy DotPlot (Sex-aware)
# ================================

# --- 0) LIBRARIES ---
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringi)
})

# --- 1) INPUT: Seurat object with required metadata ---
# Assumes `beta_cells` exists and has:
#   meta columns: beta_cluster (β1..β4), diabetes_status (ND/T2D), Sex (M/F)
#   RNA assay with gene symbols as rownames
stopifnot(exists("beta_cells"))
stopifnot(all(c("beta_cluster","diabetes_status","Sex") %in% colnames(beta_cells@meta.data)))

# --- 2) DEFINE MITOPHAGY GENE SET (with alias cleaning) ---
mitophagy_genes_raw <- c(
  "PINK1","PRKN","BNIP3","BNIP3L","NIX","FUNDC1",
  "SQSTM1","OPTN","TOMM20","MFN2","ULK1","TBK1","MAP1LC3B"
)

# Alias map: NIX -> BNIP3L (official)
alias_map <- c("NIX" = "BNIP3L")
mitophagy_genes <- unname(ifelse(mitophagy_genes_raw %in% names(alias_map),
                                 alias_map[mitophagy_genes_raw],
                                 mitophagy_genes_raw))
mitophagy_genes <- unique(mitophagy_genes)

# Filter to genes present in the object (warn if missing)
present_genes <- intersect(mitophagy_genes, rownames(beta_cells))
missing_genes <- setdiff(mitophagy_genes, present_genes)
if (length(missing_genes)) {
  message("Warning: These genes were not found and will be dropped: ",
          paste(missing_genes, collapse = ", "))
}
stopifnot(length(present_genes) > 0)

# Optional: keep a single category for faceting clarity (can add more later)
gene_categories <- list("Mitophagy Pathway" = present_genes)

# --- 3) PREPARE GROUPS (β-cluster + disease + Sex) ---
# Make sure Sex values are clean ("M"/"F")
beta_cells$Sex <- as.character(beta_cells$Sex)
beta_cells$Sex <- ifelse(beta_cells$Sex %in% c("Male","M"), "M",
                         ifelse(beta_cells$Sex %in% c("Female","F"), "F", beta_cells$Sex))

# Base group (without sex)
beta_cells$group_id_base <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, sep = "_")

# Sex-aware group used for DotPlot aggregation
beta_cells$group_id_sex  <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, beta_cells$Sex, sep = "_")

# --- 4) SET ASSAY & BUILD DOTPLOT (Sex-aware aggregation) ---
DefaultAssay(beta_cells) <- "RNA"

p_raw <- DotPlot(
  beta_cells,
  features  = present_genes,
  group.by  = "group_id_sex",
  dot.scale = 6,
  scale     = TRUE,
  col.min   = -1,
  col.max   =  1
)

dotdata <- p_raw$data

# --- 5) RESTORE GROUP PARTS & ADD ANNOTATIONS ---
# Split the sex-aware ID back into parts; keep a base group for the y-axis
dotdata <- dotdata %>%
  tidyr::separate(id, into = c("beta_cluster","diabetes_status","Sex"), sep = "_", remove = FALSE) %>%
  mutate(
    group_base = paste(beta_cluster, diabetes_status, sep = "_"),
    Sex = factor(Sex, levels = c("M","F"), labels = c("Male","Female"))
  )

# Gene → category map
gene_category_map <- stack(gene_categories)
colnames(gene_category_map) <- c("gene", "gene_category")
dotdata <- left_join(dotdata, gene_category_map, by = c("features.plot" = "gene")) %>%
  mutate(gene_category = factor(gene_category, levels = names(gene_categories)))

# Desired Y-order: β4..β1 and T2D then ND
desired_order <- c("β4_T2D","β3_T2D","β2_T2D",
                   "β4_ND", "β3_ND", "β2_ND", "β1_ND")
dotdata$group_base <- factor(dotdata$group_base, levels = desired_order)

# Ensure gene order is stable/pretty (optional: keep as provided)
dotdata$features.plot <- factor(dotdata$features.plot, levels = present_genes)

# --- 6) PLOT ---
fig5 <- ggplot(dotdata, aes(x = features.plot, y = group_base, size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, color = "black") +
  # With scale=TRUE in DotPlot, the natural midpoint is ~0
  scale_fill_gradient2(low = "white", high = "red", midpoint = -1) +
  facet_grid(. ~ gene_category + Sex, scales = "free_x", space = "free_x") +
  labs(
    title = "Mitophagy Pathway Genes in Human β-cells",
    x = "Gene",
    y = "β-cell Subtype & Diabetes Status",
    fill = "Avg Expr (scaled)",
    size = "% Expressed"
  ) +
  theme_light() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 10, face = "bold", colour = "black"),
    axis.text.y   = element_text(size = 10, face = "bold", colour = "black"),
    strip.text.x  = element_text(size = 12, face = "bold"),
    plot.title    = element_text(size = 14, face = "bold"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# Print to screen
print(fig5)

#UMAP
# --- LIBRARIES ---
library(Seurat)
library(ggplot2)
library(dplyr)

# --- 1. Set gene(s) to plot and scale range ---
genes_to_plot <- c("PINK1")  # Add more if needed
min_expr <- 0    # Minimum value for color scale (adjust as needed)
max_expr <- 0.5    # Maximum value for color scale (adjust as needed)

# --- 2. Confirm gene(s) exist in dataset ---
genes_available <- genes_to_plot[genes_to_plot %in% rownames(beta_cells)]
if (length(genes_available) == 0) stop("None of the genes found in the dataset!")

# --- 3. Loop over genes to create a plot for each ---
plot_list <- lapply(genes_available, function(gene) {
  # Fetch UMAP and expression data for plotting
  plot_df <- FetchData(beta_cells, vars = c("UMAP_1", "UMAP_2", "Sex", "diabetes_status", gene))
  colnames(plot_df)[5] <- "expr"
  
  # Factor for nice labels
  plot_df$Sex <- factor(plot_df$Sex, levels = c("M", "F"), labels = c("Male", "Female"))
  plot_df$diabetes_status <- factor(plot_df$diabetes_status, levels = c("ND", "T2D"))
  plot_df$facet_group <- factor(
    paste(plot_df$Sex, plot_df$diabetes_status, sep = " – "),
    levels = c("Male – ND", "Male – T2D", "Female – ND", "Female – T2D")
  )
  
  # ggplot2 UMAP with expression, faceted by Sex & Disease
  ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
    geom_point(size = 1.2) +
    scale_color_gradient(
      low = "lightgrey", high = "red", 
      limits = c(min_expr, max_expr), oob = scales::squish, name = paste(gene, "Expr")
    ) +
    facet_wrap(~ facet_group, nrow = 2, ncol = 2) +
    labs(
      title = paste("UMAP: Expression of", gene, "by Sex and Disease"),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_light(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 16),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
})

# --- 4. Show the plot(s) ---
if (length(plot_list) == 1) {
  print(plot_list[[1]])
} else {
  # If you plot multiple genes, show all
  for (p in plot_list) print(p)
}

# COr plots
# --- LIBRARIES ---
library(Seurat)
library(ggplot2)
library(dplyr)

# --- SETTINGS: change these if your column names differ ---
donor_col <- "Library"              # donor/sample ID column
sex_col <- "Sex"
disease_col <- "diabetes_status"
pink1_col <- "PINK1"
score_col <- "algorithme_confidence"        # change to your score column name

# --- 1. Collapse per donor × sex × disease group ---
df <- beta_cells@meta.data %>%
  filter(
    !is.na(.data[[donor_col]]),
    !is.na(.data[[sex_col]]),
    !is.na(.data[[disease_col]]),
    !is.na(.data[[score_col]])
  ) %>%
  mutate(PINK1 = FetchData(beta_cells, vars = pink1_col)[, 1]) %>%  # get per-cell PINK1
  group_by(
    Donor = .data[[donor_col]],
    Sex = .data[[sex_col]],
    Disease = .data[[disease_col]]
  ) %>%
  summarise(
    mean_PINK1 = mean(PINK1, na.rm = TRUE),
    mean_Score = mean(.data[[score_col]], na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# --- 2. Correlation stats by Sex ---
corr_tbl <- df %>%
  group_by(Sex) %>%
  summarise(
    r = cor(mean_PINK1, mean_Score, method = "pearson"),
    p = cor.test(mean_PINK1, mean_Score, method = "pearson")$p.value,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    p.label = ifelse(p < 2.2e-16, "<2.2e-16", signif(p, 2)),
    cor.label = paste0("R=", sprintf("%.2f", r), "\np=", p.label)
  )

# --- 3. Calculate y.label for annotation ---
max_y_tbl <- df %>%
  group_by(Sex) %>%
  summarise(y_max = max(mean_PINK1, na.rm = TRUE))

corr_tbl <- corr_tbl %>%
  left_join(max_y_tbl, by = "Sex") %>%
  mutate(
    y.label = y_max + 0.05 * abs(y_max),
    x.label = -Inf
  )

# --- 4. Plot ---
ggplot(df, aes(x = mean_Score, y = mean_PINK1, color = Sex, shape = Disease)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Sex, color = Sex), linetype = "solid") +
  geom_text(
    data = corr_tbl,
    aes(x = x.label, y = y.label, label = cor.label, color = Sex),
    hjust = -0.1, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("M" = "#1976D2", "F" = "#C2185B", "Male" = "#1976D2", "Female" = "#C2185B")) +
  scale_shape_manual(values = c("ND" = 16, "T2D" = 17)) +
  labs(
    x = "Mean Disease Score (per donor)",
    y = "Mean PINK1 Expression (per donor)",
    title = "Correlation of Mean PINK1 Expression with Disease Score by Donor",
    color = "Sex",
    shape = "Disease"
  ) +
  theme_minimal(base_size = 14)


# ================================
# Figure 5 — PINK1 Regulators (Sex-aware DotPlot)
# ================================

# --- 0) LIBRARIES ---
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- 1) INPUT CHECKS ---
stopifnot(exists("beta_cells"))
stopifnot(all(c("beta_cluster","diabetes_status","Sex") %in% colnames(beta_cells@meta.data)))
DefaultAssay(beta_cells) <- "RNA"

# --- 2) DEFINE PINK1 REGULATOR GENE SET ---
pink1_regulators_raw <- c(
  "NFE2L2",  # NRF2 (official symbol NFE2L2)
  "NRF1",
  "ATF4",
  "TP53",
  "CREB1",
  "PARL",
  "USP30",
  "PARK7"
)
pink1_regulators <- unique(pink1_regulators_raw)

# Filter to genes present in object (warn if missing)
present_genes <- intersect(pink1_regulators, rownames(beta_cells))
missing_genes <- setdiff(pink1_regulators, present_genes)
if (length(missing_genes)) {
  message("Warning: genes not found and dropped: ", paste(missing_genes, collapse = ", "))
}
stopifnot(length(present_genes) > 0)

# Optional: single category for faceting (extendable)
gene_categories <- list("PINK1 Regulators" = present_genes)

# --- 3) PREPARE GROUPS (β-cluster + disease + Sex) ---
beta_cells$Sex <- as.character(beta_cells$Sex)
beta_cells$Sex <- ifelse(beta_cells$Sex %in% c("Male","M"), "M",
                         ifelse(beta_cells$Sex %in% c("Female","F"), "F", beta_cells$Sex))

beta_cells$group_id_base <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, sep = "_")
beta_cells$group_id_sex  <- paste(beta_cells$beta_cluster, beta_cells$diabetes_status, beta_cells$Sex, sep = "_")

# --- 4) BUILD DOTPLOT (Sex-aware aggregation) ---
p_raw <- DotPlot(
  beta_cells,
  features  = present_genes,
  group.by  = "group_id_sex",
  dot.scale = 6,
  scale     = TRUE,
  col.min   = -1,
  col.max   =  1
)
dotdata <- p_raw$data

# --- 5) RESTORE GROUP PARTS & ANNOTATE ---
dotdata <- dotdata %>%
  tidyr::separate(id, into = c("beta_cluster","diabetes_status","Sex"), sep = "_", remove = FALSE) %>%
  mutate(
    group_base = paste(beta_cluster, diabetes_status, sep = "_"),
    Sex = factor(Sex, levels = c("M","F"), labels = c("Male","Female"))
  )

gene_category_map <- stack(gene_categories)
colnames(gene_category_map) <- c("gene","gene_category")
dotdata <- left_join(dotdata, gene_category_map, by = c("features.plot" = "gene")) %>%
  mutate(gene_category = factor(gene_category, levels = names(gene_categories)))

# Desired Y-order: β4..β1; T2D then ND
desired_order <- c("β4_T2D","β3_T2D","β2_T2D",
                   "β4_ND","β3_ND","β2_ND","β1_ND")
dotdata$group_base <- factor(dotdata$group_base, levels = desired_order)

# Lock gene order (optional: keep your input order)
dotdata$features.plot <- factor(dotdata$features.plot, levels = present_genes)

# --- 6) PLOT ---
fig5_pink1 <- ggplot(dotdata, aes(x = features.plot, y = group_base, size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradient2(low = "white", high = "red", midpoint = -1) +
  facet_grid(. ~ gene_category + Sex, scales = "free_x", space = "free_x") +
  labs(
    title = "PINK1 Regulatory Genes in Human β-cells",
    x = "Gene",
    y = "β-cell Subtype & Diabetes Status",
    fill = "Avg Expr (scaled)",
    size = "% Expressed"
  ) +
  theme_light() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 10, face = "bold", colour = "black"),
    axis.text.y   = element_text(size = 10, face = "bold", colour = "black"),
    strip.text.x  = element_text(size = 12, face = "bold"),
    plot.title    = element_text(size = 14, face = "bold"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(fig5_pink1)


# ================================
# OPTIONAL: split.by (half-dots) version
# ================================
# p_split <- DotPlot(
#   beta_cells,
#   features  = present_genes,
#   group.by  = "group_id_base",
#   split.by  = "Sex",
#   dot.scale = 6,
#   scale     = TRUE,
#   col.min   = -2,
#   col.max   =  2
# ) + RotatedAxis()
# ggsave("Figure5_PINK1_Regulators_SplitBySex.png", p_split, width = 11, height = 5.5, units = "in", dpi = 300)



# pseudo-perturbation
# --- 1. Subset to β1, ND cells ---
beta1_nd_cells <- subset(
  beta_cells, 
  subset = (beta_cluster == "β1" & diabetes_status == "ND")
)

# --- 1. Extract PINK1 raw counts and Sex for each β1 ND cell ---
# --- 1. Extract raw counts and Sex vector ---
pink1_counts <- FetchData(beta1_nd_cells, "PINK1", slot = "counts")[, 1]
sex_vec <- beta1_nd_cells$Sex

df <- data.frame(
  PINK1_count = pink1_counts,
  Sex = sex_vec
)

# --- 2a. Overlaid density plot (most visually intuitive for differences) ---
p1 <- ggplot(df, aes(x = PINK1_count, fill = Sex, color = Sex)) +
  geom_density(
    position = "identity", 
    alpha = 0.3, 
    adjust = 2 # Smoother curve, change if needed
  ) +
  scale_fill_manual(values = c("M" = "#1976D2", "F" = "#C2185B")) +
  scale_color_manual(values = c("M" = "#1976D2", "F" = "#C2185B")) +
  labs(
    x = "PINK1 Raw RNA Count",
    y = "Density",
    fill = "Sex",
    color = "Sex",
    title = "Distribution of PINK1 Raw Counts in β1 ND Cells by Sex"
  ) +
  theme_minimal(base_size = 15)

# --- 2b. Side-by-side histogram (bars side by side for each count) ---
p2 <- ggplot(df, aes(x = PINK1_count, fill = Sex)) +
  geom_histogram(
    position = "stack",
    binwidth = 1,
    color = "black",
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("M" = "#1976D2", "F" = "#C2185B")) +
  labs(
    x = "PINK1 Raw RNA Count",
    y = "Number of Cells",
    fill = "Sex",
    title = "PINK1 Expression in β1 ND Cells by Sex (Histogram)"
  ) +
  theme_minimal(base_size = 15)
(p1 / p2) + plot_annotation(
  title = "Distribution of PINK1 Expression in β1 ND Cells by Sex"
)

# --- 2c. If you want to stack bars but see color by sex, use 'position = "stack"' (default) ---

# --- 1. Create combined group label if not already present ---
beta1_nd_cells$beta_sex_diab <- paste(beta1_nd_cells$beta_cluster, beta1_nd_cells$Sex, beta1_nd_cells$diabetes_status, sep = "_")
# Example: "β1_M_ND"

# --- 2. Assign PINK1 status as before ---
pink1_counts <- FetchData(beta1_nd_cells, "PINK1", slot = "counts")[, 1]
beta1_nd_cells$PINK1_status <- "Other"
beta1_nd_cells$PINK1_status[pink1_counts >= 2] <- "PINK1_pos"
beta1_nd_cells$PINK1_status[pink1_counts == 0] <- "PINK1_neg"

# --- 3. Run DE for each sex using group.by ---
de_results <- list()

for (sex in c("M", "F")) {
  # Subset for β1 ND cells of this sex
  ident_1 <- paste("β1", sex, "ND_PINK1_pos", sep = "_")
  ident_2 <- paste("β1", sex, "ND_PINK1_neg", sep = "_")
  
  # Assign a combined identity for DE (cluster_sex_diab_PINK1status)
  beta1_nd_cells$beta_sex_diab_PINK1 <- paste(
    beta1_nd_cells$beta_cluster, 
    beta1_nd_cells$Sex, 
    beta1_nd_cells$diabetes_status, 
    beta1_nd_cells$PINK1_status,
    sep = "_"
  )
  
  # Subset only cells of this group
  cells_use <- WhichCells(beta1_nd_cells, 
                          expression = beta_sex_diab_PINK1 %in% c(ident_1, ident_2)
  )
  sub <- subset(beta1_nd_cells, cells = cells_use)
  Idents(sub) <- sub$beta_sex_diab_PINK1
  
  # Differential expression
  de <- FindMarkers(
    object = sub,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.1, # 10% change in natural log scale
    pseudocount.use = 1,
    assay = "RNA",
    group.by = NULL, # Already set identity above
    only.pos = FALSE,
    return.thresh = 0.1
  )
  
  de_results[[sex]] <- de
}

# --- 4. Access as ---
de_male <- de_results[["M"]]
de_female <- de_results[["F"]]


# For T2D
# --- 0. Load required packages ---
library(Seurat)
library(dplyr)
library(ggplot2)

# --- 1. Define sexes and initialize storage for DE results ---
sexes <- c("M", "F")  # Adjust labels if needed (e.g., "Male", "Female")
de_res_M <- NULL
de_res_F <- NULL

# --- 2. Loop over each sex ---
for (sx in sexes) {
  cat("\n-------------------\nAnalyzing Sex:", sx, "\n-------------------\n")
  
  # 2A. Subset ND β1 cells for this sex and keep only PINK1+ (>=2 counts)
  nd_b1 <- subset(beta_cells, subset = (beta_cluster == "β1" & diabetes_status == "ND" & Sex == sx))
  pink1_counts_nd_b1 <- FetchData(nd_b1, "PINK1", slot = "counts")[, 1]
  nd_b1$PINK1_status <- ifelse(pink1_counts_nd_b1 >= 2, "PINK1_pos", "Other")
  nd_b1_pos <- subset(nd_b1, subset = (PINK1_status == "PINK1_pos"))
  nd_b1_pos$Group <- "ND_b1_PINK1_pos"
  
  # 2B. Subset T2D β2/3/4 cells for this sex and keep only PINK1– (==0 counts)
  t2d_b234 <- subset(beta_cells, subset = (beta_cluster %in% c("β2", "β3", "β4") & diabetes_status == "T2D" & Sex == sx))
  pink1_counts_t2d_b234 <- FetchData(t2d_b234, "PINK1", slot = "counts")[, 1]
  t2d_b234$PINK1_status <- ifelse(pink1_counts_t2d_b234 == 0, "PINK1_neg", "Other")
  t2d_b234_neg <- subset(t2d_b234, subset = (PINK1_status == "PINK1_neg"))
  t2d_b234_neg$Group <- "T2D_b234_PINK1_neg"
  
  # 2C. Merge objects for comparison
  merged <- merge(nd_b1_pos, y = t2d_b234_neg)
  
  # 2D. Set group identities for DE analysis
  Idents(merged) <- merged$Group
  cat("Cell counts per group:\n")
  print(table(merged$Group))
  
  # 2E. Differential Expression Analysis (skip if <10 cells in any group)
  if (all(table(merged$Group) >= 10)) {
    de_res <- FindMarkers(
      object = merged,
      ident.1 = "ND_b1_PINK1_pos",
      ident.2 = "T2D_b234_PINK1_neg",
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.137504,   # ~10% change
      pseudocount.use = 1,
      assay = "RNA",
      only.pos = FALSE,
      return.thresh = 0.1
    )
    cat("Top differentially expressed genes:\n")
    print(head(de_res, 10))
    
    # 2F. Save results in appropriate dataframe
    if (sx == "M") {
      de_res_M <- de_res
    } else if (sx == "F") {
      de_res_F <- de_res
    }
    
    # 2G. Write results to CSV for external analysis
    write.csv(de_res, paste0("DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_", sx, ".csv"))
  } else {
    cat("Not enough cells in one or both groups for sex:", sx, "\n")
    # Optional: assign NA to dataframe if desired
    if (sx == "M") {
      de_res_M <- NA
    } else if (sx == "F") {
      de_res_F <- NA
    }
  }
}

# --- 3. At the end you will have two dataframes:
de_res_M  # Male DE results
de_res_F  # Female DE results
tmp <- read.csv("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/ptseqbeta1_vs_t2d/DOWN/DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_F.csv")
str(tmp)
head(tmp)

write.csv(de_res_F, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\DGE\wilcox\pseudoperterbseq\beta1_vs_t2d\DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_F.csv)")
write.csv(de_res_M, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\DGE\wilcox\pseudoperterbseq\beta1_vs_t2d\DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_M.csv)")

#ORA Perturb
# List the DGE folders to analyze
dge_folders <- list(
  ptseqbeta1_vs_t2d = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/pseudoperterbseq/beta1_vs_t2d"
)

# Set base ORA output directory
ora_base <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ORA"

for (type in names(dge_folders)) {
  input_dir <- dge_folders[[type]]
  output_up <- file.path(ora_base, type, "UP")
  output_down <- file.path(ora_base, type, "DOWN")
  
  # Create output directories if not exist
  dir.create(output_up, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_down, showWarnings = FALSE, recursive = TRUE)
  
  # List CSV files
  dge_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for (file in dge_files) {
    # Read the file
    dat <- read.csv(file, row.names = 1)
    # Standardize avg_log2FC column if needed
    if (!"avg_log2FC" %in% names(dat)) {
      logfc_col <- grep("log2FC", names(dat), value = TRUE)
      if (length(logfc_col) == 1) names(dat)[names(dat) == logfc_col] <- "avg_log2FC"
    }
    if (!"p_val_adj" %in% names(dat)) {
      adjp_col <- grep("adj", names(dat), value = TRUE)
      if (length(adjp_col) == 1) names(dat)[names(dat) == adjp_col] <- "p_val_adj"
    }
    # Get gene names
    dat$gene <- rownames(dat)
    # Genes up/down
    up <- filter(dat, p_val_adj < 0.05 & avg_log2FC > 0)$gene
    down <- filter(dat, p_val_adj < 0.05 & avg_log2FC < 0)$gene
    
    # GO for upregulated
    if (length(up) > 0) {
      go_up <- gost(up, organism = "hsapiens", significant = TRUE, user_threshold = 0.05,
                    correction_method = "fdr", domain_scope = "annotated", sources = "GO", evcodes = TRUE)
      go_up_res <- as.data.frame(go_up$result)
      # Optional: rename FDR col
      if ("p_value" %in% names(go_up_res)) names(go_up_res)[names(go_up_res) == "p_value"] <- "hypergeometric FDR"
      # Remove columns (if present)
      drop_cols <- c("parents", "source_order", "effective_domain_size", "query", "precision", "recall", "evidence_codes")
      go_up_res <- select(go_up_res, -any_of(drop_cols))
      # Write
      outfile <- file.path(output_up, paste0(basename(file)))
      write.csv(go_up_res, outfile, row.names = FALSE)
    }
    # GO for downregulated
    if (length(down) > 0) {
      go_down <- gost(down, organism = "hsapiens", significant = TRUE, user_threshold = 0.05,
                      correction_method = "fdr", domain_scope = "annotated", sources = "GO", evcodes = TRUE)
      go_down_res <- as.data.frame(go_down$result)
      if ("p_value" %in% names(go_down_res)) names(go_down_res)[names(go_down_res) == "p_value"] <- "hypergeometric FDR"
      go_down_res <- select(go_down_res, -any_of(drop_cols))
      outfile <- file.path(output_down, paste0(basename(file)))
      write.csv(go_down_res, outfile, row.names = FALSE)
    }
  }
}

cat("GO/ORA completed for all DGE folders!\n")

# MAKE REPRODUCIBLE SANKEY PLOTS
# Install if needed:
# remotes::install_github("davidsjoberg/ggsankey")
# install.packages(c("readr", "dplyr", "ggalluvial", "ggplot2", "cowplot"))
make_sankey_plot <- function(my_comparison_colors, my_pathway_colors, my_ora_dir, my_pathways) {
  # --- Libraries ---
  library(readr)
  library(dplyr)
  library(ggalluvial)
  library(ggplot2)
  library(cowplot)
  library(scales)
  
  # --- Helper: Robust pathway name cleaning ---
  nukefix <- function(x) {
    x <- gsub("β", "b", x)
    x <- trimws(x)
    x <- iconv(x, to = "ASCII//TRANSLIT")
    x <- tolower(x)
    x <- gsub("\\s+", " ", x)
    x
  }
  
  # --- DEBUG: Show what files and pathways will be used ---
  my_comp_names <- list.files(my_ora_dir, pattern = "\\.csv$", full.names = FALSE)
  cat("Found files:\n"); print(my_comp_names)
  cat("Pathways (requested):\n"); print(my_pathways)
  
  all_links <- list()
  for (fname in my_comp_names) {
    comp <- gsub("\\.csv$", "", fname)
    fpath <- file.path(my_ora_dir, fname)
    if (!file.exists(fpath)) next
    df <- read_csv(fpath, show_col_types = FALSE)
    if (!"term_name" %in% colnames(df)) next
    # Try both "hypergeometric.FDR" and "hypergeometric FDR"
    fdr_col <- if ("hypergeometric.FDR" %in% colnames(df)) "hypergeometric.FDR"
    else if ("hypergeometric FDR" %in% colnames(df)) "hypergeometric FDR"
    else NA
    if (is.na(fdr_col)) next
    
    # Only keep significant pathways
    df <- df %>% filter(!is.na(.data[[fdr_col]]), .data[[fdr_col]] < 0.05, !is.na(term_name))
    df$Pathway <- nukefix(df$term_name)
    
    # DEBUG: Show what pathways are actually present
    if (nrow(df) > 0) cat("In", fname, "found pathways:\n", paste(unique(df$Pathway), collapse=", "), "\n")
    
    for (p in my_pathways) {
      pw <- nukefix(p)
      match_row <- df[df$Pathway == pw,]
      if (nrow(match_row) == 0) next
      all_links[[length(all_links)+1]] <- data.frame(
        Comparison = comp,
        Pathway = pw,
        Value = -log10(match_row[[fdr_col]][1]) + 1,
        stringsAsFactors = FALSE
      )
    }
  }
  
  sankey_df <- dplyr::bind_rows(all_links)
  sankey_df$Comparison <- factor(sankey_df$Comparison, levels = names(my_comparison_colors))
  sankey_df$Pathway <- factor(sankey_df$Pathway, levels = nukefix(my_pathways))
  
  # --- Stop gracefully if nothing matched ---
  if (nrow(sankey_df) == 0) {
    stop("No matching pathways were found after filtering. Check pathway names and FDR thresholds.\n")
  }
  
  # --- Main Sankey plot ---
  p_main <- ggplot(sankey_df, aes(axis1 = Comparison, axis2 = Pathway, y = Value)) +
    geom_alluvium(aes(fill = Comparison), width = 1/12, alpha = 0.8) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_discrete(limits = c("Comparison", "Pathway"), expand = c(.1, .1)) +
    scale_fill_manual(
      values = c(my_comparison_colors, my_pathway_colors),
      breaks = names(my_comparison_colors)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold")
    )
  
  # --- Dot legend for FDR values ---
  fdr_values <- sankey_df$Value
  legend_vals <- pretty(fdr_values, n = 5)
  legend_vals <- sort(unique(legend_vals[legend_vals > 0]), decreasing = TRUE)
  rescale_size <- function(x, to = c(3, 10)) {
    rng <- range(x, na.rm = TRUE)
    scales::rescale(x, to = to, from = rng)
  }
  dot_sizes <- rescale_size(legend_vals, to = c(3, 10))
  legend_df <- data.frame(
    x = 1,
    y = seq(50, 38, length.out = length(dot_sizes)),
    y_label = seq(50, 38, length.out = length(dot_sizes)),
    size = dot_sizes,
    label = paste0("-log10(FDR)+1 = ", legend_vals)
  )
  
  p_legend <- ggplot(legend_df) +
    geom_point(aes(x = x, y = y, size = size), shape = 21, fill = "steelblue", color = "black", stroke = 0.25) +
    geom_text(aes(x = x + 0.4, y = y_label, label = label), hjust = 0, vjust = 0.5, size = 4) +
    theme_void() +
    coord_cartesian(clip = "off") +
    scale_size_identity() +
    scale_x_continuous(limits = c(0.9, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(plot.margin = margin(5, 20, 5, 5))
  
  # --- Combine ---
  cowplot::plot_grid(
    p_main,
    p_legend,
    rel_widths = c(4.2, 1.2),
    nrow = 1,
    axis = "none",
    align = "none"
  )
}

## PINK1 vs T2D DOWN
## Load prereqs
# --- Colors ---
my_comparison_colors <- c(
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_M" = "#274472",
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_F" = "#8C3333"
)

my_pathway_colors <- c(
  # Protein synthesis & ribosome (matte red)
  "translation"                        = "#C0392B",  # Matte muted red
  "protein metabolic process"                    = "#C0392B",
  
  # Cell Stress
  "cellular response to stress"           = "#D35400",  # Matte orange
  "stress granule assembly"          = "#D35400",
  "stress response to copper ion" = "#D35400",
  "response to endoplasmic reticulum stress"     = "#D35400",
  
  # Antioxidant/ROS (matte green)
  "cellular oxidant detoxification"    = "#27AE60",  # Muted green (same as before, but fits matte)
  "antioxidant activity"               = "#27AE60",
  
  # Insulin/hormone signaling (matte blue)
  "response to insulin"                = "#2980B9",  # Muted blue (less saturated than "#1976D2")
  "insulin receptor signaling pathway" = "#2980B9",
  "hormone secretion"                  = "#2980B9",
  "insulin processing"                 = "#2980B9",
  "insulin secretion"                  = "#2980B9",
  
  # Lipid metabolism (matte purple)
  "lipid metabolic process"            = "#7D3C98",  # Matte muted purple
  "response to fatty acid"             = "#7D3C98",
  
  # Glucose/glycemic (matte brown)
  "glucose homeostasis"                = "#8D6748",  # Matte brown
  "response to glucagon"               = "#8D6748",
  "glucose metabolic process"          = "#8D6748"
)

# --- Input ---
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/ptseqbeta1_vs_t2d/DOWN"
my_comp_names <- c(
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_M.csv",
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_F.csv"
)

my_pathways <- c(
  # Protein synthesis & ribosome (matte red)
  "translation",
  "protein metabolic process",
  
  # Cell Stress (matte orange)
  "cellular response to stress",
  "stress granule assembly",
  "stress response to copper ion",
  "response to endoplasmic reticulum stress",
  
  # Antioxidant/ROS (matte green)
  "cellular oxidant detoxification",
  "antioxidant activity",
  
  # Insulin/hormone signaling (matte blue)
  "response to insulin",
  "insulin receptor signaling pathway",
  "hormone secretion",
  "insulin processing",
  "insulin secretion",
  
  # Lipid metabolism (matte purple)
  "lipid metabolic process",
  "response to fatty acid",
  
  # Glucose/glycemic (matte brown)
  "glucose homeostasis",
  "response to glucagon",
  "glucose metabolic process"
)


# RUN plot
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/ptseqbeta1_vs_t2d/DOWN",
  my_pathways = my_pathways
)

## ## PINK1 vs T2D UP
## Load prereqs
# --- Colors ---
my_comparison_colors <- c(
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_M" = "#274472",
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_F" = "#8C3333"
)

my_pathway_colors <- c(
  # Protein synthesis & ribosome (matte red)
  "translation"                        = "#C0392B",  # Matte muted red
  "protein metabolic process"                    = "#C0392B",
  
  # Mitophagy (matte brown)
  "macroautophagy"                = "#8D6748",
  "autophagy of mitochondrion"   = "#8D6748",
  "autophagy"                    = "#8D6748",
  
  # Oxphos & mitochondrial (matte orange)
  "aerobic respiration"              = "#D35400",
  "cellular respiration"             = "#D35400",
  "oxidative phosphorylation"             = "#D35400",
  "proton transmembrane transport"        = "#D35400",
  "atp synthesis coupled electron transport"      = "#D35400",
  
  # UPR (matte green)
  "erad pathway"                         = "#27AE60",
  "perk-mediated unfolded protein response" = "#27AE60",
  "protein folding"                      = "#27AE60",
  "'de novo' protein folding"         = "#27AE60",
  "response to endoplasmic reticulum stress" = "#27AE60",
  "protein stabilization" = "#27AE60",
  
  # Protein Maturation and Transport
  "protein maturation" = "#2980B9",
  "golgi vesicle transport"          = "#2980B9",
  "intracellular protein transport"                            = "#2980B9",
  "endoplasmic reticulum to golgi vesicle-mediated transport"  = "#2980B9",
  "chaperone-mediated protein folding"           = "#2980B9"
  
  # Inflammation (matte purple)
  # "p38MAPK cascade"                                 = "#7D3C98",
  # "cellular response to interleukin-7"              = "#7D3C98",
  # "response to virus"                               = "#7D3C98",
  # "cellular response to interleukin-4"              = "#7D3C98",
  # "interleukin-11-mediated signaling pathway"       = "#7D3C98",
  # "positive regulation of NF-kappaB transcription factor activity" = "#7D3C98",
  # "inflammatory response"                           = "#7D3C98",
  # "tumor necrosis factor production"                = "#7D3C98",
  
  # Dedifferentiation / Stem-like Reprogramming (matte teal)
  # "epithelial to mesenchymal transition"            = "#148F77",  # Matte teal
  # "positive regulation of cell fate commitment"     = "#148F77",
  # "neurogenesis"                                    = "#148F77",
  # "chromatin remodeling"                            = "#148F77"
)

# --- Input ---
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/ptseqbeta1_vs_t2d/UP"
my_comp_names <- c(
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_M.csv",
  "DE_ND_b1_PINK1pos_vs_T2D_b234_PINK1neg_F.csv"
)

my_pathways <- c(
  # Protein synthesis & ribosome
  "translation",
  "protein metabolic process",
  
  # Mitophagy (matte brown)
  "macroautophagy",
  "autophagy of mitochondrion",
  "autophagy",
  
  # Oxphos & mitochondrial (matte orange)
  "aerobic respiration",
  "cellular respiration",
  "oxidative phosphorylation",
  "proton transmembrane transport",
  "atp synthesis coupled electron transport",
  
  # UPR (matte green)
  "erad pathway",
  "perk-mediated unfolded protein response",
  "protein folding",
  "'de novo' protein folding",
  "response to endoplasmic reticulum stress",
  "protein stabilization",
  
  # Protein Maturation and Transport (matte blue)
  "protein maturation",
  "golgi vesicle transport",
  "intracellular protein transport",
  "endoplasmic reticulum to golgi vesicle-mediated transport",
  "chaperone-mediated protein folding"
)

# RUN plot
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/ptseqbeta1_vs_t2d/UP",
  my_pathways = my_pathways
)



############################################################
# 0) Libraries
############################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(gprofiler2)  # ORA (gost)
  library(ggalluvial)
  library(cowplot)
  library(scales)
})

############################################################
# 1) SETTINGS / OUTPUT PATHS
############################################################
base_out <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai"
dge_base <- file.path(base_out, "DATA/DGE/wilcox/ND_beta_PINK1pos_vs_neg")
ora_base <- file.path(base_out, "DATA/ORA/ND_beta_PINK1pos_vs_neg")

dir.create(dge_base, recursive = TRUE, showWarnings = FALSE)
dir.create(ora_base, recursive = TRUE, showWarnings = FALSE)

# knobs for DE
pink1_pos_min_counts <- 2L
pink1_neg_counts     <- 0L
min_cells_per_group  <- 10L
logfc_thr            <- 0.137504   # ~10% change (ln)
min_pct_expr         <- 0.10
use_assay            <- "RNA"
sex_levels           <- c("M","F")

############################################################
# 2) HELPERS
############################################################

# Tag cells as PINK1+ / PINK1- based on counts (fallback to data if counts missing)
tag_pink1_status <- function(obj, gene = "PINK1",
                             pos_min = pink1_pos_min_counts,
                             neg_eq  = pink1_neg_counts) {
  pv <- tryCatch(FetchData(obj, gene, slot = "counts")[,1],
                 error = function(e) FetchData(obj, gene, slot = "data")[,1])
  status <- rep("Other", length(pv))
  status[pv >= pos_min] <- "PINK1_pos"
  status[pv == neg_eq]  <- "PINK1_neg"
  status
}

# Wrapper for FindMarkers with consistent args
safe_findmarkers <- function(obj, id1, id2,
                             test.use = "wilcox",
                             min.pct = min_pct_expr,
                             logfc.threshold = logfc_thr,
                             assay = use_assay) {
  FindMarkers(
    object = obj,
    ident.1 = id1,
    ident.2 = id2,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    pseudocount.use = 1,
    assay = assay,
    only.pos = FALSE,
    return.thresh = 1
  )
}

# Harmonize common DE column names (avg_log2FC / p_val_adj)
standardize_dge_cols <- function(df) {
  if (!"avg_log2FC" %in% names(df)) {
    logfc_col <- grep("log2FC|avg_log2FC|avg_logFC|avg_log|logFC",
                      names(df), value = TRUE, ignore.case = TRUE)
    if (length(logfc_col) >= 1) {
      names(df)[match(logfc_col[1], names(df))] <- "avg_log2FC"
    }
  }
  if (!"p_val_adj" %in% names(df)) {
    adjp_col <- grep("p_adj|p.val.adj|p_val_adj|padj|FDR|qvalue",
                     names(df), value = TRUE, ignore.case = TRUE)
    if (length(adjp_col) >= 1) {
      names(df)[match(adjp_col[1], names(df))] <- "p_val_adj"
    }
  }
  df
}

# Normalize term names (like your nukefix)
nukefix <- function(x) {
  x <- gsub("β", "b", x)
  x <- trimws(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("\\s+", " ", x)
  x
}

############################################################
# 3) DEFINE ND β SETS (β1 only; β1–β4 pooled)
############################################################

# β1 ND only
b1_nd <- subset(beta_cells, subset = beta_cluster == "β1" & diabetes_status == "ND")
b1_nd$PINK1_status <- tag_pink1_status(b1_nd)
b1_nd <- subset(b1_nd, subset = PINK1_status %in% c("PINK1_pos","PINK1_neg"))
Idents(b1_nd) <- b1_nd$PINK1_status
cat("β1 ND counts:\n"); print(table(Idents(b1_nd)))

# βALL (β1–β4) ND pooled
bALL_nd <- subset(beta_cells, subset = beta_cluster %in% c("β1","β2","β3","β4") & diabetes_status == "ND")
bALL_nd$PINK1_status <- tag_pink1_status(bALL_nd)
bALL_nd <- subset(bALL_nd, subset = PINK1_status %in% c("PINK1_pos","PINK1_neg"))
Idents(bALL_nd) <- bALL_nd$PINK1_status
cat("All β ND counts:\n"); print(table(Idents(bALL_nd)))

############################################################
# 4) DE ANALYSIS → write CSVs (pooled & per-sex)
############################################################
dge_files <- c()

# β1 ND pooled
if (all(table(Idents(b1_nd)) >= min_cells_per_group)) {
  de <- safe_findmarkers(b1_nd, "PINK1_pos", "PINK1_neg")
  de <- standardize_dge_cols(de)
  fp <- file.path(dge_base, "DGE_b1_ND_PINK1pos_vs_PINK1neg_pooled.csv")
  write.csv(de, fp)
  dge_files <- c(dge_files, fp)
} else {
  message("β1 ND pooled: insufficient cells; skipping DE.")
}

# βALL ND pooled
if (all(table(Idents(bALL_nd)) >= min_cells_per_group)) {
  de <- safe_findmarkers(bALL_nd, "PINK1_pos", "PINK1_neg")
  de <- standardize_dge_cols(de)
  fp <- file.path(dge_base, "DGE_bALL_ND_PINK1pos_vs_PINK1neg_pooled.csv")
  write.csv(de, fp)
  dge_files <- c(dge_files, fp)
} else {
  message("All β ND pooled: insufficient cells; skipping DE.")
}

# Per-sex
for (sx in sex_levels) {
  # β1 ND per sex
  b1_sx <- subset(b1_nd, subset = Sex == sx)
  if (ncol(b1_sx) > 0) {
    Idents(b1_sx) <- b1_sx$PINK1_status
    if (all(table(Idents(b1_sx)) >= min_cells_per_group)) {
      de <- safe_findmarkers(b1_sx, "PINK1_pos", "PINK1_neg")
      de <- standardize_dge_cols(de)
      fp <- file.path(dge_base, paste0("DGE_b1_ND_PINK1pos_vs_PINK1neg_", sx, ".csv"))
      write.csv(de, fp)
      dge_files <- c(dge_files, fp)
    } else {
      message("β1 ND ", sx, ": insufficient cells; skipping.")
    }
  }
  # βALL ND per sex
  bALL_sx <- subset(bALL_nd, subset = Sex == sx)
  if (ncol(bALL_sx) > 0) {
    Idents(bALL_sx) <- bALL_sx$PINK1_status
    if (all(table(Idents(bALL_sx)) >= min_cells_per_group)) {
      de <- safe_findmarkers(bALL_sx, "PINK1_pos", "PINK1_neg")
      de <- standardize_dge_cols(de)
      fp <- file.path(dge_base, paste0("DGE_bALL_ND_PINK1pos_vs_PINK1neg_", sx, ".csv"))
      write.csv(de, fp)
      dge_files <- c(dge_files, fp)
    } else {
      message("All β ND ", sx, ": insufficient cells; skipping.")
    }
  }
}

cat("Wrote DGE files:\n"); print(basename(dge_files))

############################################################
# 5) ORA (gost) FOR UP/DOWN → write ORA CSVs
############################################################
ora_up_dir   <- file.path(ora_base, "UP")
ora_down_dir <- file.path(ora_base, "DOWN")
dir.create(ora_up_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(ora_down_dir, recursive = TRUE, showWarnings = FALSE)

run_ora_for_file <- function(csv_path, up_dir, down_dir) {
  dat <- read.csv(csv_path, row.names = 1, check.names = FALSE)
  dat <- standardize_dge_cols(dat)
  dat$gene <- rownames(dat)
  
  up   <- dat %>% filter(!is.na(p_val_adj), p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene)
  down <- dat %>% filter(!is.na(p_val_adj), p_val_adj < 0.05, avg_log2FC < 0) %>% pull(gene)
  
  run_gost <- function(genes) {
    if (length(genes) < 3) return(NULL)
    gost(genes,
         organism = "hsapiens",
         significant = TRUE,
         user_threshold = 0.05,
         correction_method = "fdr",
         domain_scope = "annotated",
         sources = "GO",
         evcodes = TRUE)
  }
  
  # UP
  if (length(up) >= 3) {
    res <- run_gost(up)
    if (!is.null(res) && !is.null(res$result)) {
      out <- res$result
      if ("p_value" %in% names(out)) names(out)[names(out) == "p_value"] <- "hypergeometric.FDR"
      drop_cols <- intersect(c("parents","source_order","effective_domain_size","query","precision","recall","evidence_codes"), names(out))
      out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
      write.csv(out, file.path(up_dir, basename(csv_path)), row.names = FALSE)
    }
  }
  # DOWN
  if (length(down) >= 3) {
    res <- run_gost(down)
    if (!is.null(res) && !is.null(res$result)) {
      out <- res$result
      if ("p_value" %in% names(out)) names(out)[names(out) == "p_value"] <- "hypergeometric.FDR"
      drop_cols <- intersect(c("parents","source_order","effective_domain_size","query","precision","recall","evidence_codes"), names(out))
      out <- out[, setdiff(names(out), drop_cols), drop = FALSE]
      write.csv(out, file.path(down_dir, basename(csv_path)), row.names = FALSE)
    }
  }
}

invisible(lapply(dge_files, run_ora_for_file, up_dir = ora_up_dir, down_dir = ora_down_dir))
cat("ORA written to:\nUP  -> ", ora_up_dir, "\nDOWN-> ", ora_down_dir, "\n")

############################################################
# 6) SANKEY PLOT (with explicit ORDERING)
############################################################
make_sankey_plot <- function(my_comparison_colors, my_pathway_colors, my_ora_dir,
                             my_pathways, comparison_order = NULL, pathway_order = NULL) {
  comp_files <- list.files(my_ora_dir, pattern = "\\.csv$", full.names = FALSE)
  if (length(comp_files) == 0) stop("No ORA CSV files found in: ", my_ora_dir)
  message("Found ORA files: ", paste(comp_files, collapse = ", "))
  
  all_links <- list()
  for (fname in comp_files) {
    comp <- gsub("\\.csv$", "", fname)
    fpath <- file.path(my_ora_dir, fname)
    if (!file.exists(fpath)) next
    df <- read_csv(fpath, show_col_types = FALSE)
    if (!"term_name" %in% colnames(df)) next
    
    fdr_col <- if ("hypergeometric.FDR" %in% colnames(df)) "hypergeometric.FDR"
    else if ("hypergeometric FDR" %in% colnames(df)) "hypergeometric FDR"
    else NA
    if (is.na(fdr_col)) next
    
    df <- df %>%
      filter(!is.na(.data[[fdr_col]]), .data[[fdr_col]] < 0.05, !is.na(term_name)) %>%
      mutate(Pathway = nukefix(term_name))
    
    # Add links only for requested pathways (normalized)
    for (p in my_pathways) {
      pw <- nukefix(p)
      hit <- df %>% filter(Pathway == pw)
      if (nrow(hit) == 0) next
      all_links[[length(all_links)+1]] <- data.frame(
        Comparison = comp,
        Pathway    = pw,
        Value      = -log10(hit[[fdr_col]][1]) + 1,
        stringsAsFactors = FALSE
      )
    }
  }
  
  sankey_df <- dplyr::bind_rows(all_links)
  if (nrow(sankey_df) == 0) stop("No matching pathways found. Check names or FDR threshold.")
  
  # ---- ORDERING controls ----
  if (is.null(pathway_order)) pathway_order <- my_pathways
  sankey_df$Pathway <- factor(sankey_df$Pathway, levels = nukefix(pathway_order))
  
  if (is.null(comparison_order)) comparison_order <- gsub("\\.csv$", "", comp_files)
  sankey_df$Comparison <- factor(sankey_df$Comparison, levels = comparison_order)
  
  # Deterministic stacking (optional)
  sankey_df <- sankey_df[order(sankey_df$Comparison,
                               as.integer(sankey_df$Pathway),
                               -sankey_df$Value), ]
  
  # Plot
  p_main <- ggplot(sankey_df, aes(axis1 = Comparison, axis2 = Pathway, y = Value)) +
    geom_alluvium(aes(fill = Comparison), width = 1/12, alpha = 0.85) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.6) +
    scale_x_discrete(limits = c("Comparison", "Pathway"), expand = c(.1, .1)) +
    scale_fill_manual(values = c(my_comparison_colors, my_pathway_colors),
                      breaks = names(my_comparison_colors)) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      axis.text.x  = element_text(size = 14, face = "bold"),
      legend.title = element_blank()
    )
  
  # Small dot legend for -log10(FDR)+1
  fdr_values <- sankey_df$Value
  legend_vals <- pretty(fdr_values, n = 5)
  legend_vals <- sort(unique(legend_vals[legend_vals > 0]), decreasing = TRUE)
  
  rescale_size <- function(x, to = c(3, 10)) {
    rng <- range(x, na.rm = TRUE)
    scales::rescale(x, to = to, from = rng)
  }
  dot_sizes <- rescale_size(legend_vals, to = c(3, 10))
  legend_df <- data.frame(
    x = 1,
    y = seq(50, 38, length.out = length(dot_sizes)),
    y_label = seq(50, 38, length.out = length(dot_sizes)),
    size = dot_sizes,
    label = paste0("-log10(FDR)+1 = ", legend_vals)
  )
  
  p_legend <- ggplot(legend_df) +
    geom_point(aes(x = x, y = y, size = size), shape = 21, fill = "steelblue", color = "black", stroke = 0.25) +
    geom_text(aes(x = x + 0.4, y = y_label, label = label), hjust = 0, vjust = 0.5, size = 4) +
    theme_void() +
    coord_cartesian(clip = "off") +
    scale_size_identity() +
    scale_x_continuous(limits = c(0.9, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(plot.margin = margin(5, 20, 5, 5))
  
  cowplot::plot_grid(p_main, p_legend, rel_widths = c(4.2, 1.2), nrow = 1)
}

############################################################
# 7) COLOR MAPS & PATHWAYS (top-5 per theme)
############################################################
my_comparison_colors <- c(
  "DGE_b1_ND_PINK1pos_vs_PINK1neg_pooled" = "#2C3E50",
  "DGE_bALL_ND_PINK1pos_vs_PINK1neg_pooled" = "#7F8C8D",
  "DGE_b1_ND_PINK1pos_vs_PINK1neg_M" = "#274472",
  "DGE_b1_ND_PINK1pos_vs_PINK1neg_F" = "#8C3333",
  "DGE_bALL_ND_PINK1pos_vs_PINK1neg_M" = "#3B6EA5",
  "DGE_bALL_ND_PINK1pos_vs_PINK1neg_F" = "#A55D5D"
)

my_pathway_colors <- c(
  # 1) Bioenergetics & ATP/Nucleotides (matte orange)
  "oxidative phosphorylation"                      = "#D35400",
  "aerobic respiration"                            = "#D35400",
  "generation of precursor metabolites and energy" = "#D35400",
  "cellular respiration"                           = "#D35400",
  "atp metabolic process"                          = "#D35400",
  
  # 2) Proteostasis & Protein Metabolism (matte red)
  "protein metabolic process"                      = "#C0392B",
  "protein folding"                                = "#C0392B",
  "protein maturation"                             = "#C0392B",
  "translation"                                    = "#C0392B",
  "response to endoplasmic reticulum stress"       = "#C0392B",
  
  # 3) Trafficking / Localization & Secretion (matte blue)
  "localization"                                   = "#2980B9",
  "transport"                                      = "#2980B9",
  "establishment of localization"                  = "#2980B9",
  "cellular localization"                          = "#2980B9",
  "intracellular transport"                        = "#2980B9",
  
  # 4) Stress, Homeostasis & Cell Fate (matte green)
  "programmed cell death"                          = "#27AE60",
  "cell death"                                     = "#27AE60",
  "apoptotic process"                              = "#27AE60",
  "cellular response to stress"                    = "#27AE60",
  "cellular response to chemical stimulus"         = "#27AE60",
  
  # 5) Mitophagy & Mitochondrial QC (matte brown)
  "chaperone-mediated autophagy"                   = "#8D6748",
  "autophagy"                                      = "#8D6748",
  "process utilizing autophagic mechanism"         = "#8D6748",
  "macroautophagy"                                 = "#8D6748",
  "autophagy of mitochondrion"                     = "#8D6748"
)

# Pathway sets for RIGHT axis (order = display order)
up_pathways <- c(
  "oxidative phosphorylation","aerobic respiration","generation of precursor metabolites and energy",
  "cellular respiration","atp metabolic process",
  "protein metabolic process","protein folding","protein maturation","translation","response to endoplasmic reticulum stress",
  "localization","transport","establishment of localization","cellular localization","intracellular transport",
  "programmed cell death","cell death","apoptotic process","cellular response to stress","cellular response to chemical stimulus",
  "chaperone-mediated autophagy","autophagy","process utilizing autophagic mechanism","macroautophagy","autophagy of mitochondrion"
)
down_pathways <- up_pathways  # same order

# Preferred LEFT-axis order (comparisons)
my_comp_order <- c(
  "DGE_b1_ND_PINK1pos_vs_PINK1neg_pooled",
  "DGE_bALL_ND_PINK1pos_vs_PINK1neg_pooled",
  "DGE_b1_ND_PINK1pos_vs_PINK1neg_M",
  "DGE_bALL_ND_PINK1pos_vs_PINK1neg_M",
  "DGE_b1_ND_PINK1pos_vs_PINK1neg_F",
  "DGE_bALL_ND_PINK1pos_vs_PINK1neg_F"
)

############################################################
# 8) RENDER SANKEY PLOTS (returns plot objects)
############################################################

# DOWN
ora_down_dir <- file.path(ora_base, "DOWN")
message("Rendering Sankey (DOWN) from: ", ora_down_dir)
p_down <- make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors    = my_pathway_colors,
  my_ora_dir           = ora_down_dir,
  my_pathways          = down_pathways,
  comparison_order     = my_comp_order,
  pathway_order        = down_pathways
)
print(p_down)

# UP
ora_up_dir <- file.path(ora_base, "UP")
message("Rendering Sankey (UP) from: ", ora_up_dir)
p_up <- make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors    = my_pathway_colors,
  my_ora_dir           = ora_up_dir,
  my_pathways          = up_pathways,
  comparison_order     = my_comp_order,
  pathway_order        = up_pathways
)
print(p_up)




beta_cells <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)")
beta_cells_metadata <- beta_cells@meta.data
beta_cells_metadata <- beta_cells_metadata %>%
  filter(celltype_qadir == "beta")
table(beta_cells_metadata$beta_cluster)
table(is.na(beta_cells_metadata$beta_cluster))

# Replace all β with b in beta_cluster column
beta_cells_metadata$beta_cluster <- gsub("β", "b", beta_cells_metadata$beta_cluster)

# Check result
table(beta_cells_metadata$beta_cluster)

# Cross-tabulate beta cluster with diabetes status
table(beta_cells_metadata$beta_cluster, beta_cells_metadata$diabetes_status)
write.csv(beta_cells_metadata, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells_metadata.csv)")

# Your gene list
genes_to_check <- c(
  "EIF2AK3", "EIF2S1", "ATF4", "MAPK8", "MAPK9", "MAPK10", "JUN", "NFKB1", "RELA", "E2F1",
  "HIF1A", "SP1", "TP53", "CITED2", "FOXO3", "MITF", "TFEB", "TFE3", "PINK1", "PRKN",
  "BNIP3", "BNIP3L", "FUNDC1", "PHB2", "NLRX1", "BCL2L13", "FKBP8", "MUL1", "SAMM50", "MTX1",
  "PGAM5", "CSNK2A1", "CSNK2A2", "SRC", "ULK1", "OPA1", "MFN2", "TBK1", "TAX1BP1", "SQSTM1",
  "CALCOCO2", "OPTN", "NBR1", "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "BECN1", "ATG9A", "USP8",
  "USP15", "USP30", "VCP", "RAB5A", "RAB7A", "MON1A", "CCZ1", "FIS1", "BCL2L1"
)

# Get all gene names in your object (adjust assay if not "RNA")
all_genes <- rownames(beta_cells)

# Which genes are present?
present_genes <- genes_to_check[genes_to_check %in% all_genes]
missing_genes <- genes_to_check[!genes_to_check %in% all_genes]

cat("Present genes (", length(present_genes), "):\n", paste(present_genes, collapse = ", "), "\n\n")
cat("Missing genes (", length(missing_genes), "):\n", paste(missing_genes, collapse = ", "), "\n")


# Beta cell genes
beta_cells <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)")
str(beta_cells)
# These columns are named like 'feature_GENENAME_per'
mito_algo <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\PINK1_discovery\result_with_all_features.csv)")
str(mito_algo)
table(mito_algo$beta_cluster, mito_algo$diabetes_status, mito_algo$Sex)

library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(gridExtra)

# --- 1. Prepare long-format feature data ---
feature_cols <- grep("^feature_.*_per$", names(mito_algo), value = TRUE)

feature_long <- mito_algo %>%
  group_by(beta_cluster, Sex, diabetes_status) %>%
  summarise(across(all_of(feature_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "Feature",
    values_to = "MeanImportance"
  ) %>%
  mutate(
    Gene = str_remove(Feature, "feature_"),
    Gene = str_remove(Gene, "_per")
  ) %>%
  filter(!(beta_cluster == "b1" & diabetes_status == "T2D"))  # Remove group with no cells

# --- 2. Make wide matrices for each sex, order by b1/ND value ---
make_heatmat_ranked <- function(sex_value) {
  df <- feature_long %>%
    filter(Sex == sex_value) %>%
    mutate(Column = paste(beta_cluster, diabetes_status, sep = "_")) %>%
    select(Gene, Column, MeanImportance) %>%
    pivot_wider(names_from = Column, values_from = MeanImportance)
  
  # Rank genes by b1_ND value (missing => bottom)
  if ("b1_ND" %in% colnames(df)) {
    df <- df %>%
      mutate(b1_nd_sort = -ifelse(is.na(b1_ND), -Inf, b1_ND)) %>%   # Use -Inf to sort NA to bottom
      arrange(b1_nd_sort) %>%
      select(-b1_nd_sort)
  }
  
  mat <- df %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # Remove all-NA rows/columns
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  mat <- mat[, colSums(is.na(mat)) < nrow(mat), drop = FALSE]
  return(mat)
}

mat_male   <- make_heatmat_ranked("M")
mat_female <- make_heatmat_ranked("F")

# --- 3. Plot both heatmaps together, gene order set by b1_ND per sex ---
phm_male <- pheatmap(
  mat_male,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  main = "Male (Genes Ranked by b1/ND)",
  fontsize_row = 9, fontsize_col = 11, angle_col = 45,
  border_color = NA, silent = TRUE
)
phm_female <- pheatmap(
  mat_female,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  main = "Female (Genes Ranked by b1/ND)",
  fontsize_row = 9, fontsize_col = 11, angle_col = 45,
  border_color = NA, silent = TRUE
)

gridExtra::grid.arrange(phm_female$gtable, phm_male$gtable, ncol = 2)

# ... [all your previous code] ...

mat_male   <- make_heatmat_ranked("M")
mat_female <- make_heatmat_ranked("F")

# --- 2b. Limit to top 30 genes ---
get_top_n <- function(mat, n = 30) {
  top_n_genes <- rownames(mat)[1:min(n, nrow(mat))]
  mat[top_n_genes, , drop = FALSE]
}
mat_male_top30   <- get_top_n(mat_male, 30)
mat_female_top30 <- get_top_n(mat_female, 30)

# --- 3. Plot both heatmaps together, now only top 30 genes shown ---
phm_male <- pheatmap(
  mat_male_top30,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  main = "Male (Top 30 Genes by b1/ND)",
  fontsize_row = 9, fontsize_col = 11, angle_col = 45,
  border_color = NA, silent = TRUE
)
phm_female <- pheatmap(
  mat_female_top30,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  main = "Female (Top 30 Genes by b1/ND)",
  fontsize_row = 9, fontsize_col = 11, angle_col = 45,
  border_color = NA, silent = TRUE
)

gridExtra::grid.arrange(phm_female$gtable, phm_male$gtable, ncol = 2)


#GSE217775 Kidney PINK1-/-
# Define your data directory (Windows paths use forward slashes or double backslashes)
data_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/GSE217775_RAW"

# List all files including subdirectories
all_files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
cat("Total files found:", length(all_files), "\n")
print(all_files)

# Optionally, show a summary by file type/extension
file_types <- tools::file_ext(all_files)
table(file_types)

# ------------------------------------------
# STEP 1: List all uncompressed .txt files
# ------------------------------------------
data_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/GSE217775_RAW"

# This finds all .txt files but skips .txt.gz files (compressed)
txt_files <- list.files(
  path = data_dir,
  pattern = "\\.txt$",         # Only files ending in .txt
  recursive = TRUE,            # Look in subfolders
  full.names = TRUE            # Return full paths (important!)
)

print(txt_files)               # Check the file list

# ------------------------------------------
# STEP 2: Function to clean sample names
# ------------------------------------------
# Use the file name without the extension as a sample name
get_sample_name <- function(path) {
  fname <- basename(path)                  # e.g. "GSM6726594_WT_4M_1.txt"
  sub("\\.txt$", "", fname)                # Remove ".txt"
}

# ------------------------------------------
# STEP 3: Read all .txt files into a named list
# ------------------------------------------
# We'll use readr::read_tsv for tab-delimited files (typical for GEO)
library(readr)
library(purrr)

# Read all files; names will be clean sample names
data_list <- setNames(
  map(txt_files, ~ read_tsv(.x)),            # Read each file as a tibble
  nm = sapply(txt_files, get_sample_name)    # Name each element by sample name
)

# ------------------------------------------
# STEP 4: Check the structure of your loaded data
# ------------------------------------------
# This prints just the top level (names and types) of the list
str(data_list, 1)

# Optional: Peek at one sample's data
head(data_list[[1]])                        # By position
# Or, by name (replace "GSM6726594_WT_4M_1" with your sample name)
head(data_list[["GSM6726594_WT_4M_1"]])

# ------------------------------------------
# STEP 5: (Optional) Combine data if structure allows
# ------------------------------------------
# If you want to combine by gene (e.g., merge all by "Gene" column):
# combined <- purrr::reduce(data_list, dplyr::full_join, by = "Gene")

# ------------------------------------------
# NOTES:
# - Each element in data_list is a tibble of one sample.
# - Adapt 'read_tsv' to 'read_delim' or 'read_csv' if file format is different.
# - Next steps depend on what these tables contain (genes × cells, genes × sample, etc.).

# ------------------------------------------
# STEP 6: Combine all sample tables into one expression matrix
# ------------------------------------------

library(dplyr)
library(purrr)

# Start by joining on "Gene symbol"
# We'll use reduce() + full_join to merge all data by gene
# Each sample has column: 'Gene symbol' + sample_name
# We want to end up with one column 'Gene symbol', and columns for each sample

combined_counts <- purrr::reduce(
  data_list,
  .f = function(x, y) dplyr::full_join(x, y, by = "Gene symbol")
)

# Check result: genes as rows, samples as columns (first col = Gene symbol)
dim(combined_counts)    # Should be 23,183 rows × 9 cols (1 gene + 8 samples)
head(combined_counts)

# Optional: set gene symbols as rownames (if needed for downstream tools)
counts_mat <- as.data.frame(combined_counts)
rownames(counts_mat) <- counts_mat$`Gene symbol`
counts_mat$`Gene symbol` <- NULL

# Now counts_mat is a genes × samples numeric matrix
head(counts_mat)

# If you want an actual matrix (not data.frame) for some packages:
counts_matrix <- as.matrix(counts_mat)

# ------------------------------------------
# STEP 7: Check for missing values
# ------------------------------------------
# See if any NA values (shouldn't be, unless some gene missing in some sample)
sum(is.na(counts_matrix))    # Should be 0 ideally

# ------------------------------------------
# Next steps: Differential expression, visualization, normalization, etc.

# -------------------------------
# STEP 1: Setup and Load Packages
# -------------------------------
library(limma)
library(edgeR)

# Set output directory for DE files
out_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/dge"
if (!dir.exists(out_dir)) dir.create(out_dir)

# -------------------------------
# STEP 2: Prepare counts & sample info
# -------------------------------
# Subset matrix columns to correct order and name for clarity
sample_order <- c("WT_4M_1", "WT_4M_2", "WT_24M_1", "WT_24M_2",
                  "PKKO4M1", "PKKO4M2", "PKKO24M1", "PKKO24M2")
counts <- counts_matrix[, sample_order]

# Assign group labels matching your design
group <- factor(c("WT_4M", "WT_4M", "WT_24M", "WT_24M", "KO_4M", "KO_4M", "KO_24M", "KO_24M"))

# -------------------------------
# STEP 3: DGEList and normalization
# -------------------------------
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# -------------------------------
# STEP 4: Design matrix (no intercept, for easy contrasts)
# -------------------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# -------------------------------
# STEP 5: Voom transformation and model fit
# -------------------------------
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)

# -------------------------------
# STEP 6: Define contrasts
# -------------------------------
contrast.matrix <- makeContrasts(
  KO_4M_vs_WT_4M    = KO_4M - WT_4M,
  KO_24M_vs_WT_24M  = KO_24M - WT_24M,
  WT_24M_vs_WT_4M   = WT_24M - WT_4M,
  KO_24M_vs_KO_4M   = KO_24M - KO_4M,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# -------------------------------
# STEP 7: Extract and Save DE results for each comparison
# -------------------------------
# Helper function to save results
save_de <- function(fit_obj, contrast_name, out_dir) {
  res <- topTable(fit_obj, coef = contrast_name, number = Inf, adjust.method = "BH")
  out_path <- file.path(out_dir, paste0("limma_", contrast_name, ".csv"))
  write.csv(res, out_path, row.names = TRUE)
  cat("Saved:", out_path, "\n")
}

# List of contrast names
contrasts <- colnames(contrast.matrix)

# Save DE results for each contrast
for (cname in contrasts) {
  save_de(fit2, cname, out_dir)
}

# -------------------------------
# STEP 8: (Optional) Check a result in R
# -------------------------------
head(topTable(fit2, coef = "KO_4M_vs_WT_4M", number = 10))

#Plotting
# -------------------------------------------------------
# STEP 1: Load libraries
# -------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(grid)

# -------------------------------------------------------
# STEP 2: Define gene groups and filter for present genes
# -------------------------------------------------------
# Your gene groups
core_mitophagy <- c("Pink1", "Park2", "Bnip3", "Bnip3l", "Fundc1", 
                    "Bcl2l13", "Fkbp8", "Phb2", "Sqstm1", "Optn")
receptors_adapters <- c("Map1lc3a", "Map1lc3b", "Calcoco2", "Nlrp3", 
                        "Ulk1", "Tbk1", "Tax1bp1", "Nbr1")
qc_stress <- c("Mfn2", "Opa1", "Dnm1l", "Mul1")
# Check which genes are present
all_mito_genes <- unique(c(core_mitophagy, receptors_adapters, qc_stress))
present_genes <- all_mito_genes[all_mito_genes %in% rownames(counts_mat)]
missing_genes <- setdiff(all_mito_genes, present_genes)
if (length(missing_genes) > 0) {
  cat("Warning: these genes are missing and will be skipped:\n")
  print(missing_genes)
}
# Keep only present genes
core_mitophagy <- core_mitophagy[core_mitophagy %in% present_genes]
receptors_adapters <- receptors_adapters[receptors_adapters %in% present_genes]
qc_stress <- qc_stress[qc_stress %in% present_genes]

# -------------------------------------------------------
# STEP 3: Calculate Z-scored group means for each group
# -------------------------------------------------------
grouping <- c(
  rep("WT_4M", 2), rep("WT_24M", 2), rep("KO_4M", 2), rep("KO_24M", 2)
)
group_levels <- unique(grouping)
prep_zscore <- function(genes, counts_mat, grouping, group_levels) {
  expr <- counts_mat[genes, , drop = FALSE]
  expr_means <- sapply(group_levels, function(g) {
    cols <- which(grouping == g)
    rowMeans(expr[, cols, drop = FALSE])
  })
  rownames(expr_means) <- genes
  colnames(expr_means) <- group_levels
  # Row-wise Z-score normalization for heatmap (across groups)
  expr_means_z <- t(scale(t(expr_means)))
  return(expr_means_z)
}
core_expr_z <- prep_zscore(core_mitophagy, counts_mat, grouping, group_levels)
receptors_expr_z <- prep_zscore(receptors_adapters, counts_mat, grouping, group_levels)
qc_expr_z <- prep_zscore(qc_stress, counts_mat, grouping, group_levels)

library(ComplexHeatmap)
library(circlize)

# -- Prepare combined z-scored matrix with all genes, grouped by function --
all_genes <- c(core_mitophagy, receptors_adapters, qc_stress)
expr_mat <- prep_zscore(all_genes, counts_mat, grouping, group_levels)

# -- Define row annotation for functional group --
gene_groups <- c(
  rep("Core Mitophagy", length(core_mitophagy)),
  rep("Receptors/Adapters", length(receptors_adapters)),
  rep("QC/Stress", length(qc_stress))
)
row_annot <- rowAnnotation(
  FunctionalGroup = factor(gene_groups, levels = c("Core Mitophagy", "Receptors/Adapters", "QC/Stress")),
  col = list(FunctionalGroup = c("Core Mitophagy"="#b2df8a", "Receptors/Adapters"="#a6cee3", "QC/Stress"="#fb9a99")),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

# -- Cluster genes within each group independently --
# Create list of row indices per group for splitting
split_vec <- gene_groups

# -- Heatmap plotting (one column label for all, no per-group titles) --
ht <- Heatmap(
  expr_mat,
  name = "Z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("dodgerblue4", "white", "firebrick4")),
  cluster_columns = FALSE,
  cluster_rows = TRUE,              # cluster within group (split)
  split = split_vec,                # split genes by functional group for clustering
  row_title = NULL,                 # no group title on side
  show_row_names = TRUE,
  show_column_names = TRUE,         # only one column label bar at bottom
  left_annotation = row_annot,     # functional group color bar
  column_names_rot = 45,
  heatmap_legend_param = list(title = "Z-score"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  width = unit(4, "inches"),
  row_names_side = "right"
)

ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

# ---- Define your mitochondrial quality gene groups ----
fusion_fission <- c("Mfn1", "Mfn2", "Opa1", "Dnm1l", "Fis1", "Mff")
biogenesis <- c("Ppargc1a", "Nrf1", "Nfe2l2", "Tfam", "Tfb1m", "Tfb2m")
oxphos <- c("Ndufa9", "Cox4i1", "Atp5a1", "Uqcrc1", "Sdha")  # adjust/expand as needed
detox_chaperones <- c("Sod2", "Gpx1", "Prdx3", "Hspd1", "Hspa9")

# ---- Check presence in your data ----
all_quality_genes <- unique(c(fusion_fission, biogenesis, oxphos, detox_chaperones))
present_quality_genes <- all_quality_genes[all_quality_genes %in% rownames(counts_mat)]
missing_quality_genes <- setdiff(all_quality_genes, present_quality_genes)
if (length(missing_quality_genes) > 0) {
  cat("Warning: missing mitochondrial quality genes:\n")
  print(missing_quality_genes)
}
fusion_fission <- fusion_fission[fusion_fission %in% present_quality_genes]
biogenesis <- biogenesis[biogenesis %in% present_quality_genes]
oxphos <- oxphos[oxphos %in% present_quality_genes]
detox_chaperones <- detox_chaperones[detox_chaperones %in% present_quality_genes]

# ---- Make combined z-scored matrix and annotation ----
quality_genes <- c(fusion_fission, biogenesis, oxphos, detox_chaperones)
expr_mat_qual <- prep_zscore(quality_genes, counts_mat, grouping, group_levels)
gene_quality_groups <- c(
  rep("Fusion/Fission", length(fusion_fission)),
  rep("Biogenesis", length(biogenesis)),
  rep("OXPHOS", length(oxphos)),
  rep("Detox/Chaperone", length(detox_chaperones))
)
row_annot_qual <- rowAnnotation(
  QualityGroup = factor(gene_quality_groups, 
                        levels = c("Fusion/Fission", "Biogenesis", "OXPHOS", "Detox/Chaperone")),
  col = list(QualityGroup = c(
    "Fusion/Fission" = "#B3CDE3",
    "Biogenesis" = "#CCEBC5",
    "OXPHOS" = "#FBB4AE",
    "Detox/Chaperone" = "#DECBE4"
  )),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

# ---- Plot quality heatmap ----
ht_qual <- Heatmap(
  expr_mat_qual,
  name = "Z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("dodgerblue4", "white", "firebrick4")),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  split = gene_quality_groups,
  row_title = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  left_annotation = row_annot_qual,
  column_names_rot = 45,
  heatmap_legend_param = list(title = "Z-score"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  width = unit(4, "inches"),
  row_names_side = "right"
)
ComplexHeatmap::draw(ht_qual, heatmap_legend_side = "right", annotation_legend_side = "right")

# ----------------------------------------------------------
# STEP 1: Load necessary libraries
# ----------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

# ----------------------------------------------------------
# STEP 2: Define cell stress gene sets (edit as needed)
# ----------------------------------------------------------
oxidative_stress <- c("Sod2", "Gpx1", "Prdx3", "Prdx5", "Cat", "Txn2", "Gsr")
uprmt <- c("Hspd1", "Hspe1", "Clpp", "Lonp1", "Dnaja3", "Atf5", "Yme1l1", "Lars2")
inflammation <- c("Nlrp3", "Il1b", "Casp1", "Tnf", "Il6", "Ccl2")
apoptosis <- c("Bax", "Bak1", "Bcl2", "Casp3", "Casp9", "Trp53", "Bad", "Bid", "Mcl1")
dna_damage <- c("Trp53", "Gadd45a", "Cdkn1a", "Atm", "Atr")

# ----------------------------------------------------------
# STEP 3: Check gene presence in your data and filter lists
# ----------------------------------------------------------
all_stress_genes <- unique(c(oxidative_stress, uprmt, inflammation, apoptosis, dna_damage))
present_stress_genes <- all_stress_genes[all_stress_genes %in% rownames(counts_mat)]
missing_stress_genes <- setdiff(all_stress_genes, present_stress_genes)
if (length(missing_stress_genes) > 0) {
  cat("Warning: missing stress pathway genes:\n")
  print(missing_stress_genes)
}
oxidative_stress <- oxidative_stress[oxidative_stress %in% present_stress_genes]
uprmt <- uprmt[uprmt %in% present_stress_genes]
inflammation <- inflammation[inflammation %in% present_stress_genes]
apoptosis <- apoptosis[apoptosis %in% present_stress_genes]
dna_damage <- dna_damage[dna_damage %in% present_stress_genes]

# ----------------------------------------------------------
# STEP 4: Make combined z-scored group mean expression matrix
# ----------------------------------------------------------
# Combine all present genes, in pathway order
stress_genes <- c(oxidative_stress, uprmt, inflammation, apoptosis, dna_damage)
expr_mat_stress <- prep_zscore(stress_genes, counts_mat, grouping, group_levels)

# Pathway group annotation vector
stress_groups <- c(
  rep("Oxidative Stress", length(oxidative_stress)),
  rep("UPRmt", length(uprmt)),
  rep("Inflammation", length(inflammation)),
  rep("Apoptosis", length(apoptosis)),
  rep("DNA Damage", length(dna_damage))
)

# ----------------------------------------------------------
# STEP 5: Row annotation for functional pathway group
# ----------------------------------------------------------
row_annot_stress <- rowAnnotation(
  Pathway = factor(
    stress_groups, 
    levels = c("Oxidative Stress", "UPRmt", "Inflammation", "Apoptosis", "DNA Damage")
  ),
  col = list(Pathway = c(
    "Oxidative Stress" = "#a6cee3",
    "UPRmt" = "#b2df8a",
    "Inflammation" = "#fb9a99",
    "Apoptosis" = "#fdbf6f",
    "DNA Damage" = "#cab2d6"
  )),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

# ----------------------------------------------------------
# STEP 6: Plot the heatmap with ComplexHeatmap
# ----------------------------------------------------------
ht_stress <- Heatmap(
  expr_mat_stress,
  name = "Z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("dodgerblue4", "white", "firebrick4")),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  split = stress_groups,               # cluster within each pathway group only
  row_title = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  left_annotation = row_annot_stress, # pathway color bar
  column_names_rot = 45,
  heatmap_legend_param = list(title = "Z-score"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  width = unit(4, "inches"),
  row_names_side = "right"
)

# Draw the heatmap
ComplexHeatmap::draw(ht_stress, heatmap_legend_side = "right", annotation_legend_side = "right")


# --------------------------------------------------------
# STEP 1: Define a mitochondrial function gene set
# (Feel free to adjust/expand as appropriate for your study)
# --------------------------------------------------------
mito_function_genes <- c(
  # OXPHOS (one representative from each complex)
  "Ndufa9", "Ndufs1", "Sdha", "Uqcrc1", "Cox4i1", "Atp5a1",
  # Fusion/Fission
  "Mfn1", "Mfn2", "Opa1", "Dnm1l", "Fis1",
  # Biogenesis/Regulation
  "Ppargc1a", "Nrf1", "Tfam", "Tfb1m", "Tfb2m",
  # Quality Control/Detox
  "Sod2", "Gpx1", "Prdx3", "Hspd1", "Hspa9"
)

# --------------------------------------------------------
# STEP 2: Filter genes to those present in your data
# --------------------------------------------------------
mito_function_genes <- mito_function_genes[mito_function_genes %in% rownames(counts_mat)]
if (length(mito_function_genes) < 5) {
  warning("Fewer than 5 mito function genes found in your data!")
}
cat("Number of mito function genes used:", length(mito_function_genes), "\n")

# --------------------------------------------------------
# STEP 3: Calculate module score for each sample
# (Mean expression of selected genes per sample)
# --------------------------------------------------------
module_scores <- colMeans(counts_mat[mito_function_genes, , drop = FALSE], na.rm = TRUE)

# --------------------------------------------------------
# STEP 4: Create a data frame with group labels for plotting
# (Assuming you have 'grouping' as before: e.g. rep("WT_4M",2), ...)
# --------------------------------------------------------
module_score_df <- data.frame(
  Sample = colnames(counts_mat),
  Group = grouping,
  ModuleScore = module_scores
)

# Optionally, calculate group mean ± SD
library(dplyr)
module_score_summary <- module_score_df %>%
  group_by(Group) %>%
  summarise(
    MeanScore = mean(ModuleScore),
    SD = sd(ModuleScore)
  )

# --------------------------------------------------------
# STEP 5: Plot module score by group
# --------------------------------------------------------
library(ggplot2)

# Option 1: Plot all samples (jitter + mean bar)
p1 <- ggplot(module_score_df, aes(x = Group, y = ModuleScore, fill = Group)) +
  geom_jitter(width = 0.1, shape = 21, size = 3, color = "black", alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 2, color = "red") +
  labs(
    title = "Mitochondrial Function Module Score (per sample)",
    x = NULL, y = "Module Score (Mean Expression)"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# Option 2: Plot group means ± SD (barplot)
p2 <- ggplot(module_score_summary, aes(x = Group, y = MeanScore, fill = Group)) +
  geom_col(width = 0.7, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = MeanScore - SD, ymax = MeanScore + SD), width = 0.2) +
  labs(
    title = "Mitochondrial Function Module Score (group mean ± SD)",
    x = NULL, y = "Module Score (Mean Expression)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# --------------------------------------------------------
# STEP 6: Print or save plot
# --------------------------------------------------------
print(p1)
print(p2)

# Optionally, save to file:
ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/dge/mito_function_module_score.png", p1, width=6, height=4, dpi=200)
ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/dge/mito_function_module_score_bar.png", p2, width=6, height=4, dpi=200)


# -------------------------------------------------------------------------
# STEP 1: Assign correct group names to each sample column
# -------------------------------------------------------------------------
# Print the sample column names for reference
print(colnames(counts_mat))
# Example output:
# [1] "WT_4M_1"  "WT_4M_2"  "WT_24M_1" "WT_24M_2" "PKKO4M1"  "PKKO4M2"  "PKKO24M1" "PKKO24M2"

# Create a grouping vector based on the order of columns in counts_mat
grouping <- c(
  "WT_4M",    # WT_4M_1
  "WT_4M",    # WT_4M_2
  "WT_24M",   # WT_24M_1
  "WT_24M",   # WT_24M_2
  "KO_4M",    # PKKO4M1
  "KO_4M",    # PKKO4M2
  "KO_24M",   # PKKO24M1
  "KO_24M"    # PKKO24M2
)
# Give each group the sample name for easy reference (optional)
names(grouping) <- colnames(counts_mat)

# -------------------------------------------------------------------------
# STEP 2: Define gene sets for each pathway and filter for presence
# -------------------------------------------------------------------------
oxphos_genes <- c("Ndufa9", "Ndufs1", "Sdha", "Uqcrc1", "Uqcrc2",
                  "Cox4i1", "Cox5a", "Atp5a1", "Atp5b")
glycolysis_genes <- c("Hk1", "Hk2", "Gpi1", "Pfkl", "Pfkm", "Aldoa", "Aldoc",
                      "Tpi1", "Gapdh", "Pgk1", "Pgam1", "Eno1", "Eno2", "Pkm", "Ldha")
tca_genes <- c("Cs", "Aco2", "Idh3a", "Idh3b", "Idh3g", "Ogdh", "Suclg1",
               "Suclg2", "Sdha", "Sdhb", "Fh1", "Mdh2")

# Keep only genes present in your matrix
oxphos_genes <- oxphos_genes[oxphos_genes %in% rownames(counts_mat)]
glycolysis_genes <- glycolysis_genes[glycolysis_genes %in% rownames(counts_mat)]
tca_genes <- tca_genes[tca_genes %in% rownames(counts_mat)]

# Optional: print how many genes are used per module
cat("OXPHOS genes used:", length(oxphos_genes), "\n")
cat("Glycolysis genes used:", length(glycolysis_genes), "\n")
cat("TCA genes used:", length(tca_genes), "\n")

# -------------------------------------------------------------------------
# STEP 3: Calculate module scores and ratios per sample
# -------------------------------------------------------------------------
# Calculate average expression per pathway per sample (module score)
oxphos_score <- colMeans(counts_mat[oxphos_genes, , drop = FALSE], na.rm = TRUE)
glycolysis_score <- colMeans(counts_mat[glycolysis_genes, , drop = FALSE], na.rm = TRUE)
tca_score <- colMeans(counts_mat[tca_genes, , drop = FALSE], na.rm = TRUE)

# Calculate log2 ratios per sample for three pathway comparisons
log2_oxphos_gly <- log2(oxphos_score / glycolysis_score)
log2_oxphos_tca <- log2(oxphos_score / tca_score)
log2_oxphos_combo <- log2(oxphos_score / (glycolysis_score + tca_score))

# -------------------------------------------------------------------------
# STEP 4: Build a dataframe with all ratios and group information
# -------------------------------------------------------------------------
ratios_df <- data.frame(
  Sample = colnames(counts_mat),
  Group = grouping,
  Log2OXPHOS_Combo = log2_oxphos_combo,
  Log2OXPHOS_Glycolysis = log2_oxphos_gly,
  Log2OXPHOS_TCA = log2_oxphos_tca
)

# Assign Genotype and Age columns using Group
ratios_df <- ratios_df %>%
  mutate(
    Genotype = ifelse(grepl("^WT", Group), "WT", "KO"),
    Age = ifelse(grepl("24M", Group), "24M", "4M")
  )

# Check that all ages and genotypes are present
print(table(ratios_df$Age, ratios_df$Genotype))

# -------------------------------------------------------------------------
# STEP 5: Summarize by group (mean and error for each ratio)
# -------------------------------------------------------------------------
library(dplyr)
library(tidyr)

# Pivot long for easier plotting
long_ratios <- ratios_df %>%
  pivot_longer(
    cols = starts_with("Log2"),
    names_to = "Ratio",
    values_to = "Value"
  )

# Calculate mean and SEM per group, age, genotype, and ratio
plot_df <- long_ratios %>%
  group_by(Group, Genotype, Age, Ratio) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    SEM = SD / sqrt(n()),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# STEP 6: Plot results as line plot (trajectory 4M to 24M for KO/WT)
# -------------------------------------------------------------------------
library(ggplot2)

# Make Age an ordered factor for proper plot order
plot_df$Age <- factor(plot_df$Age, levels = c("4M", "24M"))

# Plot means with SEM error bars, lines for genotype
ggplot(plot_df, aes(x = Age, y = Mean, group = Genotype, color = Genotype)) +
  geom_line(aes(linetype = Genotype), size = 1.2) +   # Mean line for each genotype
  geom_point(size = 4) +                              # Mean point for each group
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.15, size = 1) +
  facet_wrap(~Ratio, scales = "free_y") +             # One panel per ratio
  labs(
    title = "Log2 Pathway Ratios: 4M vs 24M by Genotype",
    y = "Log2 Ratio (mean ± SEM)",
    x = "Age"
  ) +
  scale_color_manual(values = c("WT" = "dodgerblue4", "KO" = "firebrick4")) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

# -------------------------------------------------------------------------
# STEP 7: (Optional) Plot as heatmap of group means (summary)
# -------------------------------------------------------------------------
summary_mat <- plot_df %>%
  group_by(Group, Ratio) %>%
  summarise(Mean = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Ratio, values_from = Mean)

# Convert to matrix for heatmap
ratio_mat <- as.matrix(summary_mat[,-1])
rownames(ratio_mat) <- summary_mat$Group

library(pheatmap)
pheatmap(
  ratio_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("royalblue3", "white", "red3"))(50),
  main = "Log2 Pathway Ratios (Mean per Group)",
  fontsize_number = 12,
  angle_col = 45,
  show_rownames = TRUE
)

#ORA
# -------------------------------------------------------------------------
# STEP 1: Libraries and Directory Setup
# -------------------------------------------------------------------------
library(readr)
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(forcats)

# --- Set file paths
dge_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/dge"
ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/ora"
up_dir <- file.path(ora_dir, "UP")
down_dir <- file.path(ora_dir, "DOWN")
dir.create(up_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(down_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- c("KO_4M_vs_WT_4M", "KO_24M_vs_WT_24M")
dge_files <- file.path(dge_dir, paste0("limma_", comparisons, ".csv"))

# -------------------------------------------------------------------------
# STEP 2: Robust function to extract UP/DOWN genes
# -------------------------------------------------------------------------
get_deg_lists <- function(de_file, logfc = 0, padj = 0.1, manual_gene_col = NULL) {
  # Read file, keeping first column as gene symbol if it's unnamed
  deg <- read_csv(de_file, show_col_types = FALSE)
  cn <- colnames(deg)
  # Rename first column as 'GeneSymbol' if not already
  if (is.null(manual_gene_col)) {
    # If the first column is unnamed or ...1, treat as gene column
    if (!("GeneSymbol" %in% cn)) {
      colnames(deg)[1] <- "GeneSymbol"
    }
    gene_col <- "GeneSymbol"
  } else {
    gene_col <- manual_gene_col
    if (!(gene_col %in% colnames(deg))) stop("Gene column not found: ", gene_col)
  }
  # Standardize column names for logic below (lowercase + replace . with _)
  colnames(deg) <- tolower(gsub("[.]", "_", colnames(deg)))
  # Fix gene column (use updated name)
  gene_col <- tolower(gsub("[.]", "_", gene_col))
  # Find logFC and adjusted p-value column (robustly)
  logfc_col <- grep("^logfc$", colnames(deg), value = TRUE)
  padj_col <- grep("adj", colnames(deg), value = TRUE)
  if (length(logfc_col) == 0) stop("No logFC column detected in file: ", de_file)
  if (length(padj_col) == 0) stop("No adjusted p-value column found in file: ", de_file)
  padj_col <- padj_col[1] # Use first adjusted p-value column
  
  # Filtering
  up <- deg %>%
    filter(.data[[logfc_col]] >= logfc, .data[[padj_col]] < padj) %>%
    pull(.data[[gene_col]]) %>%
    unique()
  down <- deg %>%
    filter(.data[[logfc_col]] <= -logfc, .data[[padj_col]] < padj) %>%
    pull(.data[[gene_col]]) %>%
    unique()
  return(list(UP = up, DOWN = down))
}

# -------------------------------------------------------------------------
# STEP 3: Run ORA (gost) and Save Results for Each Comparison and Direction
# -------------------------------------------------------------------------
for (i in seq_along(comparisons)) {
  cmp <- comparisons[i]
  de_file <- dge_files[i]
  
  degs <- get_deg_lists(de_file)
  
  # --- UP genes
  if (length(degs$UP) > 0) {
    gost_up <- gost(degs$UP, organism = "mmusculus", correction_method = "fdr")
    up_table <- gost_up$result
    out_file <- file.path(up_dir, paste0("gost_", cmp, "_UP.csv"))
    write_csv(up_table, out_file)
    message("Saved ORA (UP) for ", cmp, " to: ", out_file)
  } else {
    message("No UP genes found for ", cmp)
  }
  
  # --- DOWN genes
  if (length(degs$DOWN) > 0) {
    gost_down <- gost(degs$DOWN, organism = "mmusculus", correction_method = "fdr")
    down_table <- gost_down$result
    out_file <- file.path(down_dir, paste0("gost_", cmp, "_DOWN.csv"))
    write_csv(down_table, out_file)
    message("Saved ORA (DOWN) for ", cmp, " to: ", out_file)
  } else {
    message("No DOWN genes found for ", cmp)
  }
}

# -------------------------------------------------------------------------
# STEP 4: Dotplot Function for Top Pathways (edit n_terms as needed)
# -------------------------------------------------------------------------
plot_gost_dotplot <- function(gost_result, n_terms = 15, plot_title = "") {
  plot_df <- gost_result %>%
    arrange(p_value) %>%
    head(n_terms) %>%
    mutate(term_name = fct_reorder(term_name, -log10(p_value)))
  
  ggplot(plot_df, aes(x = -log10(p_value), y = term_name, size = intersection_size, color = -log10(p_value))) +
    geom_point() +
    labs(
      title = plot_title,
      x = expression(-log[10]~"(FDR-adjusted P-value)"),
      y = "GO/Pathway Term",
      size = "Gene count",
      color = "-log10(FDR)"
    ) +
    scale_color_viridis_c(option = "A", direction = -1) +
    theme_bw(base_size = 14)
}

# -------------------------------------------------------------------------
# STEP 5: Example - Plot for KO_24M_vs_WT_24M (UP/DOWN), or loop as needed
# -------------------------------------------------------------------------
## PINK1 vs T2D DOWN
# Func
make_sankey_plot <- function(my_comparison_colors, my_pathway_colors, my_ora_dir, my_pathways) {
  # --- Libraries ---
  library(readr)
  library(dplyr)
  library(ggalluvial)
  library(ggplot2)
  library(cowplot)
  library(scales)
  
  # --- Helper: Robust pathway name cleaning ---
  nukefix <- function(x) {
    x <- gsub("β", "b", x)
    x <- trimws(x)
    x <- iconv(x, to = "ASCII//TRANSLIT")
    x <- tolower(x)
    x <- gsub("_", " ", x)    # treat underscores as spaces
    x <- gsub("\\s+", " ", x)
    x
  }
  
  # --- Find CSVs and print diagnostics ---
  my_comp_names <- list.files(my_ora_dir, pattern = "\\.csv$", full.names = FALSE)
  cat("Found files:\n"); print(my_comp_names)
  cat("Pathways (requested):\n"); print(my_pathways)
  
  nukefixed_pathways <- nukefix(my_pathways)
  all_links <- list()
  
  for (fname in my_comp_names) {
    comp <- tolower(gsub("\\.csv$", "", fname))
    fpath <- file.path(my_ora_dir, fname)
    if (!file.exists(fpath)) next
    df <- read_csv(fpath, show_col_types = FALSE)
    if (!"term_name" %in% colnames(df)) next
    
    # --- Only keep significant pathways (p_value < 0.05) ---
    if (!"p_value" %in% tolower(colnames(df))) stop("No p_value column in file: ", fname)
    pval_col <- colnames(df)[tolower(colnames(df)) == "p_value"][1]
    df <- df %>% filter(!is.na(.data[[pval_col]]), .data[[pval_col]] < 0.05, !is.na(term_name))
    df$Pathway <- nukefix(df$term_name)
    
    # --- DEBUG: Show what nukefixed pathways are present ---
    if (nrow(df) > 0) cat("In", fname, "found pathways:\n", paste(unique(df$Pathway), collapse=", "), "\n")
    
    for (i in seq_along(my_pathways)) {
      pw <- nukefixed_pathways[i]
      pretty_name <- my_pathways[i]
      match_row <- df[df$Pathway == pw,]
      if (nrow(match_row) == 0) next
      all_links[[length(all_links)+1]] <- data.frame(
        Comparison = comp,
        Pathway = pretty_name,    # preserve pretty/original
        Value = -log10(match_row[[pval_col]][1]) + 1,
        stringsAsFactors = FALSE
      )
    }
  }
  
  sankey_df <- dplyr::bind_rows(all_links)
  # Robustly match on nukefixed comparison and pathway names for factors
  sankey_df$Comparison <- factor(nukefix(sankey_df$Comparison), levels = nukefix(names(my_comparison_colors)))
  sankey_df$Pathway <- factor(nukefix(sankey_df$Pathway), levels = nukefix(my_pathways), labels = my_pathways)
  
  # --- Stop gracefully if nothing matched ---
  if (nrow(sankey_df) == 0) {
    stop("No matching pathways were found after filtering. Check pathway names and p_value thresholds.\n")
  }
  
  # --- Main Sankey plot ---
  p_main <- ggplot(sankey_df, aes(axis1 = Comparison, axis2 = Pathway, y = Value)) +
    geom_alluvium(aes(fill = Comparison), width = 1/12, alpha = 0.8) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_discrete(limits = c("Comparison", "Pathway"), expand = c(.1, .1)) +
    scale_fill_manual(
      values = c(my_comparison_colors, my_pathway_colors),
      breaks = names(my_comparison_colors)
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold")
    )
  
  # --- Dot legend for FDR values ---
  fdr_values <- sankey_df$Value
  legend_vals <- pretty(fdr_values, n = 5)
  legend_vals <- sort(unique(legend_vals[legend_vals > 0]), decreasing = TRUE)
  rescale_size <- function(x, to = c(3, 10)) {
    rng <- range(x, na.rm = TRUE)
    scales::rescale(x, to = to, from = rng)
  }
  dot_sizes <- rescale_size(legend_vals, to = c(3, 10))
  legend_df <- data.frame(
    x = 1,
    y = seq(50, 38, length.out = length(dot_sizes)),
    y_label = seq(50, 38, length.out = length(dot_sizes)),
    size = dot_sizes,
    label = paste0("-log10(p) + 1 = ", legend_vals)
  )
  
  p_legend <- ggplot(legend_df) +
    geom_point(aes(x = x, y = y, size = size), shape = 21, fill = "steelblue", color = "black", stroke = 0.25) +
    geom_text(aes(x = x + 0.4, y = y_label, label = label), hjust = 0, vjust = 0.5, size = 4) +
    theme_void() +
    coord_cartesian(clip = "off") +
    scale_size_identity() +
    scale_x_continuous(limits = c(0.9, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme(plot.margin = margin(5, 20, 5, 5))
  
  # --- Combine ---
  cowplot::plot_grid(
    p_main,
    p_legend,
    rel_widths = c(4.2, 1.2),
    nrow = 1,
    axis = "none",
    align = "none"
  )
}


## Load prereqs
# --- Colors ---
my_comparison_colors <- c(
  "gost ko 4m vs wt 4m down"  = "#274472",
  "gost ko 24m vs wt 24m down" = "#8C3333"
)

my_pathway_colors <- c(
  # --- Mitophagy (Brown) ---
  "mitophagy" = "#8D6748",
  "pink1-prkn mediated mitophagy" = "#8D6748",
  "autophagy of mitochondrion" = "#8D6748",
  "regulation of mitophagy" = "#8D6748",
  "negative regulation of mitophagy" = "#8D6748",
  
  # --- Apoptosis (Red) ---
  "intrinsic apoptotic signaling pathway in response to hydrogen peroxide" = "#C0392B",
  "programmed cell death in response to reactive oxygen species" = "#C0392B",
  "release of cytochrome c from mitochondria" = "#C0392B",
  "negative regulation of hydrogen peroxide-mediated programmed cell death" = "#C0392B",
  "cell death signalling via nrage, nrif and nade" = "#C0392B",
  
  # --- Dynamics/fission/cristae (Green) ---
  "mitochondrial fission" = "#27AE60",
  "negative regulation of mitochondrial fission" = "#27AE60",
  "positive regulation of mitochondrial fission" = "#27AE60",
  "cristae formation" = "#27AE60",
  "positive regulation of cristae formation" = "#27AE60",
  
  # --- OXPHOS/respiration (Orange) ---
  "oxidative phosphorylation" = "#D35400",
  "mitochondrial electron transport, nadh to ubiquinone" = "#D35400",
  "regulation of mitochondrial electron transport, nadh to ubiquinone" = "#D35400",
  "aerobic electron transport chain" = "#D35400",
  "positive regulation of cellular respiration" = "#D35400",
  
  # --- Stress/ROS & Cellular Response (Blue) ---
  "regulation of response to oxidative stress" = "#2980B9",
  "hydrogen peroxide metabolic process" = "#2980B9",
  "negative regulation of reactive oxygen species metabolic process" = "#2980B9",
  "response to ischemia" = "#2980B9",
  "regulation of hydrogen peroxide metabolic process" = "#2980B9",
  "response to stimulus" = "#2980B9",
  "response to stress" = "#2980B9",  # overlap #3
  "response to hypoxia" = "#2980B9",
  "positive regulation of response to stress" = "#2980B9",
  "negative regulation of cellular process" = "#2980B9",
  "negative regulation of neuron apoptotic process" = "#2980B9",
  
  # --- Development, Differentiation, Morphogenesis (Purple) ---
  "cell differentiation" = "#7B1FA2",  # overlap #1
  "developmental process" = "#7B1FA2",
  "tissue development" = "#7B1FA2",
  "neuron differentiation" = "#7B1FA2",
  "neurogenesis" = "#7B1FA2",
  "axonogenesis" = "#7B1FA2",
  "cell morphogenesis involved in neuron differentiation" = "#7B1FA2",
  "generation of neurons" = "#7B1FA2",
  "anatomical structure development" = "#7B1FA2",
  "multicellular organism development" = "#7B1FA2",
  "positive regulation of cell differentiation" = "#7B1FA2",
  "regulation of developmental process" = "#7B1FA2",
  "epithelial cell differentiation" = "#7B1FA2",
  "anatomical structure morphogenesis" = "#7B1FA2",
  
  # --- Localization, Transport, Signaling (Brown) ---
  "regulation of localization" = "#A0522D",  # overlap #2
  "localization" = "#A0522D",
  "regulation of transport" = "#A0522D",
  "transport" = "#A0522D",
  "regulation of cellular localization" = "#A0522D",
  "positive regulation of exocytosis" = "#A0522D",
  "regulation of amine transport" = "#A0522D",
  "regulation of glutamate secretion" = "#A0522D",
  "cell-cell signaling" = "#A0522D",
  "platelet-derived growth factor receptor signaling pathway" = "#A0522D",
  
  # --- ECM & Structural (Teal) ---
  "extracellular matrix" = "#008080",  # overlap #4
  
  # --- Cancer/Oncogenic Signaling (Pink) ---
  "PI3K-Akt signaling pathway" = "#E91E63"  # overlap #5
)


# --- Input ---
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/ora/DOWN"
my_comp_names <- c(
  "gost_KO_4M_vs_WT_4M_DOWN.csv",
  "gost_KO_24M_vs_WT_24M_DOWN.csv"
)

my_pathways <- names(my_pathway_colors)

# RUN plot
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/Dataset_GSE217775/ora/DOWN",
  my_pathways = my_pathways
)

####
# --- Load libraries ---
library(dplyr)
library(broom) # for tidy model output

# --- Read your dataset ---
mito_algo <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\PINK1_discovery\result_with_all_features.csv)")

# Ensure binary outcome
mito_algo$is_T2D <- as.numeric(mito_algo$is_T2D)

# --- Define mitophagy genes of interest ---
mitophagy_genes <- c(
  "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "PHB2", "NLRX1",
  "BCL2L13", "FKBP8", "MUL1", "SAMM50", "MTX1", "PGAM5",
  "CSNK2A1", "CSNK2A2", "SRC", "ULK1", "OPA1", "MFN2", "TBK1",
  "TAX1BP1", "SQSTM1", "CALCOCO2", "OPTN", "NBR1", "MAP1LC3A",
  "MAP1LC3B", "MAP1LC3C", "BECN1", "ATG9A", "USP8", "USP15",
  "USP30", "VCP", "RAB5A", "RAB7A", "MON1A", "CCZ1", "FIS1",
  "BCL2L1"
)

# --- Filter to only genes in dataset ---
mitophagy_genes <- intersect(mitophagy_genes, colnames(mito_algo))

# --- Function to fit logistic regression with interaction ---
fit_interaction <- function(gene) {
  # Skip PINK1 itself
  if (gene == "PINK1") return(NULL)
  
  formula <- as.formula(
    paste0("is_T2D ~ PINK1 * ", gene)
  )
  
  model <- glm(formula, data = mito_algo, family = binomial)
  
  # Extract interaction term coefficient
  tidy(model) %>%
    filter(term == paste0("PINK1:", gene)) %>%
    mutate(Gene = gene)
}

# --- Run for all mitophagy genes ---
results <- lapply(mitophagy_genes, fit_interaction) %>%
  bind_rows()

# --- Rank by p-value ---
results <- results %>%
  arrange(p.value)

# --- View top interactions ---
head(results, 10)

library(scales)
library(ggbreak)
#install.packages("ggbreak")
library(ggrepel)

vol_df <- results %>%
  mutate(
    FDR = p.adjust(p.value, "BH"),
    neglog10p = -log10(p.value),
    Sig = case_when(
      FDR < 0.05 ~ "FDR<0.05",
      p.value < 0.05 ~ "p<0.05",
      TRUE ~ "NS"
    )
  )

ggplot(vol_df, aes(x = estimate, y = neglog10p, color = Sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(
    data = subset(vol_df, Sig != "NS"),
    aes(label = Gene), size = 3, max.overlaps = 20
  ) +
  scale_color_manual(values = c("FDR<0.05" = "#D32F2F", "p<0.05" = "#F57C00", "NS" = "grey65")) +
  scale_x_break(c(-45, -4)) +               # creates a visible gap for the outlier
  labs(
    x = "Interaction effect (log-odds)",
    y = expression(-log[10](p)),
    title = "Volcano of PINK1 × Gene interactions\n(with axis break for MAP1LC3C outlier)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")


# By sex
# ============================================================
# Volcano plots by sex with labels for ALL significant genes
#   - Red  : FDR < 0.05
#   - Yellow: nominal p < 0.05
# Uses ggbreak::scale_x_break() to hide extreme negative outliers (optional)
# Assumes you already have:
#   - mito_algo (with column "Sex" in {"M","F"} plus columns needed by fit_interaction)
#   - mitophagy_genes (character vector)
#   - fit_interaction(gene, df) -> data.frame row with columns: Gene, estimate, p.value
# ============================================================

# ---- 0) Packages & helpers ----
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

ensure_pkg("dplyr")
ensure_pkg("ggplot2")
ensure_pkg("ggrepel")
ensure_pkg("patchwork")
ensure_pkg("ggbreak")   # for axis gap
# (Optional) reproducible label placement across runs:
set.seed(1234)

# ---- 1) Compute interaction stats per sex ----
run_interactions_df <- function(sex_value) {
  df_sex <- mito_algo %>% dplyr::filter(Sex == sex_value)
  
  results <- lapply(mitophagy_genes, fit_interaction, df = df_sex) %>%
    dplyr::bind_rows() %>%
    # Safety: drop incomplete rows to avoid plotting warnings
    dplyr::filter(!is.na(estimate), !is.na(p.value), !is.na(Gene)) %>%
    dplyr::mutate(
      FDR        = p.adjust(p.value, "BH"),
      neglog10p  = -log10(p.value),
      Sig        = dplyr::case_when(
        FDR < 0.05     ~ "FDR<0.05",
        p.value < 0.05 ~ "p<0.05",
        TRUE           ~ "NS"
      )
    ) %>%
    dplyr::arrange(FDR, p.value)
  
  results
}

# ---- 2) Reusable volcano plot function (labels ALL significant genes) ----
volcano_plot <- function(results,
                         title,
                         xlim = c(-6, 6),
                         ylim = NULL,
                         use_x_gap = FALSE,
                         x_gap = c(-45, -6),          # hidden band on x-axis (for big negatives)
                         point_size = 2,
                         text_size = 3,
                         repel_force = 0.6,
                         repel_max_overlaps = Inf) {  # allow all labels to try to place
  
  # Label ALL significant genes (red + yellow)
  lab_dat <- results %>%
    dplyr::filter(Sig != "NS", !is.na(Gene), !is.na(estimate), !is.na(neglog10p))
  
  p <- ggplot(results, aes(x = estimate, y = neglog10p, color = Sig)) +
    # Reference lines
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
    geom_vline(xintercept = 0,            linetype = 2, linewidth = 0.4) +
    # Points
    geom_point(size = point_size, alpha = 0.9, na.rm = TRUE) +
    # Labels: ALL significant
    ggrepel::geom_text_repel(
      data = lab_dat,
      aes(label = Gene),
      size = text_size,
      box.padding = 0.25,
      point.padding = 0.2,
      max.overlaps = repel_max_overlaps,
      force = repel_force,
      segment.size = 0.2,
      min.segment.length = 0,
      na.rm = TRUE
    ) +
    # Colors
    scale_color_manual(values = c("FDR<0.05" = "#D32F2F",
                                  "p<0.05"   = "#F57C00",
                                  "NS"       = "grey65")) +
    # Titles/labels
    labs(
      x = "Interaction effect (log-odds)",
      y = expression(-log[10](p)),
      title = title,
      color = NULL
    ) +
    # Theme + limits
    theme_minimal(base_size = 12) +
    theme(
      legend.position  = "top",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
  
  # Optional x-axis gap for extreme negative estimates
  if (use_x_gap) p <- p + ggbreak::scale_x_break(breaks = x_gap)
  
  p
}

# ---- 3) Run for each sex ----
res_M <- run_interactions_df("M")
res_F <- run_interactions_df("F")

# Nice rounded y-limits
ylim_M <- c(0, ceiling(max(res_M$neglog10p, na.rm = TRUE)))
ylim_F <- c(0, ceiling(max(res_F$neglog10p, na.rm = TRUE)))

# ---- 4) Build plots (label ALL significant points) ----
# Males: allow an x-gap to hide very negative outliers
plot_male <- volcano_plot(
  results   = res_M,
  title     = "PINK1 × Gene Interactions in Male Beta Cells",
  xlim      = c(-10, 5),
  ylim      = c(0, 20),
  use_x_gap = FALSE,
  #x_gap     = c(-45, -6),   # adjust if needed
  point_size = 2,
  text_size  = 3,
  repel_force = 0.7
)

# Females: tighter x-range, no gap
plot_female <- volcano_plot(
  results   = res_F,
  title     = "PINK1 × Gene Interactions in Female Beta Cells",
  xlim      = c(-50, 5),
  ylim      = c(0, 20),
  use_x_gap = TRUE,
  x_gap     = c(-47, -6),   # adjust if needed
  point_size = 2,
  text_size  = 3,
  repel_force = 0.7
)

# ---- 5) Display side by side with shared legend ----
(plot_male + plot_female) + patchwork::plot_layout(guides = "collect")

# ============================================================
# Tips:
# - If labels still collide, try:
#     * Increase repel_force (e.g., 1.0–1.5)
#     * Increase text_size slightly
#     * Widen xlim or raise ylim
#     * For very dense sets, you can set repel_max_overlaps to a large number (e.g., 500)
# - If you want to avoid ggbreak, set use_x_gap = FALSE and optionally filter outliers before plotting:
#     res_M_trim <- res_M %>% dplyr::filter(estimate >= -6)
#     plot_male <- volcano_plot(res_M_trim, "Male …", xlim = c(-6, 6), ylim = ylim_M, use_x_gap = FALSE)
# ============================================================


# ---------------------------
# DHT Analysis - β-cell mapping
# ---------------------------
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)

# --- 1) Load datasets ---
dht_dataset <- readRDS(
  r"(C:\Users\mqadir\Box\FMJ lab\2. Papers\1. Published\BARKO paper\Data_uploads\GEO Upload\R\pancreas.integrated.rds)"
)
beta_cells <- qread(
  r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)"
)

# --- 2) SCT transform for mapping ---
DefaultAssay(dht_dataset) <- "RNA"
dht_dataset <- SCTransform(dht_dataset, verbose = TRUE)
DefaultAssay(beta_cells) <- "RNA"
beta_cells <- SCTransform(beta_cells, verbose = TRUE)

# --- 3) Subset male β INS-hi/low ---
dht_male_beta <- subset(
  dht_dataset,
  subset = sex == "Male" & celltype %in% c("Beta INS-hi", "Beta INS-low")
)

# --- 4) Run UMAP on reference with model ---
beta_cells <- RunUMAP(
  beta_cells, reduction = "pca", dims = 1:8,
  return.model = TRUE
)
male_beta_cells <- subset(beta_cells, subset = Sex == "M")

DimPlot(
  beta_cells,
  reduction = "umap",
  #group.by  = "dataset",
  label     = TRUE,
  repel     = TRUE
) + theme_minimal()

# --- 5) Find anchors & MapQuery ---
anchors <- FindTransferAnchors(
  reference = beta_cells,
  query = dht_male_beta,
  normalization.method = "SCT",
  reference.reduction = "pca"
)

dht_male_beta <- MapQuery(
  anchorset = anchors,
  query = dht_male_beta,
  reference = beta_cells,
  refdata = list(
    collapsed_cluster = "collapsed_cluster",
    beta_cluster = "beta_cluster"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# --- 6) Copy ref.umap to a standard 'umap' slot in query ---
dht_male_beta[["umap"]] <- CreateDimReducObject(
  embeddings = Embeddings(dht_male_beta[["ref.umap"]]),
  key = "UMAP_",
  assay = DefaultAssay(dht_male_beta)
)

# --- 7) Merge reference + query ---
combined <- merge(
  x = beta_cells, y = dht_male_beta,
  add.cell.ids = c("REF", "DHT"),
  merge.data = TRUE
)

# --- 8) Build a unified UMAP that survives merge ---
ref_umap <- Embeddings(beta_cells[["umap"]])
rownames(ref_umap) <- paste0("REF_", rownames(ref_umap))
qry_umap <- Embeddings(dht_male_beta[["umap"]])
rownames(qry_umap) <- paste0("DHT_", rownames(qry_umap))

umap_mat <- rbind(ref_umap, qry_umap)
umap_mat <- umap_mat[Cells(combined), , drop = FALSE]
combined[["umap"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key = "UMAP_",
  assay = DefaultAssay(combined)
)

# --- 9) Fill NA collapsed_cluster with predicted values ---
combined$collapsed_cluster <- ifelse(
  is.na(combined$collapsed_cluster),
  combined$predicted.collapsed_cluster,
  combined$collapsed_cluster
)

# Replace NAs with "M"
combined$Sex <- ifelse(is.na(combined$Sex), "M", combined$Sex)

# Subset to only male cells
#combined <- subset(combined, subset = Sex == "M")

# Check
unique(combined$Sex)
table(combined$treatment, useNA = "ifany")

# --- 10) Plot unified UMAP ---
DimPlot(combined, reduction = "umap", group.by = "collapsed_cluster", label = FALSE) +
  scale_color_manual(values = c("A" = "firebrick", "B" = "darkorange",
                                "C" = "dodgerblue", "D" = "orchid"))

# Replace NAs in 'treatment' with 'Untreated'
combined$treatment <- ifelse(
  is.na(combined$treatment),
  "Untreated",
  combined$treatment
)
table(combined$treatment, useNA = "ifany")

# Now plot
DimPlot(
  combined,
  reduction = "umap",
  group.by = "treatment",
  label = FALSE,
  repel = TRUE
) +
  scale_color_manual(values = c("Untreated" = "grey65", "DHT[10nM]" = "firebrick4", "EtOH" = "dodgerblue4"))

# CHeck data
head(combined@meta.data)

DimPlot(
  combined,
  reduction = "umap",
  group.by = "beta_cluster",
  label = FALSE,
  repel = TRUE
)

# --- Dependencies ---
library(Seurat)
library(dplyr)
library(ggplot2)

# --- 1) Map collapsed_cluster -> beta_cluster labels ---
# A -> β1 (firebrick), B -> β3 (darkorange), C -> β4 (dodgerblue), D -> β2 (orchid)
map_letters_to_beta <- c(
  "A" = "β1",
  "B" = "β3",
  "C" = "β4",
  "D" = "β2"
)

# --- 2) Create/overwrite beta_cluster using the mapping (preserve NAs) ---
combined@meta.data$beta_cluster <- dplyr::recode(
  combined@meta.data$collapsed_cluster,
  !!!map_letters_to_beta,
  .default = NA_character_
)

# --- 3) Make beta_cluster an ordered factor for consistent legends & plotting ---
combined@meta.data$beta_cluster <- factor(
  combined@meta.data$beta_cluster,
  levels = c("β1", "β2", "β3", "β4")
)

# Quick sanity check
print(table(combined@meta.data$collapsed_cluster, useNA = "ifany"))
print(table(combined@meta.data$beta_cluster,      useNA = "ifany"))
print(table(combined@meta.data$treatment,      useNA = "ifany"))

# --- 4) Colors to match your figure ---
beta_cols <- c(
  "β1" = "#B22222",  # firebrick
  "β2" = "#DA70D6",  # orchid
  "β3" = "#FF8C00",  # darkorange
  "β4" = "#1E90FF",  # dodgerblue
  "NA" = "#A6A6A6"   # grey65 (for any NA if shown)
)

# --- 5) Plot by the new beta_cluster labels ---
DimPlot(
  combined,
  reduction = "umap",
  group.by  = "beta_cluster",
  label     = TRUE,
  repel     = TRUE
) +
  scale_color_manual(values = beta_cols[levels(combined@meta.data$beta_cluster)]) +
  labs(title = "UMAP colored by β-clusters (mapped from collapsed_cluster)") +
  theme_minimal()

library(dplyr)

# First, make sure treatment column has no NAs (label them as 'Untreated' or similar)
combined@meta.data <- combined@meta.data 

# Now count cells per cluster and treatment
cell_counts <- combined@meta.data %>%
  group_by(beta_cluster, treatment) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(beta_cluster, treatment)

cell_counts

# --- Prereqs ---
library(Seurat)
library(dplyr)
library(ggplot2)

# ---------- 0) Housekeeping ----------
# Replace missing treatment with "Untreated"
#combined$treatment <- ifelse(is.na(combined$treatment), "Untreated", combined$treatment)

# Ensure beta_cluster exists and is ordered
cluster_levels <- c("β1", "β2", "β3", "β4")
if (!"beta_cluster" %in% colnames(combined@meta.data)) {
  if ("collapsed_cluster" %in% colnames(combined@meta.data)) {
    combined$beta_cluster <- factor(combined$collapsed_cluster, levels = cluster_levels)
  } else {
    stop("No 'beta_cluster' or 'collapsed_cluster' column found in metadata.")
  }
} else {
  combined$beta_cluster <- factor(combined$beta_cluster, levels = cluster_levels)
}

# --- 1) Mitophagy gene set ---
mitophagy_genes <- c(
  "PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","PHB2","NLRX1","BCL2L13","FKBP8",
  "MUL1","SAMM50","MTX1","PGAM5","CSNK2A1","CSNK2A2","SRC","ULK1","OPA1",
  "MFN2","TBK1","TAX1BP1","SQSTM1","CALCOCO2","OPTN","NBR1",
  "MAP1LC3A","MAP1LC3B","MAP1LC3C","BECN1","ATG9A","USP8","USP15","USP30",
  "VCP","RAB5A","RAB7A","MON1A","CCZ1","FIS1","BCL2L1"
)

# Keep only genes present in object
genes_in <- intersect(mitophagy_genes, rownames(combined))
if (length(genes_in) == 0) {
  stop("None of the mitophagy genes are present in the Seurat object rownames.")
}

# Red palette (white → red)
red_pal <- c("white", "#fee0d2", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")

# =======================================================================================
# A) DOTPLOT: Cluster × Treatment (EtOH vs DHT[10nM] vs Untreated) across β1–β4
# =======================================================================================

# Order treatments
treat_levels <- c("EtOH", "DHT[10nM]", "Untreated")
combined$treatment <- factor(combined$treatment, levels = treat_levels)
unique(combined$treatment)
head(combined@meta.data)

# Build combined label: "βx • Treatment"
combined$cluster_treat <- with(
  combined@meta.data,
  paste0(beta_cluster, " • ", treatment)
)
unique(combined$cluster_treat)

# Order factor levels: β1_EtOH, β1_DHT, β1_Untreated, β2_EtOH, ...
combined$cluster_treat <- factor(
  combined$cluster_treat,
  levels = as.vector(outer(cluster_levels, treat_levels, 
                           function(cl, tr) paste0(cl, " • ", tr)))
)

# Create DotPlot
p_ct <- DotPlot(
  combined,
  features = genes_in,
  group.by = "cluster_treat",
  assay = DefaultAssay(combined), # Change to "SCT" if needed
  cols = red_pal
) +
  scale_colour_gradientn(colours = red_pal) +  # Average expression color
  #scale_size(range = c(1.2, 8)) +              # Dot size = % cells expressing
  labs(
    x = NULL, y = NULL,
    title = "Mitophagy genes across β-clusters by treatment",
    subtitle = "Dot size = % cells expressing; color = scaled average expression"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.y = element_blank()
  )

# Print plot
p_ct





# --- Prereqs ---
library(Seurat)
library(dplyr)
library(ggplot2)

# ---------------------------
# Setup
# ---------------------------
cluster_levels <- c("β1","β2","β3","β4")
treat_levels   <- c("EtOH", "DHT[10nM]", "Untreated")

# Ensure factors
combined$beta_cluster <- factor(combined$beta_cluster, levels = cluster_levels)
combined$treatment    <- factor(combined$treatment,    levels = treat_levels)

# Group label used for the base (detailed) DotPlot
combined$cluster_treat <- factor(
  paste0(combined$beta_cluster, " • ", combined$treatment),
  levels = as.vector(outer(cluster_levels, treat_levels, function(cl,tr) paste0(cl, " • ", tr)))
)

# Gene list
mitophagy_genes <- c(
  "PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","PHB2","NLRX1","BCL2L13","FKBP8",
  "MUL1","SAMM50","MTX1","PGAM5","CSNK2A1","CSNK2A2","SRC","ULK1","OPA1",
  "MFN2","TBK1","TAX1BP1","SQSTM1","CALCOCO2","OPTN","NBR1",
  "MAP1LC3A","MAP1LC3B","MAP1LC3C","BECN1","ATG9A","USP8","USP15","USP30",
  "VCP","RAB5A","RAB7A","MON1A","CCZ1","FIS1","BCL2L1"
)
features <- intersect(mitophagy_genes, rownames(combined))
stopifnot(length(features) > 0)

# Palette
red_pal <- c("white", "#fee0d2", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")

# ---------------------------
# 1) Build DotPlot data for the detailed groups (βx • Treatment)
# ---------------------------
p_base <- DotPlot(
  combined,
  features = features,
  group.by = "cluster_treat",
  assay = DefaultAssay(combined)
)
df_base <- p_base$data
# Ensure order of groups matches desired factor order
df_base$id <- factor(df_base$id, levels = levels(combined$cluster_treat))

# ---------------------------
# 2) Build summary rows for All EtOH and All DHT[10nM]
#    (aggregate across ALL beta clusters)
# ---------------------------
treated_obj <- subset(combined, subset = treatment %in% c("EtOH", "DHT[10nM]"))
p_overall <- DotPlot(
  treated_obj,
  features = features,
  group.by = "treatment",
  assay = DefaultAssay(combined)
)
df_overall <- p_overall$data %>%
  mutate(id = ifelse(id == "EtOH", "All EtOH", "All DHT[10nM]"))

# ---------------------------
# 3) Stack summary rows on top of detailed rows and fix factor orders
# ---------------------------
final_levels <- c("All EtOH", "All DHT[10nM]", levels(combined$cluster_treat))
df_all <- bind_rows(df_overall, df_base) %>%
  mutate(
    id = factor(id, levels = final_levels),
    features.plot = factor(features.plot, levels = features)  # keep gene order as provided
  )

# ---------------------------
# 4) Plot (replicating Seurat DotPlot style)
# ---------------------------
p <- ggplot(df_all, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_colour_gradientn(colours = red_pal, na.value = "grey90") +
  #scale_size(range = c(1.2, 8), limits = c(0, 100)) +
  labs(
    x = NULL, y = NULL,
    title = "Mitophagy genes across β-clusters by treatment",
    subtitle = "Top rows = All EtOH / All DHT[10nM] across β1–β4; others = individual β-cluster × treatment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.y = element_blank()
  )

p

print(table(combined@meta.data$cluster_treat, useNA = "ifany"))
print(table(combined@meta.data$treatment, useNA = "ifany"))


library(Seurat)
library(dplyr)
library(ggplot2)

# --- Expanded pathway gene lists (human) ---
# OXPHOS – All complexes I–V from KEGG hsa00190
oxphos_genes <- c(
  # Complex I
  "NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
  "NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11",
  "NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8",
  # Complex II
  "SDHA","SDHB","SDHC","SDHD",
  # Complex III
  "UQCRC1","UQCRC2","UQCRFS1","UQCRB","UQCRQ","UQCR10","UQCR11","UQCRH",
  # Complex IV
  "COX4I1","COX4I2","COX5A","COX5B","COX6A1","COX6A2","COX6B1","COX6B2","COX6C","COX7A1","COX7A2","COX7A2L","COX7B","COX7B2","COX7C","COX8A","COX8C",
  "MT-CO1","MT-CO2","MT-CO3",
  # Complex V
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5ME","ATP5MF","ATP5MG","ATP5PB","ATP5PD","ATP5PF","ATP5PO",
  "MT-ATP6","MT-ATP8"
)

# TCA cycle – KEGG hsa00020
tca_genes <- c(
  "CS","ACO2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","OGDH","DLST","DLD",
  "SUCLA1","SUCLA2","SUCLG1","SDHA","SDHB","SDHC","SDHD",
  "FH","MDH1","MDH2"
)

# Glycolysis – KEGG hsa00010
glyco_genes <- c(
  "HK1","HK2","HK3","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOB","ALDOC",
  "TPI1","GAPDH","PGK1","PGK2","PGAM1","PGAM2","PGAM4",
  "ENO1","ENO2","ENO3","PKLR","PKM","LDHA","LDHB","LDHC"
)

# --- Keep only genes present in dataset ---
DefaultAssay(combined) <- "RNA"
all_genes <- rownames(combined)
oxphos_genes <- intersect(oxphos_genes, all_genes)
tca_genes    <- intersect(tca_genes, all_genes)
glyco_genes  <- intersect(glyco_genes, all_genes)

cat(length(oxphos_genes), "OXPHOS genes; ",
    length(tca_genes), "TCA genes; ",
    length(glyco_genes), "Glycolysis genes\n")

# --- Add module scores ---
combined <- AddModuleScore(combined, list(oxphos_genes), name = "OXPHOS")
combined <- AddModuleScore(combined, list(tca_genes),    name = "TCA")
combined <- AddModuleScore(combined, list(glyco_genes),  name = "Glyco")

# --- Build data frame ---
df <- combined@meta.data %>%
  filter(treatment %in% c("EtOH", "DHT[10nM]")) %>%
  mutate(
    oxphos = OXPHOS1,
    tca    = TCA1,
    glyco  = Glyco1,
    oxphos_glyco = log2(oxphos / glyco),
    oxphos_tca   = log2(oxphos / tca),
    oxphos_tca_glyco = log2(oxphos / (tca + glyco))
  )

# --- Mean ± SEM ---
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

plot_df <- df %>%
  group_by(treatment) %>%
  summarise(
    OXPHOS_Glycolysis_mean = mean(oxphos_glyco, na.rm = TRUE),
    OXPHOS_Glycolysis_sem  = sem(oxphos_glyco),
    OXPHOS_TCA_mean        = mean(oxphos_tca, na.rm = TRUE),
    OXPHOS_TCA_sem         = sem(oxphos_tca),
    OXPHOS_TCA_Glyco_mean  = mean(oxphos_tca_glyco, na.rm = TRUE),
    OXPHOS_TCA_Glyco_sem   = sem(oxphos_tca_glyco)
  )

# --- Reshape for plotting ---
plot_long <- plot_df %>%
  tidyr::pivot_longer(
    cols = -treatment,
    names_to = c("metric", ".value"),
    names_pattern = "(.*)_(mean|sem)"
  )

# --- Plot ---
ggplot(plot_long, aes(x = treatment, y = mean, group = metric, color = treatment)) +
  geom_point(size = 3) +
  geom_line(aes(group = treatment), size = 0.8) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    y = expression(Log[2]~"(Module Score Ratio)"),
    x = NULL,
    title = "Pathway Ratios: OXPHOS vs TCA/Glycolysis"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

library(ggpubr)  # for stat_compare_means

# --- Stats per metric ---
# --- Helper function ---
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# --- Aggregate to treatment means ± SEM ---
plot_df <- df %>%
  group_by(treatment) %>%
  summarise(
    OXPHOS_Glycolysis_mean = mean(oxphos_glyco, na.rm = TRUE),
    OXPHOS_Glycolysis_sem  = sem(oxphos_glyco),
    OXPHOS_TCA_mean        = mean(oxphos_tca, na.rm = TRUE),
    OXPHOS_TCA_sem         = sem(oxphos_tca),
    OXPHOS_TCA_Glyco_mean  = mean(oxphos_tca_glyco, na.rm = TRUE),
    OXPHOS_TCA_Glyco_sem   = sem(oxphos_tca_glyco),
    .groups = "drop"
  )

print(plot_df)

# --- Get p-values on aggregated means ---
# Since there are only 2 groups (EtOH, DHT[10nM]), we'll just compare the two means directly
# Note: This is a simple t-test on the aggregated means (n=2), which is descriptive rather than inferential.
agg_pvals <- data.frame(
  metric = c("OXPHOS_Glycolysis", "OXPHOS_TCA", "OXPHOS_TCA_Glyco"),
  p_value = c(
    t.test(oxphos_glyco ~ treatment, data = df)$p.value,
    t.test(oxphos_tca ~ treatment, data = df)$p.value,
    t.test(oxphos_tca_glyco ~ treatment, data = df)$p.value
  )
)

print(agg_pvals)

# ======================================================================
# Mitophagy Module Score: EtOH vs DHT[10nM] — Mean ± SEM (Like your ratio code)
# ======================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Define mitophagy genes (edit if you like)
mitophagy_genes <- c(
  "PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","PHB2","NLRX1","BCL2L13","FKBP8",
  "MUL1","SAMM50","MTX1","PGAM5","CSNK2A1","CSNK2A2","SRC","ULK1","OPA1",
  "MFN2","TBK1","TAX1BP1","SQSTM1","CALCOCO2","OPTN","NBR1",
  "MAP1LC3A","MAP1LC3B","MAP1LC3C","BECN1","ATG9A","USP8","USP15","USP30",
  "VCP","RAB5A","RAB7A","MON1A","CCZ1","FIS1","BCL2L1"
)

# Keep only genes present
mitophagy_genes <- intersect(mitophagy_genes, rownames(combined))
stopifnot(length(mitophagy_genes) > 0)

# 2) Add module score (creates Mitophagy1 column; Seurat appends "1")
DefaultAssay(combined) <- DefaultAssay(combined)  # leave as-is (RNA/SCT), or set explicitly
combined <- AddModuleScore(combined, list(mitophagy_genes), name = "Mitophagy")

# 3) Build per-cell dataframe like your example, then summarise
# Define small offset to avoid log(0) or log(negative)
offset <- 1e-6

df <- combined@meta.data %>%
  filter(treatment %in% c("EtOH", "DHT[10nM]")) %>%
  mutate(Mitophagy1 = as.numeric(Mitophagy1)) %>%  # ensure numeric
  transmute(
    treatment,
    mitophagy = log2(Mitophagy1 + offset)
  )


# SEM helper
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# Mean ± SEM per treatment
plot_df <- df %>%
  group_by(treatment) %>%
  summarise(
    Mitophagy_mean = mean(mitophagy, na.rm = TRUE),
    Mitophagy_sem  = sem(mitophagy),
    n = dplyr::n(),
    .groups = "drop"
  )

print(plot_df)  # so you can see the summary numbers

# 4) Plot mean ± SEM with connecting line (EtOH -> DHT)
ggplot(plot_df, aes(x = treatment, y = Mitophagy_mean, group = 1)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = Mitophagy_mean - Mitophagy_sem,
                    ymax = Mitophagy_mean + Mitophagy_sem),
                width = 0.15) +
  labs(
    title = "Mitophagy Module Score",
    subtitle = "Mean ± SEM per treatment (EtOH vs DHT[10nM])",
    x = NULL, y = "Module Score (AddModuleScore)"
  ) +
  theme_classic(base_size = 14)

# Wilcoxon test (non-parametric)
p_val <- wilcox.test(mitophagy ~ treatment, data = df)$p.value
cat("Wilcoxon test p-value:", p_val, "\n")

# ------------------------ Optional: show per-cell distribution  ------------------------
# If you also want to see the underlying spread (like your ratio jitter panels), run:
# ggplot(df, aes(x = treatment, y = mitophagy, color = treatment)) +
#   geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
#   stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
#   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, color = "black") +
#   labs(x = NULL, y = "Mitophagy Module Score") +
#   theme_classic(base_size = 14) + theme(legend.position = "none")

# ------------------------ Optional: p-value annotation (Wilcoxon) ---------------------
# If you want a p-value in the console:
# wilcox.test(mitophagy ~ treatment, data = df)
# Or add stars/text to the jitter plot using ggpubr::stat_compare_means()


library(ggplot2)
library(dplyr)

# --- 1) Calculate mitophagy module score ---
combined <- AddModuleScore(
  combined,
  features = list(mitophagy_genes),
  name = "MitophagyScore"
)

combined$Mitophagy <- combined$MitophagyScore1

# --- 2) Filter out beta2 and prepare data ---
plot_df <- combined@meta.data %>%
  filter(treatment %in% c("EtOH", "DHT[10nM]"),
         beta_cluster %in% c("β1","β3","β4")) %>%  # remove β2
  select(beta_cluster, treatment, Mitophagy) %>%
  mutate(
    beta_cluster = factor(beta_cluster, levels = c("β1","β3","β4")),
    treatment = factor(treatment, levels = c("EtOH", "DHT[10nM]"))
  )

# --- 3) Plot non-overlapping side-by-side violins ---
dodge_val <- 0.7

ggplot(plot_df, aes(x = beta_cluster, y = Mitophagy, fill = treatment)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.6, color = "black", size = 0.3,
              position = position_dodge(width = dodge_val), width = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, color = "black",
               position = position_dodge(width = dodge_val)) +
  scale_fill_manual(values = c("EtOH" = "deeppink3", "DHT[10nM]" = "steelblue")) +
  labs(
    title = "Mitophagy Module Score",
    y = "Module Score",
    x = "Beta Cell Cluster"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90", color = NA)
  )

library(dplyr)

# Assuming `combined@meta.data` has beta_cluster, treatment, and Mitophagy columns
df_for_stats <- combined@meta.data %>%
  filter(treatment %in% c("EtOH", "DHT[10nM]"),
         beta_cluster != "β2") %>%  # remove beta2 if needed
  select(beta_cluster, treatment, Mitophagy)

# Function to compute FDR and -log10(FDR)
get_FDR_vals <- function(df) {
  p_val <- wilcox.test(Mitophagy ~ treatment, data = df)$p.value
  fdr <- p.adjust(p_val, method = "BH")
  data.frame(FDR = fdr, neg_log10FDR = -log10(fdr))
}

# Apply per beta cluster
fdr_df <- df_for_stats %>%
  group_by(beta_cluster) %>%
  do(get_FDR_vals(.)) %>%
  ungroup()

print(fdr_df)
print(table(combined@meta.data$cluster_treat,      useNA = "ifany"))


# Save combined object
#qsave(combined, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined.qs)")
combined <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined.qs)")
combined
head(combined@meta.data)

## Supplemental file F6
# --- Libraries ---
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Gene list ---
inflammation_genes <- c(
  "IL1B","IL1A","IL6","TNF",
  "CXCL8","CCL2","CCL5","CXCL10","IL18",
  "TLR4","TLR2","TLR3","TLR9","MYD88","TICAM1","IRAK1","TRAF6",
  "NFKB1","RELA","CHUK","IKBKB","IKBKG",
  "NLRP3","PYCARD","CASP1",
  "CGAS","TMEM173","IRF3","IRF7",
  "HMGB1","S100A8","S100A9","ICAM1","VCAM1","PTGS2",
  "CCL3","CCL4","CCL7","CCL8","CCL20",
  "CXCL1","CXCL2","CXCL3","CXCL5","CXCL9","CXCL11",
  "IL12A","IL12B","IL17A","IL23A","IL15","IFNG","CSF2","CSF3","IL33","IL36A","IL36G",
  "TLR1","TLR5","TLR6","TLR7","TLR8","CD14","LY96","TIRAP","TICAM2","IRAK4",
  "CASP4","CASP5","NLRC4","NLRP1","AIM2","P2RX7","NOD1","NOD2"
)

# --- Factor ordering ---
combined$beta_cluster    <- factor(combined$beta_cluster,    levels = c("β1","β2","β3","β4"))
combined$diabetes_status <- factor(combined$diabetes_status, levels = c("ND","T2D"))
combined$treatment       <- factor(combined$treatment,       levels = c("Untreated","EtOH","DHT","DHT[10nM]"))
combined$Sex             <- factor(combined$Sex,             levels = c("M","F"))

# --- Subset beta cells only ---
obj_beta <- subset(combined, cells = rownames(combined@meta.data)[!is.na(combined$beta_cluster)])

# --- Keep only genes that exist ---
genes_present <- intersect(inflammation_genes, rownames(obj_beta))
if (length(genes_present) == 0) stop("None of the inflammation genes found in dataset")

# Helper to build interleaved y-levels like:
#   for disease:  β1•ND, β1•T2D, β2•ND, β2•T2D, ...
#   for treatment: β1•DHT, β1•EtOH, β2•DHT, β2•EtOH, ...
make_levels <- function(primary_levels, secondary_levels, sep = " • ") {
  as.vector(rbind(
    sapply(primary_levels, function(p) paste0(p, sep, secondary_levels[1])),
    sapply(primary_levels, function(p) paste0(p, sep, secondary_levels[2]))
  ))
}

# ============================================================
# 1) Untreated, MALES — interleave ND/T2D for each β-subtype
# ============================================================
obj_untreated_M <- subset(obj_beta, subset = treatment == "Untreated" & Sex == "M")
obj_untreated_M$ygroup <- paste(obj_untreated_M$beta_cluster, obj_untreated_M$diabetes_status, sep = " • ")

lvl_beta <- levels(combined$beta_cluster)
lvl_dis  <- levels(combined$diabetes_status) # ND, T2D
y_levels_M <- make_levels(lvl_beta, lvl_dis, sep = " • ")
obj_untreated_M$ygroup <- factor(obj_untreated_M$ygroup, levels = y_levels_M)

dpM <- DotPlot(obj_untreated_M, features = genes_present, group.by = "ygroup")
dfM <- dpM$data %>%
  mutate(
    features.plot = factor(features.plot, levels = genes_present),
    ygroup        = factor(id, levels = y_levels_M)
  )

pM <- ggplot(dfM, aes(x = features.plot, y = ygroup)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 5), name = "% cells") +
  scale_color_gradient(low = "lightgrey", high = "red", name = "Scaled expr") +
  labs(x = "Gene", y = "β-subtype • Disease",
       title = "Untreated Male β-cells",
       subtitle = "Interleaved ND/T2D within each β-subtype") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        strip.text  = element_text(face = "bold"))

# ============================================================
# 2) Untreated, FEMALES — interleave ND/T2D for each β-subtype
# ============================================================
obj_untreated_F <- subset(obj_beta, subset = treatment == "Untreated" & Sex == "F")
obj_untreated_F$ygroup <- paste(obj_untreated_F$beta_cluster, obj_untreated_F$diabetes_status, sep = " • ")

y_levels_F <- make_levels(lvl_beta, lvl_dis, sep = " • ")
obj_untreated_F$ygroup <- factor(obj_untreated_F$ygroup, levels = y_levels_F)

dpF <- DotPlot(obj_untreated_F, features = genes_present, group.by = "ygroup")
dfF <- dpF$data %>%
  mutate(
    features.plot = factor(features.plot, levels = genes_present),
    ygroup        = factor(id, levels = y_levels_F)
  )

pF <- ggplot(dfF, aes(x = features.plot, y = ygroup)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 5), name = "% cells") +
  scale_color_gradient(low = "lightgrey", high = "red", name = "Scaled expr") +
  labs(x = "Gene", y = "β-subtype • Disease",
       title = "Untreated Female β-cells",
       subtitle = "Interleaved ND/T2D within each β-subtype") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        strip.text  = element_text(face = "bold"))

# ============================================================
# 3) Males — EtOH vs DHT — interleave DHT/EtOH for each β-subtype
# ============================================================
obj_dht <- subset(obj_beta, subset = Sex == "M" & treatment %in% c("EtOH","DHT","DHT[10nM]"))
# Collapse DHT variants to "DHT" if you want them treated the same line:
obj_dht$treat_simple <- as.character(obj_dht$treatment)
obj_dht$treat_simple[obj_dht$treat_simple == "DHT[10nM]"] <- "DHT"
obj_dht$treat_simple <- factor(obj_dht$treat_simple, levels = c("DHT","EtOH"))  # DHT first, then EtOH

obj_dht$ygroup <- paste(obj_dht$beta_cluster, obj_dht$treat_simple, sep = " • ")
y_levels_DHT <- as.vector(rbind(
  paste0(lvl_beta, " • DHT"),
  paste0(lvl_beta, " • EtOH")
))
obj_dht$ygroup <- factor(obj_dht$ygroup, levels = y_levels_DHT)

dpD <- DotPlot(obj_dht, features = genes_present, group.by = "ygroup")
dfD <- dpD$data %>%
  mutate(
    features.plot = factor(features.plot, levels = genes_present),
    ygroup        = factor(id, levels = y_levels_DHT)
  )

pDHT <- ggplot(dfD, aes(x = features.plot, y = ygroup)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 5), name = "% cells") +
  scale_color_gradient(low = "lightgrey", high = "red", name = "Scaled expr") +
  labs(x = "Gene", y = "β-subtype • Treatment",
       title = "Male β-cells: DHT vs EtOH",
       subtitle = "Interleaved DHT/EtOH within each β-subtype") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        strip.text  = element_text(face = "bold"))

# --- Show all three plots ---
print(pM)/print(pF)/print(pDHT)

# --- Optional: also show per-gene ordering by pathway or alphabetical ---
# df2$features.plot <- factor(df2$features.plot, levels = sort(genes_present))
# print(ggplot(df2, aes(features.plot, beta_cluster)) + ... same layers ...)


# ============================================================
# KAN — Train on ALL β cells (male + female), predict for all,
#       write MQI_pred to `combined`, and save PDF QC figures
# ============================================================

## 0) Packages, seed, helpers --------------------------------
suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(ggplot2)
  library(torch)
})

set.seed(42); if (requireNamespace("torch", quietly = TRUE)) torch::torch_manual_seed(42)

`%||%` <- function(a,b) if (!is.null(a)) a else b
z      <- function(x) if (all(!is.finite(x))) rep(NA_real_, length(x)) else as.numeric(scale(x))
pick_first <- function(nm, opts) { o <- opts[opts %in% nm]; if (length(o)) o[1] else NA_character_ }
AddModuleScoreSafe <- function(obj, genes, name) {
  present <- genes[genes %in% rownames(obj)]
  if (length(present) == 0) { obj[[name]] <- NA_real_; return(obj) }
  obj <- AddModuleScore(obj, features = list(present), name = name, nbin = 24, ctrl = 100, seed = 42)
  obj[[name]] <- obj[[paste0(name,"1")]]; obj[[paste0(name,"1")]] <- NULL
  obj
}
scale_safe <- function(v) { s <- sd(v, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(0, length(v)) else as.numeric(scale(v)) }

## 1) IO paths -----------------------------------------------
seu_path <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined.qs)"
outdir   <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\qc_images)"
model_pt <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\fit_all_state.pt)"
scaler_qs<- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\input_scaler.qs)"
imp_qs   <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\imp_all.qs)"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## 2) Load Seurat & standardize metadata ---------------------
stopifnot(file.exists(seu_path))
combined <- qread(seu_path); stopifnot(inherits(combined, "Seurat"))
DefaultAssay(combined) <- "RNA"
md <- combined@meta.data

# Unify β label (treat all cells as β)
combined$celltype_old       <- md$celltype
combined$celltype_qadir_old <- md$celltype_qadir
combined$celltype <- factor(rep("beta", ncol(combined)), levels = "beta")

# Normalize Sex -> Sex_std in {M,F}
combined$Sex_std <- combined$Sex
table(combined@meta.data[["Sex"]])

# Optional: pick/normalize treatment column if present
treat_col <- pick_first(colnames(md), c("treatment","Treatment","Tx","tx","condition","Condition"))

## 3) Compute MQI axes + MQI_v1 on ALL β ---------------------
uprmt_genes <- c("LONP1","CLPP","CLPX","HTRA2","HSPD1","HSPE1","HSPA9",
                 "DNAJA3","YME1L1","SPG7","AFG3L2","LRPPRC","ATF5","ATF4","DDIT3")
if (!"UPRmtScore" %in% colnames(combined@meta.data)) {
  combined <- AddModuleScoreSafe(combined, uprmt_genes, "UPRmtScore")
}

need_axes <- c("Mitophagy","OXPHOS1","percent.mt")
stopifnot(all(need_axes %in% colnames(combined@meta.data)))
if (!"TCA1"   %in% colnames(combined@meta.data)) combined$TCA1   <- NA_real_
if (!"Glyco1" %in% colnames(combined@meta.data)) combined$Glyco1 <- NA_real_

compute_MQI <- function(md,
                        weights = list(mitophagy=+1.0, oxphos=+1.0, uprmt=+0.7, tca_pen=+0.2, mtfrac=+0.5, imbalance=+0.5),
                        flips   = list(mitophagy=FALSE, oxphos=FALSE, uprmt=FALSE, tca=FALSE, glyco=FALSE, mtfrac=TRUE)) {
  a <- list(
    mitophagy = md$Mitophagy,
    oxphos    = md$OXPHOS1,
    tca       = md$TCA1 %||% NA_real_,
    glyco     = md$Glyco1 %||% NA_real_,
    mtfrac    = md$percent.mt %||% NA_real_,
    uprmt     = md$UPRmtScore %||% NA_real_
  )
  for (nm in names(a)) { if (isTRUE(flips[[nm]])) a[[nm]] <- -a[[nm]]; a[[nm]] <- z(a[[nm]]) }
  imb  <- if (all(!is.finite(a$oxphos)) || all(!is.finite(a$glyco))) rep(0, nrow(md)) else z(abs(a$oxphos - a$glyco))
  good <- weights$mitophagy*a$mitophagy + weights$oxphos*a$oxphos + weights$uprmt*a$uprmt
  pen  <- 0
  if (any(is.finite(a$tca)))    pen <- pen + weights$tca_pen * pmax(0, a$tca)
  if (any(is.finite(a$mtfrac))) pen <- pen + weights$mtfrac   * z(a$mtfrac)
  if (any(is.finite(imb)))      pen <- pen + weights$imbalance* imb
  as.numeric(good - pen)
}
combined$MQI_v1 <- compute_MQI(combined@meta.data)

# Keep only cells with finite MQI_v1
keep_cells <- colnames(combined)[is.finite(combined$MQI_v1)]
combined <- subset(combined, cells = keep_cells); DefaultAssay(combined) <- "RNA"

## 4) Regulator universe -------------------------------------
# Note: `unique(c(uprmt_genes, ...))` = include UPRmt + extras, removing duplicates
mitophagy_core <- c("PINK1","PRKN","TOMM7","TOMM20","TOMM22",
                    "SQSTM1","OPTN","CALCOCO2","TAX1BP1","NBR1","TBK1","TBKBP1",
                    "ULK1","ULK2","RB1CC1","ATG13")
mitophagy_receptors <- c("BNIP3","BNIP3L","FUNDC1","PHB2","FKBP8","BCL2L13","NLRX1","NIPSNAP1","NIPSNAP2","AMBRA1")
mito_ub_e3_dub <- c("MARCHF5","MUL1","RNF185","HUWE1","USP30","USP35","USP15")
dynamics <- c("DNM1L","MFF","MIEF1","MIEF2","FIS1","GDAP1","OPA1","MFN1","MFN2","OMA1","YME1L1","IMMT")
biogenesis <- c("PPARGC1A","PPARGC1B","TFAM","TFB2M","TFB1M","POLG","POLG2","POLRMT",
                "NRF1","NFE2L2","ESRRA","ESRRB","ESRRG","PPARA","PPARD","PPARG","CREB1")
uprmt_proteostasis <- unique(c(uprmt_genes, "SPG7","AFG3L2","YME1L1","HSPD1","HSPE1","HSPA9","DNAJA3","LRPPRC"))
antioxidant <- c("SOD2","PRDX3","PRDX5","GPX1","GPX4","TXN2","GLRX2","CAT","HMOX1","NQO1","KEAP1")
autophagy_lys <- c("BECN1","PIK3C3","PIK3R4","WIPI1","WIPI2","ATG5","ATG7","ATG3","ATG12","ATG16L1",
                   "UVRAG","TECPR1","GABARAP","GABARAPL1","GABARAPL2","MAP1LC3A","MAP1LC3B","MAP1LC3C")
er_phagy <- c("FAM134A","FAM134B","FAM134C","RTN3","ATL3","SEC62","CCPG1","TEX264","CALCOCO1","VAPA","VAPB",
              "UFM1","UFL1","DDRGK1","UFSP2","CDK5RAP3")
signaling_tfs <- c("PRKAA1","PRKAA2","PRKAB1","PRKAB2","PRKAG1","PRKAG2","PRKAG3",
                   "HIF1A","EPAS1","ARNT","ESR1","ESR2","AR","FOXO1","FOXO3","ATF2","ATF3","ATF5",
                   "TFEB","TFE3","MITF","PPARG","PPARA","PPARD","RELA","NFKB1","STAT3","JUN","FOS")

regulators_all <- unique(c(mitophagy_core, mitophagy_receptors, mito_ub_e3_dub, dynamics,
                           biogenesis, uprmt_proteostasis, antioxidant, autophagy_lys,
                           er_phagy, signaling_tfs))
present <- sort(intersect(regulators_all, rownames(combined)))

# Fallback to HVGs if too few present
if (length(present) < 12) {
  message("Few regulators present; switching to HVGs fallback.")
  if (length(VariableFeatures(combined)) == 0)
    combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  present <- VariableFeatures(combined)
}

## 5) Build X/y on ALL β (training = prediction set) ---------
cells <- colnames(combined)
df_xy <- FetchData(combined, vars = c(present, "MQI_v1"), cells = cells)
X_raw <- df_xy[, setdiff(colnames(df_xy), "MQI_v1"), drop = FALSE]
y     <- df_xy$MQI_v1; names(y) <- rownames(df_xy)

ok <- is.finite(y) & apply(X_raw, 1, function(r) all(is.finite(r)))
X_raw <- X_raw[ok, , drop = FALSE]; y <- y[ok]

# Save TRAINING scaler (mu/sd) & feature order
mu_train <- colMeans(X_raw, na.rm = TRUE)
sd_train <- apply(X_raw, 2, sd); sd_train[!is.finite(sd_train) | sd_train == 0] <- 1
feat <- names(mu_train)
X <- sweep(sweep(X_raw, 2, mu_train, "-"), 2, sd_train, "/")

stopifnot(nrow(X) == length(y), identical(rownames(X), names(y)))

## 6) Define KAN (with knots registered as buffer) -----------
KANSpline1D <- nn_module(
  "KANSpline1D",
  initialize = function(n_bins = 10, xmin = -3, xmax = 3) {
    self$n_bins <- n_bins
    self$register_buffer("knots", torch_linspace(xmin, xmax, n_bins))  # saved in state_dict
    self$coef  <- nn_parameter(torch_zeros(n_bins))
  },
  forward = function(x) {
    x_exp <- x$unsqueeze(2)
    d <- (x_exp - self$knots$unsqueeze(1))$abs()
    w <- (self$knots[2] - self$knots[1])
    phi <- (1 - d / w)$clamp(min = 0)
    (phi * self$coef)$sum(dim = 2)
  }
)

KANLayer <- nn_module(
  "KANLayer",
  initialize = function(in_features, out_features, n_bins = 10, xmin = -3, xmax = 3) {
    self$in_features <- in_features
    self$out_features <- out_features
    self$splines <- nn_module_list(lapply(1:(out_features*in_features), function(i)
      KANSpline1D(n_bins = n_bins, xmin = xmin, xmax = xmax)
    ))
    self$linear <- nn_linear(in_features, out_features)
    self$bias   <- nn_parameter(torch_zeros(out_features))
  },
  forward = function(x) {
    N <- x$size(1); P <- self$in_features
    outs <- vector("list", self$out_features); idx <- 1
    for (j in 1:self$out_features) {
      acc <- torch_zeros(N)
      for (i in 1:P) { acc <- acc + self$splines[[idx]](x[, i]); idx <- idx + 1 }
      outs[[j]] <- acc
    }
    S <- torch_stack(outs, dim = 2)
    S + self$linear(x) + self$bias$unsqueeze(1)
  }
)

KANNet <- nn_module(
  "KANNet",
  initialize = function(in_features, hidden = 48, out_features = 1, n_bins = 10) {
    self$layer1 <- KANLayer(in_features, hidden, n_bins = n_bins)
    self$act1   <- nn_relu()
    self$layer2 <- KANLayer(hidden, out_features, n_bins = n_bins)
  },
  forward = function(x) self$layer2(self$act1(self$layer1(x)))
)

train_kan <- function(X, y, epochs = 300, lr = 2e-3, n_bins = 10, hidden = 48, verbose = TRUE) {
  X_t <- torch_tensor(as.matrix(X), dtype = torch_float())
  y_t <- torch_tensor(as.numeric(y), dtype = torch_float())$unsqueeze(2)
  net <- KANNet(ncol(X), hidden = hidden, out_features = 1, n_bins = n_bins)
  opt <- optim_adam(net$parameters, lr = lr)
  loss_fn <- nn_mse_loss()
  for (e in 1:epochs) {
    net$train(); opt$zero_grad()
    pred <- net(X_t); loss <- loss_fn(pred, y_t)
    loss$backward(); opt$step()
    if (verbose && (e %% 20 == 0)) cat(sprintf("epoch %d  mse=%.4f\n", e, loss$item()))
  }
  list(model = net)
}

## 7) Train on ALL β -----------------------------------------
fit_all <- train_kan(X, y, epochs = 300, lr = 2e-3, n_bins = 10, hidden = 48, verbose = TRUE)

# Save model (torch) + scaler (qs) for future reuse
torch_save(fit_all$model$state_dict(), model_pt)
# when saving:
qsave(list(
  mu   = mu_train,
  sd   = sd_train,
  feat = feat,                 # exact column order used to train
  arch = list(hidden = 40, n_bins = 10),  # <- whatever you trained with
  info = list(seed = 42, date = as.character(Sys.time()))
), scaler_qs)

# reconstruct
library(torch); library(qs)

# 7a) Point to the folder where you saved things
outdir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN"

# 7b) Build absolute file paths
model_pt  <- file.path(outdir, "fit_all_state.pt")
scaler_qs <- file.path(outdir, "input_scaler.qs")

# 7c) Sanity checks
normalizePath(outdir)
file.exists(model_pt); file.exists(scaler_qs)
list.files(outdir, pattern = "fit_all_state|input_scaler", ignore.case = TRUE)

# 7d) Load
library(torch); library(qs)
sc    <- qread(scaler_qs)
state <- torch_load(model_pt)

# 7e) Rebuild EXACT architecture you trained with
model <- KANNet$new(
  in_features  = length(sc$feat),
  hidden       = 48,
  n_bins       = 10,
  out_features = 1
)
model$load_state_dict(state)
model$eval()

# #Paranoia
# st0 <- fit_all$model$state_dict()
# st1 <- model$state_dict()
# max_diff <- max(sapply(names(st0), function(nm)
#   max(abs(as.array(st0[[nm]]) - as.array(st1[[nm]])))))
# max_diff   # expect ~0
# 
# ## 8) Predict for ALL β and map to metadata ------------------
# fit_all$model$eval()
# pred <- as.numeric(fit_all$model$forward(torch::torch_tensor(as.matrix(X), dtype = torch::torch_float())))
# names(pred) <- rownames(X)
# r2   <- cor(pred, y)^2
# rmse <- sqrt(mean((pred - y)^2))
# cat(sprintf("ALL β: n=%d | R^2=%.3f | RMSE=%.4f\n", length(y), r2, rmse))
# 
# # Initialize / write to `combined`
# combined$MQI_obs   <- NA_real_
# combined$MQI_pred  <- NA_real_
# combined$MQI_resid <- NA_real_
# stopifnot(all(names(y)    %in% colnames(combined)),
#           all(names(pred) %in% colnames(combined)))
# combined$MQI_obs [names(y)]   <- as.numeric(y)
# combined$MQI_pred[names(pred)]<- pred
# combined$MQI_resid            <- combined$MQI_obs - combined$MQI_pred
# cat("# cells with MQI_pred in combined: ", sum(!is.na(combined$MQI_pred)), " / ", ncol(combined), "\n")

## 9) FIGURES (PDF only) -------------------------------------

# 9a) Predicted vs Observed
p_pred <- ggplot(data.frame(obs=y, pred=pred), aes(obs, pred)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.02, vjust = -0.6,
           label = sprintf("R²=%.3f  RMSE=%.3f", r2, rmse)) +
  labs(title = "Predicted vs Observed MQI (ALL β)", x = "Observed MQI_v1", y = "Predicted MQI") +
  theme_classic(base_size = 12)
ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/01_pred_vs_obs_all_beta.pdf",
       p_pred, width = 5, height = 4)

# 9b) Residuals by Sex (if available)
if (any(!is.na(combined$Sex_std))) {
  df_res_sex <- data.frame(Sex = combined$Sex_std[names(y)], resid = y - pred)
  p_res_sex <- ggplot(df_res_sex, aes(x = Sex, y = resid)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_boxplot(outlier.size = 0.6) +
    labs(title = "Residuals by Sex (ALL β)", y = "Residual (Obs - Pred)", x = NULL) +
    theme_classic(base_size = 12)
  ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/02_residuals_by_sex_all_beta.pdf",
         p_res_sex, width = 4.8, height = 4.0)
  }

# 9c) Residuals by Treatment (if available)
if (!is.na(treat_col) && !all(is.na(combined@meta.data[[treat_col]]))) {
  df_res_tx <- data.frame(Tx = combined@meta.data[names(y), treat_col], resid = y - pred)
  p_res_tx <- ggplot(df_res_tx, aes(x = Tx, y = resid)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_boxplot(outlier.size = 0.6) +
    labs(title = "Residuals by Treatment (ALL β)", y = "Residual (Obs - Pred)", x = NULL) +
    theme_classic(base_size = 12)
  ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/03_residuals_by_treatment_all_beta.pdf",
         p_res_tx, width = 6, height = 4.2)
}

# 9d) Axis densities + correlation heatmap
axis_cols <- c("Mitophagy","OXPHOS1","TCA1","Glyco1","percent.mt","UPRmtScore","MQI_v1")
axis_cols <- axis_cols[axis_cols %in% colnames(combined@meta.data)]
axis_long <- combined@meta.data[names(y), axis_cols, drop = FALSE] |>
  tibble::rownames_to_column("cell") |>
  pivot_longer(-cell, names_to = "Axis", values_to = "Value")

if (any(!is.na(combined$Sex_std))) {
  axis_long$Sex <- combined$Sex_std[axis_long$cell]
  p_den <- ggplot(axis_long, aes(Value, fill = Sex)) +
    geom_density(alpha = 0.35) +
    facet_wrap(~Axis, scales = "free", ncol = 3) +
    labs(title = "Axis distributions by Sex (ALL β)") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
} else {
  p_den <- ggplot(axis_long, aes(Value)) +
    geom_density(fill = "grey70", alpha = 0.7) +
    facet_wrap(~Axis, scales = "free", ncol = 3) +
    labs(title = "Axis distributions (ALL β)") +
    theme_classic(base_size = 11)
}
ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/04_axis_densities_all_beta.pdf",
       p_den, width = 8.5, height = 6.5)

cors <- cor(combined@meta.data[names(y), axis_cols, drop = FALSE],
            use = "pairwise.complete.obs", method = "spearman")
cors_df <- as.data.frame(as.table(cors)); colnames(cors_df) <- c("Var1","Var2","Corr")
p_corr <- ggplot(cors_df, aes(Var1, Var2, fill = Corr)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1,1)) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Axis correlation heatmap (ALL β)", x = NULL, y = NULL)
ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/05_axis_cor_heatmap_all_beta.pdf",
       p_corr, width = 5.4, height = 4.8)

# 9e) Regulator SD histogram + missing list
sd_vec <- apply(X, 2, sd)
p_sd <- ggplot(data.frame(SD = sd_vec), aes(SD)) +
  geom_histogram(bins = 40, color = "white") +
  geom_vline(xintercept = 0, linetype = 3) +
  labs(title = "Per-gene SD of inputs (z-scored, ALL β)", x = "SD", y = "Count") +
  theme_classic(base_size = 11)
ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/06_regulator_sd_hist_all_beta.pdf",
       p_sd, width = 5.2, height = 4.2)

present_genes <- colnames(X)
missing_genes <- setdiff(unique(regulators_all), present_genes)
if (length(missing_genes) > 0) {
  write.table(data.frame(missing_genes),
              file.path(outdir, "missing_regulators_all_beta.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}

# 9f) Spline gallery (Layer1) — multi-page PDF
plot_splines_paginated <- function(fit, feature_names, ncol = 4, per_page = 16,
                                   file = "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/07_splines_all_paginated_all_beta.pdf",
                                   width = 14, height = 10) {
  stopifnot(inherits(fit$model, "KANNet"))
  layer <- fit$model$layer1; P <- layer$in_features
  grDevices::pdf(file, width = width, height = height); on.exit(grDevices::dev.off(), add = TRUE)
  xs <- seq(-3, 3, length.out = 200); i <- 1
  while (i <= P) {
    n_this <- min(per_page, P - i + 1); rows <- ceiling(n_this / ncol)
    oldpar <- par(mfrow = c(rows, ncol), mar = c(2.2, 2.2, 1.6, 0.8)); on.exit(par(oldpar), add = TRUE)
    for (k in 0:(n_this - 1)) {
      idx <- i + k; vals <- matrix(0, nrow = length(xs), ncol = layer$out_features)
      for (j in 1:layer$out_features) {
        cell_idx <- (j - 1) * P + idx
        vals[, j] <- sapply(xs, function(xx)
          fit$model$layer1$splines[[cell_idx]](torch::torch_tensor(xx))$item()
        )
      }
      plot(xs, rowMeans(vals), type = "l", xlab = feature_names[idx], ylab = "spline effect",
           main = feature_names[idx]); abline(h = 0, lty = 3)
    }
    i <- i + n_this
  }
}
plot_splines_paginated(fit_all, colnames(X))

# 9g) Latent PCA (hidden layer) — 2D PDF (no plotly)
# H_t <- fit_all$model$layer1$forward(torch::torch_tensor(as.matrix(X), dtype = torch::torch_float()))
# H_t <- torch::nnf_relu(H_t); H <- as.array(H_t)
# pc <- prcomp(H, center = TRUE, scale. = TRUE)
# coords <- as.data.frame(pc$x[, 1:2]); colnames(coords) <- c("PC1","PC2")
# coords$MQI <- y
# coords$Sex <- combined$Sex_std[rownames(X)]
# 
# p_latent <- ggplot(coords, aes(PC1, PC2, color = MQI)) +
#   geom_point(size = 0.6, alpha = 0.6) +
#   labs(title = "KAN latent 2D (ALL β)", x = "PC1", y = "PC2") +
#   theme_classic(base_size = 12)
# ggsave("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/08_kan_latent_2d_all_beta.pdf",
#        p_latent, width = 5.8, height = 4.8)

# 9h) Partial dependence surface — vector PDFs (surface + mesh)
if (!requireNamespace("plot3D", quietly = TRUE)) install.packages("plot3D")
library(plot3D)
g1 <- colnames(X)[1]; g2 <- colnames(X)[2]  # or pick from importance below
grid <- seq(-2.5, 2.5, length.out = 60)
x_base <- colMeans(X, na.rm = TRUE)
G <- expand.grid(v1 = grid, v2 = grid)
Xg <- matrix(rep(x_base, each = nrow(G)), nrow = nrow(G), byrow = FALSE); colnames(Xg) <- colnames(X)
Xg[, g1] <- G$v1; Xg[, g2] <- G$v2
pred_grid <- as.numeric(fit_all$model$forward(torch::torch_tensor(Xg, dtype = torch::torch_float())))
Zmat <- matrix(pred_grid, nrow = length(grid), ncol = length(grid), byrow = FALSE)

cairo_pdf(sprintf("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/09_pdp_%s_x_%s_surface_all_beta.pdf", g1, g2),
          width = 7.5, height = 5.5, family = "Helvetica")

plot3D::persp3D(x = grid, y = grid, z = Zmat, theta = 40, phi = 25, expand = 0.8, shade = 0.3,
                border = "black", lwd = 0.1, clab = "Predicted MQI",
                colkey = list(length = 0.5, width = 0.8, dist = -0.02, cex.clab = 0.9, cex.axis = 0.8),
                ticktype = "detailed", cex.axis = 0.8,
                xlab = paste0(g1, " (z)"), ylab = paste0(g2, " (z)"), zlab = "Predicted MQI")
dev.off()

# cairo_pdf(file.path(outdir, sprintf("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/qc_images/10_pdp_%s_x_%s_mesh_all_beta.pdf", g1, g2)),
#           width = 7.5, height = 5.5, family = "Helvetica")
# plot3D::persp3D(x = grid, y = grid, z = Zmat, theta = 40, phi = 25, expand = 0.8,
#                 facets = NA, border = "black", lwd = 0.6, ticktype = "detailed", cex.axis = 0.8,
#                 xlab = paste0(g1, " (z)"), ylab = paste0(g2, " (z)"), zlab = "Predicted MQI")
# dev.off()

# 9i) Permutation importance (batched) + barplot
## 9i) Permutation importance (batched) — multi-seed stability + consensus

# If these helpers aren't already defined above in your script, uncomment:

# ================= OVERNIGHT FULL-DATA PERMUTATION IMPORTANCE =================
# Assumes:
# - fit_all : trained torch model wrapper with $model (nn_module)
# - X       : data.frame/matrix (all rows, all features)
# - y       : numeric target vector (length nrow(X))
# - optional: combined_path (to place outputs next to your dataset)
# ==============================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("qs", quietly = TRUE)) install.packages("qs")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
})

# -------- Batched prediction helper (no subsampling) --------
predict_batches <- function(model, X_mat, batch = 32768L) {
  if (!is.matrix(X_mat)) X_mat <- as.matrix(X_mat)
  n <- nrow(X_mat); out <- numeric(n); i <- 1L
  while (i <= n) {
    j <- min(i + batch - 1L, n)
    out[i:j] <- as.numeric(
      model$forward(
        torch::torch_tensor(X_mat[i:j, , drop = FALSE], dtype = torch::torch_float())
      )
    )
    i <- j + 1L
  }
  out
}

# -------- Full-data permutation importance (memory-safe; on-the-fly permute) ---
perm_importance_full <- function(
    fit, X_df, y_vec,
    n_repeats  = 3,
    seed       = 42,
    batch_rows = 32768L,
    verbose    = TRUE
) {
  if (!is.matrix(X_df)) X_df <- as.matrix(X_df)
  y_vec <- as.numeric(y_vec)
  n <- nrow(X_df); p <- ncol(X_df)
  coln <- colnames(X_df)
  
  if (verbose) message(sprintf("[PI] n=%d, p=%d, repeats=%d, seed=%d, batch=%d",
                               n, p, n_repeats, seed, batch_rows))
  
  fit$model$eval()
  set.seed(seed)
  if ("torch_manual_seed" %in% getNamespaceExports("torch")) torch::torch_manual_seed(seed)
  
  # 1) Baseline MSE (ALL rows, batched)
  if (verbose) message("[PI] Baseline predictions...")
  base_pred <- predict_batches(fit$model, X_df, batch = batch_rows)
  base_mse  <- mean((base_pred - y_vec)^2)
  
  # 2) Precompute global row permutations for each repeat
  if (verbose) message("[PI] Generating permutations...")
  perms <- replicate(n_repeats, sample.int(n, size = n), simplify = FALSE)
  
  # 3) Feature loop; permute one column per batch only
  inc <- numeric(p); names(inc) <- coln
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = p, style = 3)
  
  for (j in seq_len(p)) {
    mses <- numeric(n_repeats)
    for (r in seq_len(n_repeats)) {
      perm_r <- perms[[r]]
      se_sum <- 0
      i <- 1L
      while (i <= n) {
        k <- min(i + batch_rows - 1L, n)
        idx <- i:k
        Xb <- X_df[idx, , drop = FALSE]
        Xb[, j] <- X_df[perm_r[idx], j]               # on-the-fly col permute for this batch
        pred_b <- predict_batches(fit$model, Xb, batch = batch_rows)
        rb <- pred_b - y_vec[idx]
        se_sum <- se_sum + sum(rb * rb)
        i <- k + 1L
      }
      mses[r] <- se_sum / n
    }
    inc[j] <- mean(mses) - base_mse
    if (verbose) utils::setTxtProgressBar(pb, j)
  }
  if (verbose) close(pb)
  
  sort(inc, decreasing = TRUE)
}

# -------- Output paths --------
if (!exists("outdir")) {
  if (exists("combined_path") && !is.null(combined_path) && file.exists(combined_path)) {
    outdir <- file.path(dirname(combined_path), "qc_images")
  } else {
    outdir <- file.path(getwd(), "qc_images")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}
imp_qs <- if (exists("combined_path") && !is.null(combined_path)) {
  file.path(dirname(combined_path), "imp_all.qs")
} else {
  file.path(getwd(), "imp_all.qs")
}

# -------- Run (all rows) across multiple seeds for stability --------
seeds <- c(42, 73, 101, 202)   # edit if you want fewer/more seeds
batch_rows <- 32768L                # safe with 128GB; lower if you want gentler bursts
n_repeats  <- 3

imp_runs <- lapply(seeds, function(s) {
  perm_importance_full(
    fit        = fit_all,
    X_df       = X,
    y_vec      = y,
    n_repeats  = n_repeats,
    seed       = s,
    batch_rows = batch_rows,
    verbose    = TRUE
  )
})
names(imp_runs) <- paste0("seed", seeds)

# -------- Save outputs --------
# 1) Back-compat single vector (first seed)
imp_all <- imp_runs[[1L]]
qs::qsave(imp_all, imp_qs)

# 2) Full list of seed runs
qs::qsave(imp_runs, file.path(dirname(imp_qs), "imp_all_multiseed.qs"))

# 3) Mean ± SD across seeds (aligned features)
all_feats <- Reduce(union, lapply(imp_runs, names))
mat <- do.call(cbind, lapply(imp_runs, function(v) setNames(v[all_feats], all_feats)))
# (fill NAs if any feature missing in a run — shouldn't happen, but just in case)
for (k in seq_len(ncol(mat))) {
  miss <- is.na(mat[, k]); if (any(miss)) mat[miss, k] <- 0
}
imp_mean <- rowMeans(mat)
imp_sd   <- apply(mat, 1L, sd)
stab_df  <- data.frame(
  feature = rownames(mat),
  mean_increase_mse = as.numeric(imp_mean),
  sd_increase_mse   = as.numeric(imp_sd),
  stringsAsFactors = FALSE
)
stab_df <- stab_df[order(stab_df$mean_increase_mse, decreasing = TRUE), ]
csv_path <- file.path(dirname(imp_qs), "imp_all_multiseed_mean_sd.csv")
write.csv(stab_df, csv_path, row.names = FALSE)

# 4) Optional quick barplot of top-N features (PDF)
topN <- 30L
plot_df <- head(stab_df, topN)
pdf_path <- file.path(outdir, sprintf("perm_importance_full_top%d.pdf", topN))
suppressWarnings({
  library(ggplot2)
  gg <- ggplot(plot_df, aes(x = reorder(feature, mean_increase_mse),
                            y = mean_increase_mse)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_increase_mse - sd_increase_mse,
                      ymax = mean_increase_mse + sd_increase_mse), width = 0.2) +
    coord_flip() +
    labs(x = NULL, y = "Permutation importance (ΔMSE)",
         title = sprintf("Permutation importance (all cells) — top %d", topN),
         subtitle = sprintf("Seeds: %s; repeats per seed: %d; batch_rows=%d",
                            paste(seeds, collapse = ", "), n_repeats, batch_rows)) +
    theme_classic(base_size = 11)
  ggsave(pdf_path, gg, width = 7, height = 9, limitsize = FALSE)
})

message("[PI] DONE.")
message(sprintf("Single-seed vector: %s", imp_qs))
message(sprintf("Multi-seed list:    %s", file.path(dirname(imp_qs), "imp_all_multiseed.qs")))
message(sprintf("Mean±SD CSV:        %s", csv_path))
message(sprintf("Top-N barplot PDF:  %s", pdf_path))

#SANITY CHECKS
# Installed packages + versions
installed.packages()[, c("Package", "Version")]

# Current session info
sessionInfo()

# All objects in environment with their class + size
objs <- ls()
data.frame(
  object = objs,
  class  = sapply(objs, function(x) class(get(x))[1]),
  size   = sapply(objs, function(x) format(object.size(get(x)), units = "MB"))
)

## ==== SANITY CHECK: load + validate all artifacts (no model rebuild required) ====
base_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/Perma_analysis"

files <- c(
  state_pt   = file.path(base_dir, "fit_all_state.pt"),
  imp_vec    = file.path(base_dir, "imp_all.qs"),
  imp_list   = file.path(base_dir, "imp_all_multiseed.qs"),
  imp_csv    = file.path(base_dir, "imp_all_multiseed_mean_sd.csv"),
  scaler_qs  = file.path(base_dir, "input_scaler.qs")   # optional but recommended
)

# 1) Existence
missing <- files[!file.exists(files)]
if (length(missing)) {
  stop(sprintf("Missing files:\n- %s", paste(names(missing), missing, sep=": ", collapse = "\n")))
}

suppressPackageStartupMessages({
  library(qs)
  library(torch)
})

# 2) Load everything
imp_all        <- qs::qread(files["imp_vec"])
imp_runs       <- qs::qread(files["imp_list"])
stab_df        <- read.csv(files["imp_csv"], stringsAsFactors = FALSE)
input_scaler   <- if (file.exists(files["scaler_qs"])) qs::qread(files["scaler_qs"]) else NULL

# state dict (no model needed)
state_dict     <- torch::torch_load(files["state_pt"])

# 3) Basic structure checks
stopifnot(is.numeric(imp_all), length(imp_all) > 0)
stopifnot(is.list(imp_runs), length(imp_runs) >= 1)
stopifnot(is.data.frame(stab_df), all(c("feature","mean_increase_mse","sd_increase_mse") %in% names(stab_df)))

# Agree on feature count
p_imp_vec  <- length(imp_all)
p_imp_csv  <- length(unique(stab_df$feature))
p_imp_runs <- unique(lengths(imp_runs))  # length of each seed vector

if (length(p_imp_runs) != 1L) {
  warning("Not all seed runs have same length; check imp_all_multiseed.qs")
}
p_imp_runs <- p_imp_runs[1]

# If scaler present, infer p from it
p_scaler <- NA_integer_
if (!is.null(input_scaler)) {
  # try common field names
  cand <- c("mu","mean","center","mu_train","center_","mean_")
  found <- cand[cand %in% names(input_scaler)]
  if (length(found)) {
    p_scaler <- length(input_scaler[[found[1]]])
  }
}

# 4) State dict peek (names + a few shapes)
state_names <- names(state_dict)
some_shapes <- lapply(state_dict[head(state_names, 8)], function(t) as.integer(t$size()))

# 5) Cross-check dimensions agree (expect 141 features)
p_set <- na.omit(c(p_imp_vec, p_imp_csv, p_imp_runs, p_scaler))
agree <- length(unique(p_set)) == 1L

# 6) Report
cat("\n=== SANITY CHECK REPORT ===\n")
cat("Folder: ", base_dir, "\n")
cat("- Files present: OK\n")
cat(sprintf("- imp_all.qs:              %d features\n", p_imp_vec))
cat(sprintf("- imp_all_multiseed.qs:    %d features/run, %d seeds\n", p_imp_runs, length(imp_runs)))
cat(sprintf("- mean±sd CSV:             %d unique features\n", p_imp_csv))
if (!is.na(p_scaler)) cat(sprintf("- input_scaler.qs:         %d features (from scaler vector)\n", p_scaler)) else cat("- input_scaler.qs:         loaded (feature length not inferred)\n")
cat(sprintf("- State dict tensors:      %d params (showing first few shapes below)\n", length(state_names)))
state_names <- names(state_dict)
peek_n <- min(8L, length(state_names))
for (i in seq_len(peek_n)) {
  shp <- paste(as.integer(state_dict[[state_names[i]]]$size()), collapse = "x")
  cat(sprintf("%3d: %s  [%s]\n", i, state_names[i], shp))
}
cat(sprintf("\n- Feature counts agree across artifacts: %s\n", if (agree) "YES ✅" else "NO ❌"))

if (!agree) {
  cat("  → Check that all files were produced by the same run (same feature set/order).\n")
} else {
  cat("\nAll artifacts deserialized successfully. You’re good to proceed. 🚀\n")
}

# ==============================================================================

# --- Optional: stability across seeds (Spearman rho + Jaccard top-K) ---
feat_names <- colnames(X)
# align each vector to the same feature order
imp_aligned <- lapply(imp_runs, function(v) setNames(v[feat_names], feat_names))
S <- length(imp_aligned); seed_names <- names(imp_aligned)

rho_mat <- matrix(NA_real_, S, S, dimnames = list(seed_names, seed_names))
topK <- 30L
tops <- lapply(imp_aligned, function(v) names(sort(v, decreasing = TRUE))[1:topK])

for (i in seq_len(S)) {
  for (j in i:S) {
    r <- suppressWarnings(cor(imp_aligned[[i]], imp_aligned[[j]],
                              method = "spearman", use = "pairwise.complete.obs"))
    rho_mat[i, j] <- rho_mat[j, i] <- r
  }
}
jac <- function(a,b) length(intersect(a,b)) / length(union(a,b))
jac_mat <- matrix(NA_real_, S, S, dimnames = list(seed_names, seed_names))
for (i in seq_len(S)) for (j in i:S) {
  jj <- jac(tops[[i]], tops[[j]])
  jac_mat[i, j] <- jac_mat[j, i] <- jj
}

# (Optional) write them out
write.csv(rho_mat, file.path(dirname(imp_qs), "imp_multiseed_spearman.csv"), row.names = TRUE)
write.csv(jac_mat, file.path(dirname(imp_qs), "imp_multiseed_jaccard_top30.csv"), row.names = TRUE)


# --- 4) Convert each run to % baseline MSE on the FIXED subset
base_pred_fix <- predict_batches(fit_all$model, as.matrix(X[row_idx, , drop = FALSE]), batch = 32768L)
base_mse_fix  <- mean((base_pred_fix - as.numeric(y[row_idx]))^2)
imp_pct_runs  <- lapply(imp_runs, function(v) { vv <- 100 * as.numeric(v) / base_mse_fix; names(vv) <- names(v); vv })

# --- 5) Consensus (mean ± SD across seeds) and barplot (PDF)
imp_mat <- do.call(cbind, imp_pct_runs)   # genes x seeds
imp_mean <- rowMeans(imp_mat)
imp_sd   <- apply(imp_mat, 1, sd)
consensus <- data.frame(Gene = rownames(imp_mat), MeanPct = imp_mean, SDPct = imp_sd, row.names = NULL)
consensus <- consensus[order(-consensus$MeanPct), ]
top_cons  <- head(consensus, topK)

p_cons <- ggplot2::ggplot(top_cons, ggplot2::aes(x = reorder(Gene, -MeanPct), y = MeanPct)) +
  ggplot2::geom_col(width = 0.82) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, MeanPct - SDPct), ymax = MeanPct + SDPct), width = 0.3) +
  ggplot2::labs(x = NULL, y = "% MSE increase (perm.)",
                title = paste0("Permutation importance — consensus across seeds (top-", topK, ")")) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
ggplot2::ggsave(file.path(outdir, "11_importance_consensus_topK.pdf"), p_cons, width = 8, height = 6.4)

# --- 6) Save stability heatmaps (PDF)
as_df <- function(M, lab) { df <- as.data.frame(as.table(M)); colnames(df) <- c("A","B",lab); df }
rho_df <- as_df(rho_mat, "rho")
jac_df <- as_df(jac_mat, "jac")

# consistent fill scale from white (0) to dark blue (1)
library(ggplot2)

# Common fill scale: white → steelblue
common_scale <- scale_fill_gradient(limits = c(0.5, 1),
                                    low = "white", high = "dodgerblue4",
                                    name = "Stability (0–1)")

# Function to build heatmap with text labels
make_heatmap <- function(df, value_col, title) {
  ggplot(df, aes(A, B, fill = .data[[value_col]])) +
    geom_tile(color = "grey80") +
    geom_text(aes(label = sprintf("%.3f", .data[[value_col]])), size = 3) +
    common_scale +
    coord_equal() +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL)
}

# Build plots
p_rho <- make_heatmap(rho_df, "rho", "Permutation importance stability — Spearman ρ")
p_jac <- make_heatmap(jac_df, "jac", paste0("Permutation importance stability — Jaccard (top-", topK, ")"))

# Save
ggsave(file.path(outdir, "12_imp_seed_stability_spearman.pdf"),
       p_rho, width = 6, height = 5, device = cairo_pdf)

ggsave(file.path(outdir, "13_imp_seed_stability_jaccard_topK.pdf"),
       p_jac, width = 6, height = 5, device = cairo_pdf)

# --- 7) Optional: save consensus table
write.csv(consensus, file.path(outdir, "imp_consensus_all_genes.csv"), row.names = FALSE)

# Console summary
cat("\nPI stability (Spearman):\n"); print(round(rho_mat, 3))
cat("\nPI stability (Jaccard, top-", topK, "):\n", sep = ""); print(round(jac_mat, 3))
cat("\nSaved PDFs to: ", outdir, "\nSaved single-run PI to: ", imp_qs,
    "\nSaved multi-seed list to: ", file.path(dirname(imp_qs), "imp_all_multiseed.qs"), "\n", sep = "")

###################
#RELOAD LOCK AND LOAD
# 1) Start fresh
rm(list = ls())
gc()   # free memory

# 2) Re-load libraries
library(qs)
library(ggplot2)
library(dplyr)

# 3) Point to your base directory
base_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/KANN/Perma_analysis"
outdir   <- base_dir

# 4) Load your artifacts
imp_all       <- qread(file.path(base_dir, "imp_all.qs"))
imp_runs      <- qread(file.path(base_dir, "imp_all_multiseed.qs"))
stab_df       <- read.csv(file.path(base_dir, "imp_all_multiseed_mean_sd.csv"))
input_scaler  <- qread(file.path(base_dir, "input_scaler.qs"))

# NOTE: fit_all_state.pt is PyTorch format, not .qs — load via torch, not qread()


# 5) Recompute % importance vector
scores <- setNames(stab_df$mean_increase_mse, stab_df$feature)
scores[scores < 0] <- 0
imp_pct <- 100 * scores / sum(scores)

cat("Reload complete. Objects in memory:\n")
print(ls())


# --- Top-50 permutation-importance barplot with family coloring (+ save Seurat) ---
outdir
stopifnot(exists("imp_pct"), length(imp_pct) > 0)
stopifnot(exists("combined"))
if (!exists("outdir") || is.na(outdir) || !nzchar(outdir)) {
  outdir <- file.path(getwd(), "qc_images")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Simple family coloring for top 50 (optional and compact)
reg_sets_simple <- list(
  "Mitophagy machinery" = c("PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","FKBP8","PHB2","BCL2L13",
                            "NIPSNAP1","NIPSNAP2","SQSTM1","OPTN","TAX1BP1","CALCOCO2","NBR1",
                            "MAP1LC3A","MAP1LC3B","GABARAP","GABARAPL1","GABARAPL2","ULK1","ATG3","ATG5","TBK1"),
  "Fission/Fusion"      = c("DNM1L","FIS1","MFF","MIEF1","MIEF2","OPA1","MFN1","MFN2","YME1L1","OMA1"),
  "mtUPR/Proteostasis"  = c("HSPD1","HSPE1","HSPA9","CLPP","LONP1","HTRA2","SPG7"),
  "Import/OMM"          = c("TOMM20","TOMM70A"),
  "ISR/Transcription"   = c("ATF4","DDIT3","ATF5","ATF2")
)

# Make a tidy mapping data.frame (as *character*, not factor)
cat_df <- stack(reg_sets_simple)
colnames(cat_df) <- c("Gene","Category")
cat_df$Category <- as.character(cat_df$Category)

# Case-insensitive mapping just in case
gene2cat <- setNames(cat_df$Category, toupper(cat_df$Gene))

# Build top-k DF
topk <- 50L
df_imp <- head(
  data.frame(Gene = names(imp_pct), IncPct = as.numeric(imp_pct), stringsAsFactors = FALSE),
  topk
)

# Attach categories; fill NAs with "Other"
df_imp$Category <- gene2cat[toupper(df_imp$Gene)]
df_imp$Category[is.na(df_imp$Category)] <- "Other"
df_imp$Category <- factor(df_imp$Category, levels = c(names(reg_sets_simple), "Other"))

# Palette
pal <- c("Mitophagy machinery"="#1b9e77","Fission/Fusion"="#d95f02","mtUPR/Proteostasis"="#7570b3",
         "Import/OMM"="#e6ab02","ISR/Transcription"="#66a61e","Other"="grey60")

# Plot (ASCII-only title to avoid Windows PDF font substitution warnings)
suppressPackageStartupMessages(requireNamespace("ggplot2", quietly = TRUE))
p_imp <- ggplot2::ggplot(df_imp, ggplot2::aes(x = stats::reorder(Gene, -IncPct), y = IncPct, fill = Category)) +
  ggplot2::geom_col(width = 0.82) +
  ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
  ggplot2::labs(x = NULL, y = "% MSE increase (perm.)",
                title = "Top regulators by permutation importance (ALL beta)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "right",
                 axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot2::ggsave(file.path(outdir, "11_importance_top50_all_beta.pdf"), p_imp,
                width = 7.8, height = 6.8, limitsize = FALSE)

# --- Save updated Seurat object (with MQI columns) ---
# We ran into an issue, power outage in Tulane F'd us over. Reconstruct MQI_pred from saved artifacts.
# This is fine—we already verified this exactly reproduces previous predictions.
## === RELOAD combined, RECREATE MODEL FROM STATE, PREDICT MQI_pred, SAVE ===

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(torch)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# --- Paths (edit only if your files live elsewhere) ---
seu_path <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined.qs)"
base_dir <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\Perma_analysis)"
state_pt <- file.path(base_dir, "fit_all_state.pt")
scaler_qs<- file.path(base_dir, "input_scaler.qs")

stopifnot(file.exists(seu_path), file.exists(state_pt), file.exists(scaler_qs))

# --- 1) Load Seurat & scaler ---
combined <- qread(seu_path); stopifnot(inherits(combined, "Seurat"))
DefaultAssay(combined) <- "RNA"

sc   <- qread(scaler_qs)
feat <- sc$feat %||% names(sc$mu)            # exact feature order used in training
mu   <- as.numeric(sc$mu); names(mu)  <- feat
sdv  <- as.numeric(sc$sd); names(sdv) <- feat

# ensure we have all features
stopifnot(length(feat) > 0, length(mu) == length(sdv), all(feat %in% rownames(combined)))

# --- 2) Build X in the SAME order and scale exactly like training ---
cells <- colnames(combined)
expr  <- GetAssayData(combined, assay = DefaultAssay(combined), slot = "data")
X     <- t(as.matrix(expr[feat, cells, drop = FALSE]))   # n x p, columns exactly = feat
X     <- sweep(sweep(X, 2L, mu, "-"), 2L, sdv, "/")

# --- 3) Define KAN classes (same as training) ---
KANSpline1D <- nn_module(
  "KANSpline1D",
  initialize = function(n_bins = 10, xmin = -3, xmax = 3) {
    self$n_bins <- n_bins
    self$register_buffer("knots", torch_linspace(xmin, xmax, n_bins))
    self$coef <- nn_parameter(torch_zeros(n_bins))
  },
  forward = function(x) {
    x_exp <- x$unsqueeze(2)
    d <- (x_exp - self$knots$unsqueeze(1))$abs()
    w <- (self$knots[2] - self$knots[1])
    phi <- (1 - d / w)$clamp(min = 0)
    (phi * self$coef)$sum(dim = 2)
  }
)

KANLayer <- nn_module(
  "KANLayer",
  initialize = function(in_features, out_features, n_bins = 10, xmin = -3, xmax = 3) {
    self$in_features  <- in_features
    self$out_features <- out_features
    self$splines <- nn_module_list(lapply(1:(out_features*in_features), function(i)
      KANSpline1D(n_bins = n_bins, xmin = xmin, xmax = xmax)
    ))
    self$linear <- nn_linear(in_features, out_features)
    self$bias   <- nn_parameter(torch_zeros(out_features))
  },
  forward = function(x) {
    N <- x$size(1); P <- self$in_features
    outs <- vector("list", self$out_features); idx <- 1
    for (j in 1:self$out_features) {
      acc <- torch_zeros(N)
      for (i in 1:P) { acc <- acc + self$splines[[idx]](x[, i]); idx <- idx + 1 }
      outs[[j]] <- acc
    }
    S <- torch_stack(outs, dim = 2)
    S + self$linear(x) + self$bias$unsqueeze(1)
  }
)

KANNet <- nn_module(
  "KANNet",
  initialize = function(in_features, hidden, out_features = 1, n_bins = 10) {
    self$layer1 <- KANLayer(in_features, hidden, n_bins = n_bins)
    self$act1   <- nn_relu()
    self$layer2 <- KANLayer(hidden, out_features, n_bins = n_bins)
  },
  forward = function(x) self$layer2(self$act1(self$layer1(x)))
)

# --- 4) Load state dict and INFER the right dims (no guessing) ---
state    <- torch_load(state_pt)
hidden   <- as.integer(state[["layer1.bias"]]$size()[[1]])
n_bins   <- as.integer(state[["layer1.splines.0.coef"]]$size()[[1]])
in_feat  <- length(feat)
out_feat <- 1L

model <- KANNet$new(in_features = in_feat, hidden = hidden, out_features = out_feat, n_bins = n_bins)
model$load_state_dict(state); model$eval()

# --- 5) Batched predict ---
predict_batches <- function(model, X_mat, batch = 32768L) {
  if (!is.matrix(X_mat)) X_mat <- as.matrix(X_mat)
  n <- nrow(X_mat); out <- numeric(n); i <- 1L
  while (i <= n) {
    j <- min(i + batch - 1L, n)
    out[i:j] <- as.numeric(model$forward(torch_tensor(X_mat[i:j, , drop = FALSE], dtype = torch_float())))
    i <- j + 1L
  }
  out
}
MQI_pred <- predict_batches(model, X, batch = 32768L)
stopifnot(length(MQI_pred) == length(cells))
names(MQI_pred) <- cells

# --- 6) Attach MQI_pred (+ MQI_obs & MQI_resid if available) and save ---
# Initialize columns (safe if they already exist)
combined$MQI_obs   <- NA_real_
combined$MQI_pred  <- NA_real_
combined$MQI_resid <- NA_real_

# Predicted
combined$MQI_pred[cells] <- as.numeric(MQI_pred)

# Observed (if you still have MQI_v1 in metadata)
if ("MQI_v1" %in% colnames(combined@meta.data)) {
  combined$MQI_obs[cells]   <- as.numeric(combined@meta.data[cells, "MQI_v1"])
  combined$MQI_resid[cells] <- combined$MQI_obs[cells] - combined$MQI_pred[cells]
}

# Save updated Seurat
out_qs <- file.path(dirname(seu_path), "combined_with_MQI_all_beta.qs")
qsave(combined, out_qs)

# Also save the prediction vector (optional, for quick reuse)
qsave(MQI_pred, file.path(base_dir, "MQI_pred_allcells.qs"))
write.csv(data.frame(cell = cells, MQI_pred = MQI_pred),
          file.path(base_dir, "MQI_pred_allcells.csv"), row.names = FALSE)

message("\n✅ Done. Updated Seurat with MQI_pred (and MQI_obs/resid if present) saved to:\n", out_qs)
message("   Pred vector saved to: ", file.path(base_dir, "MQI_pred_allcells.qs"))

## ===== Fresh environment: Top-50 PI plot (Mean ± SD across seeds) with paper colors =====
rm(list = ls()); invisible(gc())

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
})

# --- Paths ---
base_dir <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\Perma_analysis)"
seu_dir  <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects)"
plot_out <- base_dir

multiseed_qs <- file.path(base_dir, "imp_all_multiseed.qs")
csv_path     <- file.path(base_dir, "imp_all_multiseed_mean_sd.csv")  # for feature order fallback
seu_with_mqi <- file.path(seu_dir, "combined_with_MQI_all_beta.qs")
seu_fallback <- file.path(seu_dir, "combined.qs")

stopifnot(file.exists(multiseed_qs))
imp_runs <- qread(multiseed_qs)
stab_df  <- if (file.exists(csv_path)) read.csv(csv_path, stringsAsFactors = FALSE) else NULL

# --- Align features & convert each seed to % of its total ΔMSE -----------------------------
feats_union <- Reduce(union, lapply(imp_runs, names))
feats <- if (!is.null(stab_df) && "feature" %in% names(stab_df)) {
  unique(c(stab_df$feature, setdiff(feats_union, stab_df$feature)))
} else feats_union

pct_list <- lapply(imp_runs, function(v) {
  v <- v[feats]; v[is.na(v)] <- 0; v[v < 0] <- 0
  s <- sum(v); if (s == 0) rep(0, length(v)) else 100 * as.numeric(v) / s
})
pct_mat <- do.call(cbind, pct_list)
rownames(pct_mat) <- feats
colnames(pct_mat) <- names(imp_runs)

Mean <- rowMeans(pct_mat)
SD   <- apply(pct_mat, 1, sd)

consensus <- data.frame(Gene = feats, Mean = Mean, SD = SD, stringsAsFactors = FALSE) |>
  arrange(desc(Mean))

# --- Load Seurat (for saving back the same object) -----------------------------------------
combined_path <- if (file.exists(seu_with_mqi)) seu_with_mqi else seu_fallback
stopifnot(file.exists(combined_path))
combined <- qread(combined_path); stopifnot(inherits(combined, "Seurat"))
DefaultAssay(combined) <- "RNA"

# --- Category map (expanded so Top-50 rarely falls to "Other") -----------------------------
reg_sets <- list(
  "Mitophagy Machinery" = c(
    # core & receptors & adaptors
    "PINK1","PRKN","SQSTM1","OPTN","CALCOCO2","TAX1BP1","NBR1",
    "BNIP3","BNIP3L","FUNDC1","PHB2","FKBP8","BCL2L13","AMBRA1",
    # LC3/GABARAPs & ULK complex
    "MAP1LC3A","MAP1LC3B","MAP1LC3C","GABARAP","GABARAPL1","GABARAPL2",
    "ULK1","ULK2","RB1CC1","ATG13","TBK1","TBKBP1",
    # ubiquitin/E3/DUB players often scoring high
    "MARCHF5","MUL1","RNF185","HUWE1","USP30","USP35","USP15",
    # TOMM often appears in your Top-50 and your legend has no Import/OMM bucket,
    # so color TOMM7/TOMM20 here to avoid grey "Other":
    "TOMM7","TOMM20"
  ),
  "Fission/Fusion" = c(
    "DNM1L","FIS1","MFF","MIEF1","MIEF2","OPA1","MFN1","MFN2",
    "YME1L1","OMA1","IMMT"
  ),
  "mtUPR/Proteostasis" = c(
    "HSPD1","HSPE1","HSPA9","DNAJA3","CLPP","LONP1","HTRA2",
    "SPG7","AFG3L2","LRPPRC","YME1L1"
  ),
  "Autophagy/Lysosome" = c(
    "BECN1","PIK3C3","PIK3R4","WIPI1","WIPI2",
    "ATG3","ATG5","ATG7","ATG12","ATG16L1",
    "UVRAG","TECPR1"
  ),
  "TFs/Signalling" = c(
    "ATF2","ATF3","ATF4","ATF5","DDIT3","STAT3",
    "FOXO1","FOXO3","ESR1","ESR2","AR",
    "EPAS1","HIF1A","NFE2L2","RELA","NFKB1","JUN","FOS",
    "TFEB","TFE3","MITF","PPARG","PPARA","PPARD"
  ),
  "ERphagy" = c(
    "FAM134A","FAM134B","FAM134C","RTN3","ATL3","SEC62",
    "CCPG1","TEX264","CALCOCO1","VAPA","VAPB", "UFM1"
  ),
  "Antioxidant" = c(
    "SOD2","PRDX3","PRDX5","GPX1","GPX4","TXN2","GLRX2",
    "CAT","HMOX1","NQO1","KEAP1"
  )
)

# tidy, case-insensitive mapping
cat_df <- stack(reg_sets); colnames(cat_df) <- c("Gene","Category")
cat_df$Category <- as.character(cat_df$Category)
gene2cat <- setNames(cat_df$Category, toupper(cat_df$Gene))

# --- Top-K table with categories -----------------------------------------------------------
topk <- 50L
df_imp <- head(consensus, topk)
df_imp$Category <- gene2cat[toupper(df_imp$Gene)]
df_imp$Category[is.na(df_imp$Category)] <- "Other"
df_imp$Category <- factor(df_imp$Category, levels = c(names(reg_sets), "Other"))

# --- Paper-style palette -------------------------------------------------------------------
pal <- c("Mitophagy Machinery"="#d62728",    # red
         "Fission/Fusion"     ="#ff7f0e",    # orange
         "mtUPR/Proteostasis" ="#9467bd",    # purple
         "Autophagy/Lysosome" ="#1f77b4",    # blue
         "TFs/Signalling"     ="#2ca02c",    # green
         "ERphagy"            ="#000000",    # black
         "Antioxidant"        ="#7f7f7f",    # gray
         "Other"              ="#bdbdbd")    # fallback

# --- Plot (Mean with SD error bars, all in % units) ----------------------------------------
outdir <- plot_out; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

p_imp <- ggplot(df_imp, aes(x = reorder(Gene, -Mean), y = Mean, fill = Category)) +
  geom_col(width = 0.82) +
  geom_errorbar(aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD), width = 0.3) +
  scale_fill_manual(values = pal, drop = FALSE) +
  labs(x = NULL, y = "% MSE increase (perm.)",
       title = "Top regulators by permutation importance (ALL beta)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p_imp
ggsave(file.path(outdir, "11_importance_top50_all_beta.pdf"),
       p_imp, width = 8.5, height = 6.6, limitsize = FALSE)

# --- Save Seurat back (no changes here; just keep the same file handy) ---------------------
qs::qsave(combined, file = file.path(seu_dir, "combined_with_MQI_all_beta.qs"))

message("\n✅ Done. Plot saved to: ", outdir, "\n")

## ============== From EMPTY env: make 3D PDP PDFs ==============

rm(list = ls()); invisible(gc())

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(torch)
  if (!requireNamespace("plot3D", quietly = TRUE)) install.packages("plot3D")
  library(plot3D)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---- Paths (edit if your files are elsewhere) ----
seu_dir  <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects)"
base_dir <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\Perma_analysis)"
out_dir  <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\KANN\qc_images)"

seu_with_mqi <- file.path(seu_dir, "combined_with_MQI_all_beta.qs")
seu_fallback <- file.path(seu_dir, "combined.qs")
state_pt     <- file.path(base_dir, "fit_all_state.pt")
scaler_qs    <- file.path(base_dir, "input_scaler.qs")

stopifnot(file.exists(state_pt), file.exists(scaler_qs))
seu_path <- if (file.exists(seu_with_mqi)) seu_with_mqi else seu_fallback
stopifnot(file.exists(seu_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load Seurat & scaler; build X exactly as in training ----
combined <- qread(seu_path); stopifnot(inherits(combined, "Seurat"))
DefaultAssay(combined) <- "RNA"

sc   <- qread(scaler_qs)
feat <- sc$feat %||% names(sc$mu)                  # exact feature order used in training
mu   <- as.numeric(sc$mu); names(mu) <- feat
sdv  <- as.numeric(sc$sd); names(sdv) <- feat

expr <- GetAssayData(combined, assay = DefaultAssay(combined), slot = "data")
stopifnot(all(feat %in% rownames(expr)))
cells <- colnames(combined)

X <- t(as.matrix(expr[feat, cells, drop = FALSE]))  # n x p (already log-normalized in Seurat)
X <- sweep(sweep(X, 2L, mu, "-"), 2L, sdv, "/")     # z-score with training scaler

# ---- KAN architecture (must match training), then load state ----
KANSpline1D <- nn_module(
  "KANSpline1D",
  initialize = function(n_bins = 10, xmin = -3, xmax = 3) {
    self$n_bins <- n_bins
    self$register_buffer("knots", torch_linspace(xmin, xmax, n_bins))
    self$coef <- nn_parameter(torch_zeros(n_bins))
  },
  forward = function(x) {
    x_exp <- x$unsqueeze(2)
    w <- (self$knots[2] - self$knots[1])
    phi <- (1 - (x_exp - self$knots$unsqueeze(1))$abs() / w)$clamp(min = 0)
    (phi * self$coef)$sum(dim = 2)
  }
)

KANLayer <- nn_module(
  "KANLayer",
  initialize = function(in_features, out_features, n_bins = 10, xmin = -3, xmax = 3) {
    self$in_features <- in_features
    self$out_features <- out_features
    self$splines <- nn_module_list(lapply(1:(out_features*in_features), function(i)
      KANSpline1D(n_bins = n_bins, xmin = xmin, xmax = xmax)
    ))
    self$linear <- nn_linear(in_features, out_features)
    self$bias   <- nn_parameter(torch_zeros(out_features))
  },
  forward = function(x) {
    N <- x$size(1); P <- self$in_features
    outs <- vector("list", self$out_features); idx <- 1L
    for (j in 1:self$out_features) {
      acc <- torch_zeros(N)
      for (i in 1:P) { acc <- acc + self$splines[[idx]](x[, i]); idx <- idx + 1L }
      outs[[j]] <- acc
    }
    S <- torch_stack(outs, dim = 2)
    S + self$linear(x) + self$bias$unsqueeze(1)
  }
)

KANNet <- nn_module(
  "KANNet",
  initialize = function(in_features, hidden, out_features = 1, n_bins = 10) {
    self$layer1 <- KANLayer(in_features, hidden, n_bins = n_bins)
    self$act1   <- nn_relu()
    self$layer2 <- KANLayer(hidden, out_features, n_bins = n_bins)
  },
  forward = function(x) self$layer2(self$act1(self$layer1(x)))
)

# infer dims from state dict
state   <- torch_load(state_pt)
hidden  <- as.integer(state[["layer1.bias"]]$size()[[1]])
n_bins  <- as.integer(state[["layer1.splines.0.coef"]]$size()[[1]])
model   <- KANNet$new(in_features = ncol(X), hidden = hidden, out_features = 1, n_bins = n_bins)
model$load_state_dict(state); model$eval()

# ---- Minimal PDP helpers (no retraining, no extra deps) ----
.pdp_predict <- function(mod, X_mat, batch = 32768L) {
  if (!is.matrix(X_mat)) X_mat <- as.matrix(X_mat)
  n <- nrow(X_mat); out <- numeric(n); i <- 1L
  while (i <= n) {
    j <- min(i + batch - 1L, n)
    out[i:j] <- as.numeric(mod$forward(torch_tensor(X_mat[i:j, , drop = FALSE], dtype = torch_float())))
    i <- j + 1L
  }
  out
}

pdp2d <- function(g1, g2, out_dir,
                  qrange = c(0.02, 0.98), grid_n = 70,
                  center_baseline = TRUE, symmetric_z = TRUE) {
  
  stopifnot(g1 %in% colnames(X), g2 %in% colnames(X))
  
  # use observed quantile range to avoid “ramp”
  xr <- as.numeric(quantile(X[, g1], qrange, na.rm = TRUE))
  yr <- as.numeric(quantile(X[, g2], qrange, na.rm = TRUE))
  gx <- seq(xr[1], xr[2], length.out = grid_n)
  gy <- seq(yr[1], yr[2], length.out = grid_n)
  
  base <- colMeans(X, na.rm = TRUE)
  G <- expand.grid(v1 = gx, v2 = gy)
  Xg <- matrix(rep(base, each = nrow(G)), nrow = nrow(G), byrow = FALSE)
  colnames(Xg) <- colnames(X)
  Xg[, g1] <- G$v1; Xg[, g2] <- G$v2
  
  z <- .pdp_predict(model, Xg)
  Z <- matrix(z, nrow = length(gx), ncol = length(gy), byrow = FALSE)
  
  if (center_baseline) {
    z0 <- .pdp_predict(model, t(base))
    Z  <- Z - z0
  }
  
  zlim <- range(Z, finite = TRUE)
  if (symmetric_z) { m <- max(abs(zlim)); zlim <- c(-m, m) }
  
  f_surface <- file.path(out_dir, sprintf("pdp_%s_x_%s_surface.pdf", g1, g2))
  f_mesh    <- file.path(out_dir, sprintf("pdp_%s_x_%s_mesh.pdf", g1, g2))
  
  grDevices::cairo_pdf(f_surface, width = 7.2, height = 5.4, family = "sans")
  persp3D(x = gx, y = gy, z = Z, theta = 40, phi = 25, expand = 0.85,
          shade = 0.35, border = "black", lwd = 0.2,
          colkey = list(length = 0.5, width = 0.8, dist = -0.02,
                        cex.clab = 0.9, cex.axis = 0.8),
          ticktype = "detailed", cex.axis = 0.9,
          xlab = paste0(g1, " (z)"), ylab = paste0(g2, " (z)"), zlab = "Predicted MQI",
          zlim = zlim)
  dev.off()
  
  grDevices::cairo_pdf(f_mesh, width = 7.2, height = 5.4, family = "sans")
  persp3D(x = gx, y = gy, z = Z, theta = 40, phi = 25, expand = 0.85,
          facets = NA, border = "black", lwd = 0.6,
          ticktype = "detailed", cex.axis = 0.9,
          xlab = paste0(g1, " (z)"), ylab = paste0(g2, " (z)"), zlab = "Predicted MQI",
          zlim = zlim)
  dev.off()
  
  message("Saved:\n- ", f_surface, "\n- ", f_mesh)
}

# ---- Make your two plots (edit genes if you want different pairs) ----
pdp2d("SQSTM1", "HSPE1", out_dir)
pdp2d("PINK1",  "BNIP3", out_dir)

## =============================================================
# Load combined file:
combined <- qread(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined_with_MQI_all_beta.qs)")

library(Seurat)
library(ggplot2)
library(patchwork)

# 0) Ensure pooled treatment + collapsed diabetes status exist
md <- combined@meta.data
head(md)

# pooled treatment: Untreated + EtOH -> "Untreated"; any DHT* -> "DHT"
if (!"treat_simple" %in% names(md)) {
  tr <- as.character(md$treatment)
  tr <- ifelse(is.na(tr), NA, tr)
  tr_simple <- ifelse(tr %in% c("Untreated","EtOH"), "Untreated",
                      ifelse(grepl("^DHT", tr), "DHT", NA))
  md$treat_simple <- factor(tr_simple, levels = c("Untreated","DHT"))
}

# diabetes_status2: NA -> ND
if (!"diabetes_status2" %in% names(md)) {
  ds <- as.character(md$diabetes_status)
  ds[is.na(ds)] <- "ND"
  md$diabetes_status2 <- factor(ds, levels = c("ND","T2D"))
}

combined@meta.data <- md

# 1) Helper: EXCLUDE DHT by using treat_simple == "Untreated"
plot_umap_split <- function(obj, sex_val, feature = "MQI_pred", limits = c(0,1)) {
  keep_cells <- WhichCells(obj, expression = Sex == sex_val & treat_simple == "Untreated")
  sub <- subset(obj, cells = keep_cells)
  
  p_list <- FeaturePlot(
    sub,
    features = feature,
    reduction = "umap",
    split.by  = "diabetes_status2",
    combine   = FALSE
  )
  p_list <- lapply(p_list, function(p)
    p + scale_color_gradientn(
      colours = c("dodgerblue4","white","firebrick3"),
      limits  = limits,
      oob     = scales::squish
    ) + labs(title = paste(feature, "-", sex_val))
  )
  wrap_plots(p_list, ncol = 2)  # vertical: ND on top, T2D below
}

# 2) Build plots (Untreated pool only; DHT excluded)
p_female <- plot_umap_split(combined, sex_val = "F", feature = "MQI_pred", limits = c(0,1))
p_male   <- plot_umap_split(combined, sex_val = "M", feature = "MQI_pred", limits = c(0,1))

# 3) Stack Female over Male and save
p_all <- p_female / p_male
#ggsave("MQI_pred_UntreatedOnly_ND_vs_T2D_Female_over_Male.pdf", p_all, width = 9, height = 10)
p_all


# Deps
library(Seurat)
library(ggplot2)

# Make sure you have these metadata prepared once:
# - treat_simple: "Untreated" (pools Untreated+EtOH) vs "DHT"
# - diabetes_status2: "ND" (NA pooled to ND) vs "T2D"

# ---------- Function: ND vs T2D violins, Untreated only, with p + FC ----------
plot_nd_t2d_violins <- function(
    obj,
    feature = "MQI_pred",
    beta_col = "beta_cluster",
    sex_filter = NULL,                  # e.g., "F", "M", or NULL for all sexes
    outfile = "MQI_pred_Untreated_ND_vs_T2D_violins.pdf",
    fill_colors = c(ND = "#1976D2", T2D = "#C2185B"),  # customize as you like
    alpha_violin = 0.85,
    add_box = TRUE
) {
  stopifnot(all(c(feature, beta_col, "treat_simple", "diabetes_status2") %in% colnames(obj@meta.data)))
  
  # Filter: Untreated pool only (EtOH + Untreated), optionally by sex
  keep <- WhichCells(obj, expression = treat_simple == "Untreated")
  sub  <- subset(obj, cells = keep)
  if (!is.null(sex_filter)) {
    keep <- WhichCells(sub, expression = Sex == sex_filter)
    sub  <- subset(sub, cells = keep)
  }
  
  # Pull tidy data
  df <- FetchData(sub, vars = c(feature, beta_col, "diabetes_status2"))
  names(df) <- c("val", "beta", "dx")
  df <- df[complete.cases(df), , drop = FALSE]
  
  # Keep only ND/T2D present; order β clusters (drop empties)
  dx_levels <- c("ND","T2D")
  df$dx   <- factor(df$dx, levels = intersect(dx_levels, unique(df$dx)))
  beta_lv <- c("β1","β2","β3","β4")
  df$beta <- factor(df$beta, levels = intersect(beta_lv, unique(df$beta)))
  
  if (nrow(df) == 0L || nlevels(df$dx) < 2L) {
    message("No usable cells for ND vs T2D after filtering.")
    return(invisible(NULL))
  }
  
  # For boxplots, require >=2 per (beta, dx)
  tab <- as.data.frame(table(df$beta, df$dx)); names(tab) <- c("beta","dx","n")
  ok  <- tab[tab$n >= 2, c("beta","dx")]
  df_box <- if (nrow(ok)) merge(df, ok, by = c("beta","dx")) else df[0, , drop = FALSE]
  
  # Compute per-β stats (Wilcoxon p, median FC and Δ)
  by_beta <- split(df, df$beta, drop = TRUE)
  ymax <- tapply(df$val, df$beta, max, na.rm = TRUE)
  yrng <- range(df$val, na.rm = TRUE); ypad <- diff(yrng) * 0.06
  
  ann <- lapply(names(by_beta), function(b) {
    d <- by_beta[[b]]
    n_nd  <- sum(d$dx == "ND");  n_t2d <- sum(d$dx == "T2D")
    if (n_nd >= 2 && n_t2d >= 2) {
      med_nd  <- median(d$val[d$dx == "ND"],  na.rm = TRUE)
      med_t2d <- median(d$val[d$dx == "T2D"], na.rm = TRUE)
      fc      <- if (is.finite(med_nd) && med_nd != 0) med_t2d / med_nd else NA_real_
      diffm   <- med_t2d - med_nd
      pval <- tryCatch(
        wilcox.test(d$val[d$dx == "ND"], d$val[d$dx == "T2D"])$p.value,
        error = function(e) NA_real_
      )
      data.frame(
        beta = b,
        y    = ymax[[b]] + ypad,
        label = sprintf("p=%.2e | FC=%.2f | Δ=%.3f", pval, fc, diffm),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
  ann <- do.call(rbind, ann)
  
  # Plot
  pd <- position_dodge(width = 0.75)
  p <- ggplot(df, aes(x = beta, y = val,
                      fill = dx,
                      group = interaction(beta, dx))) +
    geom_violin(width = 0.7, trim = TRUE, alpha = alpha_violin, position = pd) +
    {if (add_box) geom_boxplot(data = df_box, width = 0.18, position = pd,
                               outlier.shape = NA, color = "black")} +
    scale_fill_manual(values = fill_colors[levels(df$dx)], drop = FALSE) +
    theme_classic(base_size = 12) +
    labs(title = paste(feature, "— Untreated (ND vs T2D)"),
         x = "β-cluster", y = feature, fill = "Status")
  
  if (!is.null(ann) && nrow(ann) > 0) {
    p <- p + annotate("text", x = ann$beta, y = ann$y, label = ann$label, size = 3)
  }
  
  # Save
  outdir <- "figures_violin"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outpath <- file.path(outdir, outfile)
  if (capabilities("cairo")) {
    ggsave(outpath, p, width = 9, height = 5.5, device = grDevices::cairo_pdf)
  } else {
    ggsave(outpath, p, width = 9, height = 5.5)
  }
  p
}

# All sexes pooled (Untreated pool only)
p_all <- plot_nd_t2d_violins(
  combined,
  feature = "MQI_pred",
  outfile = "MQI_pred_Untreated_ND_vs_T2D_allSex.pdf",
  fill_colors = c(ND = "#1976D2", T2D = "#C2185B")
)

# Females only
p_f <- plot_nd_t2d_violins(
  combined,
  feature = "MQI_pred",
  sex_filter = "F",
  outfile = "MQI_pred_Untreated_ND_vs_T2D_Female.pdf",
  fill_colors = c(ND = "#1976D2", T2D = "#C2185B")
)

# Males only
p_m <- plot_nd_t2d_violins(
  combined,
  feature = "MQI_pred",
  sex_filter = "M",
  outfile = "MQI_pred_Untreated_ND_vs_T2D_Male.pdf",
  fill_colors = c(ND = "#1976D2", T2D = "#C2185B")
)

p_all / p_f / p_m
p_f
p_m


# 0) Libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# 1) Parameters
feature   <- "MQI_pred"
beta_col  <- "beta_cluster"
dx_col    <- "diabetes_status2"
sex_col   <- "Sex"
threshold <- 1
fill_cols <- c("High (>1)" = "#C2185B", "Low (<1)" = "#1976D2")

# 2) Filter: Untreated only (EtOH+Untreated if desired)
keep <- WhichCells(combined, expression = treat_simple == "Untreated")
sub  <- subset(combined, cells = keep)

# 3) Extract data
df <- FetchData(sub, vars = c(feature, beta_col, dx_col, sex_col))
names(df) <- c("val", "beta", "dx", "sex")
df <- df[complete.cases(df), , drop = FALSE]

# 4) Keep ND/T2D only, order factors
df <- df[df$dx %in% c("ND","T2D"), ]
df$beta <- factor(df$beta, levels = c("β1","β2","β3","β4"))
df$dx   <- factor(df$dx,   levels = c("ND","T2D"))
df$sex  <- factor(df$sex,  levels = c("M","F"))

# 5) Categorize High vs Low
df <- df[df$val != threshold, ]
df$cat <- ifelse(df$val > threshold, "High (>1)", "Low (<1)")

# 6) Tabulate counts → percentages  (force dplyr::count)
tab <- df %>%
  dplyr::count(beta, dx, sex, cat, name = "n") %>%
  dplyr::group_by(beta, dx, sex) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()

# 7) Plot stacked % barplot, faceted by sex × dx
p <- ggplot(tab, aes(x = beta, y = pct, fill = cat)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  facet_grid(sex ~ dx) +
  scale_fill_manual(values = fill_cols) +
  labs(title = "Percentage of High vs Low MQI_pred cells (Untreated)",
       x = "β-cluster", y = "% of cells", fill = "MQI_pred") +
  theme_classic(base_size = 13)

# 8) Save
ggsave("figures_pctbar/MQI_pred_highlow_pct_Untreated_byBeta_NDvsT2D_bySex.pdf",
       p, width = 10, height = 6)

p

# DGE
# ---- Minimal DGE: β3 MQI High vs Low ----

# 1) Subset β3 females
cells_F <- WhichCells(combined, expression = beta_cluster == "β3" & Sex == "F" & treat_simple == "Untreated")
sub_F <- subset(combined, cells = cells_F)
sub_F$MQI_group <- ifelse(sub_F$MQI_pred > 1, "High", "Low")
Idents(sub_F) <- "MQI_group"

# Female DGE
dge_F <- FindMarkers(sub_F, ident.1 = "High", ident.2 = "Low",
                     test.use = "wilcox", logfc.threshold = 0.263034,
                     min.pct = 0.1, only.pos = FALSE)
dge_F <- tibble::rownames_to_column(dge_F, "gene")
head(dge_F)

# 2) Subset β3 males
cells_M <- WhichCells(combined, expression = beta_cluster == "β3" & Sex == "M" & treat_simple == "Untreated")
sub_M <- subset(combined, cells = cells_M)
sub_M$MQI_group <- ifelse(sub_M$MQI_pred > 1, "High", "Low")
Idents(sub_M) <- "MQI_group"

# Male DGE
dge_M <- FindMarkers(sub_M, ident.1 = "High", ident.2 = "Low",
                     test.use = "wilcox", logfc.threshold = 0.263034,
                     min.pct = 0.1, only.pos = FALSE)
dge_M <- tibble::rownames_to_column(dge_M, "gene")
head(dge_M)

# 3) Save results
#dir.create("DGE_beta3_MQI", showWarnings = FALSE)
readr::write_csv(dge_F, "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/MFI/beta3_F_MQI_High_vs_Low.csv")
readr::write_csv(dge_M, "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/MFI/beta3_M_MQI_High_vs_Low.csv")

#ORA Anlalysis
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(gprofiler2)
})

# Helper to flatten list columns from gprofiler
.flatten_lists <- function(df) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        where(is.list),
        ~ vapply(., function(x) paste(as.character(x), collapse = ";"), character(1))
      )
    )
}

# Output directory
outdir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/MFI/ORA_beta3_MQI"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Wrapper function to run ORA on one DEG file
run_ORA <- function(csv_file, tag){
  deg <- read_csv(csv_file)
  sig <- deg %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.263034)
  
  genes_up   <- sig %>% filter(avg_log2FC > 0) %>% pull(gene)
  genes_down <- sig %>% filter(avg_log2FC < 0) %>% pull(gene)
  
  # Up in High
  if (length(genes_up) >= 5) {
    ora_up <- tryCatch(
      gost(query = genes_up, organism = "hsapiens",
           correction_method = "fdr", user_threshold = 0.05,
           significant = TRUE, sources = "GO", evcodes = TRUE),
      error = function(e) { message("gprofiler2 up-error: ", e$message); NULL }
    )
    if (!is.null(ora_up) && !is.null(ora_up$result)) {
      ora_up_out <- .flatten_lists(ora_up$result)
      write.csv(ora_up_out,
                file.path(outdir, paste0("ORA_", tag, "_Up_in_High.csv")),
                row.names = FALSE)
    } else {
      message(tag, ": no significant terms for Up-in-High.")
    }
  } else {
    message(tag, ": <5 up-genes; skipping ORA_up.")
  }
  
  # Down in High
  if (length(genes_down) >= 5) {
    ora_down <- tryCatch(
      gost(query = genes_down, organism = "hsapiens",
           correction_method = "fdr", user_threshold = 0.05,
           significant = TRUE, sources = "GO", evcodes = TRUE),
      error = function(e) { message("gprofiler2 down-error: ", e$message); NULL }
    )
    if (!is.null(ora_down) && !is.null(ora_down$result)) {
      ora_down_out <- .flatten_lists(ora_down$result)
      write.csv(ora_down_out,
                file.path(outdir, paste0("ORA_", tag, "_Down_in_High.csv")),
                row.names = FALSE)
    } else {
      message(tag, ": no significant terms for Down-in-High.")
    }
  } else {
    message(tag, ": <5 down-genes; skipping ORA_down.")
  }
}

# Run for Female and Male β3
run_ORA("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/MFI/beta3_F_MQI_High_vs_Low.csv", "beta3_F")
run_ORA("C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/MFI/beta3_M_MQI_High_vs_Low.csv", "beta3_M")

# ---------- Sankey util (your function, unchanged) ----------
make_sankey_plot <- function(my_comparison_colors, my_pathway_colors, my_ora_dir, my_pathways) {
  library(readr); library(dplyr); library(ggalluvial); library(ggplot2); library(cowplot); library(scales)
  nukefix <- function(x) {
    x <- gsub("β", "b", x)
    x <- trimws(x)
    x <- iconv(x, to = "ASCII//TRANSLIT")
    x <- tolower(x)
    x <- gsub("_", " ", x)
    x <- gsub("\\s+", " ", x)
    x
  }
  my_comp_names <- list.files(my_ora_dir, pattern = "\\.csv$", full.names = FALSE)
  cat("Found files:\n"); print(my_comp_names)
  cat("Pathways (requested):\n"); print(my_pathways)
  
  nukefixed_pathways <- nukefix(my_pathways)
  all_links <- list()
  
  for (fname in my_comp_names) {
    comp <- tolower(gsub("\\.csv$", "", fname))
    fpath <- file.path(my_ora_dir, fname)
    if (!file.exists(fpath)) next
    df <- read_csv(fpath, show_col_types = FALSE)
    if (!"term_name" %in% colnames(df)) next
    if (!"p_value" %in% tolower(colnames(df))) stop("No p_value column in file: ", fname)
    pval_col <- colnames(df)[tolower(colnames(df)) == "p_value"][1]
    df <- df %>% filter(!is.na(.data[[pval_col]]), .data[[pval_col]] < 0.05, !is.na(term_name))
    df$Pathway <- nukefix(df$term_name)
    
    if (nrow(df) > 0) cat("In", fname, "found pathways:\n", paste(unique(df$Pathway), collapse=", "), "\n")
    
    for (i in seq_along(my_pathways)) {
      pw <- nukefixed_pathways[i]
      pretty_name <- my_pathways[i]
      match_row <- df[df$Pathway == pw,]
      if (nrow(match_row) == 0) next
      all_links[[length(all_links)+1]] <- data.frame(
        Comparison = comp,
        Pathway = pretty_name,
        Value = -log10(match_row[[pval_col]][1]) + 1,
        stringsAsFactors = FALSE
      )
    }
  }
  
  sankey_df <- dplyr::bind_rows(all_links)
  sankey_df$Comparison <- factor(nukefix(sankey_df$Comparison), levels = nukefix(names(my_comparison_colors)))
  sankey_df$Pathway <- factor(nukefix(sankey_df$Pathway), levels = nukefix(my_pathways), labels = my_pathways)
  if (nrow(sankey_df) == 0) stop("No matching pathways were found after filtering. Check pathway names and p_value thresholds.\n")
  
  p_main <- ggplot(sankey_df, aes(axis1 = Comparison, axis2 = Pathway, y = Value)) +
    geom_alluvium(aes(fill = Comparison), width = 1/12, alpha = 0.8) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_discrete(limits = c("Comparison", "Pathway"), expand = c(.1, .1)) +
    scale_fill_manual(values = c(my_comparison_colors, my_pathway_colors), breaks = names(my_comparison_colors)) +
    theme_minimal(base_size = 16) +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_text(size = 14, face = "bold"))
  
  fdr_values <- sankey_df$Value
  legend_vals <- pretty(fdr_values, n = 5); legend_vals <- sort(unique(legend_vals[legend_vals > 0]), decreasing = TRUE)
  rescale_size <- function(x, to = c(3, 10)) { rng <- range(x, na.rm = TRUE); scales::rescale(x, to = to, from = rng) }
  dot_sizes <- rescale_size(legend_vals, to = c(3, 10))
  legend_df <- data.frame(x = 1, y = seq(50, 38, length.out = length(dot_sizes)), y_label = seq(50, 38, length.out = length(dot_sizes)),
                          size = dot_sizes, label = paste0("-log10(p) + 1 = ", legend_vals))
  
  p_legend <- ggplot(legend_df) +
    geom_point(aes(x = x, y = y, size = size), shape = 21, fill = "steelblue", color = "black", stroke = 0.25) +
    geom_text(aes(x = x + 0.4, y = y_label, label = label), hjust = 0, vjust = 0.5, size = 4) +
    theme_void() + coord_cartesian(clip = "off") +
    scale_size_identity() + scale_x_continuous(limits = c(0.9, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) + theme(plot.margin = margin(5, 20, 5, 5))
  
  cowplot::plot_grid(p_main, p_legend, rel_widths = c(4.2, 1.2), nrow = 1, axis = "none", align = "none")
}

# ---------- Inputs (directory + two files) ----------
my_ora_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/ora/MFI/ORA_beta3_MQI"
# Expecting:
#   ORA_beta3_F_Up_in_High.csv
#   ORA_beta3_M_Up_in_High.csv

# Map comparisons by the file basenames (lowercased in the function)
my_comparison_colors <- c(
  "ora beta3 f up in high" = "#C2185B",  # Female (pink/red)
  "ora beta3 m up in high" = "#1976D2"   # Male (blue)
)

# ---------- Pathways (exact names you gave) ----------
# OXPHOS
oxphos <- c(
  "oxidative phosphorylation",
  "respiratory electron transport chain",
  "ATP synthesis coupled electron transport",
  "mitochondrial ATP synthesis coupled electron transport",
  "aerobic electron transport chain",
  "mitochondrial electron transport, NADH to ubiquinone",
  "mitochondrial electron transport, ubiquinol to cytochrome c",
  "mitochondrial respiratory chain complex I assembly",
  "mitochondrial respirasome assembly"
)
# Translation
translation <- c(
  "cytoplasmic translation",
  "translation",
  "translational initiation",
  "cytoplasmic translational initiation",
  "translational elongation",
  "ribosome biogenesis",
  "ribosome assembly",
  "ribosomal small subunit biogenesis",
  "ribosomal large subunit biogenesis"
)
# Proteostasis
proteostasis <- c(
  "protein folding",
  "protein refolding",
  "chaperone-mediated protein folding",
  "response to unfolded protein",
  "endoplasmic reticulum unfolded protein response",
  "ubiquitin-dependent protein catabolic process",
  "proteasomal protein catabolic process",
  "regulation of protein ubiquitination",
  "regulation of proteasomal ubiquitin-dependent protein catabolic process"
)
# Mitophagy
mitophagy <- c(
  "mitophagy",
  "autophagy of mitochondrion",
  "macroautophagy",
  "process utilizing autophagic mechanism",
  "positive regulation of autophagy of mitochondrion",
  "regulation of autophagy of mitochondrion",
  "mitochondrial membrane organization",
  "regulation of mitochondrion organization"
)
# Ion transport
ion_transport <- c(
  "proton transmembrane transport",
  "electron transport coupled proton transport",
  "inorganic cation transmembrane transport",
  "monoatomic cation transmembrane transport",
  "monoatomic ion transmembrane transport",
  "inorganic ion transmembrane transport",
  "lysosomal lumen acidification",
  "vacuolar acidification",
  "calcium ion export across plasma membrane"
)

my_pathways <- c(oxphos, translation, proteostasis, mitophagy, ion_transport)

# ---------- Pathway colors by category ----------
my_pathway_colors <- c(
  # OXPHOS (Orange)
  setNames(rep("#D35400", length(oxphos)), oxphos),
  # Translation (Purple)
  setNames(rep("#7B1FA2", length(translation)), translation),
  # Proteostasis / UPR (Teal-ish green to stand out from comparison colors)
  setNames(rep("#2E8B57", length(proteostasis)), proteostasis),
  # Mitophagy (Brown)
  setNames(rep("#8D6748", length(mitophagy)), mitophagy),
  # Ion transport (Blue)
  setNames(rep("#2980B9", length(ion_transport)), ion_transport)
)

# ---------- RUN ----------
make_sankey_plot(
  my_comparison_colors = my_comparison_colors,
  my_pathway_colors = my_pathway_colors,
  my_ora_dir = my_ora_dir,
  my_pathways = my_pathways
)


## ================================
## β1 (Untreated): HbA1c vs % cells
## 4 lines = (M/F) × (High/Low MQI)
## With R and p annotations per line
## ================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(rlang)
})

## Avoid plyr/dplyr conflicts
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

## -------- CONFIG (edit paths/cols if needed) --------
csv_path     <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\Donor_Summary_186.csv)"
hba1c_col_in <- "hba1c"          # HbA1c column name in CSV
obj          <- combined         # your Seurat object
beta_col     <- "beta_cluster"   # β1–β4
sex_col      <- "Sex"            # "M"/"F"
treat_col    <- "treat_simple"   # "Untreated"
donor_col    <- "Library"
feature      <- "MQI_pred"       # score in meta.data
beta_target  <- "β1"
threshold    <- 1                # split at 1, exclude exactly 1

# Colors: your current blue/red, lighter shade for "Low"
col_M_high <- "#1976D2"; col_M_low <- "#90CAF9"
col_F_high <- "#C2185B"; col_F_low <- "#F48FB1"

## 1) HbA1c lookup (CSV + manual overrides)
clinical_df <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(
    Library = toupper(str_trim(donor_ID)),
    hba1c   = as.numeric(.data[[hba1c_col_in]])
  ) %>%
  select(Library, hba1c) %>%
  filter(!is.na(Library), !is.na(hba1c)) %>%
  mutate(.src = "csv")

hba1c_missing <- tibble::tribble(
  ~Library,        ~hba1c,
  "HP2022801",     5.5,
  "HP2024001",     5.4,
  "SAMN15877725",  5.6,
  "HP2031401",     5.4,
  "HP2105501",     5.6,
  "HP2106201",     5.3,
  "HP2107001",     5.1,
  "HP2107901",     5.2,
  "HP2108601",     5.1,
  "HP2108901",     5.9,
  "HP2110001",     5.5,
  "HP2123201",     5.3,
  "HP2132801",     5.5,
  "HP2202101",     5.5,
  "HP2121601",     5.8
) %>%
  mutate(
    Library = toupper(str_trim(Library)),
    hba1c   = as.numeric(hba1c),
    .src    = "manual"
  )

hba1c_full <- bind_rows(clinical_df, hba1c_missing) %>%
  arrange(Library, .src) %>%
  group_by(Library) %>%
  slice_tail(n = 1) %>%            # manual overrides csv if duplicated
  ungroup() %>%
  select(Library, hba1c)

## 2) β1 Untreated cells → label High/Low per cell
stopifnot(all(c(feature, beta_col, sex_col, treat_col, donor_col) %in% colnames(obj@meta.data)))
sub_beta1 <- subset(obj, subset = !!sym(treat_col) == "Untreated" & !!sym(beta_col) == beta_target)

df_cells <- sub_beta1@meta.data %>%
  as.data.frame() %>%
  transmute(
    Library = toupper(str_trim(.data[[donor_col]])),
    Sex     = .data[[sex_col]],
    val     = .data[[feature]]
  ) %>%
  filter(!is.na(Library), !is.na(Sex), !is.na(val), val != threshold) %>%  # exclude exactly 1.0
  mutate(
    Sex   = factor(Sex, levels = c("M","F")),
    Group = ifelse(val > threshold, "High (>1)", "Low (<1)")
  )

if (nrow(df_cells) == 0L) stop("No usable β1 Untreated cells after filtering.")

## 3) Per-donor × Sex × Group: % of β1 cells in that group, join HbA1c
by_donor <- df_cells %>%
  dplyr::count(Library, Sex, Group, name = "n") %>%
  dplyr::group_by(Library, Sex) %>%
  dplyr::mutate(
    n_total   = sum(n),
    pct_group = 100 * n / n_total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(hba1c_full, by = "Library") %>%
  dplyr::filter(!is.na(hba1c))

if (nrow(by_donor) == 0L) stop("No donors with HbA1c after join.")

## 4) Correlation (R, p) per line = (Sex × Group)
# Build a single key that matches legend labels
by_donor$SexGroup <- interaction(by_donor$Sex, by_donor$Group, sep = " · ")

corr_tbl <- by_donor %>%
  dplyr::group_by(Sex, Group, SexGroup) %>%
  dplyr::summarise(
    n = dplyr::n(),
    r = suppressWarnings(if (n >= 3) cor(hba1c, pct_group, method = "pearson") else NA_real_),
    p = tryCatch(if (n >= 3) cor.test(hba1c, pct_group, method = "pearson")$p.value else NA_real_,
                 error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p.label   = dplyr::case_when(
      is.na(p)            ~ "NA",
      p < 2.2e-16         ~ "<2.2e-16",
      TRUE                ~ as.character(signif(p, 2))
    ),
    r.label   = ifelse(is.na(r), "NA", sprintf("%.2f", r)),
    cor.label = paste0("R=", r.label, "\np=", p.label)
  )

# Place each label slightly above its group's max Y, at left side of plot
pos_tbl <- by_donor %>%
  dplyr::group_by(SexGroup) %>%
  dplyr::summarise(
    y_max = max(pct_group, na.rm = TRUE),
    .groups = "drop"
  )
x_min <- min(by_donor$hba1c, na.rm = TRUE)

corr_tbl <- corr_tbl %>%
  dplyr::left_join(pos_tbl, by = "SexGroup") %>%
  dplyr::mutate(
    y.label = y_max + 0.05 * abs(y_max),
    x.label = x_min
  )

## 5) Plot: 4 lines (M/F × High/Low), with R/p annotations
pal <- c(
  "M · High (>1)" = col_M_high,
  "M · Low (<1)"  = col_M_low,
  "F · High (>1)" = col_F_high,
  "F · Low (<1)"  = col_F_low
)

p <- ggplot(by_donor, aes(x = hba1c, y = pct_group, color = SexGroup)) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, size = 1.05) +
  # correlation labels (colored to match the line)
  geom_text(
    data = corr_tbl,
    aes(x = x.label, y = y.label, label = cor.label, color = SexGroup),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1, size = 4, fontface = "bold"
  ) +
  scale_color_manual(values = pal, name = "Sex × MQI group") +
  labs(
    title = "β1 (Untreated): HbA1c vs % of cells in MQI_pred High/Low",
    x = "HbA1c (%)",
    y = "% of β1 cells per donor"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
p
# Save (optional)
dir.create("figures_regression", showWarnings = FALSE, recursive = TRUE)
ggsave("figures_regression/HbA1c_vs_pctHighLow_beta1_Untreated_4lines_with_stats.pdf",
       p, width = 8.5, height = 6.2)

p


###############################################
## DGE + ORA for β1 (Untreated)
## Compare β1 cells with MQI_pred > 1 (High) vs < 1 (Low)
## - Subset + group
## - DGE with Seurat::FindMarkers (Wilcoxon)
## - Volcano plot (optional labels)
## - ORA via gprofiler2 (handles list-columns)
## - ORA via clusterProfiler (offline option)
## - All outputs saved to an output folder
###############################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

## ------------------------ CONFIG ------------------------
obj           <- combined        # your Seurat object
treat_col     <- "treat_simple"  # column that marks 'Untreated'
beta_col      <- "beta_cluster"  # β1–β4 annotation
feature_col   <- "MQI_pred"      # metric to split High/Low
beta_target   <- "β1"            # focus on β1 only
logfc_min     <- 0.25            # ~1.2x fold-change
min_pct       <- 0.10            # min pct cells expressing a gene in a group
p_adj_cutoff  <- 0.05            # FDR cutoff
# Use forward slashes on Windows to avoid \U unicode escapes
outdir        <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/MFI"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------------- (Optional) pick assay explicitly -------
# DefaultAssay(obj) <- "RNA"

## ------------------ 1) SUBSET & GROUPING -----------------
# Keep β1 + Untreated only
beta1_untreated <- subset(
  obj,
  subset = !!rlang::sym(treat_col) == "Untreated" & !!rlang::sym(beta_col) == beta_target
)

# Create a High/Low label by MQI_pred, exclude exactly 1.0 and NAs
vals <- beta1_untreated@meta.data[[feature_col]]
grp  <- ifelse(vals > 1, "High", ifelse(vals < 1, "Low", NA))

# Add to meta.data safely (ensures cell-order alignment)
beta1_untreated <- AddMetaData(
  beta1_untreated,
  metadata = data.frame(MQI_group = grp, row.names = colnames(beta1_untreated))
)

# Drop cells with NA group explicitly (avoids expression parsing issues)
keep_cells <- rownames(beta1_untreated@meta.data)[!is.na(beta1_untreated@meta.data$MQI_group)]
beta1_untreated <- subset(beta1_untreated, cells = keep_cells)
beta1_untreated$MQI_group <- factor(beta1_untreated$MQI_group, levels = c("Low","High"))

# Quick sanity: counts per group; stop early if too few cells
cat("Cell counts — High/Low:\n"); print(table(beta1_untreated$MQI_group))
if (any(table(beta1_untreated$MQI_group) < 5)) {
  stop("One of the groups has <5 cells; consider relaxing filters or pooling.")
}

# Set identities for DGE
Idents(beta1_untreated) <- "MQI_group"

## ------------------ 2) DGE: High vs Low ------------------
# Wilcoxon test (Seurat default). Consider MAST if you need covariates.
dge_beta1 <- FindMarkers(
  beta1_untreated,
  ident.1 = "High",
  ident.2 = "Low",
  test.use = "wilcox",
  min.pct = min_pct,
  logfc.threshold = logfc_min
)

# Tidy table: gene symbol column + sort by FDR then by effect size
dge_beta1 <- dge_beta1 %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj, desc(avg_log2FC))

# Save DGE
write.csv(dge_beta1, file.path(outdir, "DGE_beta1_Untreated_MQI_High_vs_Low.csv"), row.names = FALSE)

# Print a quick peek
cat("\nTop 10 DGE genes:\n"); print(head(dge_beta1, 10))

## --------------- 3) SPLIT UP/DOWN GENE SETS --------------
dge_sig <- dge_beta1 %>%
  filter(!is.na(p_val_adj), p_val_adj < p_adj_cutoff)

genes_up   <- dge_sig %>% filter(avg_log2FC > 0) %>% pull(gene)
genes_down <- dge_sig %>% filter(avg_log2FC < 0) %>% pull(gene)

cat("\nSig genes at FDR <", p_adj_cutoff, ":\n",
    "Up:", length(genes_up), " | Down:", length(genes_down), "\n", sep = "")

writeLines(genes_up,   file.path(outdir, "genes_up_in_High.txt"))
writeLines(genes_down, file.path(outdir, "genes_down_in_High.txt"))

## ---------------- 4) VOLCANO PLOT (quick) ----------------
volc_df <- dge_beta1 %>%
  mutate(
    neglog10FDR = -log10(p_val_adj + 1e-300),
    sig = case_when(
      p_val_adj < p_adj_cutoff & avg_log2FC > 0 ~ "Up in High",
      p_val_adj < p_adj_cutoff & avg_log2FC < 0 ~ "Down in High",
      TRUE ~ "NS"
    )
  )

p_volcano <- ggplot(volc_df, aes(x = avg_log2FC, y = neglog10FDR, color = sig)) +
  geom_point(alpha = 0.7, size = 1.3) +
  scale_color_manual(values = c("Up in High" = "#C2185B", "Down in High" = "#1976D2", "NS" = "grey70")) +
  geom_vline(xintercept = c(-logfc_min, logfc_min), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_adj_cutoff), linetype = "dashed") +
  labs(
    title = "β1 (Untreated): High (>1) vs Low (<1) — Volcano",
    x = "avg_log2FC (High vs Low)",
    y = "-log10(FDR)",
    color = ""
  ) +
  theme_classic(base_size = 13)
p_volcano +
       geom_text_repel(
         data = subset(volc_df, sig != "NS"),
         aes(label = gene),
         size = 2.5,
         max.overlaps = Inf
       )
# Save volcano
ggsave(file.path(outdir, "Volcano_beta1_Untreated_MQI_High_vs_Low.pdf"),
       p_volcano, width = 6.5, height = 5.3)

## (Optional) label points with ggrepel
# if (requireNamespace("ggrepel", quietly = TRUE)) {
#   library(ggrepel)
#   p_volcano_lbl <- p_volcano +
#     geom_text_repel(
#       data = subset(volc_df, sig != "NS"),
#       aes(label = gene),
#       size = 2.5,
#       max.overlaps = Inf
#     )
#   ggsave(file.path(outdir, "Volcano_beta1_Untreated_MQI_High_vs_Low_labeled.pdf"),
#          p_volcano_lbl, width = 6.5, height = 5.3)
# }

## -------------------- 5) ORA: gprofiler2 ------------------
# Uses g:Profiler web API (internet required).
# NOTE: gprofiler2 returns list-columns; we flatten them before writing.

# helper: flatten any list-columns to semicolon-joined strings
.flatten_lists <- function(df) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        where(is.list),
        ~ vapply(., function(x) paste(as.character(x), collapse = ";"), character(1))
      )
    )
}

if (requireNamespace("gprofiler2", quietly = TRUE)) {
  library(gprofiler2)
  
  # Up in High
  if (length(genes_up) >= 5) {
    ora_up <- tryCatch(
      gost(query = genes_up, organism = "hsapiens",
           correction_method = "fdr", user_threshold = 0.05, significant = TRUE, sources = "GO", evcodes = TRUE),
      error = function(e) { message("gprofiler2 up-error: ", e$message); NULL }
    )
    if (!is.null(ora_up) && !is.null(ora_up$result)) {
      ora_up_out <- .flatten_lists(ora_up$result)
      write.csv(ora_up_out, file.path(outdir, "ORA_gprof_Up_in_High.csv"), row.names = FALSE)
      # Optional static plot
      # gp <- gprofiler2::gostplot(ora_up, capped = TRUE, interactive = FALSE)
      # ggsave(file.path(outdir, "ORA_gprof_Up_in_High_plot.pdf"), gp, width = 7, height = 5)
    } else {
      message("gprofiler2: no significant terms for Up-in-High.")
    }
  } else {
    message("gprofiler2: <5 up-genes; skipping ORA_up.")
  }
  
  # Down in High
  if (length(genes_down) >= 5) {
    ora_down <- tryCatch(
      gost(query = genes_down, organism = "hsapiens",
           correction_method = "fdr", user_threshold = 0.05, significant = TRUE, sources = "GO", evcodes = TRUE),
      error = function(e) { message("gprofiler2 down-error: ", e$message); NULL }
    )
    if (!is.null(ora_down) && !is.null(ora_down$result)) {
      ora_down_out <- .flatten_lists(ora_down$result)
      write.csv(ora_down_out, file.path(outdir, "ORA_gprof_Down_in_High.csv"), row.names = FALSE)
      # Optional static plot
      # gp2 <- gprofiler2::gostplot(ora_down, capped = TRUE, interactive = FALSE)
      # ggsave(file.path(outdir, "ORA_gprof_Down_in_High_plot.pdf"), gp2, width = 7, height = 5)
    } else {
      message("gprofiler2: no significant terms for Down-in-High.")
    }
  } else {
    message("gprofiler2: <5 down-genes; skipping ORA_down.")
  }
} else {
  message("Package 'gprofiler2' not installed; skipping g:Profiler ORA.")
}

## ------------- 6) ORA: clusterProfiler (offline) ---------
# Works without internet; uses local org.Hs.eg.db.
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  if (requireNamespace("enrichplot", quietly = TRUE)) library(enrichplot)
  
  # GO BP enrichment — Up in High
  if (length(genes_up) >= 5) {
    ego_up <- enrichGO(
      gene          = genes_up,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    saveRDS(ego_up, file.path(outdir, "enrichGO_Up_in_High.rds"))
    write.csv(as.data.frame(ego_up), file.path(outdir, "enrichGO_Up_in_High.csv"), row.names = FALSE)
    if (exists("dotplot")) {
      pdf(file.path(outdir, "enrichGO_Up_in_High_dotplot.pdf"), width = 7, height = 5)
      print(dotplot(ego_up, showCategory = 15) + ggtitle("GO BP — Up in High (β1 Untreated)"))
      dev.off()
    }
  } else {
    message("clusterProfiler: <5 up-genes; skipping GO_up.")
  }
  
  # GO BP enrichment — Down in High
  if (length(genes_down) >= 5) {
    ego_down <- enrichGO(
      gene          = genes_down,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    saveRDS(ego_down, file.path(outdir, "enrichGO_Down_in_High.rds"))
    write.csv(as.data.frame(ego_down), file.path(outdir, "enrichGO_Down_in_High.csv"), row.names = FALSE)
    if (exists("dotplot")) {
      pdf(file.path(outdir, "enrichGO_Down_in_High_dotplot.pdf"), width = 7, height = 5)
      print(dotplot(ego_down, showCategory = 15) + ggtitle("GO BP — Down in High (β1 Untreated)"))
      dev.off()
    }
  } else {
    message("clusterProfiler: <5 down-genes; skipping GO_down.")
  }
  
  # (Optional) KEGG enrichment if you want:
  # up_entrez   <- bitr(genes_up,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  # down_entrez <- bitr(genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  # if (length(up_entrez) >= 10) {
  #   ek_up <- enrichKEGG(gene = up_entrez, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)
  #   write.csv(as.data.frame(ek_up), file.path(outdir, "enrichKEGG_Up_in_High.csv"), row.names = FALSE)
  # }
  # if (length(down_entrez) >= 10) {
  #   ek_down <- enrichKEGG(gene = down_entrez, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)
  #   write.csv(as.data.frame(ek_down), file.path(outdir, "enrichKEGG_Down_in_High.csv"), row.names = FALSE)
  # }
  
} else {
  message("clusterProfiler and/or org.Hs.eg.db not installed; skipping offline ORA.")
}

## ----------------------- DONE ---------------------------
cat("\nDone. Outputs written to:\n", normalizePath(outdir), "\n")
# Output files:
# - DGE_beta1_Untreated_MQI_High_vs_Low.csv
# - Volcano_beta1_Untreated_MQI_High_vs_Low.pdf
# - ORA_gprof_Up_in_High.csv (if any) / ORA_gprof_Down_in_High.csv (if any)
# - enrichGO_Up_in_High.csv / enrichGO_Down_in_High.csv (+ dotplots if enrichplot is available)


###############################################################
## Circos chord diagrams (custom selection) with GENE COLORS
## colored by log2 fold-change from your DGE results.
##
## What this script does:
## - Reads your existing g:Profiler CSVs:
##     ORA_gprof_Up_in_High.csv
##     ORA_gprof_Down_in_High.csv
## - Lets you PICK which pathways to show (by exact name)
## - Reads your DGE table (for avg_log2FC) and your gene lists:
##     genes_up_in_High.txt, genes_down_in_High.txt
## - Builds Pathway ↔ Gene edges from ORA intersections
## - Colors the **gene sectors** (and links) by avg_log2FC:
##       blue  = negative (down in High)
##       grey  = ~0
##       red   = positive (up in High)
## - Pathway sectors keep distinct colors.
## - Writes two PDFs: Up and Down selections (separately)
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(circlize)  # for chordDiagram() and colorRamp2()
})

## ---------------------- CONFIG ----------------------
base_dir <- "C:/Users/mqadir/Box/Lab 2301/1. R_Coding Scripts/Algo.Ai/DATA/DGE/wilcox/MFI"

# ORA CSVs (already flattened earlier; may contain 'intersection' or 'intersection_symbols')
file_up_ora   <- file.path(base_dir, "ORA_gprof_Up_in_High.csv")
file_down_ora <- file.path(base_dir, "ORA_gprof_Down_in_High.csv")

# DGE results CSV (from your DGE script; must contain columns: gene, avg_log2FC)
file_dge <- file.path(base_dir, "DGE_beta1_Untreated_MQI_High_vs_Low.csv")

# Optional gene lists (used for info/consistency checks; not strictly required)
file_genes_up   <- file.path(base_dir, "genes_up_in_High.txt")
file_genes_down <- file.path(base_dir, "genes_down_in_High.txt")

# >>> Provide the pathway names you want to plot (as they appear in your CSVs)
chosen_up_paths <- c(
  "oxidative phosphorylation",
  "aerobic electron transport chain",
  "ATP biosynthetic process",
  "mitochondrion organization",
  "de novo' protein folding",
  "amylin receptor signaling pathway",
  "mitochondrial unfolded protein response"
)

chosen_down_paths <- c(
  "hormone secretion",
  "response to calcium ion",
  "insulin-like growth factor binding protein complex",
  "regulation of programmed cell death",
  "lipid droplet"
)

# Matching options for pathway names
exact_names   <- TRUE    # match names exactly
ignore_case   <- TRUE    # case-insensitive comparison

# Output PDFs
circos_up_pdf   <- file.path(base_dir, "Circos_Custom_Up_in_High_FCcolored.pdf")
circos_down_pdf <- file.path(base_dir, "Circos_Custom_Down_in_High_FCcolored.pdf")

## ----------------- Helpers -----------------

# Safely return the first column that exists in df among 'candidates'
first_existing_col <- function(df, candidates, default = NA) {
  for (nm in candidates) {
    if (nm %in% names(df)) return(df[[nm]])
  }
  rep(default, nrow(df))
}

# Read g:Profiler CSV and standardize key columns (robust to varying headers)
read_gprof_csv <- function(path) {
  stopifnot(file.exists(path))
  df <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  names(df) <- tolower(names(df))  # normalize headers to lowercase
  
  # Prefer already-encoded symbols if present; otherwise intersections (IDs)
  term_vec  <- first_existing_col(df, c("term_name","term","name"))
  genes_vec <- first_existing_col(df, c("intersection_symbols","intersection"))
  padj_vec  <- first_existing_col(df, c("p_value_adjusted","padj","p_adjust","p_value"))
  
  out <- df %>%
    mutate(
      term  = as.character(term_vec),
      genes = as.character(genes_vec),
      padj  = suppressWarnings(as.numeric(padj_vec))
    ) %>%
    filter(!is.na(term), term != "")
  out
}

# Filter to chosen pathways (exact or exact+case-insensitive)
filter_terms <- function(df, wanted, exact = TRUE, ignore_case = TRUE) {
  if (length(wanted) == 0) stop("Please provide at least one pathway name in `chosen_*_paths`.")
  if (exact) {
    if (ignore_case) {
      df %>% filter(tolower(term) %in% tolower(wanted))
    } else {
      df %>% filter(term %in% wanted)
    }
  } else {
    # partial matching variant could be added if needed
    stop("Currently only exact matching is implemented. Set exact_names <- TRUE.")
  }
}

# Build Pathway↔Gene edge list from the 'genes' column (semicolon- or comma-separated)
edges_from_intersection <- function(df) {
  if (!nrow(df)) return(tibble(Pathway = character(), Gene = character()))
  df %>%
    transmute(term,
              genes = if_else(is.na(genes), "", genes)) %>%
    separate_rows(genes, sep = "\\s*;\\s*|\\s*,\\s*") %>%
    mutate(
      Gene = toupper(str_trim(genes)),
      Pathway = term
    ) %>%
    filter(Gene != "") %>%
    select(Pathway, Gene) %>%
    distinct()
}

# Build a color function for log2FC values
#    low (neg) -> blue, mid (0) -> light grey, high (pos) -> red
fc_color_fun <- function(fc_values) {
  rng <- range(fc_values, na.rm = TRUE)
  # Ensure symmetric-ish range if possible; avoids all red/blue when unbalanced
  lim <- max(abs(rng), na.rm = TRUE)
  # Fallback in case all NAs or all zeros
  if (!is.finite(lim) || lim == 0) lim <- 1
  circlize::colorRamp2(c(-lim, 0, lim), c("#1976D2", "#D9D9D9", "#C2185B"))
}

# Draw circos chord diagram with gene colors from log2FC
plot_circos_fc <- function(edges, fc_tbl, title = "Circos", pdf_file = NULL) {
  if (!nrow(edges)) {
    message("No edges to plot for: ", title)
    return(invisible(NULL))
  }
  
  # Join fold changes onto edges (by Gene symbol)
  fc_tbl <- fc_tbl %>%
    mutate(Gene = toupper(gene)) %>%
    distinct(Gene, avg_log2FC, .keep_all = TRUE)
  
  edges_fc <- edges %>%
    left_join(fc_tbl %>% select(Gene, avg_log2FC), by = "Gene")
  
  # Build color function from present FCs
  col_fun <- fc_color_fun(edges_fc$avg_log2FC)
  
  pathways <- unique(edges_fc$Pathway)
  genes    <- unique(edges_fc$Gene)
  
  # Colors: pathways get distinct hues; genes colored by FC (NA -> grey)
  set.seed(42)
  pal_path <- grDevices::rainbow(length(pathways), s = 0.6, v = 0.9)
  names(pal_path) <- pathways
  
  gene_cols <- vapply(edges_fc$Gene, function(g) {
    val <- edges_fc$avg_log2FC[match(g, edges_fc$Gene)]
    if (is.na(val)) "#BBBBBB" else col_fun(val)
  }, FUN.VALUE = character(1))
  
  # grid.col expects a named vector for ALL sectors (unique)
  grid_gene_cols <- gene_cols[match(genes, edges_fc$Gene)]
  names(grid_gene_cols) <- genes
  
  grid_cols <- c(pal_path, grid_gene_cols)
  
  # Edge (link) colors: use the gene color for each edge row (with some transparency)
  edge_cols <- vapply(edges_fc$Gene, function(g) grid_gene_cols[g], FUN.VALUE = character(1))
  edge_cols <- adjustcolor(edge_cols, alpha.f = 0.35)
  
  if (!is.null(pdf_file)) pdf(pdf_file, width = 10, height = 7)
  
  circlize::circos.clear()
  circlize::circos.par(gap.after = c(rep(6, length(pathways)-1), 12, rep(2, length(genes)-1), 12))
  
  circlize::chordDiagram(
    x = edges_fc,
    grid.col = grid_cols,
    col = edge_cols,
    transparency = 0,                 # we controlled alpha in edge_cols
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.08)
  )
  
  # Add sector labels
  circlize::circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector = circlize::get.cell.meta.data("sector.index")
      xlim   = circlize::get.cell.meta.data("xlim")
      ylim   = circlize::get.cell.meta.data("ylim")
      circlize::circos.text(
        x = mean(xlim), y = ylim[1] + 0.1, labels = sector,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6, col = "black"
      )
    },
    bg.border = NA
  )
  
  title(title, line = -1)
  
  # (Optional) add a small legend for FC colors
  # Draw in base plotting region:
  par(xpd = NA)
  lg_x <- grconvertX(0.9, from = "npc", to = "user")
  lg_y <- grconvertY(0.9, from = "npc", to = "user")
  # Color bar
  z <- seq(-1, 1, length.out = 100)
  cols <- col_fun(z)
  rasterImage(as.raster(matrix(cols, ncol = 1)), lg_x, lg_y - 0.3, lg_x + 0.02, lg_y)
  text(lg_x + 0.04, lg_y,        labels = "log2FC +", cex = 0.7, adj = 0)
  text(lg_x + 0.04, lg_y - 0.15, labels = "0",        cex = 0.7, adj = 0)
  text(lg_x + 0.04, lg_y - 0.3,  labels = "log2FC −", cex = 0.7, adj = 0)
  
  if (!is.null(pdf_file)) dev.off()
  invisible(NULL)
}

## -------------------- Load inputs --------------------
# ORA CSVs
stopifnot(file.exists(file_up_ora) || file.exists(file_down_ora))
up_df_ora   <- if (file.exists(file_up_ora))   read_gprof_csv(file_up_ora)   else tibble()
down_df_ora <- if (file.exists(file_down_ora)) read_gprof_csv(file_down_ora) else tibble()

message(sprintf("Loaded ORA: Up terms = %d, Down terms = %d", nrow(up_df_ora), nrow(down_df_ora)))

# DGE CSV (for log2FC)
stopifnot(file.exists(file_dge))
dge_tbl <- readr::read_csv(file_dge, show_col_types = FALSE)
stopifnot(all(c("gene","avg_log2FC") %in% names(dge_tbl)))

# Optional gene lists (informational)
genes_up   <- if (file.exists(file_genes_up))   readLines(file_genes_up)   else character(0)
genes_down <- if (file.exists(file_genes_down)) readLines(file_genes_down) else character(0)

## -------------------- Select pathways --------------------
if (nrow(up_df_ora) > 0) {
  up_sel <- filter_terms(up_df_ora, chosen_up_paths, exact = exact_names, ignore_case = ignore_case)
  if (nrow(up_sel) == 0) warning("No Up pathways matched your selections.")
  edges_up <- edges_from_intersection(up_sel)
} else {
  up_sel <- tibble(); edges_up <- tibble(Pathway = character(), Gene = character())
}

if (nrow(down_df_ora) > 0) {
  down_sel <- filter_terms(down_df_ora, chosen_down_paths, exact = exact_names, ignore_case = ignore_case)
  if (nrow(down_sel) == 0) warning("No Down pathways matched your selections.")
  edges_down <- edges_from_intersection(down_sel)
} else {
  down_sel <- tibble(); edges_down <- tibble(Pathway = character(), Gene = character())
}

## -------------------- Plot circos with FC colors --------------------
if (nrow(edges_up) > 0) {
  plot_circos_fc(
    edges = edges_up,
    fc_tbl = dge_tbl,  # uses avg_log2FC from this table
    title = "Custom (Up in High) — Pathways ↔ Genes (gene colors = log2FC)",
    pdf_file = circos_up_pdf
  )
  message("Wrote: ", circos_up_pdf)
}

if (nrow(edges_down) > 0) {
  plot_circos_fc(
    edges = edges_down,
    fc_tbl = dge_tbl,
    title = "Custom (Down in High) — Pathways ↔ Genes (gene colors = log2FC)",
    pdf_file = circos_down_pdf
  )
  message("Wrote: ", circos_down_pdf)
}

## -------------------- Console feedback --------------------
if (exists("up_sel") && nrow(up_sel) > 0) {
  cat("\nSelected Up pathways:\n"); print(up_sel$term)
}
if (exists("down_sel") && nrow(down_sel) > 0) {
  cat("\nSelected Down pathways:\n"); print(down_sel$term)
}

# Notes:
# - If your ORA CSVs contain only IDs (not symbols) in 'intersection', re-save them with
#   a symbol column (e.g., 'intersection_symbols') or replace 'genes' with symbols before running.
# - Any gene in the edges that lacks an avg_log2FC in the DGE CSV is colored grey.
# - You can adjust the blue/grey/red hex codes in fc_color_fun() to match your palette.

###############################################################
## β1 Untreated — Per-donor robustness (meta-analysis + forest)
## Builds df_cells safely, then runs metafor and plots a forest.
###############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(metafor)
})

## ---- CONFIG: edit if your column names differ ----
obj         <- combined
treat_col   <- "treat_simple"   # "Untreated"
beta_col    <- "beta_cluster"   # β1–β4
beta_target <- "β1"
feature_col <- "MQI_pred"
donor_col   <- "Library"
sex_col     <- "Sex"
threshold   <- 1

## ---- Safety: detach plyr to avoid count()/mutate masking ----
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

## ---- 1) Subset β1 Untreated and build per-cell table ----
stopifnot(all(c(treat_col, beta_col, feature_col, donor_col, sex_col) %in% colnames(obj@meta.data)))

sub_beta1 <- subset(
  obj,
  subset = !!rlang::sym(treat_col) == "Untreated" & !!rlang::sym(beta_col) == beta_target
)

df_cells <- sub_beta1@meta.data %>%
  as.data.frame() %>%
  transmute(
    Library = .data[[donor_col]],
    Sex     = .data[[sex_col]],
    val     = .data[[feature_col]]
  ) %>%
  # keep valid rows and exclude exactly threshold
  filter(!is.na(Library), !is.na(Sex), !is.na(val), val != threshold) %>%
  mutate(
    Library = as.character(Library),
    Sex     = factor(Sex, levels = c("M","F")),
    Group   = factor(ifelse(val > threshold, "High (>1)", "Low (<1)"),
                     levels = c("Low (<1)", "High (>1)"))
  )

if (nrow(df_cells) == 0) stop("No usable β1 Untreated cells after filtering.")

## ---- 2) Per-donor counts and proportions ----
tab_donor <- df_cells %>%
  group_by(Library, Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  mutate(
    High  = `High (>1)`,
    Low   = `Low (<1)`,
    total = High + Low,
    prop_high = High / total
  ) %>%
  arrange(desc(prop_high)) %>%
  as_tibble()

print(tab_donor)

## ---- 3) Logit-transformed proportions + variances ----
escalc_dat <- metafor::escalc(
  measure = "PLO",            # logit of proportion
  xi = tab_donor$High,
  ni = tab_donor$total,
  slab = tab_donor$Library
)

## ---- 4) Random-effects meta-analysis ----
res <- metafor::rma(yi, vi, data = escalc_dat, method = "REML")
print(res)

## ---- 5) Forest plot (back-transformed to proportion) ----
out_pdf <- file.path(getwd(), "Forest_Beta1_Untreated_MQI_High_vs_Low.pdf")
pdf(out_pdf, width = 7.2, height = 6.4)

metafor::forest(
  res,
  xlab = "Proportion of β1 cells with MQI_pred > 1",
  transf = metafor::transf.ilogit,                     # logit → proportion
  at = qlogis(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)),
  ilab = cbind(tab_donor$High, tab_donor$total),
  ilab.xpos = c(-2.7, -1.9),
  ilab.pos = 4,
  order = order(escalc_dat$yi),
  cex = 0.85
)
text(-2.7,  length(escalc_dat$yi) + 2, "High",  cex = 0.85)
text(-1.90, length(escalc_dat$yi) + 2, "Total", cex = 0.85)

dev.off()
message("Forest plot written to: ", out_pdf)

## ---- Subgroup meta-analysis by Sex (robust) ----
## ---- Subgroup meta-analysis by Sex (robust, fixed) ----
# Majority Sex per donor from df_cells
donor_sex <- df_cells %>%
  dplyr::group_by(Library, Sex) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(Library, Sex)

# escalc_dat rows are in the same order as the input vectors (tab_donor)
# so we can attach the Library IDs from tab_donor directly.
escalc_sex <- escalc_dat %>%
  as.data.frame() %>%
  dplyr::mutate(Library = tab_donor$Library) %>%
  dplyr::left_join(donor_sex, by = "Library")

# Helper to fit subgroup model if it has enough donors
fit_if <- function(d, label) {
  d <- subset(d, !is.na(Sex) & Sex == label)
  if (nrow(d) >= 2) metafor::rma(yi, vi, data = d, method = "REML") else NULL
}

res_M <- fit_if(escalc_sex, "M")
res_F <- fit_if(escalc_sex, "F")

if (!is.null(res_M)) { cat("\nMale subgroup:\n"); print(res_M) }
if (!is.null(res_F)) { cat("\nFemale subgroup:\n"); print(res_F) }

# (Optional) write subgroup forest plots
plot_subset <- function(fit, title, outfile) {
  if (is.null(fit)) return(invisible(NULL))
  pdf(outfile, width = 7.2, height = 6.0)
  metafor::forest(
    fit,
    xlab   = paste0(title, " — Proportion MQI_pred > 1"),
    transf = metafor::transf.ilogit,
    at     = qlogis(c(0.1, 0.25, 0.5, 0.75, 0.9)),
    cex    = 0.9
  )
  dev.off()
  message("Wrote: ", outfile)
}

plot_subset(res_M,
            title = "β1 Untreated (Male donors)",
            outfile = file.path(getwd(), "Forest_Beta1_Untreated_Male.pdf"))

plot_subset(res_F,
            title = "β1 Untreated (Female donors)",
            outfile = file.path(getwd(), "Forest_Beta1_Untreated_Female.pdf"))


###############################################################
## Panel K (fixed) — Up × Down Pathway Overlap Heatmap
## Uses FULL GO:BP gene sets from MSigDB via msigdbr (direction-agnostic),
## so overlaps between Up and Down are no longer trivially zero.
###############################################################
install.packages("msigdbr")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(readr); library(ggplot2); library(msigdbr)
})

k_out_pdf <- file.path(base_dir, "K_Pathway_Overlap_Heatmap_Up_vs_Down_FULL_GO.pdf")

# --- Helper: normalize term names for fuzzy matching ---
.norm <- function(x) {
  x %>% tolower() %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_trim() %>%
    str_squish()
}

# --- Load MSigDB GO:BP gene sets (species = human) ---
m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  transmute(term = gs_name, name = .norm(gs_name), gene = toupper(gene_symbol))

# Also keep a lookup by plain English term if available
# (MSigDB BP names look like "GOBP_OXIDATIVE_PHOSPHORYLATION")
m_alias <- m_df %>%
  mutate(alias = .norm(str_remove(gs_name, "^GOBP_") %>% str_replace_all("_", " "))) %>%
  select(term, alias) %>% distinct()

m_terms <- m_df %>% select(term) %>% distinct()

# Build full gene sets as list: TERM -> unique genes
full_sets <- m_df %>%
  group_by(term) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  { setNames(.$genes, .$term) }

# --- Map your selected terms to MSigDB terms ---
# We’ll try exact alias, then contains match as fallback.
map_to_msig <- function(user_terms) {
  out <- vector("list", length(user_terms))
  names(out) <- user_terms
  for (t in user_terms) {
    q <- .norm(t)
    # exact alias match first
    hits <- m_alias$term[m_alias$alias == q]
    if (length(hits) == 0) {
      # fallback: partial alias contains
      hits <- m_alias$term[str_detect(m_alias$alias, fixed(q))]
    }
    if (length(hits) == 0) {
      out[[t]] <- character(0)  # not found
    } else {
      out[[t]] <- unique(hits)
    }
  }
  out
}

# Use the SAME selections you already made
#   up_sel / down_sel computed earlier; they hold term names from g:Profiler
up_names   <- if (exists("up_sel") && nrow(up_sel))   unique(up_sel$term)   else character(0)
down_names <- if (exists("down_sel") && nrow(down_sel)) unique(down_sel$term) else character(0)

stopifnot(length(up_names) > 0, length(down_names) > 0)

map_up   <- map_to_msig(up_names)
map_down <- map_to_msig(down_names)

# Convert to TERM->genes lists using MSigDB (fallback to ORA intersections if unfound)
term_gene_list_from_full <- function(names_vec, map_list, df_sel) {
  # fallback from your ORA intersections (direction-specific)
  fallback <- term_gene_list(df_sel)  # already defined earlier in your script
  out <- vector("list", length(names_vec)); names(out) <- names_vec
  for (nm in names_vec) {
    msig_terms <- map_list[[nm]]
    if (length(msig_terms) > 0) {
      genes <- unique(unlist(full_sets[msig_terms], use.names = FALSE))
      out[[nm]] <- genes
    } else {
      # fallback: still show something rather than drop the term
      out[[nm]] <- fallback[[nm]] %||% character(0)
    }
  }
  out
}

up_full_sets   <- term_gene_list_from_full(up_names,   map_up,   up_sel)
down_full_sets <- term_gene_list_from_full(down_names, map_down, down_sel)

# --- Build all Up × Down pairs with FULL gene sets ---
pairs_df <- tidyr::expand_grid(
  up   = names(up_full_sets),
  down = names(down_full_sets)
) %>%
  rowwise() %>%
  mutate(
    shared  = length(intersect(up_full_sets[[up]], down_full_sets[[down]])),
    unionn  = length(union(up_full_sets[[up]],   down_full_sets[[down]])),
    jaccard = ifelse(unionn > 0, shared / unionn, 0)
  ) %>%
  ungroup()

# Order rows/cols by signal
row_order <- pairs_df %>%
  group_by(up) %>%
  summarise(sum_j = sum(jaccard, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(sum_j)) %>%
  pull(up)

col_order <- pairs_df %>%
  group_by(down) %>%
  summarise(sum_j = sum(jaccard, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(sum_j)) %>%
  pull(down)

pairs_df <- pairs_df %>%
  mutate(up = factor(up, levels = row_order),
         down = factor(down, levels = col_order))

# Pretty labels (original user names) on axes
lab_up   <- setNames(names(up_full_sets),   names(up_full_sets))
lab_down <- setNames(names(down_full_sets), names(down_full_sets))

# Plot
p_heat <- ggplot(pairs_df, aes(x = down, y = up, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(shared > 0, shared, "")), size = 3) +
  scale_fill_gradient(
    name = "Jaccard",
    low = "#F7FBFF", high = "#08306B",
    limits = c(0, max(1e-8, max(pairs_df$jaccard, na.rm = TRUE)))
  ) +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5)) +
  labs(
    x = "Downregulated pathways (full GO sets)",
    y = "Upregulated pathways (full GO sets)",
    title = "Panel K — Up × Down Pathway Overlap (direction-agnostic GO gene sets)",
    subtitle = "Fill = Jaccard index (shared/union); Numbers = shared genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(k_out_pdf, p_heat, width = 9.5, height = 6.8, units = "in")
message("Panel K heatmap written: ", k_out_pdf)

# Inspect strongest crosstalk
pairs_df %>%
  arrange(desc(shared), desc(jaccard)) %>%
  head(15) %>%
  { cat("\nTop overlaps (Up ↔ Down, FULL GO sets): shared genes, Jaccard\n"); print(.) }


## ============================================================
## Panel K — β1 (Untreated): HbA1c vs % SQSTM1+ cells per donor
## 4 lines = (M/F) × (SQSTM1+ / SQSTM1−) with R and p shown
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(rlang)
  library(Seurat)
})

## --- Avoid plyr/dplyr conflicts (if plyr is attached accidentally) ---
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

## -------- CONFIG (edit paths/cols if needed) --------
csv_path        <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\Donor_Summary_186.csv)"
hba1c_col_in    <- "hba1c"        # HbA1c column name in the CSV
obj             <- combined       # your Seurat object already in memory
assay_name      <- "RNA"          # assay holding SQSTM1 expression (log1p-normalized typical)
gene_symbol     <- "SQSTM1"       # target gene symbol
beta_col        <- "beta_cluster" # meta.data col with β1–β4
beta_target     <- "β1"           # value of β1 in beta_col
sex_col         <- "Sex"          # "M"/"F"
treat_col       <- "treat_simple" # "Untreated"
donor_col       <- "Library"      # donor/sample ID in meta.data

## Threshold for SQSTM1 positivity (log1p scale). 0 = any detected counts.
## For a stricter cut (≈ raw count ≥1 before log1p), use 0.25.
sq_pos_threshold <- 0

## Colors: blue/red for positive; lighter tints for negative
col_M_pos <- "#1976D2"; col_M_neg <- "#90CAF9"
col_F_pos <- "#C2185B"; col_F_neg <- "#F48FB1"

label_pos <- "SQSTM1+"
label_neg <- "SQSTM1−"

## -------- 1) HbA1c lookup (CSV + manual overrides) --------
clinical_df <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(
    Library = toupper(str_trim(donor_ID)),
    hba1c   = as.numeric(.data[[hba1c_col_in]])
  ) %>%
  select(Library, hba1c) %>%
  filter(!is.na(Library), !is.na(hba1c)) %>%
  mutate(.src = "csv")

## Manual fill-ins for donors missing HbA1c in the CSV
hba1c_missing <- tibble::tribble(
  ~Library,        ~hba1c,
  "HP2022801",     5.5,
  "HP2024001",     5.4,
  "SAMN15877725",  5.6,
  "HP2031401",     5.4,
  "HP2105501",     5.6,
  "HP2106201",     5.3,
  "HP2107001",     5.1,
  "HP2107901",     5.2,
  "HP2108601",     5.1,
  "HP2108901",     5.9,
  "HP2110001",     5.5,
  "HP2123201",     5.3,
  "HP2132801",     5.5,
  "HP2202101",     5.5,
  "HP2121601",     5.8
) %>%
  mutate(
    Library = toupper(str_trim(Library)),
    hba1c   = as.numeric(hba1c),
    .src    = "manual"
  )

## Prefer manual values when duplicates exist
hba1c_full <- bind_rows(clinical_df, hba1c_missing) %>%
  arrange(Library, .src) %>%
  group_by(Library) %>%
  slice_tail(n = 1) %>%  # keep last (manual overrides csv)
  ungroup() %>%
  select(Library, hba1c)

## -------- 2) Subset β1 Untreated & extract SQSTM1 expression --------

## Safety checks on Seurat object and assay
if (!inherits(obj, "Seurat")) stop("`obj` is not a Seurat object.")
assays_available <- names(obj@assays)
if (!(assay_name %in% assays_available)) {
  message("Assay '", assay_name, "' not found; using DefaultAssay(obj) = ", DefaultAssay(obj))
  assay_name <- DefaultAssay(obj)
}

## Ensure needed meta columns exist
stopifnot(all(c(beta_col, sex_col, treat_col, donor_col) %in% colnames(obj@meta.data)))

## Subset β1 Untreated cells
sub_beta1 <- subset(
  obj,
  subset = !!rlang::sym(treat_col) == "Untreated" & !!rlang::sym(beta_col) == beta_target
)
if (ncol(sub_beta1) == 0) stop("Subset has 0 cells: check treat_col/beta_col values.")

## Extract SQSTM1 expression from assay
expr_vec <- tryCatch({
  mat <- GetAssayData(sub_beta1, assay = assay_name, slot = "data") # log1p-normalized
  if (!(gene_symbol %in% rownames(mat))) stop(sprintf("Gene %s not found in assay %s", gene_symbol, assay_name))
  as.numeric(mat[gene_symbol, ])
}, error = function(e) {
  ## Fallback to counts if 'data' missing (less common)
  mat <- GetAssayData(sub_beta1, assay = assay_name, slot = "counts")
  if (!(gene_symbol %in% rownames(mat))) stop(sprintf("Gene %s not found in assay %s (counts slot fallback)", gene_symbol, assay_name))
  as.numeric(mat[gene_symbol, ])
})

## Build per-cell frame with SQSTM1 positivity
df_cells <- sub_beta1@meta.data %>%
  as.data.frame() %>%
  mutate(
    Library = toupper(str_trim(.data[[donor_col]])),
    Sex     = .data[[sex_col]],
    SQSTM1  = expr_vec
  ) %>%
  filter(!is.na(Library), !is.na(Sex), !is.na(SQSTM1)) %>%
  mutate(
    Sex   = factor(Sex, levels = c("M","F")),
    Group = ifelse(SQSTM1 > sq_pos_threshold, label_pos, label_neg)
  )

if (nrow(df_cells) == 0L) stop("No usable β1 Untreated cells after filtering.")

## (Optional) Adaptive per-donor thresholding (median split) instead of fixed cut:
# df_cells <- df_cells %>%
#   group_by(Library) %>%
#   mutate(Group = ifelse(SQSTM1 > median(SQSTM1, na.rm = TRUE), label_pos, label_neg)) %>%
#   ungroup()

## -------- 3) Per-donor × Sex × Group: % SQSTM1+/− β1 cells; join HbA1c --------
by_donor <- df_cells %>%
  dplyr::count(Library, Sex, Group, name = "n") %>%
  dplyr::group_by(Library, Sex) %>%
  dplyr::mutate(
    n_total   = sum(n),
    pct_group = 100 * n / n_total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(hba1c_full, by = "Library") %>%
  dplyr::filter(!is.na(hba1c))

if (nrow(by_donor) == 0L) stop("No donors with HbA1c after join.")

## -------- 4) Correlation per line = (Sex × SQSTM1 group) --------
## Build legend key that matches plotted groups
by_donor$SexGroup <- interaction(by_donor$Sex, by_donor$Group, sep = " · ")

## Compute Pearson r and p for each line (needs n>=3 points)
corr_tbl <- by_donor %>%
  dplyr::group_by(Sex, Group, SexGroup) %>%
  dplyr::summarise(
    n = dplyr::n(),
    r = suppressWarnings(if (n >= 3) cor(hba1c, pct_group, method = "pearson") else NA_real_),
    p = tryCatch(if (n >= 3) cor.test(hba1c, pct_group, method = "pearson")$p.value else NA_real_,
                 error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p.label   = dplyr::case_when(
      is.na(p)    ~ "NA",
      p < 2.2e-16 ~ "<2.2e-16",
      TRUE        ~ as.character(signif(p, 2))
    ),
    r.label   = ifelse(is.na(r), "NA", sprintf("%.2f", r)),
    cor.label = paste0("R=", r.label, "\np=", p.label)
  )

## Place each label slightly above that line’s max Y, left side of the x-axis
pos_tbl <- by_donor %>%
  dplyr::group_by(SexGroup) %>%
  dplyr::summarise(y_max = max(pct_group, na.rm = TRUE), .groups = "drop")
x_min <- min(by_donor$hba1c, na.rm = TRUE)

corr_tbl <- corr_tbl %>%
  dplyr::left_join(pos_tbl, by = "SexGroup") %>%
  dplyr::mutate(
    y.label = y_max + 0.05 * abs(y_max),
    x.label = x_min
  )

## -------- 5) Plot: 4 lines (M/F × SQSTM1+ / SQSTM1−) + R/p labels --------
pal <- setNames(
  c(col_M_pos,               col_M_neg,               col_F_pos,               col_F_neg),
  c(paste0("M · ", label_pos), paste0("M · ", label_neg), paste0("F · ", label_pos), paste0("F · ", label_neg))
)

p <- ggplot(by_donor, aes(x = hba1c, y = pct_group, color = SexGroup)) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, size = 1.05) +
  ## correlation labels (colored to match their line)
  geom_text(
    data = corr_tbl,
    aes(x = x.label, y = y.label, label = cor.label, color = SexGroup),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1, size = 4, fontface = "bold"
  ) +
  scale_color_manual(values = pal, name = "Sex × SQSTM1 group") +
  labs(
    title = "β1 (Untreated): HbA1c vs % of β1 cells that are SQSTM1-positive/negative",
    x = "HbA1c (%)",
    y = "% of β1 cells per donor"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
p
## Save
dir.create("figures_regression", showWarnings = FALSE, recursive = TRUE)
ggsave("figures_regression/HbA1c_vs_pct_SQSTM1pos_beta1_Untreated_4lines_with_stats.pdf",
       p, width = 8.5, height = 6.2)

p

# ---------------------------
# Patch-seq Analysis - β-cell mapping
# ---------------------------
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)

# --- 1) Load datasets ---
# Reference (β1–β4 clusters, already SCT-normalized)
beta_cells <- qread(
  r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\beta_cells.qs)"
)

# Query (Patrick Macdonald’s Patch-seq wide-format table)
patch_file <- r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\macdonald\18082025 For Dr Qadir.csv)"
patch_raw  <- read.csv(patch_file, row.names = 1, check.names = FALSE)

# Split Patch-seq into genes vs metadata
ref_genes    <- rownames(beta_cells)
common_genes <- intersect(rownames(patch_raw), ref_genes)
patch_expr   <- as.matrix(patch_raw[common_genes, , drop = FALSE]); mode(patch_expr) <- "numeric"

meta_rows    <- setdiff(rownames(patch_raw), common_genes)
patch_meta   <- as.data.frame(t(patch_raw[meta_rows, , drop = FALSE]), stringsAsFactors = FALSE)
patch_meta   <- tibble::rownames_to_column(patch_meta, "cell_id") |>
  mutate(across(-cell_id, ~ type.convert(as.character(.), as.is = TRUE)))
rownames(patch_meta) <- patch_meta$cell_id; patch_meta$cell_id <- NULL

# Make Seurat object
patch_obj <- CreateSeuratObject(counts = patch_expr, meta.data = patch_meta)

# Keep only β cells (adjust col name if needed)
ctype_col <- c("cell_type","celltype","ctype","Cell.Type","Cell_type")
ctype_col <- ctype_col[ctype_col %in% colnames(patch_obj@meta.data)][1]
stopifnot(!is.na(ctype_col))
patch_obj$cell_type_tmp <- tolower(as.character(patch_obj@meta.data[[ctype_col]]))
patch_beta <- subset(patch_obj, subset = cell_type_tmp %in% c("beta","β","beta cell"))

# --- 2) SCT transform for mapping ---
DefaultAssay(patch_beta) <- "RNA"
patch_beta <- SCTransform(patch_beta, verbose = TRUE)

DefaultAssay(beta_cells) <- "RNA"
beta_cells <- SCTransform(beta_cells, verbose = TRUE)

# --- 3) Subset male β INS-hi/low ---
# (skip for Patch-seq: you already subset to β cells)

# --- 4) Run UMAP on reference with model ---
beta_cells <- RunUMAP(
  beta_cells, reduction = "pca", dims = 1:8,
  return.model = TRUE
)

DimPlot(beta_cells, reduction = "umap", label = TRUE, repel = TRUE) + theme_minimal()

# --- 5) Find anchors & MapQuery ---
anchors <- FindTransferAnchors(
  reference = beta_cells,
  query = patch_beta,
  normalization.method = "SCT",
  reference.reduction = "pca"
)

patch_beta <- MapQuery(
  anchorset = anchors,
  query = patch_beta,
  reference = beta_cells,
  refdata = list(
    collapsed_cluster = "collapsed_cluster",
    beta_cluster = "beta_cluster"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# --- 6) Copy ref.umap to a standard 'umap' slot in query ---
patch_beta[["umap"]] <- CreateDimReducObject(
  embeddings = Embeddings(patch_beta[["ref.umap"]]),
  key = "UMAP_",
  assay = DefaultAssay(patch_beta)
)

# --- 7) Merge reference + query ---
combined <- merge(
  x = beta_cells, y = patch_beta,
  add.cell.ids = c("REF", "PATCH"),
  merge.data = TRUE
)

# --- 8) Build a unified UMAP that survives merge ---
ref_umap <- Embeddings(beta_cells[["umap"]]); rownames(ref_umap) <- paste0("REF_", rownames(ref_umap))
qry_umap <- Embeddings(patch_beta[["umap"]]); rownames(qry_umap) <- paste0("PATCH_", rownames(qry_umap))
umap_mat <- rbind(ref_umap, qry_umap); umap_mat <- umap_mat[Cells(combined), , drop = FALSE]

combined[["umap"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key = "UMAP_",
  assay = DefaultAssay(combined)
)

# --- 9) Fill NA collapsed_cluster with predicted values ---
combined$collapsed_cluster <- ifelse(
  is.na(combined$collapsed_cluster),
  combined$predicted.collapsed_cluster,
  combined$collapsed_cluster
)
combined$beta_cluster <- ifelse(
  is.na(combined$beta_cluster),
  combined$predicted.beta_cluster,
  combined$beta_cluster
)

# --- 10) Plot unified UMAP ---
DimPlot(combined, reduction = "umap", group.by = "collapsed_cluster", label = FALSE) +
  scale_color_manual(values = c("A" = "firebrick", "B" = "darkorange",
                                "C" = "dodgerblue", "D" = "orchid"))

DimPlot(combined, reduction = "umap", group.by = "beta_cluster", label = FALSE, repel = TRUE)

# Complete dataset
#qsave(combined, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined_mcdonald.qs)")

list_meta_cols <- function(seu, n = 100) {
  stopifnot(inherits(seu, "Seurat"))
  cols <- colnames(seu@meta.data)
  head(cols, n)
}

## usage
list_meta_cols(combined, 100)

cut <- "NormalizedLateCaCurrentAmplitude_pA.pF"
cols <- colnames(combined@meta.data)
i <- match(cut, cols)
if (is.na(i)) stop("Column not found: ", cut)
combined@meta.data <- combined@meta.data[, 1:i, drop = FALSE]

# quick check
tail(colnames(combined@meta.data), 10)
#qsave(combined, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined_mcdonald_correctmetadat.qs)")


combined <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined.qs)")
table(combined@meta.data[["treatment"]])

library(dplyr)
library(ggplot2)

## 1) metadata for EtOH + DHT only (coerce to plain df)
md <- as.data.frame(combined@meta.data, stringsAsFactors = FALSE) |>
  subset(treatment %in% c("EtOH","DHT[10nM]")) |>
  transform(
    treatment_plot = factor(ifelse(treatment == "DHT[10nM]", "DHT", "EtOH"),
                            levels = c("EtOH","DHT"))
  )

## 2) counts -> percentages (force dplyr::count explicitly)
tab <- dplyr::count(md, treatment_plot, beta_cluster, name = "n") |>
  dplyr::group_by(treatment_plot) |>
  dplyr::mutate(pct = 100 * n / sum(n)) |>
  dplyr::ungroup()

## quick sanity check
head(tab)

## 3) 100% stacked bar with % labels
ggplot(tab, aes(x = treatment_plot, y = pct, fill = beta_cluster)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = ifelse(pct >= 4, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  labs(x = NULL, y = "Cell fraction (%)", fill = "β-cluster",
       title = "β-cluster composition by treatment (EtOH vs DHT)") +
  theme_classic(base_size = 12)

# Contine McDonald Analysis
combined_md_pc <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined_mcdonald_correctmetadat.qs)")
head(combined_md_pc@meta.data)

table(combined_md_pc@meta.data$'Cell Type')
table(combined_md_pc@meta.data$cell_type)
table(combined_md_pc@meta.data$celltype_qadir)

md <- combined_md_pc@meta.data

# 1) Create unified column
#    - For McDonald cells, keep their "cell_type"
#    - For all other cells, keep "celltype_qadir"
md$celltype_unified <- md$celltype_qadir
md$celltype_unified[!is.na(md$cell_type)] <- md$cell_type[!is.na(md$cell_type)]

# 2) Drop the old `cell_type` column
md$cell_type <- NULL

# 3) Write back into Seurat object
combined_md_pc@meta.data <- md

# 4) Quick sanity check
table(combined_md_pc$celltype_unified, useNA = "ifany")

md <- combined_md_pc@meta.data

# Check which sex columns exist
sex_cols <- grep("^sex$|^Sex$", colnames(md), value = TRUE)
print(sex_cols)
# e.g. might return c("Sex", "sex")

# 1) Create unified Sex column
# Prefer the one that is more complete, but overwrite with values if available
md$Sex_unified <- md[[sex_cols[1]]]   # start from first column
if (length(sex_cols) > 1) {
  # Overwrite with second column values when present
  overwrite_vals <- !is.na(md[[sex_cols[2]]]) & md[[sex_cols[2]]] != ""
  md$Sex_unified[overwrite_vals] <- md[[sex_cols[2]]][overwrite_vals]
}

# 2) Standardize to "M" and "F"
md$Sex_unified <- toupper(substr(md$Sex_unified, 1, 1))  # first letter, uppercase
md$Sex_unified[md$Sex_unified %in% c("MALE","M")] <- "M"
md$Sex_unified[md$Sex_unified %in% c("FEMALE","F")] <- "F"

# 3) Drop old duplicates, keep only one column
md <- md %>% select(-any_of(sex_cols))
colnames(md)[colnames(md) == "Sex_unified"] <- "Sex"

# 4) Write back
combined_md_pc@meta.data <- md

# 5) Quick check
table(combined_md_pc$Sex, useNA = "ifany")

table(combined_md_pc$Diabetes.Status, useNA = "ifany")
table(combined_md_pc$diabetes_status, useNA = "ifany")

md <- combined_md_pc@meta.data

# Source columns
src_cols <- c("Diabetes.Status", "diabetes_status")
src_cols <- intersect(src_cols, colnames(md))
if (length(src_cols) == 0) stop("No diabetes status columns found.")

# Start from scRNA-seq column, then overwrite with patch-seq where available
ds <- md[["diabetes_status"]]
if ("Diabetes.Status" %in% src_cols) {
  overwrite <- !is.na(md[["Diabetes.Status"]]) & md[["Diabetes.Status"]] != ""
  ds[overwrite] <- md[["Diabetes.Status"]][overwrite]
}

# Normalize to "ND" / "T2D"
norm <- toupper(trimws(as.character(ds)))
norm <- gsub("[[:punct:]]", "", norm)          # remove punctuation like "." or "-"
norm[norm %in% c("ND","NONDIABETIC","NONDIABETES","CTRL","CONTROL")] <- "ND"
norm[grepl("^T2", norm) | norm %in% c("T2D","T2DM","TYPE 2 DIABETES")] <- "T2D"

# Assign unified column with requested name
md[["Diabetes.status"]] <- ifelse(norm %in% c("ND","T2D"), norm, NA)

# Drop old columns
md$Diabetes.Status <- NULL
md$diabetes_status <- NULL

# Write back
combined_md_pc@meta.data <- md

# Sanity check
table(combined_md_pc$Diabetes.status, useNA = "ifany")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# 0) Sanity checks
stopifnot("beta_cluster" %in% colnames(combined_md_pc@meta.data))
stopifnot("umap" %in% names(combined_md_pc@reductions))

# 1) Keep only cells with a beta_cluster label
obj_b <- subset(combined_md_pc, subset = !is.na(beta_cluster) & beta_cluster != "")

# Optional: order clusters nicely if they’re labeled β1–β4
lvl <- c("β1","β2","β3","β4")
if (all(lvl %in% unique(obj_b$beta_cluster))) {
  obj_b$beta_cluster <- factor(obj_b$beta_cluster, levels = lvl)
}

# 2) Quick counts
cat("Cells per β-cluster:\n")
print(sort(table(obj_b$beta_cluster), decreasing = TRUE))

# 3) UMAP colored by beta_cluster
p_umap <- DimPlot(
  obj_b,
  reduction = "umap",
  group.by  = "beta_cluster",
  label     = TRUE,
  repel     = TRUE,
  raster    = TRUE,
  pt.size   = 0.5
) + ggtitle("UMAP — β-clusters") +
  theme(plot.title = element_text(hjust = 0.5))

print(p_umap)

md <- combined_md_pc@meta.data

# Add a column marking whether the cell is from McDonald patch-seq or Reference
md$dataset_origin <- ifelse(!is.na(md$Donor.ID) & md$Donor.ID != "",
                            "McDonald_patchseq", "Reference_scRNA")

# Write back into object
combined_md_pc@meta.data <- md

# Quick counts
table(combined_md_pc$dataset_origin)

# 3) UMAP colored by dataset_origin
DimPlot(
  combined_md_pc,
  reduction = "umap",
  group.by  = "dataset_origin",
  label     = FALSE,
  repel     = TRUE,
  #raster    = TRUE,
  pt.size   = 0.2
) + ggtitle("UMAP — dataset_origin") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(
    values = c(
      "Reference_scRNA"     = "dodgerblue3",
      "McDonald_patchseq"   = "firebrick3"
    )
  )

library(dplyr)
library(ggplot2)

md_tbl <- tibble::as_tibble(combined_md_pc@meta.data)

# 1) Keep only McDonald patch-seq cells
df_md <- md_tbl %>%
  filter(dataset_origin == "McDonald_patchseq") %>%
  select(beta_cluster, Diabetes.status)

# 2) Count by beta_cluster and Diabetes.status
tab_md <- df_md %>%
  dplyr::count(beta_cluster, Diabetes.status) %>%
  dplyr::group_by(beta_cluster) %>%
  dplyr::mutate(percent = 100 * n / sum(n))

print(tab_md)

# 3) Plot stacked percentage bars
p <- ggplot(tab_md, aes(x = beta_cluster, y = percent, fill = Diabetes.status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percent, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3.5) +
  labs(x = "β-cluster", y = "Percent distribution",
       title = "McDonald patch-seq cells by Diabetes.status") +
  scale_fill_manual(values = c("ND" = "dodgerblue3",
                               "T2D" = "firebrick3")) +
  theme_classic(base_size = 13)

print(p)


#PLOT data McDonald
# ============================================================
#   Heatmap of electrophysiology features in McDonald β cells
#   Panels: All β combined vs β1 only
#   Splits: Diabetes status (ND vs T2D) × Sex (F/M)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

md <- combined_md_pc@meta.data

# --- 1) All β combined ---
df_md_all <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster %in% c("β1","β3","β4")) %>%
  mutate(beta_cluster = "All β") %>%
  select(beta_cluster, Sex, Diabetes.status,
         CellSize_pF,
         NormalizedTotalCapacitance_fF.pF,
         NormalizedLateDepolarizationCapacitance,
         CalciumIntegralNormalizedtoCellSize_pC.pF,
         NormalizedPeakSodiumCurrentAmplitude_pA.pF,
         HalfInactivationSodiumCurrent_mV,
         NormalizedEarlyPeakCaCurrentAmplitude_pA.pF,
         NormalizedLateCaCurrentAmplitude_pA.pF)

# --- 2) β1 only ---
df_md_beta1 <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster == "β1") %>%
  select(beta_cluster, Sex, Diabetes.status,
         CellSize_pF,
         NormalizedTotalCapacitance_fF.pF,
         NormalizedLateDepolarizationCapacitance,
         CalciumIntegralNormalizedtoCellSize_pC.pF,
         NormalizedPeakSodiumCurrentAmplitude_pA.pF,
         HalfInactivationSodiumCurrent_mV,
         NormalizedEarlyPeakCaCurrentAmplitude_pA.pF,
         NormalizedLateCaCurrentAmplitude_pA.pF)

# --- 3) Combine both datasets ---
df_md <- bind_rows(df_md_all, df_md_beta1)

# --- 4) Reshape to long format ---
df_long <- df_md %>%
  pivot_longer(-c(beta_cluster, Sex, Diabetes.status),
               names_to = "variable", values_to = "value")

# --- 5) Compute means ---
# --- Run Wilcoxon tests (ND vs T2D) ---
stats_table <- df_long %>%
  group_by(beta_cluster, variable, Sex) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(value ~ Diabetes.status)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

# Inspect results
stats_table

# Scale within each feature
heatmap_data <- summary_stats %>%
  group_by(variable) %>%
  mutate(scaled = scale(mean)) %>%
  ungroup() %>%
  mutate(variable_sex = paste0(variable, "_", Sex))

# --- 6) Plot heatmap (side by side: All β vs β1) ---
ggplot(heatmap_data,
       aes(x = Diabetes.status, y = variable_sex, fill = scaled)) +
  geom_tile(color = "white") +
  facet_wrap(~beta_cluster, nrow = 1) +  # Panels: All β | β1
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Electrophysiology features — McDonald patch-seq",
       x = "Diabetes status", y = "Feature × Sex",
       fill = "Scaled mean") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ============================================================
#   Heatmap of electrophysiology features in McDonald β cells
#   Panels: All β combined vs β1 only
#   Splits: PINK1+ vs PINK1- × Sex × Diabetes status (ND, T2D)
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)

md <- combined_md_pc@meta.data

# --- 1) Ensure PINK1 status column ---
if (!"PINK1_status" %in% colnames(md)) {
  md$PINK1_status <- ifelse(combined_md_pc@assays$RNA@data["PINK1", ] > 0, "PINK1+", "PINK1-")
}

# --- 2) Subset McDonald β1 + All β (collapse β1+β3+β4) ---
df_md_all <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster %in% c("β1","β3","β4")) %>%
  mutate(beta_cluster = "All β")

df_md_beta1 <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster == "β1")

# Merge them
df_md <- bind_rows(df_md_all, df_md_beta1) %>%
  select(beta_cluster, Sex, Diabetes.status, PINK1_status,
         CellSize_pF,
         NormalizedTotalCapacitance_fF.pF,
         NormalizedLateDepolarizationCapacitance,
         CalciumIntegralNormalizedtoCellSize_pC.pF,
         NormalizedPeakSodiumCurrentAmplitude_pA.pF,
         HalfInactivationSodiumCurrent_mV,
         NormalizedEarlyPeakCaCurrentAmplitude_pA.pF,
         NormalizedLateCaCurrentAmplitude_pA.pF)

# --- 3) Reshape long ---
df_long <- df_md %>%
  pivot_longer(-c(beta_cluster, Sex, Diabetes.status, PINK1_status),
               names_to = "variable", values_to = "value")

# --- 4) Compute group means ---
summary_stats <- df_long %>%
  group_by(beta_cluster, variable, Sex, Diabetes.status, PINK1_status) %>%
  summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")

# Scale by feature
heatmap_data <- summary_stats %>%
  group_by(variable) %>%
  mutate(scaled = scale(mean)) %>%
  ungroup() %>%
  mutate(variable_sex = paste0(variable, "_", Sex))

# --- 5) Run Wilcoxon tests ---
stats_table <- df_long %>%
  group_by(beta_cluster, variable, Sex, Diabetes.status) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(value ~ PINK1_status)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

print(stats_table)

# --- 6) Plot heatmap ---
ggplot(heatmap_data,
       aes(x = PINK1_status, y = variable_sex, fill = scaled)) +
  geom_tile(color = "white") +
  facet_grid(~beta_cluster + Diabetes.status, scales = "free_x") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = -0.5) +
  labs(title = "Electrophysiology features — PINK1+ vs PINK1-",
       x = "PINK1 status", y = "Feature × Sex",
       fill = "Scaled mean") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ============================================================
# Donor-level correlations: MQI_v1 (mitochondrial index) vs
# McDonald patch-seq electrophysiology features
# Outputs:
#   1) 'corr_results' table (Spearman rho, p, FDR)
#   2) faceted scatterplot (one panel per feature)
# ============================================================

# ============================================================
# MQI_v1 (no percent.mt) + donor×sex correlations (McDonald β)
# - Axes: Mitophagy, OXPHOS, UPRmt, TCA, Glycolysis
# - MQI_v1 excludes percent.mt completely
# - Correlate donor-mean MQI_v1 vs electrophys per Sex (F/M)
# - Plot with shaded CI and triangle shape for T2D
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# ---------- 0) Pick Seurat object ----------
obj <- if (exists("combined")) combined else if (exists("combined_md_pc")) combined_md_pc else stop("No Seurat object named 'combined' or 'combined_md_pc' found.")
stopifnot(inherits(obj, "Seurat"))
DefaultAssay(obj) <- "RNA"

# ---------- 1) dataset_origin flag (McDonald vs reference) ----------
md <- obj@meta.data
if (!"dataset_origin" %in% colnames(md)) {
  stopifnot("Donor.ID" %in% colnames(md))
  obj$dataset_origin <- ifelse(!is.na(md$Donor.ID) & md$Donor.ID != "", "McDonald_patchseq", "Reference_scRNA")
  md <- obj@meta.data
}

# ---------- 2) Helpers ----------
set.seed(42)
`%||%` <- function(a,b) if (!is.null(a)) a else b
z <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
# NA-robust z: replaces NA with 0 after z-scoring so missing components don't nuke the sum
z0 <- function(x) { out <- z(x); out[!is.finite(out)] <- 0; out }

# Normalize Sex to {F, M}
if (!"Sex" %in% colnames(obj@meta.data)) stop("No 'Sex' column in meta.data.")
obj$Sex <- as.character(obj$Sex)
obj$Sex <- dplyr::case_when(
  tolower(obj$Sex) %in% c("f","female","woman","women") ~ "F",
  tolower(obj$Sex) %in% c("m","male","man","men")       ~ "M",
  TRUE                                                 ~ NA_character_
)

AddModuleScoreSafe <- function(obj_, genes, name) {
  present <- genes[genes %in% rownames(obj_)]
  if (length(present) == 0) {
    obj_[[name]] <- NA_real_
    message(sprintf("'%s': 0/%d genes present — set to NA.", name, length(genes)))
    return(obj_)
  }
  obj_ <- AddModuleScore(obj_, features = list(present), name = name,
                         nbin = 24, ctrl = 100, seed = 42)
  obj_[[name]] <- obj_[[paste0(name,"1")]]; obj_[[paste0(name,"1")]] <- NULL
  message(sprintf("'%s': %d/%d genes present.", name, length(present), length(genes)))
  obj_
}

# ---------- 3) Gene sets ----------
uprmt_genes <- c("LONP1","CLPP","CLPX","HTRA2","HSPD1","HSPE1","HSPA9",
                 "DNAJA3","YME1L1","SPG7","AFG3L2","LRPPRC","ATF5","ATF4","DDIT3")

mitophagy_core <- c("PINK1","PRKN","SQSTM1","OPTN","CALCOCO2","TAX1BP1","NBR1","TBK1","TBKBP1")
mitophagy_receptors <- c("BNIP3","BNIP3L","FUNDC1","PHB2","FKBP8","BCL2L13","NLRX1","NIPSNAP1","NIPSNAP2","AMBRA1")
mito_ub_e3_dub <- c("MARCHF5","MUL1","RNF185","HUWE1","USP30","USP35","USP15")
mitophagy_genes <- unique(c(mitophagy_core, mitophagy_receptors, mito_ub_e3_dub))

oxphos_genes <- unique(c(
  paste0("NDUFA", 1:13), "NDUFAB1", paste0("NDUFB",1:11), "NDUFC1","NDUFC2",
  paste0("NDUFS",1:8), paste0("NDUFV",1:3),
  "SDHA","SDHB","SDHC","SDHD",
  "UQCRC1","UQCRC2","UQCRB","UQCRH","UQCR10","UQCR11","CYC1","UQCRQ",
  "COX4I1","COX4I2","COX5A","COX5B","COX6A1","COX6A2","COX6B1","COX6B2","COX6C",
  "COX7A1","COX7A2","COX7A2L","COX7B","COX7C","COX8A",
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
  "ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5MF","ATP5PO","ATP5PB","ATP5PD","ATP5PF","ATP5MD"
))

tca_genes <- c("PDHA1","PDHB","PDHX","DLAT","DLD","PC","CS","ACO2","IDH2","IDH3A","IDH3B","IDH3G",
               "OGDH","DLST","SUCLG1","SUCLG2","SUCLA2","SDHA","SDHB","SDHC","SDHD","FH","MDH2")

glyco_genes <- c("SLC2A1","SLC2A3","SLC2A2","HK1","HK2","HK3","GPI","PFKM","PFKL","PFKP",
                 "ALDOA","ALDOB","ALDOC","TPI1","GAPDH","PGK1","PGK2","PGAM1","PGAM2",
                 "ENO1","ENO2","ENO3","PKM","PKLR","LDHA","LDHB")

# ---------- 4) Axis scores (no percent.mt) ----------
DefaultAssay(obj) <- "RNA"
obj <- AddModuleScoreSafe(obj, uprmt_genes,    "UPRmtScore")
obj <- AddModuleScoreSafe(obj, mitophagy_genes,"Mitophagy")
obj <- AddModuleScoreSafe(obj, oxphos_genes,   "OXPHOS1")
obj <- AddModuleScoreSafe(obj, tca_genes,      "TCA1")
obj <- AddModuleScoreSafe(obj, glyco_genes,    "Glyco1")

# ---------- 5) MQI (exclude percent.mt entirely) ----------
compute_MQI_noMT <- function(md,
                             weights = list(mitophagy=+1.0, oxphos=+1.0, uprmt=+0.7,
                                            tca_pen=+0.2, imbalance=+0.5)) {
  z_mito <- z0(md$Mitophagy)
  z_ox   <- z0(md$OXPHOS1)
  z_upr  <- z0(md$UPRmtScore)
  z_tca  <- z0(md$TCA1)
  z_gly  <- z0(md$Glyco1)
  # imbalance penalty between OXPHOS and Glycolysis
  imb    <- z0(abs(z_ox - z_gly))
  
  good <- weights$mitophagy * z_mito +
    weights$oxphos    * z_ox   +
    weights$uprmt     * z_upr
  
  pen  <- 0
  pen  <- pen + weights$tca_pen  * pmax(0, z_tca)   # only penalize positive TCA z-scores
  pen  <- pen + weights$imbalance* imb
  
  as.numeric(good - pen)
}

mq <- compute_MQI_noMT(obj@meta.data)
names(mq) <- colnames(obj)
obj <- AddMetaData(obj, mq, col.name = "MQI_v1")  # create/overwrite MQI_v1

# ---------- 6) Donor×Sex×Disease aggregation (McDonald β1/β3/β4) ----------
# Electrophysiology variables to analyze
electro_vars <- c(
  "CellSize_pF",
  "NormalizedTotalCapacitance_fF.pF",
  "NormalizedLateDepolarizationCapacitance",
  "CalciumIntegralNormalizedtoCellSize_pC.pF",
  "NormalizedPeakSodiumCurrentAmplitude_pA.pF",
  "HalfInactivationSodiumCurrent_mV",
  "NormalizedEarlyPeakCaCurrentAmplitude_pA.pF",
  "NormalizedLateCaCurrentAmplitude_pA.pF"
)
missing_ephys <- setdiff(electro_vars, colnames(obj@meta.data))
if (length(missing_ephys)) stop("Missing electrophysiology columns: ", paste(missing_ephys, collapse = ", "))

md <- obj@meta.data

# Carry Diabetes.status so we can map shapes in plots
df_donors <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster %in% c("β1","β3","β4"),
         !is.na(Sex)) %>%
  select(Donor.ID, Sex, Diabetes.status, MQI_v1, all_of(electro_vars)) %>%
  mutate(Diabetes.status = dplyr::case_when(
    is.na(Diabetes.status) ~ "Unknown",
    TRUE ~ as.character(Diabetes.status)
  ))

# Donor-mean per Sex × Diabetes.status (if a donor appears once per sex, this is 1 row per donor×sex)
df_corr <- df_donors %>%
  group_by(Donor.ID, Sex, Diabetes.status) %>%
  summarise(across(c(MQI_v1, all_of(electro_vars)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(is.finite(MQI_v1)) %>%
  mutate(Diabetes.status = factor(Diabetes.status, levels = c("ND","T2D","Unknown")))

# ---------- 7) Spearman correlations (per Sex, pooling disease) ----------
spearman_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(tibble(rho = NA_real_, p_value = NA_real_))
  out <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(rho = unname(out$estimate), p_value = out$p.value)
}

# sanity: donors per sex
print(table(df_corr$Sex, useNA = "ifany"))

corr_results <- map_dfr(electro_vars, function(v) {
  # F
  res_F <- spearman_one(
    df_corr$MQI_v1[df_corr$Sex == "F"],
    df_corr[[v]][df_corr$Sex == "F"]
  ) %>% mutate(variable = v, Sex = "F")
  
  # M
  res_M <- spearman_one(
    df_corr$MQI_v1[df_corr$Sex == "M"],
    df_corr[[v]][df_corr$Sex == "M"]
  ) %>% mutate(variable = v, Sex = "M")
  
  bind_rows(res_F, res_M)
}) %>%
  group_by(Sex) %>%                        # FDR within each sex (change if you prefer global BH)
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(sig = case_when(
    is.na(p_value)  ~ "",
    p_value < 1e-3  ~ "***",
    p_value < 1e-2  ~ "**",
    p_value < 5e-2  ~ "*",
    TRUE            ~ ""
  ))

print(corr_results)

# ---------- 8) Plot: shaded CI + triangle shape for T2D ----------
df_long <- df_corr %>%
  pivot_longer(-c(Donor.ID, Sex, Diabetes.status, MQI_v1), names_to = "Feature", values_to = "Value")

# per-panel labels (ρ, p) per Feature × Sex
lab_df <- corr_results %>%
  transmute(Feature = variable, Sex,
            label = paste0(Sex, ": ρ=", ifelse(is.na(rho), "NA", sprintf("%.2f", rho)),
                           ", p=", ifelse(is.na(p_value), "NA", formatC(p_value, format="e", digits=2)),
                           ifelse(sig=="","", paste0(" ", sig))))

# label positions
pos_df <- df_long %>%
  group_by(Feature, Sex) %>%
  summarise(x_min = min(MQI_v1, na.rm = TRUE),
            y_max = max(Value,  na.rm = TRUE), .groups = "drop") %>%
  left_join(lab_df, by = c("Feature","Sex"))

p <- ggplot(df_long, aes(x = MQI_v1, y = Value,
                         color = Sex, shape = Diabetes.status)) +
  geom_point(size = 2.6, alpha = 0.9) +
  # shaded 95% CI for sex-specific lines (pooled over disease)
  geom_smooth(aes(group = Sex), method = "lm", se = TRUE, alpha = 0.15) +
  # on-panel stats labels
  geom_text(
    data = pos_df,
    aes(x = x_min, y = y_max, label = label, color = Sex),
    hjust = 0, vjust = 1, size = 3.1, show.legend = FALSE,
    inherit.aes = FALSE
  )+
  facet_wrap(~ Feature, scales = "free_y") +
  scale_shape_manual(values = c("ND" = 16, "T2D" = 17, "Unknown" = 1), drop = FALSE) +
  labs(title = "Donor×Sex correlations: MQI_v1 (no %mt) vs electrophysiology (McDonald β)",
       x = "MQI_v1 (donor mean)",
       y = "Electrophysiology (donor mean)",
       shape = "Diabetes", color = "Sex") +
  theme_classic(base_size = 12)

print(p)
df_long$Diabetes.status

qsave(obj, r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined_mcdonald_correctmetadat.qs)")
# ---------- 9) (Optional) Save ----------
# ggsave("donor_sexsplit_MQI_noMT_vs_ephys_facets.pdf", p, width = 11, height = 8.5, useDingbats = FALSE)

# ---------- 10) (Optional) If you want per-Sex×Disease correlations too ----------
# Uncomment to compute separate ρ/p for each Feature within each Sex×Disease stratum.
# strat_results <- map_dfr(electro_vars, function(v) {
#   df_tmp <- df_corr %>%
#     select(Sex, Diabetes.status, MQI_v1, v = all_of(v))
#   df_tmp %>%
#     group_by(Sex, Diabetes.status) %>%
#     group_modify(~ {
#       ok <- is.finite(.x$MQI_v1) & is.finite(.x$v)
#       if (sum(ok) < 3) return(tibble(variable = v, rho = NA_real_, p_value = NA_real_))
#       out <- suppressWarnings(cor.test(.x$MQI_v1[ok], .x$v[ok], method = "spearman", exact = FALSE))
#       tibble(variable = v, rho = unname(out$estimate), p_value = out$p.value)
#     }) %>% ungroup()
# }) %>%
#   group_by(Sex, Diabetes.status) %>%
#   mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
#   ungroup()
# print(strat_results)

# Load Objec
combined_ps <- qread(file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Algo.Ai\DATA\seurat_objects\combined_mcdonald_correctmetadat.qs)")
sum(!is.na(combined_ps@meta.data[["MQI_pred"]]))
sum(!is.na(combined_ps@meta.data[["MQI_v1"]]))
table(combined_ps@meta.data[["dataset_origin"]])
table(combined_ps@meta.data[["Glucose_mM"]])

# A) Electrophysiology heatmap — All β vs β1, split by Sex × Diabetes × Glucose
# ============================================================
#   Heatmap of electrophysiology features in McDonald β cells
#   Panels: All β combined vs β1 only
#   Splits: Sex × Diabetes (ND/T2D) × Glucose (1/5/10 mM)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

md <- combined_ps@meta.data

# Order factors (optional but makes panels consistent)
md <- md %>%
  mutate(
    Diabetes.status = factor(Diabetes.status, levels = c("ND","T2D")),
    Sex             = factor(Sex, levels = c("F","M")),
    Glucose_mM      = factor(Glucose_mM, levels = c(1,5,10), labels = c("1 mM","5 mM","10 mM"))
  )

# --- 1) All β combined ---
df_md_all <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster %in% c("β1","β3","β4")) %>%
  mutate(beta_cluster = "All β") %>%
  select(beta_cluster, Sex, Diabetes.status, Glucose_mM,
         CellSize_pF,
         NormalizedTotalCapacitance_fF.pF,
         NormalizedLateDepolarizationCapacitance,
         CalciumIntegralNormalizedtoCellSize_pC.pF,
         NormalizedPeakSodiumCurrentAmplitude_pA.pF,
         HalfInactivationSodiumCurrent_mV,
         NormalizedEarlyPeakCaCurrentAmplitude_pA.pF,
         NormalizedLateCaCurrentAmplitude_pA.pF)

# --- 2) β1 only ---
df_md_beta1 <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster == "β1") %>%
  select(beta_cluster, Sex, Diabetes.status, Glucose_mM,
         CellSize_pF,
         NormalizedTotalCapacitance_fF.pF,
         NormalizedLateDepolarizationCapacitance,
         CalciumIntegralNormalizedtoCellSize_pC.pF,
         NormalizedPeakSodiumCurrentAmplitude_pA.pF,
         HalfInactivationSodiumCurrent_mV,
         NormalizedEarlyPeakCaCurrentAmplitude_pA.pF,
         NormalizedLateCaCurrentAmplitude_pA.pF)

# --- 3) Combine both datasets ---
df_md <- bind_rows(df_md_all, df_md_beta1)

# --- 4) Reshape to long format ---
df_long <- df_md %>%
  pivot_longer(-c(beta_cluster, Sex, Diabetes.status, Glucose_mM),
               names_to = "variable", values_to = "value")

# --- 5) Compute means (by Glucose_mM as well) ---
summary_stats <- df_long %>%
  group_by(beta_cluster, variable, Sex, Diabetes.status, Glucose_mM) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            n = dplyr::n(), .groups = "drop")

# Scale within each feature across *all* panels (keeps comparability)
heatmap_data <- summary_stats %>%
  group_by(variable) %>%
  mutate(scaled = as.numeric(scale(mean))) %>%
  ungroup() %>%
  mutate(variable_sex = paste0(variable, "_", Sex))

# --- Optional: ND vs T2D Wilcoxon within Sex × Glucose × beta_cluster ---
min_n <- 3
stats_table <- df_long %>%
  group_by(beta_cluster, variable, Sex, Glucose_mM) %>%
  summarise(
    n_ND  = sum(Diabetes.status == "ND"  & is.finite(value)),
    n_T2D = sum(Diabetes.status == "T2D" & is.finite(value)),
    p_value = if (n_ND >= min_n && n_T2D >= min_n)
      tryCatch(wilcox.test(value ~ Diabetes.status)$p.value, error = function(e) NA_real_)
    else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(sig = case_when(
    is.finite(p_value) & p_value < 0.001 ~ "***",
    is.finite(p_value) & p_value < 0.01  ~ "**",
    is.finite(p_value) & p_value < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  # place the star on the T2D column tile
  mutate(Diabetes.status = factor("T2D", levels = c("ND","T2D")))

# --- 6) Plot heatmap (facet by Glucose and beta_cluster) ---
ggplot(heatmap_data,
       aes(x = Diabetes.status, y = variable_sex, fill = scaled)) +
  geom_tile(color = "white") +
  facet_grid(Glucose_mM ~ beta_cluster, scales = "free_x") +
  scale_fill_gradient2(
    low = "dodgerblue4", mid = "white", high = "firebrick4", midpoint = 0,
    limits = c(-1, 1), oob = scales::squish,
    breaks = c(-1, -0.5, 0, 0.5, 1), name = "Scaled mean"
  ) +
  labs(title = "Electrophysiology features — McDonald patch-seq (by glucose)",
       x = "Diabetes status", y = "Feature × Sex",
       fill = "Scaled mean") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # significance stars on T2D tiles (optional)
  geom_text(data = stats_table,
            aes(x = Diabetes.status, y = paste0(variable, "_", Sex), label = sig),
            inherit.aes = FALSE, vjust = 0.5, hjust = 0.5, size = 3)

# B) Electrophysiology heatmap — PINK1+ vs PINK1-, split by Sex × Diabetes × Glucose
# ============================================================
#   Electrophysiology features — PINK1+ vs PINK1- (sex pooled)
#   Panels: All β combined vs β1 only
#   Splits: Diabetes (ND/T2D) × Glucose (1/5/10 mM)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

md <- combined_ps@meta.data

# Ensure PINK1 status column (RNA counts > 0)
if (!"PINK1_status" %in% colnames(md)) {
  md$PINK1_status <- ifelse(combined_ps@assays$RNA@data["PINK1", ] > 0, "PINK1+", "PINK1-")
}

# Order factors
md <- md %>%
  mutate(
    Diabetes.status = factor(Diabetes.status, levels = c("ND","T2D")),
    Glucose_mM      = factor(Glucose_mM, levels = c(1,5,10), labels = c("1 mM","5 mM","10 mM")),
    PINK1_status    = factor(PINK1_status, levels = c("PINK1-","PINK1+"))
  )

# --- 1) Subset McDonald β1 + All β (collapse β1+β3+β4) ---
df_md_all <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster %in% c("β1","β3","β4")) %>%
  mutate(beta_cluster = "All β")

df_md_beta1 <- md %>%
  filter(dataset_origin == "McDonald_patchseq",
         beta_cluster == "β1")

# Merge them (sex is intentionally NOT kept; pooled downstream)
df_md <- bind_rows(df_md_all, df_md_beta1) %>%
  select(beta_cluster, Diabetes.status, Glucose_mM, PINK1_status,
         CellSize_pF,
         NormalizedTotalCapacitance_fF.pF,
         NormalizedLateDepolarizationCapacitance,
         CalciumIntegralNormalizedtoCellSize_pC.pF,
         NormalizedPeakSodiumCurrentAmplitude_pA.pF,
         HalfInactivationSodiumCurrent_mV,
         NormalizedEarlyPeakCaCurrentAmplitude_pA.pF,
         NormalizedLateCaCurrentAmplitude_pA.pF)

# --- 2) Long format ---
df_long <- df_md %>%
  pivot_longer(-c(beta_cluster, Diabetes.status, Glucose_mM, PINK1_status),
               names_to = "variable", values_to = "value") %>%
  filter(is.finite(value))

# --- 3) Group means (sex pooled) ---
summary_stats <- df_long %>%
  group_by(beta_cluster, variable, Diabetes.status, Glucose_mM, PINK1_status) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            n = dplyr::n(), .groups = "drop")

# Scale per feature across all panels (keeps colors comparable)
heatmap_data <- summary_stats %>%
  group_by(variable) %>%
  mutate(scaled = as.numeric(scale(mean))) %>%
  ungroup()

# --- 4) PINK1+ vs PINK1- Wilcoxon within each panel (sex pooled) ---
min_n <- 3
stats_table <- df_long %>%
  group_by(beta_cluster, variable, Diabetes.status, Glucose_mM) %>%
  summarise(
    n_pos = sum(PINK1_status == "PINK1+"),
    n_neg = sum(PINK1_status == "PINK1-"),
    p_value = if (n_pos >= min_n && n_neg >= min_n)
      tryCatch(wilcox.test(value ~ PINK1_status)$p.value, error = function(e) NA_real_)
    else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(sig = case_when(
    is.finite(p_value) & p_value < 0.001 ~ "***",
    is.finite(p_value) & p_value < 0.01  ~ "**",
    is.finite(p_value) & p_value < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  # place star on the PINK1+ tile
  mutate(PINK1_status = factor("PINK1+", levels = c("PINK1-","PINK1+")))

# --- 5) Plot ---
ggplot(heatmap_data,
       aes(x = PINK1_status, y = variable, fill = scaled)) +
  geom_tile(color = "white") +
  facet_grid(Glucose_mM ~ beta_cluster + Diabetes.status, scales = "free_x") +
  scale_fill_gradient2(
    low = "dodgerblue4", mid = "white", high = "firebrick4", midpoint = 0,
    limits = c(-1, 1), oob = scales::squish,
    breaks = c(-1, -0.5, 0, 0.5, 1), name = "Scaled mean"
  ) +
  labs(title = "Electrophysiology features — PINK1+ vs PINK1- (sex pooled, by glucose)",
       x = "PINK1 status", y = "Feature",
       fill = "Scaled mean") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # significance stars on PINK1+ tiles (optional)
  geom_text(data = stats_table,
            aes(x = PINK1_status, y = variable, label = sig),
            inherit.aes = FALSE, vjust = 0.5, hjust = 0.5, size = 3)


# ============================================================
#   E-phys heatmap: HighMQI & PINK1+  vs  LowMQI & PINK1−
#   Sex pooled; faceted by Glucose (rows) × Diabetes (cols)
#   Fill = z-scored group means (for contrast)
#   Overlay = asterisks from RAW Wilcoxon p-values (least stringent)
#   Star key: * p<0.10, ** p<0.05, *** p<0.01, **** p<0.001
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(forcats); library(tibble); library(Seurat)
})

# ------------------ Data: McDonald patch-seq only ------------------
obj  <- combined_ps
meta <- obj@meta.data %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::filter(dataset_origin == "McDonald_patchseq")

stopifnot(all(c("Diabetes.status","Glucose_mM","MQI_v1") %in% names(meta)))
mqi_var <- "MQI_v1"

# --- Build PINK1_status aligned to the McDonald cells ---
if (!"PINK1_status" %in% names(meta)) {
  assay_name <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
  stopifnot("PINK1 not found in the chosen assay" = "PINK1" %in% rownames(obj[[assay_name]]@data))
  pink_df <- FetchData(obj, vars = "PINK1", cells = meta$cell, slot = "data")
  meta$PINK1_status <- ifelse(pink_df$PINK1 > 0, "PINK1+", "PINK1-")
}

# --- Factors (sex pooled) ---
meta <- meta %>%
  mutate(
    Diabetes.status = factor(Diabetes.status, levels = c("ND","T2D")),
    Glucose_mM      = factor(Glucose_mM, levels = c(1,5,10),
                             labels = c("1 mM","5 mM","10 mM")),
    PINK1_status    = factor(PINK1_status, levels = c("PINK1-","PINK1+"))
  )

# ---- Electrophysiology features ----
pattern <- "(CellSize_pF|Calcium|Sodium|CaCurrent|Inactivation|Capacitance|Integral|HalfInactivation)"
num_cols   <- names(meta)[vapply(meta, is.numeric, logical(1))]
ephys_cols <- grep(pattern, num_cols, value = TRUE)
ephys_cols <- setdiff(ephys_cols, mqi_var)
non_na_counts <- vapply(ephys_cols, function(cn) sum(is.finite(meta[[cn]])), integer(1))
ephys_cols    <- ephys_cols[non_na_counts >= 20]
stopifnot("No e-phys columns left after filtering!" = length(ephys_cols) > 0)

# ---- Long format ----
long <- meta %>%
  select(cell, Diabetes.status, Glucose_mM, PINK1_status,
         all_of(mqi_var), all_of(ephys_cols)) %>%
  pivot_longer(all_of(ephys_cols), names_to = "feature", values_to = "value") %>%
  filter(is.finite(value), is.finite(.data[[mqi_var]]))

# ---- MQI quartiles within Diabetes × Glucose (sex pooled) ----
long <- long %>%
  group_by(Diabetes.status, Glucose_mM) %>%
  mutate(
    Q25 = quantile(.data[[mqi_var]], 0.45, na.rm = TRUE),
    Q75 = quantile(.data[[mqi_var]], 0.55, na.rm = TRUE),
    MQI_bin = case_when(
      .data[[mqi_var]] <= Q25 ~ "Low",
      .data[[mqi_var]] >= Q75 ~ "High",
      TRUE ~ "Mid"
    )
  ) %>%
  ungroup()

# ---- Keep only the two target groups ----
target <- long %>%
  filter(
    (MQI_bin == "High" & PINK1_status == "PINK1+") |
      (MQI_bin == "Low"  & PINK1_status == "PINK1-")
  ) %>%
  mutate(
    Group = ifelse(MQI_bin == "High", "HighMQI&PINK1+", "LowMQI&PINK1-"),
    Group = factor(Group, levels = c("LowMQI&PINK1-","HighMQI&PINK1+"))
  )

# ---- Compute means (for fill) ----
summary_stats <- target %>%
  group_by(feature, Diabetes.status, Glucose_mM, Group) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            n    = dplyr::n(),
            .groups = "drop")

# z-score the means per feature across all panels (color comparable)
heatmap_data <- summary_stats %>%
  group_by(feature) %>%
  mutate(scaled = as.numeric(scale(mean))) %>%
  ungroup()

# ---- Wilcoxon tests: (HighMQI&PINK1+) vs (LowMQI&PINK1-) within each Diabetes × Glucose × feature ----
min_n <- 5  # per group
p_df <- target %>%
  group_by(feature, Diabetes.status, Glucose_mM) %>%
  summarise(
    n_hi = sum(Group == "HighMQI&PINK1+"),
    n_lo = sum(Group == "LowMQI&PINK1-"),
    p    = {
      if (n_hi >= min_n && n_lo >= min_n) {
        x <- value[Group == "HighMQI&PINK1+"]
        y <- value[Group == "LowMQI&PINK1-"]
        wilcox.test(x, y)$p.value
      } else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    # "least stringent": include p<0.10 as *
    stars = case_when(
      is.finite(p) & p < 0.001 ~ "****",
      is.finite(p) & p < 0.01  ~ "***",
      is.finite(p) & p < 0.05  ~ "**",
      is.finite(p) & p < 0.10  ~ "*",
      TRUE ~ ""
    ),
    # place star on the HighMQI&PINK1+ tile
    Group = factor("HighMQI&PINK1+", levels = c("LowMQI&PINK1-","HighMQI&PINK1+"))
  )

# ---- Plot: two columns per panel (Low vs High), stars on High tile ----
ggplot(heatmap_data,
       aes(x = Group, y = feature, fill = scaled)) +
  geom_tile(color = "white") +
  facet_grid(Glucose_mM ~ Diabetes.status, scales = "free_x") +
  scale_fill_gradient2(
    low = "dodgerblue4", mid = "white", high = "firebrick4", midpoint = 0,
    limits = c(-1, 1), oob = scales::squish,
    breaks = c(-1, -0.5, 0, 0.5, 1), name = "Scaled mean"
  ) +
  labs(
    title = "E-phys: HighMQI & PINK1+  vs  LowMQI & PINK1−\nSex pooled, split by glucose and diabetes",
    x = NULL, y = "Electrophysiology feature",
    caption = "Stars = raw Wilcoxon p-values (High vs Low): * p<0.10, ** p<0.05, *** p<0.01, **** p<0.001"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold")) +
  geom_text(data = p_df,
            aes(x = Group, y = feature, label = stars),
            inherit.aes = FALSE, vjust = 0.5, hjust = 0.5, size = 3)



# ============================================================
#   E-phys heatmap: High-MQI vs Low-MQI (no PINK1), sex pooled
#   Facets: Glucose (rows) × Diabetes (cols)
#   Fill = z-scored group means; Stars = RAW Wilcoxon p-values
#   Star key: * p<0.10, ** p<0.05, *** p<0.01, **** p<0.001
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(forcats); library(tibble); library(Seurat)
})

# ------------------ Data: McDonald patch-seq only ------------------
obj  <- combined_ps
meta <- obj@meta.data %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::filter(dataset_origin == "McDonald_patchseq")

stopifnot(all(c("Diabetes.status","Glucose_mM","MQI_v1") %in% names(meta)))
mqi_var <- "MQI_v1"

# Factors (sex pooled)
meta <- meta %>%
  mutate(
    Diabetes.status = factor(Diabetes.status, levels = c("ND","T2D")),
    Glucose_mM      = factor(Glucose_mM, levels = c(1,5,10),
                             labels = c("1 mM","5 mM","10 mM"))
  )

# ---- Electrophysiology features ----
pattern <- "(CellSize_pF|Calcium|Sodium|CaCurrent|Inactivation|Capacitance|Integral|HalfInactivation)"
num_cols   <- names(meta)[vapply(meta, is.numeric, logical(1))]
ephys_cols <- grep(pattern, num_cols, value = TRUE)
ephys_cols <- setdiff(ephys_cols, mqi_var)
non_na_counts <- vapply(ephys_cols, function(cn) sum(is.finite(meta[[cn]])), integer(1))
ephys_cols    <- ephys_cols[non_na_counts >= 20]
stopifnot("No e-phys columns left after filtering!" = length(ephys_cols) > 0)

# ---- Long format (sex pooled) ----
long <- meta %>%
  select(cell, Diabetes.status, Glucose_mM, all_of(mqi_var), all_of(ephys_cols)) %>%
  pivot_longer(all_of(ephys_cols), names_to = "feature", values_to = "value") %>%
  filter(is.finite(value), is.finite(.data[[mqi_var]]))

# ---- MQI quartiles within Diabetes × Glucose ----
# Set cutpoints here (default 0.25/0.75). If you want narrower bands, change to 0.45/0.55.
q_low  <- 0.25
q_high <- 0.75

long <- long %>%
  group_by(Diabetes.status, Glucose_mM) %>%
  mutate(
    QL = quantile(.data[[mqi_var]], q_low,  na.rm = TRUE),
    QH = quantile(.data[[mqi_var]], q_high, na.rm = TRUE),
    MQI_bin = case_when(
      .data[[mqi_var]] <= QL ~ "LowMQI",
      .data[[mqi_var]] >= QH ~ "HighMQI",
      TRUE ~ "Mid"
    )
  ) %>%
  ungroup()

# ---- Keep only High vs Low MQI groups ----
target <- long %>%
  filter(MQI_bin %in% c("HighMQI","LowMQI")) %>%
  mutate(
    Group = factor(MQI_bin, levels = c("LowMQI","HighMQI"))
  )

# ---- Compute means (for fill) ----
summary_stats <- target %>%
  group_by(feature, Diabetes.status, Glucose_mM, Group) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            n    = dplyr::n(),
            .groups = "drop")

# z-score the means per feature across all panels (for comparable colors)
heatmap_data <- summary_stats %>%
  group_by(feature) %>%
  mutate(scaled = as.numeric(scale(mean))) %>%
  ungroup()

# ---- Wilcoxon tests: HighMQI vs LowMQI within each Diabetes × Glucose × feature ----
min_n <- 5  # per group
p_df <- target %>%
  group_by(feature, Diabetes.status, Glucose_mM) %>%
  summarise(
    n_hi = sum(Group == "HighMQI"),
    n_lo = sum(Group == "LowMQI"),
    p    = {
      if (n_hi >= min_n && n_lo >= min_n) {
        x <- value[Group == "HighMQI"]
        y <- value[Group == "LowMQI"]
        wilcox.test(x, y)$p.value
      } else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    stars = case_when(
      is.finite(p) & p < 0.001 ~ "****",
      is.finite(p) & p < 0.01  ~ "***",
      is.finite(p) & p < 0.05  ~ "**",
      is.finite(p) & p < 0.10  ~ "*",
      TRUE ~ ""
    ),
    # place star on the HighMQI tile
    Group = factor("HighMQI", levels = c("LowMQI","HighMQI"))
  )

# ---- Plot: two columns per panel (LowMQI vs HighMQI), stars on High tile ----
ggplot(heatmap_data,
       aes(x = Group, y = feature, fill = scaled)) +
  geom_tile(color = "white") +
  facet_grid(Glucose_mM ~ Diabetes.status, scales = "free_x") +
  scale_fill_gradient2(
    low = "dodgerblue4", mid = "white", high = "firebrick4", midpoint = 0,
    limits = c(-1, 1), oob = scales::squish,
    breaks = c(-1, -0.5, 0, 0.5, 1), name = "Scaled mean"
  ) +
  labs(
    title = "E-phys features: High-MQI vs Low-MQI (sex pooled)\nStratified by bath glucose and diabetes status",
    x = NULL, y = "Electrophysiology feature",
    caption = "Stars = raw Wilcoxon p-values (High vs Low): * p<0.10, ** p<0.05, *** p<0.01, **** p<0.001"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold")) +
  geom_text(data = p_df,
            aes(x = Group, y = feature, label = stars),
            inherit.aes = FALSE, vjust = 0.5, hjust = 0.5, size = 3)
# ggsave("ephys_means_High_vs_Low_MQI_by_glucose_diabetes.pdf", width = 12, height = 9)
