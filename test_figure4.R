
# --- Load libraries ---
library(tidyverse)
library(pheatmap)
library(tidyr)
library(broom)
library(grid)

# --- Load data ---
my_data <- read_csv("../data/02_dat_clean.csv")
# view(my_data)  # optional

# --- Encode Grade as numeric (malignancy scale) ---
# G1 / Low Grade = 1 (low malignancy)
# G2             = 2 (medium malignancy)
# G3 / High Grade= 3 (high malignancy)
# G4             = 4 (highest malignancy)
changed_data <- my_data %>%
  mutate(
    Grade_in_numbers = case_when(
      Grade == "G1" ~ 1,
      Grade == "Low Grade" ~ 1,
      Grade == "G2" ~ 2,
      Grade == "G3" ~ 3,
      Grade == "High Grade" ~ 3,
      Grade == "G4" ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  drop_na(Grade_in_numbers)

# --- Recode Cancer study IDs to short labels (to match figure) ---
changed_data$Cancer <- dplyr::recode(changed_data$Cancer,
                                     "blca_tcga_pan_can_atlas_2018"     = "BLCA",
                                     "brca_tcga_pan_can_atlas_2018"     = "BRCA",
                                     "cesc_tcga_pan_can_atlas_2018"     = "CESC",
                                     "coadread_tcga_pan_can_atlas_2018" = "COAD",
                                     "esca_tcga_pan_can_atlas_2018"     = "ESCA",
                                     "gbm_tcga_pan_can_atlas_2018"      = "GBM",
                                     "hnsc_tcga_pan_can_atlas_2018"     = "HNSC",
                                     "kirc_tcga_pan_can_atlas_2018"     = "KIRC",
                                     "kirp_tcga_pan_can_atlas_2018"     = "KIRP",
                                     "lgg_tcga_pan_can_atlas_2018"      = "LGG",
                                     "lihc_tcga_pan_can_atlas_2018"     = "LIHC",
                                     "luad_tcga_pan_can_atlas_2018"     = "LUAD",
                                     "lusc_tcga_pan_can_atlas_2018"     = "LUSC",
                                     "ov_tcga_pan_can_atlas_2018"       = "OV",
                                     "paad_tcga_pan_can_atlas_2018"     = "PAAD",
                                     "prad_tcga_pan_can_atlas_2018"     = "PRAD",
                                     "sarc_tcga_pan_can_atlas_2018"     = "SARC",
                                     "skcm_tcga_pan_can_atlas_2018"     = "SKCM",
                                     "stad_tcga_pan_can_atlas_2018"     = "STAD",
                                     "tgct_tcga_pan_can_atlas_2018"     = "TGCT",
                                     "thca_tcga_pan_can_atlas_2018"     = "THCA",
                                     "ucec_tcga_pan_can_atlas_2018"     = "UCEC"
)

# --- Fix the cancer order to match the figure (no clustering) ---
cancer_order <- c(
  "LIHC", "UCEC", "KIRC", "LGG",
  "OV", "PAAD",
  "HNSC", "ESCA", "CESC", "STAD",
  "COAD", "BLCA"
)
cancers_use <- intersect(cancer_order, unique(changed_data$Cancer))

changed_data <- changed_data %>%
  filter(Cancer %in% cancers_use) %>%
  mutate(Cancer = factor(Cancer, levels = cancers_use))

# --- Define genes in the plotting order (HDAC/SIRT classes) ---
genes_to_test <- c(
  # class I
  "HDAC3", "HDAC2", "HDAC1", "HDAC8",
  # class IIA
  "HDAC7", "HDAC9", "HDAC5", "HDAC4",
  # class III
  "SIRT1", "SIRT7", "SIRT6", "SIRT4", "SIRT3", "SIRT2", "SIRT5",
  # class IIB
  "HDAC6", "HDAC10",
  # class IV
  "HDAC11"
)

# --- Compute Spearman correlations and p-values (per gene x cancer) ---
s_correlation <- changed_data %>%
  pivot_longer(
    cols = any_of(genes_to_test),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(Cancer, gene) %>%
  summarize(
    summary_test = broom::tidy(cor.test(expression, Grade_in_numbers, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) %>%
  unnest(summary_test) %>%
  select(Cancer, gene, rho = estimate, p.value)

# --- Build matrices (rho and p) in the fixed order ---
matrix_rho <- s_correlation %>%
  select(gene, Cancer, rho) %>%
  pivot_wider(names_from = Cancer, values_from = rho) %>%
  column_to_rownames("gene") %>%
  as.matrix()

matrix_p <- s_correlation %>%
  select(gene, Cancer, p.value) %>%
  pivot_wider(names_from = Cancer, values_from = p.value) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# --- Filter and sort matrices to match desired gene/cancer order ---
final_genes   <- genes_to_test[genes_to_test %in% rownames(matrix_rho)]
final_cancers <- cancers_use[cancers_use %in% colnames(matrix_rho)]

matrix_rho    <- matrix_rho[final_genes, final_cancers, drop = FALSE]
matrix_p      <- matrix_p[final_genes, final_cancers, drop = FALSE]

# --- Replace NA correlations with 0 for plotting (optional) ---
matrix_rho_plot <- matrix_rho
matrix_rho_plot[is.na(matrix_rho_plot)] <- 0

# --- BH correction on p-values to get q-values (matrix form) ---
matrix_q <- matrix(
  p.adjust(as.vector(matrix_p), method = "BH"),
  nrow = nrow(matrix_p),
  ncol = ncol(matrix_p),
  byrow = FALSE,
  dimnames = dimnames(matrix_p)
)

# --- Build display annotations (symbols) based on q-values ---
# ■ for q < 0.001, ▲ for q < 0.01, • for q < 0.05, blank otherwise
display_matrix <- ifelse(matrix_q < 0.001, "■",
                         ifelse(matrix_q < 0.01,  "▲",
                                ifelse(matrix_q < 0.05,  "•", "")))

# --- Row annotation (HDAC classes) ---
gene_classes <- c(
  HDAC3="class I", HDAC2="class I", HDAC1="class I", HDAC8="class I",
  HDAC7="class IIA", HDAC9="class IIA", HDAC5="class IIA", HDAC4="class IIA",
  SIRT1="class III", SIRT7="class III", SIRT6="class III", SIRT4="class III",
  SIRT3="class III", SIRT2="class III", SIRT5="class III",
  HDAC6="class IIB", HDAC10="class IIB",
  HDAC11="class IV"
)

annot_row <- data.frame(
  HDAC_family = as.character(gene_classes[final_genes]),
  row.names = final_genes,
  stringsAsFactors = FALSE
)

annot_row$HDAC_family <- factor(
  annot_row$HDAC_family,
  levels = c("class I", "class IIA", "class IIB", "class III", "class IV")
)

ann_colors <- list(
  HDAC_family = c(
    "class I"   = "salmon",
    "class IIA" = "#E69F52",
    "class IIB" = "#BDB600",
    "class III" = "#009E60",
    "class IV"  = "#008ECE"
  )
)

# --- Heatmap color scale and breaks ---
my_breaks  <- seq(-0.3, 0.3, length.out = 101)
my_colours <- colorRampPalette(c("blue","white","red"))(100)

# --- Plot heatmap (no clustering, fixed order) ---
pheatmap(
  matrix_rho_plot,
  color = my_colours,
  breaks = my_breaks,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annot_row,
  annotation_colors = ann_colors,
  display_numbers = display_matrix,   # significance symbols (q-based)
  fontsize_number = 6,
  main = " ",
  number_color = "black",
  show_colnames = TRUE,
  show_rownames = TRUE,
  annotation_legend = TRUE,
  fontsize = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 13,
  cellheight = 10
)
