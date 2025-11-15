
if (!dir.exists("_raw")) dir.create("_raw", showWarnings = FALSE)

if (!require("cBioPortalData", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("cBioPortalData")
}
library(cBioPortalData)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Genes HDAC/SIRT
genes <- c(
  "HDAC1","HDAC2","HDAC3","HDAC8",
  "HDAC4","HDAC5","HDAC7","HDAC9",
  "HDAC6","HDAC10",
  "SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7",
  "HDAC11"
)

# selected studies
studies <- c(
  "brca_tcga_pan_can_atlas_2018",
  "luad_tcga_pan_can_atlas_2018",
  "coadread_tcga_pan_can_atlas_2018"
)

study_tbl <- tibble(
  studyId   = studies,
  profileId = paste0(studies, "_rna_seq_v2_mrna")
)

cbio <- cBioPortal()

# download genical expression
study_tbl %>%
  mutate(
    expr = map2(studyId, profileId, ~ {
      message("Scarico espressione per ", .x)
      se <- getDataByGenes(
        api = cbio,
        studyId = .x,
        genes = genes,
        molecularProfileIds = .y,
        by = "hugoGeneSymbol"
      )[[1]]
      df <- as.data.frame(se)
      df %>% select(patientId, hugoGeneSymbol, value)
    }),
    file = paste0("_raw/expr_", studyId, ".csv"),
    written = map2(expr, file, write_csv)
  ) %>% invisible()

#download clinical data
study_tbl %>%
  mutate(
    clin = map(studyId, ~ {
      message("Scarico clinici per ", .x)
      as.data.frame(clinicalData(cbio, .x))
    }),
    file = paste0("_raw/clinical_", studyId, ".csv"),
    written = map2(clin, file, write_csv)
  ) %>% invisible()

# load data
load_study_data <- function(study) {
  expr <- read_csv(paste0("_raw/expr_", study, ".csv"), show_col_types = FALSE)
  clin <- read_csv(paste0("_raw/clinical_", study, ".csv"), show_col_types = FALSE)

  mr_col <- grep("mRNA|stem|Stemness", names(clin), value = TRUE, ignore.case = TRUE)
  if (length(mr_col) == 0) {
    message("⚠️ mRNAsi non trovato per ", study)
    return(NULL)
  }

  clin <- clin %>% rename(mRNAsi = all_of(mr_col[1]))
  expr_wide <- expr %>% pivot_wider(names_from = hugoGeneSymbol, values_from = value)
  expr_wide %>% inner_join(clin, by = "patientId")
}

# Correlations
compute_correlations <- function(df, study) {
  map_dfr(genes, function(g) {
    if (!g %in% names(df)) return(NULL)
    test <- suppressWarnings(cor.test(df[[g]], df$mRNAsi, method = "spearman"))
    tibble(study = study, gene = g, rho = test$estimate, p = test$p.value)
  })
}

# Loop on valid studies
all_results <- map_dfr(studies, ~{
  df <- load_study_data(.x)
  if (is.null(df)) return(NULL)
  compute_correlations(df, .x)
})

# Correction
all_results <- all_results %>%
  group_by(gene) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

# Matrixes for heatmap
cor_mat <- all_results %>%
  select(gene, study, rho) %>%
  pivot_wider(names_from = study, values_from = rho) %>%
  column_to_rownames("gene") %>%
  as.matrix()

q_mat <- all_results %>%
  select(gene, study, q) %>%
  pivot_wider(names_from = study, values_from = q) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Symbols
pch_mat <- matrix(NA, nrow = nrow(q_mat), ncol = ncol(q_mat))
for(i in seq_len(nrow(q_mat))) {
  for(j in seq_len(ncol(q_mat))) {
    pch_mat[i,j] <- case_when(
      q_mat[i,j] < 0.001 ~ 16,  # circle
      q_mat[i,j] < 0.01  ~ 17,  # triangle
      q_mat[i,j] < 0.05  ~ 15,  # square
      TRUE               ~ NA_real_
    )
  }
}

# Classes HDAC
hdac_class <- tibble(
  gene = rownames(cor_mat),
  class = case_when(
    gene %in% c("HDAC1","HDAC2","HDAC3","HDAC8") ~ "class I",
    gene %in% c("HDAC4","HDAC5","HDAC7","HDAC9") ~ "class IIA",
    gene %in% c("HDAC6","HDAC10") ~ "class IIB",
    gene %in% c("SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7") ~ "class III",
    gene == "HDAC11" ~ "class IV",
    TRUE ~ "other"
  )
)

class_colors <- c(
  "class I"   = "#F8766D",
  "class IIA" = "#E68613",
  "class IIB" = "#00BA38",
  "class III" = "#619CFF",
  "class IV"  = "#C77CFF"
)

row_ha <- rowAnnotation(
  HDAC_class = hdac_class$class,
  col = list(HDAC_class = class_colors),
  width = unit(5, "mm")
)

# Figure 4A – Heatmap
col_fun <- colorRamp2(c(-0.7, 0, 0.8), c("#313695", "white", "#A50026"))

Heatmap(
  cor_mat,
  name = "rho",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  left_annotation = row_ha,
  cell_fun = function(j, i, x, y, w, h, fill) {
    pch <- pch_mat[i,j]
    if(!is.na(pch)) {
      grid.points(x, y, pch = pch, size = unit(3.5, "mm"))
    }
  },
  heatmap_legend_param = list(
    title = "Spearman rho",
    legend_height = unit(4, "cm")
  )
)

# Figure 4B – Barplot
figB <- all_results %>%
  mutate(
    sign = case_when(
      q > 0.05 ~ "ns",
      rho > 0  ~ "positive",
      TRUE     ~ "negative"
    )
  ) %>%
  group_by(gene, sign) %>%
  summarise(n = n(), .groups = "drop")

ggplot(figB, aes(x = n, y = gene, fill = sign)) +
  geom_col() +
  scale_fill_manual(
    values = c(positive = "red", negative = "blue", ns = "grey80")
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Number of TCGA cohorts",
    y = "",
    title = "Number of TCGA cohorts with positive or negative correlation"
  )
