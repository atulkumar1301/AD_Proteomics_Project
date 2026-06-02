# 1. Environment Setup 
library(tidyverse)    # Includes dplyr, stringr, ggplot2, tidyr
library(matrixStats)
library(preprocessCore)
library(impute)

options(scipen = 999)

# Load Data
raw_data <- read.delim("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Mali_Proteomics_Data.txt", sep = "\t", check.names = FALSE)
raw_data[is.na(raw_data)] <- NA
anno_cols <- c("PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions")

# 1A & 1B. Isolate and Clean Matrices
quant_df <- raw_data %>%
  select(all_of(anno_cols), contains("PG.Quantity")) %>%
  rename_with(~ paste0("Sample_6", str_match(.x, "_(\\d+)\\.raw")[, 2]), contains("PG.Quantity"))

qual_df <- raw_data %>%
  select(all_of(anno_cols), contains("PG.NrOfModifiedSequencesIdentified")) %>%
  rename_with(~ paste0("Sample_6", str_match(.x, "_(\\d+)\\.raw")[, 2]), contains("PG.NrOfModifiedSequencesIdentified"))

# Convert to matrices
quant_mat <- as.matrix(quant_df[, -c(1:3)])
qual_mat  <- as.matrix(qual_df[, -c(1:3)])

# CRITICAL FIX: Use NA instead of NaN for downstream impute.knn compatibility
quant_mat[qual_mat < 2] <- NA 

# Re-assemble
filtered_quant_df <- bind_cols(quant_df[, 1:3], as.data.frame(quant_mat))

# --- Step 3: Group-Aware Missing Value Filtering ---
metadata <- data.frame(
  SampleID = c("Sample_6017", "Sample_6091", "Sample_6092", "Sample_6093", "Sample_6095", 
               "Sample_6096", "Sample_6097", "Sample_6098", "Sample_6099", "Sample_6100", 
               "Sample_6101", "Sample_6102", "Sample_6103", "Sample_6104", "Sample_6105"),
  Group = c("BC2_noMG_NP", "BC2_withMG_CTRL", "BC2_noMG_CTRL", "BC2_withMG_AD", "BC2_withMG_AD_NP",
            "BC2_withMG_AD_NP", "AWT_withMG_AD_NP", "AWT_noMG_AD_NP", "BC2_withMG_NAD_NP", "AWT_withMG_NAD_NP",
            "BC2_noMG_NAD_NP", "BC2_noMG_NAD_NP", "Tissue_AD", "Tissue_NAD", "BC2_withMG_NP"),
  stringsAsFactors = FALSE
)

numeric_mat <- as.matrix(filtered_quant_df[, metadata$SampleID])

# Global valid counts filter
keep_global <- rowSums(!is.na(numeric_mat)) >= 3

# Group-aware 70% rule (Handles single-replicate groups cleanly)
keep_group <- rep(FALSE, nrow(numeric_mat))
for (g in unique(metadata$Group)) {
  group_samples <- metadata$SampleID[metadata$Group == g]
  sub_mat <- numeric_mat[, group_samples, drop = FALSE]
  
  # If group has 1 sample, it must be valid. If >1 sample, must meet 70% threshold.
  if (ncol(sub_mat) == 1) {
    valid_group <- !is.na(sub_mat[, 1])
  } else {
    valid_group <- rowMeans(!is.na(sub_mat)) >= 0.70
  }
  keep_group <- keep_group | valid_group
}

# Apply filters
final_keep_indices <- which(keep_group & keep_global)
filtered_quant_df <- filtered_quant_df[final_keep_indices, ]
cat("Proteins remaining after stricter filtering:", nrow(filtered_quant_df), "\n")


# --- Step 4: Log2 Transformation & Quantile Normalization ---
log_mat <- log2(as.matrix(filtered_quant_df[, metadata$SampleID]))

# Log2 of 0 yields -Inf. Convert -Inf to NA, not NaN.
log_mat[is.infinite(log_mat)] <- NA 

# Quantile normalization
norm_mat <- preprocessCore::normalize.quantiles(log_mat)
dimnames(norm_mat) <- dimnames(log_mat)


# --- Step 5: Imputation (k-NN) ---
# impute.knn breaks if a row is completely missing. 
# We filter out rows that are entirely NA after normalization just in case.
keep_non_empty_rows <- rowSums(!is.na(norm_mat)) > 0
norm_mat_clean <- norm_mat[keep_non_empty_rows, , drop = FALSE]
filtered_quant_df_clean <- filtered_quant_df[keep_non_empty_rows, ]

imputed_output <- impute.knn(norm_mat_clean, k = 10, rowmax = 0.8)
final_numeric_mat <- imputed_output$data

final_clean_dataset <- bind_cols(filtered_quant_df_clean[, 1:3], as.data.frame(final_numeric_mat))

# Save output
write.table(final_clean_dataset, "~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/fully_processed_proteomics_matrix_v2.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


# --- Step 6: Re-Plot the Distributions ---
plot_df <- final_clean_dataset %>%
  pivot_longer(cols = all_of(metadata$SampleID), names_to = "SampleID", values_to = "Log2_Intensity") %>%
  left_join(metadata, by = "SampleID")

ggplot(plot_df, aes(x = Log2_Intensity, group = SampleID, color = Group)) +
  geom_density(linewidth = 0.7, alpha = 0.6) +
  theme_minimal(base_family = "serif") + 
  labs(
    title = "Log2 Intensity Distribution (Optimized Missing Value Filter)",
    x = "Normalized Log2 Intensity",
    y = "Density",
    color = "Experimental Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9)
  )
