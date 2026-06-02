library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel) # Smooth label layout without overlaps

#Step 1: Define the Design Matrix in Limma
# final_clean_dataset comes from the optimized filter step
final_clean_dataset <- read.delim("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Fully_processed_proteomics_matrix.txt", sep = "\t", check.names = FALSE)
metadata <- data.frame(
  SampleID = c("Sample_6017", "Sample_6091", "Sample_6092", "Sample_6093", "Sample_6095", 
               "Sample_6096", "Sample_6097", "Sample_6098", "Sample_6099", "Sample_6100", 
               "Sample_6101", "Sample_6102", "Sample_6103", "Sample_6104", "Sample_6105"),
  Group = c("BC2_noMG_NP", "BC2_withMG_CTRL", "BC2_noMG_CTRL", "BC2_withMG_AD", "BC2_withMG_AD_NP",
            "BC2_withMG_AD_NP", "AWT_withMG_AD_NP", "AWT_noMG_AD_NP", "BC2_withMG_NAD_NP", "AWT_withMG_NAD_NP",
            "BC2_noMG_NAD_NP", "BC2_noMG_NAD_NP", "Tissue_AD", "Tissue_NAD", "BC2_withMG_NP"),
  stringsAsFactors = FALSE
)

expression_matrix <- as.matrix(final_clean_dataset[, metadata$SampleID])
rownames(expression_matrix) <- final_clean_dataset$PG.ProteinGroups

# 1. Create a factor of the experimental groups
group_factor <- factor(metadata$Group)

# 2. Build the cell-means design matrix (omitting the intercept)
design <- model.matrix(~ 0 + group_factor)
# Clean up column names so they match your contrast strings exactly ===
colnames(design) <- levels(group_factor) 
# This removes the automatic "group_factor" prefix added by R

# 3. Fit the initial linear model
fit <- lmFit(expression_matrix, design)

#Step 2: Define Custom Contrasts for Your Research Questions

contrast_matrix <- makeContrasts(
  # --- RQ 1: AD Brain Extract Exposure Effects ---
  AD_vs_CTRL = BC2_withMG_AD - BC2_withMG_CTRL,
  AD_NP_vs_CTRL = BC2_withMG_AD_NP - BC2_withMG_CTRL,
  
  # --- RQ 2: Microglia Status & Interactions ---
  Microglia_Effect_AD = BC2_withMG_AD_NP - BC2_noMG_NP,
  Microglia_Interaction = (BC2_withMG_AD_NP - BC2_noMG_NP) - (BC2_withMG_CTRL - BC2_noMG_CTRL),
  
  # --- RQ 3: Cell Line Differences (BC2 vs A-WT) ---
  CellLine_Diff_AD = BC2_withMG_AD_NP - AWT_withMG_AD_NP,
  CellLine_Diff_NAD = BC2_withMG_NAD_NP - AWT_withMG_NAD_NP,
  CellLine_Interaction = (BC2_withMG_AD_NP - AWT_withMG_AD_NP) - (BC2_withMG_NAD_NP - AWT_withMG_NAD_NP),
  
  # --- RQ 5: Human Tissue Baseline Pathology ---
  Tissue_AD_vs_NAD = Tissue_AD - Tissue_NAD,
  
  levels = design
)

# Compute contrasts and apply Empirical Bayes moderation
fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_bayes <- eBayes(fit_contrasts)

#Step 3: Extract and Export Results for Each Comparison

# Ensure directory exists for organizing outputs
dir.create("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Proteomics_Results", showWarnings = FALSE)

# Get the names of all defined contrasts
contrast_names <- colnames(contrast_matrix)

# Loop through each contrast and export the results
for (contrast in contrast_names) {
  
  # Extract all data rows without a p-value filter to preserve the background proteome
  res_table <- topTable(fit_bayes, coef = contrast, number = Inf, adjust.method = "BH")
  
  # Merge metadata back to the result rows
  res_table_clean <- res_table %>%
    rownames_to_column("PG.ProteinGroups") %>%
    left_join(final_clean_dataset[, c("PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions")], by = "PG.ProteinGroups") %>%
    select(PG.ProteinGroups, PG.Genes, PG.ProteinDescriptions, logFC, AveExpr, t, P.Value, adj.P.Val, B)
  
  # Save file
  filename <- paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Proteomics_Results/", contrast, "_Full_Results.txt")
  write.table(res_table_clean, filename, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Log progress
  cat("Exported successfully:", filename, "\n")
}

## Volcano Plot

# Create a directory to house the plots
dir.create("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Proteomics_Plots", showWarnings = FALSE)

for (contrast in contrast_names) {
  
  # Read back the saved table data
  plot_data <- read.delim(paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Proteomics_Results/", contrast, "_Full_Results.txt"), sep = "\t", check.names = FALSE)
  
  # Calculate -log10 Raw P-Value explicitly for data limits and plotting
  plot_data$log10_P <- -log10(plot_data$P.Value)
  
  # Dynamic Axis Range and Break Generation (Interval of 1)
  x_min <- floor(min(plot_data$logFC, na.rm = TRUE))
  x_max <- ceiling(max(plot_data$logFC, na.rm = TRUE))
  x_breaks <- seq(x_min, x_max, by = 1)
  
  y_min <- 0  # P-values log transformed start at 0
  y_max <- ceiling(max(plot_data$log10_P, na.rm = TRUE))
  y_breaks <- seq(y_min, y_max, by = 1)
  
  # Fail-safe check for adjusted p-value threshold line
  sig_genes <- plot_data %>% filter(adj.P.Val < 0.05)
  if(nrow(sig_genes) > 0) {
    raw_cutoff_for_adjusted <- max(sig_genes$P.Value, na.rm = TRUE)
    fdr_line <- geom_hline(yintercept = -log10(raw_cutoff_for_adjusted), linetype = "dashed", color = "blue", alpha = 0.7)
    subtitle_text <- "Solid Line: Raw P < 0.05 | Dashed Line: FDR < 0.05"
  } else {
    fdr_line <- geom_blank() # Does not crash the plot loop if 0 hits found
    subtitle_text = "Solid Line: Raw P < 0.05 | No genes passed FDR < 0.05"
  }
  
  # 1. Establish custom coloring categories (4 groups based on adj.P.Val and logFC)
  # Note: Adjust the 'logFC' thresholds if you want non-significant points to strictly respect the 0.5 boundary as well.
  plot_data <- plot_data %>%
    mutate(Significance = case_when(
      adj.P.Val < 0.05 & logFC >= 0    ~ "Significant upregulated",
      adj.P.Val < 0.05 & logFC < 0     ~ "Significant downregulated",
      adj.P.Val >= 0.05 & logFC >= 0   ~ "Non-significant upregulated",
      adj.P.Val >= 0.05 & logFC < 0    ~ "Non-significant downregulated"
    ))
  
  # 2. Select the top 5 most significant from each directional category
  top_up <- plot_data %>%
    filter(Significance == "Significant upregulated") %>%
    arrange(adj.P.Val) %>%
    head(5)
  
  top_down <- plot_data %>%
    filter(Significance == "Significant downregulated") %>%
    arrange(adj.P.Val) %>%
    head(5)
  
  # Combine them into a single labeling dataset
  top_labels <- bind_rows(top_up, top_down)
  
  # 3. Render the Volcano Plot
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value))) +
    # Background points mapped to the 4 new categories
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    
    # Text annotations for top genes (Font forced to serif)
    geom_text_repel(data = top_labels, aes(label = PG.Genes), 
                    family = "serif", fontface = "bold", size = 3, 
                    max.overlaps = 15, box.padding = 0.3) +
    
    # Custom color palette for the 4 distinct categories
    scale_color_manual(values = c("Significant upregulated" = "#009E73", 
                                  "Significant downregulated" = "#FC4E2A", 
                                  "Non-significant upregulated" = "#0072B2",
                                  "Non-significant downregulated" = "#CC79A7")) +
    
    # Apply custom axis limits and breaks of 1
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_breaks) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks) +
    
    # Vertical signposts for fold change reference
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", alpha = 0.7) +
    
    # Two horizontal lines: One for raw p-value and one for adjusted p-value
    geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "darkred", alpha = 0.7) + 
    geom_hline(yintercept = -log10(raw_cutoff_for_adjusted), linetype = "dashed", color = "blue", alpha = 0.7) + 
    
    # Formatting text architecture (Enforcing serif globally)
    theme_minimal(base_family = "serif") +
    labs(
      title = paste("Volcano Plot:", gsub("_", " ", contrast)),
      subtitle = "Solid Line: Raw P-Value < 0.05 | Dashed Line: Adjusted P-Value < 0.05",
      x = "Log2 Fold Change",
      y = "-Log10 Raw P-Value",
      color = "Expression Status"
    ) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10)
    )
  
  # 4. Save high-resolution vector graphics in A4 Landscape format (11.7 x 8.3 inches)
  plot_filename <- paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Proteomics_Plots/", contrast, "_Volcano_Plot.pdf")
  ggsave(plot_filename, plot = p, width = 11.7, height = 8.3, device = cairo_pdf)
}

cat("\nAll Volcano Plots successfully generated in the 'Proteomics_Plots' folder!\n")
