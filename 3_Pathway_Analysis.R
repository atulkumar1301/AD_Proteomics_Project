if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"))

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# Create an output directory for pathway tables and figures
dir.create("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Pathway_Analysis_Results", showWarnings = FALSE)

# 1. Define the target comparisons to analyze
target_contrasts <- c(
  "AD_vs_CTRL", 
  "AD_NP_vs_CTRL", 
  "Microglia_Effect_AD_NP", 
  "Microglia_Effect_CTRL",
  "Microglia_Interaction",
  "CellLine_Diff_AD", 
  "CellLine_Diff_NAD", 
  "CellLine_Interaction", 
  "Tissue_AD_vs_NAD"
)

# 2. Keywords to automatically flag and prioritize your requested pathways
target_keywords <- "inflam|immun|neuro|synap|axon|mitochon|metabol|oxida|respirat|matrix|collagen|myelin|cell adhesion"

# Loop through each target comparison
for (contrast in target_contrasts) {
  cat("\n========================================\n")
  cat("Processing Pathway Enrichment for:", contrast, "\n")
  cat("========================================\n")
  
  # Load the full results file generated in the previous step
  file_path <- paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Proteomics_Results/", contrast, "_Full_Results.txt")
  if (!file.exists(file_path)) {
    cat("Warning: Target file not found:", file_path, "\nSkipping...\n")
    next
  }
  
  res_table <- read.delim(file_path, sep = "\t", , check.names = FALSE)
  
  # Filter for raw P-Value < 0.05 and an absolute Log2FC > 0.5
  sig_data <- res_table %>%
    filter(P.Value < 0.05 & abs(logFC) > 0.5) %>%
    filter(PG.Genes != "" & !is.na(PG.Genes))
  
  # Extract and clean gene symbols (split cases like "GENE1;GENE2" into individual elements)
  sig_genes <- unlist(strsplit(as.character(sig_data$PG.Genes), ";")) %>%
    trimws() %>%
    unique()
  
  cat("Number of unique significant genes found:", length(sig_genes), "\n")
  if (length(sig_genes) < 5) {
    cat("Too few genes for robust enrichment. Skipping to next contrast.\n")
    next
  }
  
  # ------------------------------------------------------------
  # A. GENE ONTOLOGY (GO) BIOLOGICAL PROCESS ENRICHMENT
  # ------------------------------------------------------------
  cat("Running GO Biological Process enrichment...\n")
  go_res <- enrichGO(
    gene          = sig_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP", # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20
  )
  
  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
    go_df <- as.data.frame(go_res)
    
    # Label your specifically requested pathway domains for easy report grouping
    go_df <- go_df %>%
      mutate(Requested_Domain = ifelse(grepl(target_keywords, Description, ignore.case = TRUE), "Target Pathway", "Other"))
    
    write.table(go_df, paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Pathway_Analysis_Results/", contrast, "_GO_BP_Enrichment.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save a publication barplot featuring serif typography
    p_go <- barplot(go_res, showCategory = 15, title = paste("GO Biological Process:", gsub("_", " ", contrast))) +
      theme_minimal(base_family = "serif") +
      theme(plot.title = element_text(face = "bold", size = 12))
    
    ggsave(paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Pathway_Analysis_Results/", contrast, "_GO_Barplot.pdf"), plot = p_go, width = 11.7, height = 8.3, device = cairo_pdf)
  } else {
    cat("No significant GO terms found.\n")
  }
  
  # ------------------------------------------------------------
  # B. REACTOME PATHWAY ENRICHMENT 
  # ------------------------------------------------------------
  cat("Running Reactome Pathway enrichment...\n")
  gene_mapping <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_mapping) > 0) {
    react_res <- enrichPathway(
      gene         = gene_mapping$ENTREZID,
      organism     = "human",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20
    )
    
    if (!is.null(react_res) && nrow(as.data.frame(react_res)) > 0) {
      react_df <- as.data.frame(react_res) %>%
        mutate(Requested_Domain = ifelse(grepl(target_keywords, Description, ignore.case = TRUE), "Target Pathway", "Other"))
      
      write.table(react_df, paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Pathway_Analysis_Results/", contrast, "_Reactome_Enrichment.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Save a publication dotplot featuring serif typography
      p_react <- dotplot(react_res, showCategory = 15, title = paste("Reactome Pathways:", gsub("_", " ", contrast))) +
        theme_minimal(base_family = "serif") +
        theme(plot.title = element_text(face = "bold", size = 12))
      
      ggsave(paste0("~/OneDrive - University of Eastern Finland/Projects/Aakash_Mali/Results/Pathway_Analysis_Results/", contrast, "_Reactome_Dotplot.pdf"), plot = p_react, width = 11.7, height = 8.3, device = cairo_pdf)
    } else {
      cat("No significant KEGG terms found.\n")
    }
  } else {
    cat("Failed to map gene symbols to Entrez IDs for KEGG.\n")
  }
}

cat("\n========================================\n")
cat("All functional pathway analyses complete!\n")
cat("Check the 'Pathway_Analysis_Results' folder.\n")
