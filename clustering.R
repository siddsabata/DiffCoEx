# Jason Hyun (jasonhyu)
# Siddharth Sabata (ssabata)
# Darrick Lo (ddlo)
# Katie Wang (kcw2)

# Dec 1, 2024

# NOTE: Generative AI used to produce following code:

library(coXpress)
library(WGCNA)          ###used for topological overlap calculation and clustering steps
library(RColorBrewer)   ###used to create nicer colour palettes
library(preprocessCore) ###used by the quantile normalization function
library(flashClust)
library(ggplot2)
library(tidyverse)
theme_set(theme_classic())

performClustering <- function(data_paths) {
  # Read the preprocessed data from CSV files for DiffCoEx
  condition1_diffcoex <- read.csv(data_paths$diffcoex$condition1, row.names = 1)
  condition2_diffcoex <- read.csv(data_paths$diffcoex$condition2, row.names = 1)
  
  # Read the preprocessed data from CSV files for coXpress
  condition1_coxpress <- read.csv(data_paths$coxpress$condition1, row.names = 1)
  condition2_coxpress <- read.csv(data_paths$coxpress$condition2, row.names = 1)
  
  # Print initial dimensions
  cat("DiffCoEx dimensions:\n")
  cat("condition1:", dim(condition1_diffcoex), "\n")
  cat("condition2:", dim(condition2_diffcoex), "\n")
  
  cat("\ncoXpress dimensions:\n")
  cat("condition1:", dim(condition1_coxpress), "\n")
  cat("condition2:", dim(condition2_coxpress), "\n")
  
  # Convert to matrices while preserving structure
  condition1_diffcoex_matrix <- as.matrix(condition1_diffcoex)
  condition2_diffcoex_matrix <- as.matrix(condition2_diffcoex)
  condition1_coxpress_matrix <- as.matrix(condition1_coxpress)
  condition2_coxpress_matrix <- as.matrix(condition2_coxpress)
  
  # Ensure all values are numeric
  mode(condition1_diffcoex_matrix) <- "numeric"
  mode(condition2_diffcoex_matrix) <- "numeric"
  mode(condition1_coxpress_matrix) <- "numeric"
  mode(condition2_coxpress_matrix) <- "numeric"
  
  # Replace NAs with row means for each matrix
  for(i in 1:nrow(condition1_diffcoex_matrix)) {
    row_mean <- mean(condition1_diffcoex_matrix[i,], na.rm = TRUE)
    condition1_diffcoex_matrix[i, is.na(condition1_diffcoex_matrix[i,])] <- row_mean
  }
  
  for(i in 1:nrow(condition2_diffcoex_matrix)) {
    row_mean <- mean(condition2_diffcoex_matrix[i,], na.rm = TRUE)
    condition2_diffcoex_matrix[i, is.na(condition2_diffcoex_matrix[i,])] <- row_mean
  }
  
  for(i in 1:nrow(condition1_coxpress_matrix)) {
    row_mean <- mean(condition1_coxpress_matrix[i,], na.rm = TRUE)
    condition1_coxpress_matrix[i, is.na(condition1_coxpress_matrix[i,])] <- row_mean
  }
  
  for(i in 1:nrow(condition2_coxpress_matrix)) {
    row_mean <- mean(condition2_coxpress_matrix[i,], na.rm = TRUE)
    condition2_coxpress_matrix[i, is.na(condition2_coxpress_matrix[i,])] <- row_mean
  }
  
  # Remove rows where all values are the same
  valid_rows_diffcoex <- apply(condition1_diffcoex_matrix, 1, function(x) length(unique(x)) > 1) &
                        apply(condition2_diffcoex_matrix, 1, function(x) length(unique(x)) > 1)
  
  valid_rows_coxpress <- apply(condition1_coxpress_matrix, 1, function(x) length(unique(x)) > 1) &
                        apply(condition2_coxpress_matrix, 1, function(x) length(unique(x)) > 1)
  
  condition1_diffcoex_matrix <- condition1_diffcoex_matrix[valid_rows_diffcoex, ]
  condition2_diffcoex_matrix <- condition2_diffcoex_matrix[valid_rows_diffcoex, ]
  condition1_coxpress_matrix <- condition1_coxpress_matrix[valid_rows_coxpress, ]
  condition2_coxpress_matrix <- condition2_coxpress_matrix[valid_rows_coxpress, ]
  
  # DiffCoEx analysis
  datC1 <- t(condition1_diffcoex_matrix) # First condition
  datC2 <- t(condition2_diffcoex_matrix) # Second condition
  
  beta1=6 #user defined parameter for soft thresholding
  AdjMatC1<-sign(cor(datC1,method="spearman"))*(cor(datC1,method="spearman"))^2
  AdjMatC2<-sign(cor(datC2,method="spearman"))*(cor(datC2,method="spearman"))^2
  diag(AdjMatC1)<-0
  diag(AdjMatC2)<-0
  collectGarbage()
  
  dissTOMC1C2=TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))
  collectGarbage()
  
  geneTreeC1C2 = flashClust(as.dist(dissTOMC1C2), method = "average")
  
  dynamicModsHybridC1C2 = cutreeDynamic(dendro = geneTreeC1C2, distM = dissTOMC1C2,
                                       method="hybrid",cutHeight=.996,deepSplit = T,
                                       pamRespectsDendro = FALSE,minClusterSize = 20)
  
  dynamicColorsHybridC1C2 = labels2colors(dynamicModsHybridC1C2)
  
  mergedColorC1C2<-mergeCloseModules(rbind(datC1,datC2),dynamicColorsHybridC1C2,cutHeight=.2)$color
  colorh1C1C2<-mergedColorC1C2
  
  # coXpress analysis
  hc.gene <- cluster.gene(condition1_coxpress_matrix, s="pearson", m="average")
  g <- cutree(hc.gene, h=0.4)
  
  # Create results dataframes
  diffcoex_results <- data.frame(
    Gene = rownames(condition1_diffcoex_matrix),
    Cluster = colorh1C1C2
  )
  
  coxpress_results <- data.frame(
    Gene = rownames(condition1_coxpress_matrix),
    Cluster = g
  )
  
  # Count genes per cluster
  genes_per_cluster_diffcoex <- as.numeric(table(diffcoex_results$Cluster))
  genes_per_cluster_coxpress <- as.numeric(table(coxpress_results$Cluster))
  
  # Compute summary statistics
  diffcoex_summary <- list(
    Total_Genes = nrow(diffcoex_results),
    Total_Clusters = length(unique(diffcoex_results$Cluster)),
    Genes_Per_Cluster = summary(genes_per_cluster_diffcoex)
  )
  
  coxpress_summary <- list(
    Total_Genes = nrow(coxpress_results),
    Total_Clusters = length(unique(coxpress_results$Cluster)),
    Genes_Per_Cluster = summary(genes_per_cluster_coxpress)
  )
  
  # Create visualizations
  # Number of clusters comparison
  num_clusters <- data.frame(
    Method = c("DiffCoEx", "CoXpress"),
    Clusters = c(length(unique(diffcoex_results$Cluster)), 
                length(unique(coxpress_results$Cluster)))
  )
  
  num_clusters_plot <- ggplot(num_clusters, aes(x = Method, y = Clusters, fill = Method)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "Number of Clusters by Method", x = "Method", y = "Number of Clusters") +
    theme_minimal() +
    theme(text = element_text(size = 14))
  
  # Genes per cluster distributions
  genes_per_cluster_df_diffcoex <- as.data.frame(table(diffcoex_results$Cluster))
  colnames(genes_per_cluster_df_diffcoex) <- c("Cluster", "Gene_Count")
  
  genes_per_cluster_df_coxpress <- as.data.frame(table(coxpress_results$Cluster))
  colnames(genes_per_cluster_df_coxpress) <- c("Cluster", "Gene_Count")
  
  # Create cluster overlap heatmap data
  merged_data <- merge(diffcoex_results, coxpress_results, by = "Gene", suffixes = c("_DiffCoEx", "_CoXpress"))
  shared_genes_table <- table(merged_data$Cluster_DiffCoEx, merged_data$Cluster_CoXpress)
  
  # Convert to data frame and filter for significant overlaps
  shared_genes_df <- as.data.frame(as.table(shared_genes_table))
  colnames(shared_genes_df) <- c("DiffCoEx_Module", "CoXpress_Module", "Gene_Count")
  
  # Filter for overlaps > 5 genes
  shared_genes_df <- shared_genes_df[shared_genes_df$Gene_Count > 5, ]
  
  # Get unique modules that have significant overlaps
  significant_diffcoex_modules <- unique(shared_genes_df$DiffCoEx_Module)
  significant_coxpress_modules <- unique(shared_genes_df$CoXpress_Module)
  
  # Filter the data to only include these modules
  shared_genes_df <- shared_genes_df[
    shared_genes_df$DiffCoEx_Module %in% significant_diffcoex_modules & 
    shared_genes_df$CoXpress_Module %in% significant_coxpress_modules, 
  ]
  
  # Create the filtered overlap plot
  cluster_overlap_plot <- ggplot(shared_genes_df, 
    aes(x = CoXpress_Module, y = DiffCoEx_Module, fill = Gene_Count)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = "Cluster Overlap Between Methods\n(>5 genes)", 
         x = "CoXpress Clusters", 
         y = "DiffCoEx Clusters") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  # Save DiffCoEx module assignments to CSV
  diffcoex_module_map <- data.frame(
    Gene = rownames(condition1_diffcoex),
    Module = colorh1C1C2
  )
  
  # Create directories if they don't exist
  dir.create("data/golub", recursive = TRUE, showWarnings = FALSE)
  dir.create("data/rat", recursive = TRUE, showWarnings = FALSE)
  
  # Save module map based on dataset type
  if (grepl("golub", data_paths$diffcoex$condition1)) {
    write.csv(diffcoex_module_map, "data/golub/golub_diffcoex.csv", row.names = FALSE)
  } else {
    write.csv(diffcoex_module_map, "data/rat/rat_diffcoex.csv", row.names = FALSE)
  }
  
  # Return results
  return(list(
    diffcoex_results = diffcoex_results,
    coxpress_results = coxpress_results,
    diffcoex_summary = diffcoex_summary,
    coxpress_summary = coxpress_summary,
    plots = list(
      cluster_comparison = num_clusters_plot,
      diffcoex_distribution = genes_per_cluster_df_diffcoex,
      coxpress_distribution = genes_per_cluster_df_coxpress,
      cluster_overlap = cluster_overlap_plot
    )
  ))
}