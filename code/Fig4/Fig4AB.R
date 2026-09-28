## PCA and Volcano plot figure 4A-B
# Date: 15 Oct 2024
# Programmer: Cyriel Huijer
# The reference proteome that was used is UP000001075

# libraries
library(tidyverse)
library(ggrepel)

# clear environment
rm(list=ls())

# set working environment to /code
setwd("")

# load data
complete_table <- read.csv("data/proteomics/output/fmol_ug_table.csv")

# Figure 4A
data_pca <- complete_table[, c(13:18)] # This line removes the outliers

# Calculate row variances
row_variances <- apply(data_pca, 1, var)

# select the top 500 most variable rows
top_500_indices <- order(row_variances, decreasing = TRUE)[1:500]

# subset the data to keep only the top 500 rows
data_pca_top500 <- data_pca[top_500_indices, ]

# transpose the data
data_pca_transposed <- t(data_pca_top500)

# Perform PCA
pca_result <- prcomp(data_pca_transposed, center = TRUE, scale. = TRUE)

explained_variance <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
explained_variance_df <- data.frame(
  PC = paste0("PC", 1:length(explained_variance)),
  Variance = explained_variance
)
explained_variance_PC1 <- round(explained_variance[1], 0)
explained_variance_PC2 <- round(explained_variance[2], 0)

# Create a dataframe with the PCA results
pca_data <- as.data.frame(pca_result$x)

# Add column names as a new variable to the PCA dataframe
original_colnames <- colnames(data_pca)
extracted_labels <- gsub("fmol_per_ug_", "", original_colnames)


pca_data$Label <- extracted_labels
pca_data$condition <- c("2E2","WT","2E2","WT","WT","2E2")

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,color = condition)) +
  geom_point(size = 8) +  # Plot points for each column
  geom_text(aes(label = condition), vjust = -0.5, hjust = 0.5, size = 8) +  # Add labels (column names)
  scale_color_manual(
    values = c(
      "2E2" = "#0070C0",
      "WT"  = "#FF7043"
    )
  )+
  theme_bw(base_size = 15) +
  theme(legend.position = "none") +
  labs(
    x = paste0("PC1 (", explained_variance_PC1, "% Variance)"),
    y = paste0("PC2 (", explained_variance_PC2, "% Variance)")
  )

# Figure 4A
pca_plot

# Figure 4B

for (i in 1:nrow(complete_table)) {
  # Extract WT and KO replicates
  WT <- as.numeric(complete_table[i, grep("fmol_per_ug_WT_", names(complete_table))])
  KO <- as.numeric(complete_table[i, grep("fmol_per_ug_X2E2_", names(complete_table))])
  # Perform t-test
  t_test <- t.test(KO, WT)
  
  # Calculate log2 fold change
  log2FC <- mean(log2(KO)) - mean(log2(WT))
  
  if (i == 1420){
    t_test_result <- t_test
    print(WT)
    print(KO)
  }
  
  # Store results
  complete_table$mean_WT[i] <- mean(WT)
  complete_table$mean_2E2[i] <- mean(KO)
  complete_table$log2FC[i] <- log2FC
  complete_table$p_value[i] <- t_test$p.value
  
}

complete_table$adjusted_p_value <- p.adjust(complete_table$p_value, method = "BH")
significant <- subset(complete_table, adjusted_p_value < 0.05 & abs(log2FC) > 0.5)

significant$short_description <- gsub(" OS=.*", "", significant$Description)

# Figure 4B
p <- ggplot(complete_table, aes(x = log2FC, y = -log10(adjusted_p_value))) +
  geom_point(alpha = 0.3, color = "gray") +  # Plot all points in gray
  geom_point(data = significant, aes(x = log2FC, y = -log10(adjusted_p_value)), color = "red", alpha = 0.6) +  # Highlight significant points
  theme_bw(base_size = 15) + 
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted p-value)") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +   # Significance threshold line
  geom_vline(xintercept = c(-0.5, 0.5), color = "blue", linetype = "dashed") +  # Fold change threshold lines
  geom_text_repel(data = significant, aes(label = short_description), size = 3, box.padding = 1, max.overlaps = Inf)+
  xlim(-2.5,2.5)# Repel labels to avoid overlap
p



