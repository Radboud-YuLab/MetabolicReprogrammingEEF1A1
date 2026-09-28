## CHO 2E2 proteomics analysis - DEP package
# Date: 15 Oct 2024
# Programmer: Cyriel Huijer
# The reference proteome that was used is UP000001075

# packages
library(DEP)
library(dplyr)
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(UniProt.ws)
library(OmnipathR)
library(ape)
library(ggrepel)
library(reshape2)
library(VennDiagram)

# clear environment
rm(list=ls())

# set wd to /code
setwd("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/WesternOntario_collab/9_manuscript/code/")

# Load data
data <- read_excel("data/proteomics/202407231_Exploris1_SA_CH_CHO_tmt_50ug_50percent.xlsx")
df_master <- subset(data, Master == "Master Protein" & Contaminant == "FALSE")
data <- df_master
data$protein_name <- sub(" OS=.*", "", data$Description)
data_unique <- make_unique(data, "Accession", "protein_name")
#data_unique <- make_unique(data, "Accession","protein_name")
# Metadata

sample_info <- data.frame(
  label = c("TMT.WT_1", "TMT.2E2_1", "TMT.WT_2", "TMT.2E2_2", "TMT.WT_3", "TMT.2E2_3", "TMT.WT_4", "TMT.2E2_4","pooled"),  # Add all sample names
  condition = c("WT", "2E2", "WT", "2E2", "WT", "2E2", "WT", "2E2","pooled"),              # Define conditions
  replicate = c(1,1,2,2,3,3,4,4,5) 
)
# Create a SummarizedExperiment object

colnames(data_unique)[36:44] <- sample_info$label

se <- make_se(data_unique, columns = c(36:44), expdesign = sample_info)
TMT_columns <- grep("TMT.", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- data_se_parsed <- make_se_parse(data_unique, TMT_columns)

# Plot a barplot of the protein identification overlap between samples
plot_frequency(se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(se, thr = 0)
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)


# Normalize the data
data_norm <- normalize_vsn(data_filt)
data_for_iBAQ <- data_norm
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(se, data_norm)

# Test every sample versus control
data_diff <- test_diff(data_norm, type = "control", control = "pooled")
# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_norm, type = "all")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 0.5)

# Generate a results table
data_results <- get_results(dep)

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)



# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
volcano_plot<- plot_volcano(dep, contrast = "WT_vs_X2E2", 
             label_size = 3, 
             add_names = TRUE)
volcano_plot

## Rerun analysis with some outliers deleted: #############################################################################################

data <- read_excel("data/proteomics/202407231_Exploris1_SA_CH_CHO_tmt_50ug_50percent.xlsx")
df_master <- subset(data, Master == "Master Protein")
data <- df_master
data$protein_name <- sub(" OS=.*", "", data$Description)
data_unique <- make_unique(data, "protein_name", "Accession")
# Metadata

sample_info <- data.frame(
  label = c("TMT.2E2_1", "TMT.WT_2", "TMT.2E2_2", "TMT.WT_3", "TMT.2E2_3", "TMT.WT_4"),  # Add all  sample names
  condition = c("2E2", "WT", "2E2", "WT", "2E2", "WT"),              # Define conditions
  replicate = c(1,2,2,3,3,4) 
)

data <- read_excel("6_results/proteomics/tmt/202407231_Exploris1_SA_CH_CHO_tmt_50ug_50percent.xlsx")
df_master <- subset(data, Master == "Master Protein")
data <- df_master
data$protein_name <- sub(" OS=.*", "", data$Description)
data_unique <- make_unique(data, "protein_name", "Accession")

sample_info <- data.frame(
  label = c("TMT.2E2_1", "TMT.WT_2", "TMT.2E2_2", "TMT.WT_3", "TMT.2E2_3", "TMT.WT_4"),  # Add all sample names
  condition = c("2E2", "WT", "2E2", "WT", "2E2", "WT"),              # Define conditions
  replicate = c(1,2,2,3,3,4) 
)

colnames(data_unique)[37:42] <- sample_info$label

se <- make_se(data_unique, columns = c(37:42), expdesign = sample_info)
TMT_columns <- grep("TMT.", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- data_se_parsed <- make_se_parse(data_unique, TMT_columns)

# Plot a barplot of the protein identification overlap between samples
plot_frequency(se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(se, thr = 0)
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(se, data_norm)

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_norm, type = "control", control = "pooled")
# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_norm, type = "all")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 0.5)

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
volcano_plot<- plot_volcano(dep, contrast = "WT_vs_X2E2", 
                            label_size = 3, 
                            add_names = TRUE,
                            adjusted = TRUE)
volcano_plot

## Integration of UPS2 data (absolute quantities derived from iBAQ) ##########################

absolute_quant <- read.csv("data/proteomics/output/UPS2_quantification_results.csv", header = TRUE)
data_corr <- data.frame(data_for_iBAQ@assays@data@listData[[1]]) # Take this column from the data normalization in the first DEP analysis
data_corr$Accession <- rownames(data_corr)

# Merge data with UPS2 data based on Accession ID

merged_data <- merge(data_corr, absolute_quant, by = "Accession", all.x = TRUE) # Merge the dataframes based on Accession numbers
filtered_data_no_na <- na.omit(merged_data) # Remove NAs


for(i in 2:9) {
  # Create column names for 
  column_name <- paste("fmol_per_ug", colnames(filtered_data_no_na)[i], sep = "_")
  
  # Calculate femtomoles per condition by multiplying the log2 ratio by the femtomoles of the pooled condition 
  #filtered_data_no_na[[column_name]] <- (2^(filtered_data_no_na[[i]] - filtered_data_no_na[[10]]) * filtered_data_no_na$femtomoles)/15
  filtered_data_no_na[[column_name]] <- ((2^(filtered_data_no_na[[i]]) / 2^(filtered_data_no_na[[10]])) * filtered_data_no_na$femtomoles)/15 # Divide by 15 to calculate per µg protein
}

# This is the final table, now rerun the analysis using the DEP package:

complete_table <- merge(filtered_data_no_na, data, by = "Accession", all.x = TRUE)

complete_table$protein_name <- sub(" OS=.*", "", complete_table$Description)
complete_table <- complete_table[ , -c(12, 17)] # Remove outliers
#write.csv(complete_table, "data/proteomics/output/fmol_ug_table.csv")

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
tableS2 <- complete_table %>% dplyr::select(Accession,protein_name,mean_WT,mean_2E2,log2FC,p_value,adjusted_p_value)
#write.csv(tableS2, "9_manuscript/supplementary_tables/TableS2.csv",row.names = F)
