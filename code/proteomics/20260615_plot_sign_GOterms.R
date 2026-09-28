## CHO 2E2 proteomics analysis - Plot significant GO terms
# Date: 09 Jan 2024
# Programmer: Cyriel Huijer
# The reference proteome that was used is UP000001075

root <- "C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/WesternOntario_collab/"
setwd(root)
rm(list=ls())


# Load in packages:
library(tidyverse)
library(broom)
library(ggpubr)
library(ggrepel)

# Load in data:
data <- read.csv("5_processedData/Proteomics/go_term_analysis/complete_table_GOs_updated.csv")
data <- data %>%
  filter(WT_mean > 1e-7)
# Normalize the columns prot_fractions according to the sum:
data$norm_protfract_WT_1 <- data$protfract_WT_1 / sum(data$protfract_WT_1)
data$norm_protfract_WT_2 <- data$protfract_WT_2 / sum(data$protfract_WT_2)
data$norm_protfract_WT_3 <- data$protfract_WT_3 / sum(data$protfract_WT_3)
data$norm_protfract_2E2_1 <- data$protfract_2E2_1 / sum(data$protfract_2E2_1)
data$norm_protfract_2E2_2 <- data$protfract_2E2_2 / sum(data$protfract_2E2_2)
data$norm_protfract_2E2_3 <- data$protfract_2E2_3 / sum(data$protfract_2E2_3)
# Extract GO terms:
go_terms <- unique(unlist(strsplit(data$Gene.Ontology..biological.process., ";")))
go_terms <- gsub("^\\s+", "", go_terms)
go_terms <- unique(go_terms)
# Format square brackets in order to be able to use grepl
go_terms <- gsub("\\[", "\\\\\\[", go_terms)  # Escape '['
go_terms <- gsub("\\]", "\\\\\\]", go_terms)  # Escape ']'

# For testing:
# go_terms <- c("translation \\[GO:0006412\\]","glycolytic process \\[GO:0006096\\]")
# gsub("\\\\", "", "translation \\[GO:0006412\\]")

# Define an empty list to store significant terms
significant_terms <- c()
find_term <- c()
mean_WT <- c()
mean_2E2 <- c()
pvalue <- c()
n_proteins_associated <- c()
WT_sum_1 <- c()
WT_sum_2 <- c()
WT_sum_3 <- c()
X2E2_sum_1 <- c()
X2E2_sum_2 <- c()
X2E2_sum_3 <- c()

# Loop over each go_term
for (term in go_terms) {
  # Create a logical vector for rows where the go_term is present
  rows_with_term <- grepl(term, data$Gene.Ontology..biological.process.)
  
  # Check if there are more than 4 rows with the term
  if (sum(rows_with_term) > 4) {
    clean_term <- gsub("\\\\", "", term)
    
    # Filter the dataframe to keep only rows where the term is present
    temp_dataframe <- data[rows_with_term, ]
    
    # Pivot to a longer format
    temp_df_long <- temp_dataframe %>%
      pivot_longer(
        cols = starts_with("norm_protfract_WT_") | starts_with("norm_protfract_2E2_"), # Select columns for replicates
        names_to = "replicate",
        values_to = "prot_fraction"
      )
    
    temp_df_long <- temp_df_long %>%
      mutate(
        condition = ifelse(grepl("WT", replicate), "WT", "2E2")
      )
    # Sum protein fractions per replicate, unit, and condition
    go_sums <- temp_df_long %>%
      group_by(condition, replicate) %>%
      summarise(sum_fraction = sum(prot_fraction, na.rm = TRUE), .groups = "drop")
    
    go_sums$condition <- factor(go_sums$condition, levels = c("WT", "2E2"))
    go_sums_long <- go_sums %>%
      dplyr::select(-condition) %>%  # Drop the 'condition' column
      pivot_wider(
        names_from = replicate,    # Use 'replicate' for new column names
        values_from = sum_fraction # Fill values from 'sum_fraction'
      )
    

    
    # Perform t-test for each unit, comparing conditions
    t_test_results <- go_sums %>%
      do({
        # Perform t-test for each unit between conditions
        t_test <- t.test(sum_fraction ~ condition, data = .)
        # Tidy the result of the t-test
        tidy(t_test) 
      })
    significant_terms <- append(significant_terms, clean_term)
    find_term <- append(find_term,term)
    mean_WT <- append(mean_WT,(t_test_results$estimate1*100))
    mean_2E2 <- append(mean_2E2,(t_test_results$estimate2*100))
    pvalue <- append(pvalue,t_test_results$p.value)
    n_proteins_associated <- append(n_proteins_associated,sum(rows_with_term))
    WT_sum_1 <- append(WT_sum_1,go_sums_long$norm_protfract_WT_1)
    WT_sum_2 <- append(WT_sum_2,go_sums_long$norm_protfract_WT_2)
    WT_sum_3 <- append(WT_sum_3,go_sums_long$norm_protfract_WT_3)
    X2E2_sum_1 <- append(X2E2_sum_1,go_sums_long$norm_protfract_2E2_1)
    X2E2_sum_2 <- append(X2E2_sum_2,go_sums_long$norm_protfract_2E2_2)
    X2E2_sum_3 <- append(X2E2_sum_3,go_sums_long$norm_protfract_2E2_3)
    
    # If p-value is significant
    if (t_test_results$p.value[1] < 0.05) {
      # Append the clean_term to the list
      
      term_for_save <- gsub(" ", "_", clean_term)  # Replace spaces with underscores
      term_for_save <- gsub("/", "_", term_for_save)
      term_for_save <- gsub(":", "_", term_for_save)
      term_for_save <- gsub("\\.", "_", term_for_save)
      term_for_save <- substr(term_for_save, 1, 60)  # Limit to the first 50 characters
      # Plot
      plot <- ggplot(go_sums, aes(x = condition, y = 100 * sum_fraction, color = condition)) +
        geom_bar(stat = "summary", fun = "mean", alpha = 0.7, position = position_dodge(width = 0.8), aes(fill = condition)) + # Bar for mean
        geom_point(width = 0.2, size = 2) +  # Jitter for points
        labs(
          title = paste0("Proteome Fraction (", clean_term, ")"),
          x = "",
          y = "Proteome Fraction (%)"
        ) +
        theme_bw() +
        theme(legend.position = 'none')
      
      # Print the plot
      #print(plot)
      # Save plots
      #ggsave(paste0(root,"6_results/proteomics/proteome_fraction/sign_GO_terms/pf_", term_for_save, ".pdf"), plot, width = 6, height = 4, dpi = 300)
    }
  }
}

results_GO_terms <- data.frame(significant_terms,find_term,mean_WT,mean_2E2,pvalue,
                               n_proteins_associated,WT_sum_1,WT_sum_2,WT_sum_3,X2E2_sum_1,
                               X2E2_sum_2,X2E2_sum_3)
results_GO_terms <- results_GO_terms %>%
  filter(mean_WT > 1e-6)
ggplot(results_GO_terms, aes(x = log10(mean_WT), y = log10(mean_2E2))) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "GO term proteome allocation: WT vs 2E2",
    x = "Mean proteome fraction (WT)",
    y = "Mean proteome fraction (2E2)"
  ) +
  theme_bw(base_size = 14)
results_GO_terms$ratio <- results_GO_terms$mean_WT / results_GO_terms$mean_2E2

results_GO_terms$adjusted_pvalue <- p.adjust(results_GO_terms$pvalue, method = "BH")

TableS3 <- results_GO_terms %>% dplyr::select(significant_terms,mean_WT,mean_2E2,n_proteins_associated,pvalue,adjusted_pvalue)
colnames(TableS3)[4] <- "n_proteins"
#write.csv(TableS3,"9_manuscript/supplementary_tables/TableS3.csv",row.names = F)
# Extract GO terms  from significnat_terms
# results_GO_terms$GO_term <- sub(".*\\[(GO:[0-9]+)\\].*", "\\1", results_GO_terms$significant_terms)
# 
# length(intersect(results_GO_terms$GO_term,children_terms))
# metabolic_results <- results_GO_terms[results_GO_terms$GO_term %in% children_terms, ]
# metabolic_results <- metabolic_results %>%
#   filter(mean_WT > 0.1 | mean_2E2 > 0.1) %>%
#   filter(mean_WT <  4)
# 
# ggplot(metabolic_results, aes(x = mean_WT, y=mean_2E2))+
#   geom_point()+
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")


#write.csv(results_GO_terms, file = paste0(root,"6_results/proteomics/proteome_fraction/results_GO_terms_all.csv"), row.names = FALSE)


sign_go_terms <- results_GO_terms #%>%
  filter(adjusted_pvalue < 0.10)

sign_go_terms_high <-sign_go_terms %>%
  filter(mean_WT > 0.5 | mean_2E2 > 0.5)
# Plot:
# Pivot data to long format
sign_go_terms_long <- sign_go_terms %>%
  pivot_longer(
    cols = c(mean_WT, mean_2E2),
    names_to = "condition",
    values_to = "mean_value"
  )
sign_go_terms_long <- sign_go_terms_long %>%
  mutate(significant_terms = fct_reorder(significant_terms, mean_value))
# Create the bar plot
ggplot(sign_go_terms_long, aes(x = significant_terms, y = mean_value, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), aes(color = condition)) +
  labs(
    title = "Proteome Fraction per GO-term (WT vs 2E2)",
    x = "GO Terms",
    y = "Proteome Fraction (%)")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels 45 degrees
    legend.position = "right"
  )


# Plot terms: Glycolysis, TCA cycle, 

up <- c("response to endoplasmic reticulum stress [GO:0034976]","positive regulation of gene expression [GO:0010628]","regulation of mRNA stability [GO:0043488]","negative regulation of apoptotic process [GO:0043066]","acetyl-CoA biosynthetic process [GO:0006085]")
down <- c("glycolytic process [GO:0006096]","cristae formation [GO:0042407]","ribosomal small subunit export from nucleus [GO:0000056]","ribosomal large subunit export from nucleus [GO:0000055]","chromatin remodeling [GO:0006338]")


# GO terms of interest
terms_to_plot <- c(
  "response to endoplasmic reticulum stress [GO:0034976]",
  "positive regulation of gene expression [GO:0010628]",
  #"regulation of mRNA stability [GO:0043488]",
  "negative regulation of apoptotic process [GO:0043066]",
  "chromatin remodeling [GO:0006338]",
  "ribosomal small subunit biogenesis [GO:0042274]",
  "ribosomal large subunit export from nucleus [GO:0000055]",
  #"cristae formation [GO:0042407]",
  "glycolytic process [GO:0006096]",
  "tricarboxylic acid cycle [GO:0006099]",
  #"fatty acid beta-oxidation [GO:0006635]",
  "acetyl-CoA biosynthetic process [GO:0006085]"
  
)

# sign_go_terms <- sign_go_terms %>%
#   filter(ratio > 1.1 | ratio < 0.9)

label_terms <- sign_go_terms_high |>
  dplyr::slice_min(ratio, n = 10) |>      # 10 lowest
  dplyr::bind_rows(
    sign_go_terms |> dplyr::slice_max(ratio, n = 10)  # 10 highest
  )

ggplot(sign_go_terms_high, aes(x = mean_WT, y = mean_2E2)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "GO term proteome allocation: WT vs 2E2",
    x = "Mean proteome fraction (WT)",
    y = "Mean proteome fraction (2E2)"
  ) +
  geom_text_repel(
    data = label_terms,
    aes(label = significant_terms),        # adjust if your column name differs
    size = 4,
    min.segment.length = 0
  ) +
  theme_bw(base_size = 14)


plots <- list()

# Loop over each GO
for (i in seq_along(terms_to_plot)) {
  term <- terms_to_plot[i]
  
  # Filter and reshape for that term
  df <- sign_go_terms %>%
    filter(significant_terms == term) %>%
    dplyr::select(significant_terms, starts_with("WT_sum"), starts_with("X2E2_sum")) %>%
    pivot_longer(
      cols = starts_with("WT_sum") | starts_with("X2E2_sum"),
      names_to = c("condition", "rep"),
      names_pattern = "(WT|X2E2)_sum_(\\d+)",
      values_to = "fraction"
    ) %>%
    mutate(
      condition = recode(condition, WT = "WT", X2E2 = "2E2"),
      condition = factor(condition, levels = c("WT", "2E2")) # ensure WT first
    )
  
  # Compute per-term y max
  max_y <- 100*max(df$fraction, na.rm = TRUE)
  
  # Plot
  p <- ggplot(df, aes(x = condition, y = 100*fraction, color = condition, fill = condition)) +
    geom_bar(
      stat = "summary", fun = "mean", position = "dodge",
      width = 0.7, alpha = 0.7, color = "black"
    ) +
    geom_point(size = 3, alpha = 0.9, position = position_jitter(width = 0.05)) +
    labs(
      title = term,
      x = "",
      y = "Proteome fraction (%)"
    ) +
    scale_fill_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0")) +
    scale_color_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0")) +
    coord_cartesian(ylim = c(0, max_y)) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  plots[[i]] <- p
}
# Combine first 3 plots in one row, Fig 4I-K
up <- ggarrange(plotlist = plots[1:3], nrow = 1, ncol = 3)
up
ggsave("6_results/final_figures/figure_2/20260204_upregulated_2E2.pdf", up, width = 10, height = 3, dpi = 300)

# Combine next 3 plots in one row, FigF-H
down <- ggarrange(plotlist = plots[4:6], nrow = 1, ncol = 3)
down
ggsave("6_results/final_figures/figure_2/20260204_downregulated_2E2.pdf", down, width = 10, height = 3, dpi = 300)

# Combine next 3 plots in one row, FigC-E
metabolic <- ggarrange(plotlist = plots[7:9], nrow = 1, ncol = 3)
metabolic
ggsave("6_results/final_figures/figure_2/20260204_metabolic_2E2.pdf", metabolic, width = 10, height = 3, dpi = 300)

# Define a function to recursively get all child terms of a GO term
getAllChildren <- function(go_id) {
  # Get direct children of the GO term
  children <- as.list(GOBPOFFSPRING)[[go_id]]
  if (is.null(children)) {
    return(character(0))
  }
  return(children)
}

# Get all children for GO:0008152
children_terms <- getAllChildren("GO:0044238")

# Print the results
print(children_terms)











# Now add go_terms together:
# Summing up the values for selected GO terms
selected_terms <- c("ribosomal small subunit export from nucleus [GO:0000056]", "ribosomal large subunit export from nucleus [GO:0000055]", 
                    "ribosomal small subunit biogenesis [GO:0042274]")  
summed_data <- results_GO_terms %>%
  filter(grepl(paste(selected_terms, collapse = "|"), significant_terms)) %>%
  summarize(
    total_mean_WT = sum(mean_WT),
    total_mean_2E2 = sum(mean_2E2),
    total_WT_sum_1 = sum(WT_sum_1),
    total_WT_sum_2 = sum(WT_sum_2),
    total_WT_sum_3 = sum(WT_sum_3),
    total_X2E2_sum_1 = sum(X2E2_sum_1),
    total_X2E2_sum_2 = sum(X2E2_sum_2),
    total_X2E2_sum_3 = sum(X2E2_sum_3)
  )

# Reshape data for plotting
plot_data <- summed_data %>%
  pivot_longer(
    cols = everything(),
    names_to = "category",
    values_to = "value"
  ) %>%
  mutate(
    condition = case_when(
      grepl("mean_WT|WT_sum", category) ~ "WT",
      grepl("mean_2E2|X2E2_sum", category) ~ "2E2"
    ),
    variable = case_when(
      grepl("mean", category) ~ "Mean",
      grepl("sum", category) ~ "Sum"
    )
  )

# Plotting
ggplot(plot_data, aes(x = variable, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Summed Mean and Sums of Selected GO Terms",
    x = "Metric",
    y = "Value",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 15) +
  scale_fill_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0"))

# Plot up and downregulated GO terms:
upregulated <- read.csv("6_results/proteomics/proteome_fraction/upregulated_GO_terms.csv")
upregulated_long <- upregulated %>%
  pivot_longer(
    cols = starts_with("WT_") | starts_with("KO_"),
    names_to = c("Condition", "Replicate"),
    names_pattern = "(WT|KO)_(\\d+)",
    values_to = "Value"
  ) %>% drop_na()

# Replace "KO" with "2E2" and reorder 'Condition' levels
upregulated_long <- upregulated_long %>%
  mutate(
    Condition = recode(Condition, "KO" = "2E2"),
    Condition = factor(Condition, levels = c("WT", "2E2")) # Explicitly set order
  )
  

p <- ggplot(upregulated_long, aes(x = GO_term, y = Value, fill = Condition)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.7,color = "black",alpha = 0.7) +
  geom_point(aes(color = Condition), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
             size = 3) +
  labs(
    title = "2E2 upregulated GO-terms",
    x = "GO Term",
    y = "Proteome Fraction (%)",
    fill = "Condition",
    color = "Condition"
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0")) +
  scale_color_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0"))
p
ggsave("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/WesternOntario_collab/6_results/final_figures/figure_2/upreg_goterm_6_6.pdf", plot = p, width = 6, height = 6, units = "in", dpi = 300)


# Downregulated:
downregulated <- read.csv("6_results/proteomics/proteome_fraction/downregulated_GO_terms.csv")  
downregulated_long <- downregulated %>%
  pivot_longer(
    cols = starts_with("WT_") | starts_with("KO_"),
    names_to = c("Condition", "Replicate"),
    names_pattern = "(WT|KO)_(\\d+)",
    values_to = "Value"
  ) %>% drop_na()

# Replace "KO" with "2E2" and reorder 'Condition' levels
downregulated_long <- downregulated_long %>%
  mutate(
    Condition = recode(Condition, "KO" = "2E2"),
    Condition = factor(Condition, levels = c("WT", "2E2")) # Explicitly set order
  )


p <- ggplot(downregulated_long, aes(x = GO_term, y = Value, fill = Condition)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.7,color = "black",alpha = 0.7) +
  geom_point(aes(color = Condition), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
             size = 3) +
  labs(
    title = "2E2 downregulated GO-terms",
    x = "GO Term",
    y = "Proteome Fraction (%)",
    fill = "Condition",
    color = "Condition"
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0")) +
  scale_color_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0"))
p 
ggsave("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/WesternOntario_collab/6_results/final_figures/figure_2/downreg_goterm_6_6.pdf", plot = p, width = 6, height = 6, units = "in", dpi = 300)


#### For printing single go terms: ---------------------------------------------

# go_terms <- c("translation \\[GO:0006412\\]","glycolytic process \\[GO:0006096\\]")
# gsub("\\\\", "", "translation \\[GO:0006412\\]")
# 
# for (term in go_terms) {
#   # Create a logical vector for rows where the go_term is present
#   rows_with_term <- grepl(term, data$Gene.Ontology..biological.process.)
#   print(sum(rows_with_term))
#   
# 
#     clean_term <- gsub("\\\\", "", term)
#     
#     # Filter the dataframe to keep only rows where the term is present
#     temp_dataframe <- data[rows_with_term, ]
#     
#     # Pivot to a longer format for manipulation
#     temp_df_long <- temp_dataframe %>%
#       pivot_longer(
#         cols = starts_with("norm_protfract_WT_") | starts_with("norm_protfract_2E2_"), # Select columns for replicates
#         names_to = "replicate",
#         values_to = "prot_fraction"
#       )
#     
#     temp_df_long <- temp_df_long %>%
#       mutate(
#         condition = ifelse(grepl("WT", replicate), "WT", "2E2")
#       )
#     
#     # Sum protein fractions per replicate, unit, and condition
#     go_sums <- temp_df_long %>%
#       group_by(condition, replicate) %>%
#       summarise(sum_fraction = sum(prot_fraction, na.rm = TRUE), .groups = "drop")
#     
#     go_sums$condition <- factor(go_sums$condition, levels = c("WT", "2E2"))
#     
#     # Create the plot
#     plot <- ggplot(go_sums, aes(x = condition, y = 100 * sum_fraction, color = condition)) +
#       geom_bar(stat = "summary", fun = "mean", alpha = 0.7, position = position_dodge(width = 0.8), aes(fill = condition)) + # Bar for mean
#       geom_jitter(width = 0.2, size = 2) + # Jitter for points
#       labs(
#         title = paste0("Proteome Fraction (", clean_term, ")"),
#         x = "",
#         y = "Proteome Fraction (%)"
#       ) +
#       theme_bw() +
#       theme(legend.position = 'none')
#     
#     # Print the plot
#     print(plot)
# }  # Close the 'for' loop
