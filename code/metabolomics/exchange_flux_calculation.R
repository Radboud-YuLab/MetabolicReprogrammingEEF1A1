# Exchange flux calculation UWO Collaboration
# Programmer: Cyriel Huijer
# Last updated: 22 June 2026
library(dplyr)
library(tidyr)
library(ggrepel)

# set wd to /code
setwd("")

cell_number_data <- read.csv("data/growth_rate/cell_numbers_24h_doublingTime.csv")
amino_acids <- read.csv("data/metabolomics/amino_acids_HisTrpupdated.csv")
gluc_lact <- read.csv("data/metabolomics/glucose_lactate.csv")
manual_inspection <- read.csv(file = "data/metabolomics/manual_inspection_curves_keep_ala.csv")

# Format amino acids & lactate/glucose so it can be fitted in a single dataframe;

for (i in 1:nrow(gluc_lact)){
  if (gluc_lact$condition[i] == "KO"){
    gluc_lact$condition[i] <- "2E2"
  }
}
for (i in 1:nrow(amino_acids)){
  if (amino_acids$condition[i] == "KO"){
    amino_acids$condition[i] <- "2E2"
  }
}

for (i in 1:nrow(cell_number_data)){
  if (cell_number_data$condition[i] == "KO"){
    cell_number_data$condition[i] <- "2E2"
  }
}


amino_acids <- amino_acids %>%
  dplyr::select(condition, replicate, timepoint, amino_acid, concentration_sample_uM)
colnames(amino_acids)[4] <- "compound"

gluc_lact_long <- gluc_lact %>%
  pivot_longer(cols = c(glucose, lactate), 
               names_to = "compound", 
               values_to = "concentration_sample")

# Calculate molarity of glucose and lactate
gluc_lact_long$concentration_sample_uM <- NA
for (i in 1:nrow(gluc_lact_long)){
  if (gluc_lact_long$compound[i] == "glucose"){
    gluc_lact_long$concentration_sample_uM[i] <- (gluc_lact_long$concentration_sample[i]/180.156)*1e6
  }
  else{
    gluc_lact_long$concentration_sample_uM[i] <- (gluc_lact_long$concentration_sample[i]/90.08)*1e6
  }
}
gluc_lact_long <- as.data.frame(gluc_lact_long)

# Ensure column names match and remove unnecessary columns
gluc_lact_long <- gluc_lact_long %>%
  dplyr::select(condition, replicate, timepoint, compound, concentration_sample_uM)

# Convert data types if necessary
amino_acids$condition <- as.character(amino_acids$condition)
gluc_lact_long$condition <- as.character(gluc_lact_long$condition)
amino_acids$compound <- as.character(amino_acids$compound)
gluc_lact_long$compound <- as.character(gluc_lact_long$compound)

metabolite_concentrations <- rbind(gluc_lact_long,amino_acids)

metabolite_concentrations$concentration_sample_mM <- metabolite_concentrations$concentration_sample_uM /1e3

# Fix cell number as cells were grown in 96 well plates instead of 10 cm plates:

#cell_number_data$cells_adjusted <- cell_number_data$cells * (56.7/0.32)
colnames(cell_number_data)[7] <- "timepoint"
merged_df <- left_join(metabolite_concentrations, cell_number_data, by = c("condition", "timepoint"))
merged_df$volume_corr <- merged_df$concentration_sample_mM * 0.01 # Cells were cultured in 10 mL
merged_df$mmol_cell <- merged_df$volume_corr / merged_df$cells
merged_df$condition <- factor(merged_df$condition, levels = c("WT", "2E2"))

merged_df$mmol_gDW <- merged_df$mmol_cell *1/264e-12 # cell dry weight of 260 pg per cell taken

tableS1 <- merged_df %>%
  dplyr::select(condition, replicate, timepoint, compound, concentration_sample_mM,
                cells, mmol_cell, mmol_gDW)
colnames(tableS1)[4] <- "metabolite"
tableS1 <- tableS1 %>% filter(metabolite %in% c('glucose','lactate','Asn','Gln','Asp','Glu','Met','Phe','Pro','Ser','Trp','Tyr','Val'))

#write.csv(tableS1,"9_manuscript/supplementary_tables/TableS1.csv",row.names = F)


## Average cell number calculation + plots ################################################# 

# Plot exchange flux calculations per metabolite
AA_list <- unique(merged_df$compound)
p = list()
for (i in 1:length(AA_list)){
  AA <- AA_list[i]
  filtered_data <- merged_df %>%
    filter(compound == AA)
  # <- subset(mean_concentration, amino_acid == AA)
  p[[i]] <- ggplot(filtered_data,aes(x=timepoint,y=mmol_gDW))+
    geom_point(aes(shape=as.factor(replicate)))+
    ylim(0,max(filtered_data$mmol_gDW))+
    labs(title = filtered_data$compound)+
    ylab(paste("mmol/gDW of ",AA))+
    xlab("Timepoint (in hours)")+
    facet_wrap(~condition)+
    theme_bw(base_size = 15)+
    geom_smooth(method="lm",se = FALSE)
  print(p[[i]])
  
}

# Plot without the legend

AA_list <- unique(merged_df$compound)
p = list()

for (i in 1:length(AA_list)){
  AA <- AA_list[i]
  filtered_data <- merged_df %>%
    filter(compound == AA)
  
  # Create the plot and remove the legend
  p[[i]] <- ggplot(filtered_data, aes(x = timepoint, y = mmol_gDW)) +
    geom_point(aes(shape = as.factor(replicate))) +
    ylim(0, max(filtered_data$mmol_gDW)) +
    labs(title = AA) +
    ylab(paste("mmol/cell of ", AA)) +
    xlab("Timepoint (in hours)") +
    facet_wrap(~condition) +
    theme_bw(base_size = 15) +
    geom_smooth(method = "lm", se = FALSE) +
    theme(legend.position = "none")  # Hides the legend
  
  print(p[[i]])
}

# # Save the plots to a PDF with 5 plots per page
# pdf("amino_acid_plots_no_legend.pdf", width = 10, height = 12)  # Adjust width and height as needed
# marrangeGrob(grobs = p, nrow = 2, ncol = 2)  # 2 rows, 3 columns for a 5-plot layout
# dev.off()

# ------ updated 02 January 2025 -----------------------------------------------
# In merged_df, we have per amino acid and conditon the variable timepoint and mmol/cell, we want to calculate the slopes per single replicate.
# Files will be created and manual QC of single replicates is performed.  

merged_df_qc <- merge(merged_df, manual_inspection, by = c("replicate", "compound"))
metab_concentration <- subset(merged_df_qc, keep_manual_inspection != "n")

output_dir <- ""
AA_list <- unique(metab_concentration$compound)
p = list()
for (i in 1:length(AA_list)){
  AA <- AA_list[i]
  filtered_data <- metab_concentration %>%
    filter(compound == AA)
  p[[i]] <- ggplot(filtered_data,aes(x=timepoint,y=mmol_gDW,shape = as.factor(replicate)))+
    geom_point()+
    labs(title = filtered_data$compound)+
    ylab(paste("mmol/gDW ",AA))+
    xlab("Timepoint (in hours)")+
    facet_wrap(~condition)+
    theme_bw(base_size = 15)+
    geom_smooth(method = "lm", se = F) +
    ylim(0,max(filtered_data$mmol_gDW)) +
    theme(legend.position = "none")
  file_path <- file.path(output_dir, paste0(AA, "_concentration_plot.pdf"))
  #ggsave(file_path, plot = p[[i]], width = 6, height = 6, units = "in", dpi = 300)
  print(p[[i]])
}

calculate_slope <- function(data) {
  fit <- lm(mmol_gDW ~ timepoint, data = data)
  slope <- coef(fit)["timepoint"]
  R2 <- summary(fit)$r.squared
  p_value <- summary(fit)$coefficients["timepoint", "Pr(>|t|)"]
  return(list(slope = slope, R2 = R2, p_value = p_value))
}

# No cutoff for R2 was taken, since all metabolites are above R2 > 0.5, except for lactate in 2E2s,
# this is acceptable, since there seems to be no lactate production in 2E2s!

slopes_df <- metab_concentration %>%
  group_by(condition, compound, replicate) %>%
  summarize(
    slope = calculate_slope(cur_data())$slope,
    R2 = calculate_slope(cur_data())$R2,
    p_value = calculate_slope(cur_data())$p_value,
    .groups = 'drop'
  ) 

# reverse slope for lactate plot
slopes_df$slope_reverse <- slopes_df$slope

for (i in 1:nrow(slopes_df)){
  if (slopes_df$compound[i] == 'lactate'){
    slopes_df$slope_reverse[i] <- slopes_df$slope_reverse[i]
  }
  else {
    slopes_df$slope_reverse[i] <- slopes_df$slope[i] * -1
  }
}

# For fig 1E-F
#write.csv(slopes_df, "data/exchange_fluxes/exch_fluxes_positive.csv", row.names = F)

unique_aminos <- unique(slopes_df$compound)
# Plot uptake/release plots:
for (amino in unique_aminos) {
  # Filter the data fo2r the current amino acid
  amino_data <- slopes_df %>%
    filter(compound == amino)
  
  # Calculate min and max values for y-axis scaling
  y_min <- min(amino_data$slope_reverse, na.rm = TRUE)
  y_max <- max(amino_data$slope_reverse, na.rm = TRUE)
  
  # Ensure 'replicate' is a factor for shape mapping
  amino_data$replicate <- as.factor(amino_data$replicate)
  
  # Create the plot
  p <- ggplot(amino_data, aes(x = condition, y = slope_reverse, color = condition,fill=condition)) +
    geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.7,color = "black",alpha = 0.7) +
    geom_point(aes(shape = replicate),size = 3) +  # Adjust jitter to prevent overlap
    labs(title = paste("Uptake/release of", amino),
         x = "Condition",
         y = paste("Exchange flux of ",amino,"(in mmol/gDW/hr)")) +
    scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043"))+
    scale_fill_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0"))+
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none") +  # Rotate x-axis labels for better readability
    scale_y_continuous(limits = c(min(y_min, 0), max(y_max, 0)))  # Set y-axis limits
  
  # Print the plot
  print(p)

  # ggsave(filename = paste("6_results/metabolomics/uptake_release_plots_after_QC_gDW/uptake_release_slope_", gsub(" ", "_", amino), ".pdf", sep = ""),
  #        plot = p,
  #        width = 6,
  #        height = 6,
  #        dpi = 300)
}

## Here the plot shows the differences in exch. fluxes and whether there is a sign. difference. 

# For fig 1G: 
#write.csv(metab_concentration,"data/exchange_fluxes/metab_concentration.csv")

slopes_per_rep_metab <- metab_concentration %>%
  group_by(condition, compound, replicate) %>%
  summarize(
    slope = calculate_slope(cur_data())$slope,
    R2 = calculate_slope(cur_data())$R2,
    .groups = 'drop'
  ) 
slopes_per_rep_metab_filtered <- slopes_per_rep_metab %>%
  group_by(condition, compound) %>%
  filter(n() > 1) %>%  # Keep only groups with more than 1 row
  ungroup()
results_per_metab <- slopes_per_rep_metab_filtered %>%
  group_by(compound) %>%
  summarise(
    mean_slope_2E2 = mean(slope[condition == "2E2"], na.rm = TRUE),
    mean_slope_WT = mean(slope[condition == "WT"], na.rm = TRUE),
    p_value = tryCatch(
      t.test(slope ~ condition)$p.value,
      error = function(e) NA
    )
  )
plot_data <- results_per_metab %>%
  mutate(
    log10_mean_slope_2E2 = log10(abs(mean_slope_2E2)),
    log10_mean_slope_WT = log10(abs(mean_slope_WT))
  )

# Only Asp and lactate are significant, therefore, exchange fluxes will be calculated from all slopes for the other metabolites:
slopes_df_ns <- slopes_df[!slopes_df$compound %in% c("Asp", "lactate"), ]
slopes_df_lb_ub_ns <- slopes_df_ns %>%
  group_by(compound) %>%
  summarize(upper_bound = mean(slope) + sd(slope),
            lower_bound = mean(slope) - sd(slope)) %>%
  ungroup()
upper_bound_ns <- slopes_df_lb_ub_ns %>%
  filter(!is.na(upper_bound)) %>%
  dplyr::select(compound, upper_bound)
lower_bound_ns <- slopes_df_lb_ub_ns %>%
  filter(!is.na(lower_bound)) %>%
  dplyr::select(compound, lower_bound)
upper_bound_ns$WT <-upper_bound_ns$upper_bound
lower_bound_ns$WT <-lower_bound_ns$lower_bound
colnames(upper_bound_ns) <- c("metabolite","2E2", "WT")
colnames(lower_bound_ns) <- c("metabolite","2E2", "WT")  
  
upper_bound <- upper_bound[upper_bound$compound %in% c("Asp", "lactate"), ]
lower_bound <- lower_bound[lower_bound$compound %in% c("Asp", "lactate"), ]
colnames(upper_bound)[1] <- "metabolite"
colnames(lower_bound)[1] <- "metabolite"

ub_final <- rbind(upper_bound,upper_bound_ns)
lb_final <- rbind(lower_bound,lower_bound_ns)

#write.csv(lb_final, "data/exchange_fluxes/for_modeling/exch_fluxes_lb_lm_sign_diff.csv", row.names = FALSE)
#write.csv(ub_final, "data/exchange_fluxes/for_modeling/exch_fluxes_ub_lm_sign_diff.csv", row.names = FALSE)

