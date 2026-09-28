## CHO UPS2 analysis 
# Programmer: Cyriel Huijer
# Date: 25 Sept 2024
# Data from: CH00028 (UPS2 run)
# Filename: 20240802_Exploris1_SA_CH_CHO_UPS2_18ug_50percent.raw

library(ggplot2)
library(tidyverse)
library(ggpmisc)
setwd("C:/Users/chuijer/surfdrive/PhD/Projects/WesternOntario_collab/")

ups2_meta <- read.csv("4_supplementalData/UPS_proteome_set.csv",header = TRUE)
raw_data <- read.csv("6_results/proteomics/ups2/MaxQuant_UPS2_output/proteinGroups.txt", header = TRUE, sep = "\t")

# Filter out contaminants etc:
filtered_data <- raw_data[raw_data$Only.identified.by.site != "+" &
                            raw_data$Reverse != "+" &
                            raw_data$Potential.contaminant != "+", ]
filtered_data_with_ups <- filtered_data %>%
  filter(grepl("ups", Protein.IDs)&
  !grepl("^G", Protein.IDs))%>%
  mutate(Protein.IDs_6 = substr(Protein.IDs, 1, 6))  # Extract the first 6 characters

merged_data <- merge(filtered_data_with_ups, ups2_meta, 
                     by.x = "Protein.IDs_6", 
                     by.y = "UniProt.Accession.Number", 
                     all.x = TRUE)  # Keep all rows from filtered_data_with_ups_no_G

merged_data <- merged_data %>%
  mutate(
    iBAQ_log10 = log10(iBAQ + 1),  # Log10 transformation of iBAQ
    # The UPS2 amount in fmol is multiplied by 3/10.6, as only 3 ug out of 10.6 ug was used. 
    UPS2_log10 = log10(UPS2.Amount..fmol.*(3/10.6) + 1)  # Log10 transformation of UPS2Amount..fmol.
  )

filtered_merged_data <- merged_data %>%
  filter(Protein.IDs_6 != "P08758" &
           Protein.IDs_6 != "P02787")

x = filtered_merged_data$UPS2_log10
y = filtered_merged_data$iBAQ_log10
formula <- y ~ x
ggplot(filtered_merged_data, aes(x = UPS2_log10, y = iBAQ_log10)) +
  geom_point(alpha = 0.6,size = 3) +  # Add points with some transparency
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Add a linear regression line
  labs(
    title = "UPS2 Spike-in iBAQ signal",
    x = "UPS2 (log10-transformed)",
    y = "iBAQ (log10-transformed)"
  ) +
  stat_poly_eq(use_label(c("eq","adj.R2","P")),formula=formula)+
  theme_bw(base_size = 15) # Use a minimal theme

model <- lm(iBAQ_log10~UPS2_log10, data=filtered_merged_data)
summary(model)# Adjusted R2 = 0.92
model$coefficients["(Intercept)"]

# Formula for calculating the protein amounts for other proteins: 
# iBAQ = 5.99 + 0.967[femtomoles] 
# To retrieve the femtomoles per protein, calculate as follows:
# [femtomoles] = (iBAQ - 5.99)/0.967


filtered_data <- filtered_data %>%
  mutate(
    iBAQ_log10 = log10(iBAQ + 1),  # Log10 transformation of iBAQ
  )
filtered_data$femtomoles <- 10^((filtered_data$iBAQ_log10 - model$coefficients["(Intercept)"])/model$coefficients["UPS2_log10"])
filtered_data <- filtered_data[!grepl("ups", filtered_data$Protein.IDs), ]

filtered_data$Accession <- substr(filtered_data$Protein.IDs, 1, 6)


quant_results <- filtered_data[, c("Accession", "femtomoles")]
write.csv(quant_results, file = "5_processedData/UPS2_quantification_results.csv", row.names = FALSE)






