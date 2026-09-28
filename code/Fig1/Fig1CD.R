## Volcano plot figure 1CD
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


# Plot EEF1A1/EEF1A2 expression levels (Fig 1C-D)

# Filter the row where Accession column contains (EEF1A1 = "G3HH39", EEF1A2 = "G3HGY8")
g3hh39_row <- complete_table[complete_table$Accession == "G3HH39", ]

# Select the columns that start with 'fmol_per_ug_WT_' and 'fmol_per_ug_X2E2_'
wt_columns <- grep("fmol_per_ug_WT_", names(g3hh39_row), value = TRUE)
ko_columns <- grep("fmol_per_ug_X2E2_", names(g3hh39_row), value = TRUE)

# Extract the WT and KO values
wt_values <- as.numeric(g3hh39_row[, wt_columns])
ko_values <- as.numeric(g3hh39_row[, ko_columns])

# Create a dataframe for plotting with KO labeled as 2E2
g3hh39_data <- data.frame(
  Condition = factor(c(rep("WT", length(wt_values)), rep("2E2", length(ko_values))),
                     levels = c("WT", "2E2")),  # Set KO as 2E2 and WT after
  Value = c(wt_values, ko_values)
)

# Plot the bar plot with individual points
p <- ggplot(g3hh39_data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black", alpha = 0.7) +  
  geom_point(width = 0.2, size = 3, aes(color = Condition), alpha = 0.9) +  
  theme_bw(base_size = 15) +
  labs(title = "Absolute Protein Quantity (EEF1A1)", 
       x = "Condition", 
       y = "fmol/µg total protein") +
  scale_fill_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +  
  scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +  
  theme(legend.position = "none") 
# Fig 1C
p

# EEF1A2 (Fig 1D)

# Filter the row where Accession column contains (EEF1A1 = "G3HH39", EEF1A2 = "G3HGY8")
g3hh39_row <- complete_table[complete_table$Accession == "G3HGY8", ]

# Select the columns that start with 'fmol_per_ug_WT_' and 'fmol_per_ug_X2E2_'
wt_columns <- grep("fmol_per_ug_WT_", names(g3hh39_row), value = TRUE)
ko_columns <- grep("fmol_per_ug_X2E2_", names(g3hh39_row), value = TRUE)

# Extract the WT and KO values
wt_values <- as.numeric(g3hh39_row[, wt_columns])
ko_values <- as.numeric(g3hh39_row[, ko_columns])

# Create a dataframe for plotting with KO labeled as 2E2
g3hh39_data <- data.frame(
  Condition = factor(c(rep("WT", length(wt_values)), rep("2E2", length(ko_values))),
                     levels = c("WT", "2E2")),  # Set KO as 2E2 and WT after
  Value = c(wt_values, ko_values)
)

# Fig 1D
p <- ggplot(g3hh39_data, aes(x = Condition, y = Value, fill = Condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black", alpha = 0.7) +  
  geom_point(width = 0.2, size = 3, aes(color = Condition), alpha = 0.9) +
  theme_bw(base_size = 15) +
  labs(title = "Absolute Protein Quantity (EEF1A2)", 
       x = "Condition", 
       y = "fmol/µg total protein") +
  scale_fill_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +
  scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) + 
  theme(legend.position = "none") 
# Fig 1D
p




