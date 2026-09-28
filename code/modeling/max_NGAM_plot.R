# Script for plotting max(NGAM):
# Cyriel Huijer - 24 Jan 2025

library(tidyverse)
root <- "C:/Users/chuijer/surfdrive/PhD/Projects/WesternOntario_collab/"
setwd(root)


# Load in csv
file <- read.csv("6_results/metabolic_modeling/20250120_rsWithRosemary/max_ngam.csv")
file$condition <- c("WT","2E2")
file$condition <- factor(file$condition, levels = c("WT", "2E2"))


# Plot the bar plot with individual points
p <- ggplot(file, aes(x = condition, y = max_NGAM, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black", alpha = 0.7) +  # Bars showing the mean
  theme_bw(base_size = 15) +
  labs(title = "max(ATP) production", 
       x = "Condition", 
       y = "max(NGAM) mmol/gDW/hr") +
  scale_fill_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +  # Custom colors for bars
  scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +  # Custom colors for points
  theme(legend.position = "none")  # Hide legend for a cleaner look
# Fig 2A
p
