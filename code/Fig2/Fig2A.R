# Script for Fig2A max(NGAM):
# Cyriel Huijer - 24 Jan 2025
# last updated: 22 June 2026

library(tidyverse)

# clear environment
rm(list=ls())

# working directory to /code
setwd("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/WesternOntario_collab/9_manuscript/code/")

# Load in csv
file <- read.csv("data/modeling/max_ngam.csv")
file$condition <- c("WT","2E2")
file$condition <- factor(file$condition, levels = c("WT", "2E2"))

# Plot the bar plot with individual points
p <- ggplot(file, aes(x = condition, y = max_NGAM, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black", alpha = 0.7) + 
  theme_bw(base_size = 15) +
  labs(title = "max(ATP) production", 
       x = "Condition", 
       y = "max(NGAM) mmol/gDW/hr") +
  scale_fill_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +  
  scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) + 
  theme(legend.position = "none")
# Fig 2A
p
