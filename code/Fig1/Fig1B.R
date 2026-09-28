# Script for calculation of growth curves between WT and 2E2 curves
# Updated: 24 Jan 2025
# Programmer: Cyriel Huijer

# Load in packages:
library(tidyverse)

# clear environment
rm(list=ls())

# set working directory to /code
setwd("")

# Load in cell count/protein data
file <- read_csv("data/proteomics/protein_per_cell.csv")
file$sample <- c("WT","2E2","WT","2E2","WT","2E2","WT","2E2")
file$sample <- factor(file$sample, levels = c("WT", "2E2"))


# Figure 1B 
p <- ggplot(file, aes(x = sample, y = pg_cell, fill = sample)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black", alpha = 0.7) + #Bars showing the mean
  geom_point(width = 0.2, size = 3, aes(color = sample), alpha = 0.9) + # Individual data points
  theme_bw(base_size = 15) +
  labs(title = "Protein Content/cell in picograms", 
       x = "Condition", 
       y = "pg total protein/cell") +
  scale_fill_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +# Custom colors for bars
  scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043")) +#Custom colors for points
  theme(legend.position = "none")  # Hide legend for a cleaner look
p
#ggsave("6_results/final_figures/figure_2/protein_content_cell.pdf", plot = p, width = 6, height = 6, units = "in", dpi = 300)

# do t-test:
t_test_result <- t.test(pg_cell ~ sample, data = file)
print(t_test_result)

