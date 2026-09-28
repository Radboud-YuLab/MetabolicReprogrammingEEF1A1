# Exchange flux plots 1E-G
# Programmer: Cyriel Huijer
# Last updated: 22 Oct 2025
#   - Added data export for JT

library(tidyverse)

# clear environment
rm(list=ls())

# set to /code
setwd("")

# Function:
calculate_slope <- function(data) {
  fit <- lm(mmol_gDW ~ timepoint, data = data)
  slope <- coef(fit)["timepoint"]
  R2 <- summary(fit)$r.squared
  p_value <- summary(fit)$coefficients["timepoint", "Pr(>|t|)"]
  return(list(slope = slope, R2 = R2, p_value = p_value))
}

slopes_df <- read.csv("data/exchange_fluxes/exch_fluxes_positive.csv")

# factor WT first, 2E2 second
slopes_df$condition <- factor(slopes_df$condition,
                            levels = c("WT", "2E2"))

# Fig 1E
lactate <- slopes_df %>%
  filter(compound == "lactate")

# Ensure 'replicate is a factor for shape mapping
lactate$replicate <- as.factor(lactate$replicate)

# Create the plot
p <- ggplot(lactate, aes(x = condition, y = slope_reverse, color = condition,fill=condition)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.7,color = "black",alpha = 0.7) +
  geom_point(aes(shape = replicate),size = 3) +
  labs(title = paste("Uptake/release of lactate"),
       x = "Condition",
       y = paste("Exchange flux of lactate in mmol/gDW/hr)")) +
  scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043"))+
  scale_fill_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0"))+
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
# Print the plot
print(p)

# Fig 1F
Asp <- slopes_df %>%
  filter(compound == "Asp")

# Ensure 'replicate is a factor for shape mapping
Asp$replicate <- as.factor(lactate$replicate)

# Create the plot
p <- ggplot(Asp, aes(x = condition, y = slope_reverse, color = condition,fill=condition)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.7,color = "black",alpha = 0.7) +
  geom_point(aes(shape = replicate),size = 3) +
  labs(title = paste("Uptake/release of Asp"),
       x = "Condition",
       y = paste("Exchange flux of Asp in mmol/gDW/hr)")) +
  scale_color_manual(values = c("2E2" = "#0070C0", "WT" = "#FF7043"))+
  scale_fill_manual(values = c("WT" = "#FF7043", "2E2" = "#0070C0"))+
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
# Fig 1F
print(p)


