## Volcano plot figure 1CD
# Programmer: Cyriel Huijer
# Last updated: 22 June 2026

# clear environment
rm(list=ls())

# set wd to /code
setwd("")

# Function to calculate slope
calculate_slope <- function(data) {
  fit <- lm(mmol_gDW ~ timepoint, data = data)
  slope <- coef(fit)["timepoint"]
  R2 <- summary(fit)$r.squared
  p_value <- summary(fit)$coefficients["timepoint", "Pr(>|t|)"]
  return(list(slope = slope, R2 = R2, p_value = p_value))
}

metab_concentration <- read.csv("data/exchange_fluxes/metab_concentration.csv")

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

# Plot the data
ggplot(plot_data, aes(x = log10_mean_slope_2E2, y = log10_mean_slope_WT)) +
  geom_point(color = "blue", size = 3) +
  labs(
    x = "Log10(Absolute Mean Slope 2E2)",
    y = "Log10(Absolute Mean Slope WT)",
    title = "Log10(Absolute Mean Slopes for 2E2 vs WT)"
  ) +
  theme_bw(base_size = 15) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")

results_per_metab$adjusted_pvalue <- p.adjust(results_per_metab$p_value, method = "BH")


# Create a new dataframe for plotting
plot_data <- results_per_metab %>%
  mutate(
    log10_mean_slope_2E2 = log10(abs(mean_slope_2E2)),
    log10_mean_slope_WT = log10(abs(mean_slope_WT)),
    significant = ifelse(adjusted_pvalue < 0.05, "Yes", "No")  # "Yes" for p_value < 0.05, else "No"
  )

plot_data <- plot_data[!plot_data$compound %in% c("Ala"), ]

# Fig 1G
p <- ggplot(plot_data, aes(x = log10_mean_slope_WT, y = log10_mean_slope_2E2)) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  
  # Points with color and shape based on significance
  geom_point(aes(shape = significant, fill = significant), size = 3.5, color = "black") +
  
  # Labels for selected compounds
  geom_text_repel(
    aes(label = compound),
    size = 4,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.5,
    min.segment.length = 0,
    segment.color = "grey60",
    segment.size = 0.3
  ) +
  labs(
    title = "All exchange fluxes",
    x = expression(Log[10]*"(Absolute Mean Exchange Flux, WT)"),
    y = expression(Log[10]*"(Absolute Mean Exchange Flux, 2E2)")
  ) +
  scale_shape_manual(values = c("Yes" = 21, "No" = 21)) +
  scale_fill_manual(values = c("Yes" = "#E64B35FF", "No" = "#4DBBD5FF")) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")
p
