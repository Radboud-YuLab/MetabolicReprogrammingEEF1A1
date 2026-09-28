# Carbon flux in cells analysis
# written by: CYriel Huijer
# last updated: 22 June 2026
library(tidyverse)
library(reshape2)
#clear env
rm(list=ls())

# working directory to /code
setwd("")

# Load in file:
data <- read.csv(file = "data/modeling/exchange_fluxes_solution.csv")

# Number of carbons per metabolite from here:
# https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/

carbons <- read.csv(file = "data/modeling/carbon_counts.csv")

data <- data %>% filter(WT_mean < -2e-2 | WT_mean > 2e-5)

# merge so carbons match
data <- data %>%
  left_join(carbons, by = c("met" = "metabolite"))

# calculate total carbon input/release per metabolite
data$WT_carbon <- data$WT_mean*data$carbons
data$KO_carbon <- data$KO_mean*data$carbons

# Calculate total carbon input
total_WT_carbon_input <- sum(data$WT_carbon[data$WT_carbon < 0])
total_KO_carbon_input <- sum(data$KO_carbon[data$KO_carbon < 0])

# Select rows for lactate, CO2 and alanine:
release_metabolites <- c("ala__L_e", "lac__L_e", "co2_e")
release_data <- data %>% 
  filter(met %in% release_metabolites)
release_data$WT_fraction <- abs(release_data$WT_carbon / total_WT_carbon_input)
release_data$KO_fraction <- abs(release_data$KO_carbon / total_KO_carbon_input)

# Add biomass:
biomass_WT_fraction <- 1 - sum(release_data$WT_fraction)
biomass_KO_fraction <- 1 - sum(release_data$KO_fraction)

# Add a biomass row
biomass_row <- data.frame(
  met = "biomass",
  WT_carbon = NA,
  KO_carbon = NA,
  carbons = NA,
  WT_fraction = biomass_WT_fraction,
  KO_fraction = biomass_KO_fraction
)

release_data <- bind_rows(release_data, biomass_row)

# Melt data to long format
release_data_melt <- melt(release_data,
                  id.vars = "met",
                  measure.vars = c("WT_fraction", "KO_fraction"),
                  variable.name = "condition",
                  value.name = "fraction")
# Change labels and factor:
release_data_melt$condition <- factor(release_data_melt$condition,
                              levels = c("WT_fraction", "KO_fraction"),
                              labels = c("WT", "2E2"))


# Plot stacked barplot
ggplot(release_data_melt, aes(x = condition, y = fraction, fill = met)) +
  geom_bar(stat = "identity") +
  labs(x = "Condition", y = "Fraction of Carbon Input", fill = "Metabolite") +
  theme_bw()

# Rename metabolites to friendly names for colors
release_data_melt$met <- recode(release_data_melt$met,
                                "ala__L_e" = "alanine",
                                "lac__L_e" = "lactate",
                                "co2_e" = "CO2",      # Capitalize CO2 for style
                                "biomass" = "biomass")

conc_colors <- c(
  "alanine" = "#666666",
  "biomass" = "#1b9e77",
  "lactate" = "#d95f02",
  "CO2"     = "#7570b3"
)

# Set factor levels in order bottom to top
release_data_melt$met <- factor(release_data_melt$met,
                                levels = c("lactate", "CO2", "alanine", "biomass"))

p <- ggplot(release_data_melt, aes(x = condition, y = fraction, fill = met)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = conc_colors) +
  labs(x = "Condition", y = "Carbon release (%)", fill = "Metabolite") +
  theme_bw(base_size = 15)

#Fig 3C
p
ggsave("../../../../6_results/final_figures/carbon_release.pdf", plot = p, width = 6, height = 6, units = "in")
