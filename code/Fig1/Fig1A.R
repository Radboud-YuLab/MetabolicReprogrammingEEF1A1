# Script for calculation of growth curves between WT and 2E2 curves
# Updated: 07 Jan 2025
# Programmer: Cyriel Huijer

library(tidyverse)
library(broom)

# set working directory to /code
setwd("")

# Load data

data <- read.csv("data/growth_rate/cell_count_overview.csv")

# Functions:

# Function to calculate slope and R^2
calc_slope_R2 <- function(data) {
  # Fit linear regression model
  model <- lm(cells ~ timepoint, data = data)
  
  # Extract slope and R^2
  slope <- coef(model)[2]
  R2 <- summary(model)$r.squared
  intercept <- coef(model)[1]
  
  return(c(slope = slope, R2 = R2,intercept = intercept))
}

for (i in 1:nrow(data)){
  if (data$condition[i] == "KO"){
    data$condition[i] <- "2E2"
  }
}

# Plotting the raw data:

# ggplot(data,aes(x = timepoint,y = log(cells),color = as.factor(experiment_id))) +
#   geom_point() +
#   facet_wrap(~condition) +
#   theme_bw(base_size = 15)
# ggplot(data,aes(x = timepoint,y = cells,color=as.factor(experiment_id))) +
#   geom_point() +
#   facet_wrap(~condition)+
#   ylim(0,max(data$cells))+
#   theme_bw(base_size = 15)

# Between t=12 and t=36, cells appear to be in exponential growth phase, 
# therefore, only datapoints between these timepoints are included to calculate growth rate. 

exp_data <- data %>%
  subset(timepoint >= 12 & timepoint <= 36)
ggplot(exp_data,aes(x=timepoint,y=cells)) +
  geom_point() +
  facet_wrap(~condition)+
  ylim(0,max(exp_data$cells))
# ggplot(exp_data,aes(x = timepoint, y = log(cells),color=as.factor(experiment_id))) +
#   geom_point() +
#   facet_wrap(~condition) +
#   geom_smooth(method = "lm",se = FALSE)

exp_data$log_cells <- log(exp_data$cells)

slopes <- exp_data %>%
  group_by(experiment_id, condition) %>%
  do(tidy(lm(log_cells ~ timepoint, data = .))) %>%
  filter(term == "timepoint") %>%
  dplyr::select(experiment_id, condition, slope = estimate)

#slopes <- slopes %>% filter(slope < 0.0875)

colors <- c("WT" = "#FF7043", "2E2" = "#0070C0")

# Figure 1A
p <- ggplot(slopes, aes(x = condition, y = slope, color = condition, fill = condition)) +
  stat_summary(fun = mean, geom = "bar", alpha = 0.7, color = "black", width = 0.6) +  # bars with black border
  scale_color_manual(values = colors) +  
  geom_point(size = 3, shape = 21, stroke = 1) +  # points with black border
  scale_fill_manual(values = colors) +
  scale_x_discrete(limits = c("WT", "2E2")) +
  labs(x = "Condition", y = "Growth rate", color = "Condition", fill = "Condition") +
  theme_bw(base_size = 15)
p
#ggsave("../../../../6_results/final_figures/growth_rate.pdf", plot = p, width = 6, height = 6, units = "in")


# T test, result not significant
t_test_result <- t.test(slope ~ condition, data = slopes)
print(t_test_result)





