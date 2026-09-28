root <- "C:/Users/chuijer/surfdrive/PhD/Projects/WesternOntario_collab/"
setwd(root)

# Load libraries:
library(ggplot2)

# Clear workspace:
rm(list=ls())

# working directory to /code
setwd("")

alt_subs <- read.csv("data/modeling/alt_subs.csv")
for (i in 1:nrow(alt_subs)){
  if (alt_subs$subsytem[i] == "NUCLEOTIDES"){
    alt_subs$subsytem[i] <- "NUCLEOTIDE METABOLISM"
  }
}

alt_subs <- alt_subs %>% filter(!subsytem %in% c("UNASSIGNED","TRANSPORT, NUCLEAR","TRANSPORT, MITOCHONDRIAL","TRANSPORT, PEROXISOMAL", "TRANSPORT, EXTRACELLULAR"))

p <- ggplot(alt_subs, aes(
  x = reorder(subsytem, pct_rxn_change),   # order by percentage change
  y = pct_rxn_change,
  size = counts_altered
)) +
  geom_point(color = "#0070C0", alpha = 0.7) +
  coord_flip() +  # flips to make top bubble the highest pct_rxn_change
  scale_size(range = c(3, 12)) +
  labs(
    x = "Subsystem",
    y = "% Reactions Changed",
    size = "Altered Reactions",
    title = "Subsystem Reaction Changes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p
#ggsave("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/WesternOntario_collab/6_results/final_figures/figure_4/20260118_alt_subs_12_6.pdf", plot = p, width = 12, height = 6, units = "in", dpi = 300)
