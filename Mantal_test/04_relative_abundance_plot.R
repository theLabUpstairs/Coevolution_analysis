# Clear workplace
#------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
#------------------------------------------------------------------------------
library(tidyverse)

# Load data
#------------------------------------------------------------------------------

df <- readRDS("data/03_amend.csv")

# Plot
#------------------------------------------------------------------------------

df <- df %>%
  mutate(condition = recode(
    condition,
    "anoxic_control" = "A_ctl",
    "microaerophilic_control" = "M_ctl",
    "anoxic_kno3" = "A_kno3",
    "microaerophilic_kno3" = "M_kno3"
  ))

phylum_colors <- c(
  "Bacteroidota" = "#39FF14",
  "Firmicutes" = "#FF2079",
  "Campilobacterota" = "#FFD300",
  "Desulfobacterota" = "#FF6F61",
  "Proteobacteria" = "#7F00FF",
  "Spirochaetota" = "#FF5F1F",
  "Verrucomicrobiota" = "#0ABAB5",
  "other" = "#808080"
)


ggplot(df, aes(x = condition, y = genus)) +
  # Shading for conditions
  annotate("rect", xmin = 0.4, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.5) +
  #annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = "lightpink", alpha = 0.5) +
  geom_point(aes(size = rel_abundance, fill = phylum), shape = 21, stroke = 1.1) +
  facet_grid(~ strain) +
  scale_fill_manual(values = phylum_colors,
                    guide = guide_legend(override.aes = list(size = 5))) + 
  # Custom size legend labels for relative abundance
  scale_size_continuous(name = "Relative Abundance", 
                        breaks = c(0.1,0.3, 0.5),
                        labels = c("10%",  "30%",  "50%")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.grid.major = element_line(color = "lightgray", size = 0.1),  # Light gray major grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  labs(y = "Genus", x = "Condition", fill = "Phylum")



# Save figures
#------------------------------------------------------------------------------
ggsave("figures/relative_abundance_bubble_plot.png", width = 8, height = 6)
