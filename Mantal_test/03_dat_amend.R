# Clear workplace
#------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
#------------------------------------------------------------------------------
library(tidyverse)

# Load data
#------------------------------------------------------------------------------

df <- read_csv("data/02_wrangle.csv")

# Amend
#------------------------------------------------------------------------------
df <- df %>%
  mutate(condition = str_c(oxygen.treat,culture.kno3, sep = '_'))

df_meta <- df %>%  
  group_by(index,strain, condition) %>% 
  summarise(n(), .groups = "drop")

df <- df %>%
  select(index, domain, phylum,class,order, family, genus, abundance) %>% 
  pivot_wider(names_from = index, values_from = abundance)

# relative abundance
df <- df %>% mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE)))
#df %>% summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

df_long <- df %>%
  pivot_longer(cols = where(is.numeric), names_to = 'sample', values_to = 'rel_abundance') %>% 
  select('genus', 'phylum', sample, rel_abundance)

df_long <- df_long %>%
  left_join(df_meta, by=c("sample" = "index"))

df_longer <- df_long %>% group_by(strain, condition, phylum, genus) %>% 
  summarise(rel_abundance = mean(rel_abundance), .groups = "drop")
#df_longer %>%  group_by(strain, condition) %>% summarise(sum(rel_abundance), .groups = "drop")

threshold <- 0.01
df_longer <- df_longer %>%
  group_by(strain, condition) %>%
  mutate(
    phylum = if_else(rel_abundance < threshold, "Other", phylum),
    genus = if_else(rel_abundance < threshold, "Other", genus)
  ) %>%
  group_by(strain, condition, phylum, genus) %>%
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop")



######################## Factorize variables ###################################
# Factorize variables to plot them in desired order

# Extract all unique phyla and sort them in alphabetic order
phylum_order <- df_longer %>% pull(phylum) %>% unique() %>% sort()
# "Other" is plasted in the last position of the list
phylum_order<- c(setdiff(phylum_order, "Other"), "Other")

# Reorder the phylum factor
df_longer <- df_longer %>%
  mutate(phylum = factor(phylum, levels = phylum_order))

# Arrange data by phylum and genus to define the order
# then use unique() to capture genera in that order
genus_order <- df_longer %>%
  arrange(phylum, genus) %>%
  pull(genus) %>%
  unique() %>% 
  rev()

# Set genus as a factor with the custom order
df_longer <- df_longer %>%
  mutate(genus = factor(genus, levels = genus_order))

# Reorder the 'strain' factor to control facet order
strain_order <- c("LRM1b", "PCB", "LRM2N6", "PCE", "FB10N2")
df_longer <- df_longer %>%
  mutate(strain = factor(strain, levels = strain_order))


# Write data
#------------------------------------------------------------------------------

df_longer %>% saveRDS("data/03_amend.csv")