
# Clear environment and load libraries
#------------------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(ape)
library(vegan)


# Load input files
#------------------------------------------------------------------------------

# Copy Number Variation (CNV)–normalized abundance counts (from Qiime2 pipeline)
df <- read_csv("_raw/all_taxonomy_abundance_CNV.csv")

# Tree file, the 5 Briviate strains (unrooted)
tree <- read.tree("_raw/18S_align_mask.fasta.MFPuf.treefile")


# Script seting
#------------------------------------------------------------------------------

set.seed(42)

# Process data
#------------------------------------------------------------------------------

# Rename tree tip labels to match strain names
tree$tip.label <- c("PCE", "FB10N2", "LRM1b", "PCB", "LRM2N6")

# Calculate phylogenetic distance matrix
phylo_dist <- cophenetic(tree)


# Bray-Curtis (Arcobacter-only) + PCoA Visualization
#------------------------------------------------------------------------------

# Filter the abundance data to only inlcude the arcobacter species
taxa_cols <- grep("f__Arcobacteraceae", colnames(df), value = TRUE)
print(taxa_cols)
taxa_table <- df %>%
  select(index, strain, all_of(taxa_cols))


# Calculate Bray-Curtis distances between samples
taxa_dist <- taxa_table %>%
  column_to_rownames("index") %>%
  select(-strain) %>%
  vegdist(method = "bray")

# PCoA
pcoa_result <- cmdscale(taxa_dist, k = 2, eig = TRUE)
eig <- pcoa_result$eig
var_explained <- eig / sum(eig) * 100

pcoa_df <- as.data.frame(pcoa_result$points) %>%
  rownames_to_column("index") %>%
  left_join(select(taxa_table, index, strain), by = "index")

# Plot
ggplot(pcoa_df, aes(x = V1, y = V2, color = strain)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCoA of Arcobacter (Bray-Curtis)",
    x = paste0("PCoA 1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(var_explained[2], 1), "%)"),
    color = "Strain"
  )
ggsave("figures/pcoa_arcobacter.svg", width = 6, height = 5)

#------------------------------------------------------------------------------
# Calculate mean Bray-Curtis / strain - to use in Mantel test
#------------------------------------------------------------------------------

df_strain <- taxa_table %>%
  group_by(strain) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("strain")

# Ensure no extra spaces or factors
rownames(df_strain) <- str_trim(rownames(df_strain))
rownames(phylo_dist) <- str_trim(rownames(phylo_dist))

# Bray-Curtis matrix per strain
bray_strain <- vegdist(df_strain, method = "bray")

#------------------------------------------------------------------------------
# Mantel Test
#------------------------------------------------------------------------------

# Final check before running test
mantel_result <- mantel(phylo_dist, bray_strain, method = "pearson", permutations = 999)
mantel_result

#Mantel statistic based on Pearson's product-moment correlation 

#Mantel statistic r: 0.2502 
#      Significance: 0.3 
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.463 0.529 0.573 0.586 
#Permutation: free
#Number of permutations: 119



