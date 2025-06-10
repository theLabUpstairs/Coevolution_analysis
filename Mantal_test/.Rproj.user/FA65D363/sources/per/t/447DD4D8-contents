
# Clear environment and load libraries
#------------------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(ape)
library(vegan)


# Load input files
#------------------------------------------------------------------------------

#Copy Number Variation (CNV)–normalized abundance counts, Karla's Qiime2 pipeline
df <- read_csv("_raw/all_taxonomy_abundance_CNV.csv")

# Tree file by Courtney, the 5 Briviate strains (unrooted)
tree <- read.tree("_raw/18S_align_mask.fasta.MFPuf.treefile")


# Process data
#------------------------------§------------------------------------------------

# Rename tree tip labels to match strain names
tree$tip.label <- c("PCE", "FB10N2", "LRM1b", "PCB", "LRM2N6")

# Calculate phylogenetic distance matrix
phylo_dist <- cophenetic(tree)


# Bray-Curtis (Arcobacter-only) + PCoA Visualization
#------------------------------------------------------------------------------

"d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio"
"d__Bacteria;p__Campilobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;"
"d__Bacteria;p__Campilobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Malaciobacter"
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Fusibacteraceae;g__Fusibacter"
"d__Bacteria;p__Desulfobacterota;c__Desulfovibrionia;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Desulfovibrio"

taxa_cols <- grep("f__Arcobacteraceae", colnames(df), value = TRUE)
print(taxa_cols)
taxa_table <- df %>%
  select(index, strain, all_of(taxa_cols))

taxa_table %>%
  rowwise() %>%
  mutate(total_abundance = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(index, strain, total_abundance) %>%
  filter(total_abundance > 0) %>% 
  count(strain)

# 2. Calculate Bray-Curtis distances between samples
taxa_dist <- taxa_table %>%
  column_to_rownames("index") %>%
  select(-strain) %>%
  vegdist(method = "bray")

# 3. PCoA
pcoa_result <- cmdscale(taxa_dist, k = 2, eig = TRUE)
eig <- pcoa_result$eig
var_explained <- eig / sum(eig) * 100

pcoa_df <- as.data.frame(pcoa_result$points) %>%
  rownames_to_column("index") %>%
  left_join(select(taxa_table, index, strain), by = "index")

# 4. Plot
ggplot(pcoa_df, aes(x = V1, y = V2, color = strain)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCoA of Arcobacter (Bray-Curtis)",
    x = paste0("PCoA 1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(var_explained[2], 1), "%)"),
    color = "Strain"
  )

#------------------------------------------------------------------------------
# Summarize data per strain (mean) + Bray-Curtis between strains
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

  


