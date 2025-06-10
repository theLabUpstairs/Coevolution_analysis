# Clear workplace
#------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
#------------------------------------------------------------------------------
library(tidyverse)
library(ape)

library(phytools)

# Load data
#------------------------------------------------------------------------------

association_df <- read_csv("data/01_association_df.csv")
symbiont_tree <- read.tree("data/01_symbiont_tree.newick")
symbiont_tree_desulfo <-  read.tree("data/01_symbiont_tree_desulfo.newick")
host_tree <- read.tree("data/01_host_tree.newick")


# Prune trees (arcobacteraea)
#------------------------------------------------------------------------------

# Pull tips to inlcude in the trees
symbiont_names <- association_df %>% pull(Symbiont) %>% unique()
host_names <- association_df %>% pull(Host) %>% unique()

# Make strain names the same as in the association matrix
symbiont_tree$tip.label <- gsub("'", "", symbiont_tree$tip.label)
host_tree$tip.label <- gsub("'", "", host_tree$tip.label)

# Prune tree to only include strains in assocciation matrix
symbiont_tree <- drop.tip(symbiont_tree, setdiff(symbiont_tree$tip.label, symbiont_names))
host_tree <- drop.tip(host_tree, setdiff(host_tree$tip.label, host_names))

# Prune: remove Lenisia and EP1
symbiont_tree <- drop.tip(symbiont_tree, "Arcobacter_EP1_GCA_001655195.1_ASM165519v1_genomic")
host_tree <- drop.tip(host_tree, "nrpr2neighbors_KT023596.1.1885_U|Obazoa|Breviatea|Breviatea_XX|Breviatea_XXX|Breviatea_XXXX|Lenisia|Lenisia_limosa")
association_df <- association_df %>% filter(Symbiont!="Arcobacter_EP1_GCA_001655195.1_ASM165519v1_genomic")


# Rename strains names (arcobacteraea)
#------------------------------------------------------------------------------

# Rename host labels
rename_host <- c(
  "18S_undescribed_breviate_PCE" = "PCE",
  "18S_undescribed_breviate_FB10N2_clone" = "FB10N2",
  "18S_undescribed_breviate_LRM1b" = "LRM1b",
  "LRM2N6_metagenomeRC" = "LRM2N6",
  "nrpr2neighbors_KC433554.1_Pygsuia_biforma_strain_PCbi66_18S_ribosomal_RNA_gene_partial_sequence" = "Pygsuia",
  "nrpr2neighbors_KT023596.1.1885_U|Obazoa|Breviatea|Breviatea_XX|Breviatea_XXX|Breviatea_XXXX|Lenisia|Lenisia_limosa" = "Lenisia_limosa"
)

# Rename symbiont labels
rename_symbiont <- c(
  "1-PCE.1-PCE.Arco_Marinarcus_sp002382325_trycycler_medaka_polypolish" = "Arco_PCE",
  "1-PCE.Halarcobacter_sp004116455_trycycler_medaka_polypolish" = "Halarco_PCE",
  "1-PCE.Halarcobacter_sp003252105_trycycler_medaka_polypolish" = "Halarco2_PCE",
  "2-FB10N2.Malaciobacter_marinus_medaka_polypolish" = "Malacio_FB10N2",
  "2-FB10N2.Halarcobacter_sp003252105_medaka_polypolish" = "Halarco_FB10N2",
  "4-LRM1b.Arcobacter_sp002869535_medaka_polypolish" = "Arco_LRM1b",
  "5-LRM2N6.Halarcobacter_sp004116455_medaka_polypolish" = "Halarco_LRM2N6",
  "5-LRM2N6.Malaciobacter_marinus_medaka_polypolish" = "Malacio_LRM2N6",
  "7-Pygsuia.Arcobacter_medaka" = "Arco_Pygsuia",
  "Arcobacter_EP1_GCA_001655195.1_ASM165519v1_genomic" = "Arco_EP1"
)

host_tree$tip.label <- rename_host[host_tree$tip.label]
symbiont_tree$tip.label <- rename_symbiont[symbiont_tree$tip.label]
association_df <- association_df %>%
  mutate(
    Symbiont = recode(Symbiont, !!!rename_symbiont),
    Host = recode(Host, !!!rename_host)
  )


# Transform into distance matrices and binary association matrix (arcobacteraea)
#------------------------------------------------------------------------------
symbiont_dist <- cophenetic(symbiont_tree)
host_dist <- cophenetic(host_tree)
assoc_matrix <- as.matrix(table(association_df$Symbiont, association_df$Host))



# Process the Desulfo data
#------------------------------------------------------------------------------

# Extract names of tree tips that includes "desulfo"
desulfo_labels <- grep("desulfo", symbiont_tree_desulfo$tip.label, ignore.case = TRUE, value = TRUE)
# Prune tree to only inlcude desulfo tips
symbiont_tree_desulfo <- drop.tip(symbiont_tree_desulfo, setdiff(symbiont_tree_desulfo$tip.label, desulfo_labels))
# Prune tree to remove "DesulfovibrioDMSS1_2576861818", lack information about breviate strain
symbiont_tree_desulfo <- drop.tip(symbiont_tree_desulfo, "DesulfovibrioDMSS1_2576861818")

# Remove "'" flanking the strain name in tree
symbiont_tree_desulfo$tip.label <- gsub("^'|'$", "", symbiont_tree_desulfo$tip.label)

# Rename tips 
rename_desulfo <- c(
  "7-Pygsuia.11108931849_Desulfovibrio7M6.assembly" = "Desulfo7M6_Pygsuia",
  "1-PCE.Desulfovibrionaceae_trycycler_medaka_polypolish" = "Desulfovibrionaceae_PCE",
  "1-PCE.Maridesulfovibrio_sp_trycycler_medaka_polypolish" = "Maridesulfovibrio_PCE",
  "5-LRM2N6.11108969682_Maridesulfovibrio5SRBS1.assembly" = "Maridesulfovibrio5SRBS1_LRM2N6",
  "5-LRM2N6.Maridesulfovibrio_medaka_polypolish" = "Maridesulfovibrio_LRM2N6",
  "2-FB10N2.11108969682_Maridesulfovibrio2SRBS4_assembly" = "Maridesulfovibrio2SRBS4_FB10N2",
  "7-Pygsuia.11108969682_Desulfovibrio7SRBS1_assembly" = "Desulfo7SRBS1_Pygsuia")
symbiont_tree_desulfo$tip.label <- rename_desulfo[symbiont_tree_desulfo$tip.label]

# Bulid an association matrix
# Strains to use in the matrix (LRM1b lacks desulfovibrio)
breviate_ids <- c("PCE", "FB10N2", "LRM2N6", "Pygsuia")
desulfo_labels <- symbiont_tree_desulfo$tip.label
# An empty matrix
assoc_matrix_desulfo <- matrix(0, nrow = length(breviate_ids), ncol = length(desulfo_labels),
                       dimnames = list(breviate_ids, desulfo_labels))
# Populate the empty matrix
for (i in seq_along(breviate_ids)) {
  for (j in seq_along(desulfo_labels)) {
    if (grepl(breviate_ids[i], desulfo_labels[j])) {
      assoc_matrix_desulfo[i, j] <- 1
    }
  }
}

# Build phylogeny distance trees 
symbiont_dist_desulfo <- cophenetic(symbiont_tree_desulfo)
# LMR1b is dropped, it lacks desulfo strains
host_tree_desulfo <- drop.tip(host_tree, "LRM1b")
host_dist_desulfo <- cophenetic(host_tree_desulfo)


# Save data
#------------------------------------------------------------------------------

saveRDS(symbiont_dist, "data/02_symbiont_dist.rds")
saveRDS(host_dist, "data/02_host_dist.rds")
saveRDS(assoc_matrix, "data/02_assoc_matrix.rds")
write_csv(association_df, "data/02_assoc_df.csv")
write.tree(host_tree, file = "data/02_host_tree.newick")
write.tree(symbiont_tree, file = "data/02_symbiont_tree.newick")

# Desulfo related
saveRDS(symbiont_dist_desulfo, "data/02_symbiont_dist_desulfo.rds")
saveRDS(host_dist_desulfo, "data/02_host_dist_desulfo.rds")
saveRDS(assoc_matrix_desulfo, "data/02_assoc_matrix_desulfo.rds")
write.tree(host_tree_desulfo, file = "data/02_host_tree_desulfo.newick")
write.tree(symbiont_tree_desulfo, file = "data/02_symbiont_tree_desulfo.newick")