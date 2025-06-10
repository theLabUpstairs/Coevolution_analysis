# Clear workplace
#------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
#------------------------------------------------------------------------------
library(tidyverse)
library(ape)

# Load data
#------------------------------------------------------------------------------

# Loading breviate - arcobacter association matrix
association_file <- "data/_raw/parasite_host_associations_filtered_using_option_0.txt"
association_df <- read_table(association_file, col_names = FALSE) %>% 
  filter(!X1 == "Parasite") %>% 
  set_names(c("Symbiont", "Host"))

# Loading tree files
symbiont_tree_file <- "data/_raw/de_novo_wf_outgroup_f__Campylobacteraceae_gtdbtk.bac120.decorated-itol_withTax.newick"
symbiont_tree_file_dsulfo <- "data/_raw/de_novo_wf_outgroup_f__Desulfohalobiaceae_gtdbtk.bac120.decorated-itol.tree"
host_tree_file <- "data/_raw/SSU150_SSU-align_regions_aligned.MFPbs.treefile.newick"

symbiont_tree <- read.tree(symbiont_tree_file)
symbiont_tree_desulfo <- read.tree(symbiont_tree_file_dsulfo)
host_tree <- read.tree(host_tree_file)


# Write data
#------------------------------------------------------------------------------

write_csv(association_df,"data/01_association_df.csv")
write.tree(symbiont_tree, file = "data/01_symbiont_tree.newick")
write.tree(symbiont_tree_desulfo, file = "data/01_symbiont_tree_desulfo.newick")
write.tree(host_tree, file = "data/01_host_tree.newick")