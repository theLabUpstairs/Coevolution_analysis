# Clear workplace
#------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
#------------------------------------------------------------------------------

library(tidyverse)
library(ape)
library(paco)
library(phytools)

# Load data
#------------------------------------------------------------------------------

symbiont_dist <- readRDS("data/02_symbiont_dist.rds")
host_dist <- readRDS("data/02_host_dist.rds")
assoc_matrix <- readRDS("data/02_assoc_matrix.rds")
association_df <- read_csv("data/02_assoc_df.csv")
symbiont_tree <- read.tree("data/02_symbiont_tree.newick")
host_tree <- read.tree("data/02_host_tree.newick")

# Run PACo
#------------------------------------------------------------------------------

# Preparing the input data
paco_data <- prepare_paco_data(H = host_dist, P = symbiont_dist, HP = assoc_matrix)

# Dimensional reduction using PCoA
paco_data <- add_pcoord(paco_data)

# Run PACo
paco_result <- PACo(paco_data, nperm = 1000, method = "r2", symmetric = TRUE)

# Results
# - p-value (0.279)
paco_result$gof$p 
# - Goodness-of-fit, r2 (0.7528783)
paco_result$gof$ss
# - Number of permutations (1000)
paco_result$gof$n


# Create tanglegram
#------------------------------------------------------------------------------

# Get matching  labels in the correct order
matching_hosts <- association_df$Host
matching_symbionts <- association_df$Symbiont

# Root the tree to reseble the tree in the manuscript
host_tree <- root(host_tree, outgroup = "FB10N2", resolve.root = TRUE)
symbiont_tree <- root(symbiont_tree, outgroup = "Arco_PCE", resolve.root = TRUE)

cophylo_obj <- cophylo(host_tree, symbiont_tree, assoc = cbind(matching_hosts, matching_symbionts))

# Save as PNG
png("figures/cophylo_plot.png", width = 2000, height = 2000, res = 300)
plot(cophylo_obj, fsize = 1.2, link.type = "curved", link.lwd = 2, link.col = "blue")
dev.off()

# Save as SVG
svg("figures/cophylo_plot.svg", width = 10, height = 10)
plot(cophylo_obj, fsize = 1.2, link.type = "curved", link.lwd = 2, link.col = "blue")
dev.off()
