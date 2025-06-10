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

symbiont_dist_desulfo <- readRDS("data/02_symbiont_dist_desulfo.rds")
host_dist_desulfo <- readRDS("data/02_host_dist_desulfo.rds")
assoc_matrix_desulfo <- readRDS("data/02_assoc_matrix_desulfo.rds")
host_tree_desulfo <- read.tree("data/02_host_tree_desulfo.newick")
symbiont_tree_desulfo <- read.tree("data/02_symbiont_tree_desulfo.newick")

# Script settings
#------------------------------------------------------------------------------

# set seed
set.seed(42)

png_width <- 3000
png_hight <- 2000
png_res <- 300

# Run PACo (Arcobacteraea)
#------------------------------------------------------------------------------

# Preparing the input data
paco_data <- prepare_paco_data(H = host_dist, P = symbiont_dist, HP = assoc_matrix)

# Dimensional reduction using PCoA
paco_data <- add_pcoord(paco_data)

# Run PACo
paco_result <- PACo(paco_data, nperm = 1000, method = "r2", symmetric = TRUE)

# Results
# - p-value (0.266)
paco_result$gof$p 
# - Goodness-of-fit, r2 (0.7528783)
paco_result$gof$ss


# Create tanglegram (Arcobacteraea)
#------------------------------------------------------------------------------

# Get matching  labels in the correct order
matching_hosts <- association_df$Host
matching_symbionts <- association_df$Symbiont

# Root the tree to reseble the tree in the manuscript
host_tree <- root(host_tree, outgroup = "FB10N2", resolve.root = TRUE)
symbiont_tree <- root(symbiont_tree, outgroup = "Arco_PCE", resolve.root = TRUE)

cophylo_obj <- cophylo(host_tree, symbiont_tree, assoc = cbind(matching_hosts, matching_symbionts))

# Save as PNG
png("figures/cophylo_plot_arcobacteraea.png", width = png_width, height = png_hight, res = png_res)
plot(cophylo_obj, fsize = 1.2, link.type = "curved", link.lwd = 2, link.col = "blue")
dev.off()

# Save as SVG
#svg("figures/cophylo_plot.svg", width = 10, height = 10)
#plot(cophylo_obj, fsize = 1.2, link.type = "curved", link.lwd = 2, link.col = "blue")
#dev.off()


# Run PACo (Desulfovibrionaceae)
#------------------------------------------------------------------------------

# Preparing the input data
paco_data <- prepare_paco_data(H = host_dist_desulfo, P = symbiont_dist_desulfo, HP = assoc_matrix_desulfo)
# Dimensional reduction using PCoA
paco_data <- add_pcoord(paco_data)
# Run PACo
paco_result <- PACo(paco_data, nperm = 1000, method = "r2", symmetric = TRUE)

# Results
# - p-value (0.283)
paco_result$gof$p 
# - Goodness-of-fit, r2 (0.7273486)
paco_result$gof$ss


# Create tanglegram (Desulfovibrionaceae)
#------------------------------------------------------------------------------

# Create an association dataframe for desulfo
association_df <- which(assoc_matrix_desulfo == 1, arr.ind = TRUE)
association_df <- data.frame(
  Host = rownames(assoc_matrix_desulfo)[association_df[, 1]],
  Symbiont = colnames(assoc_matrix_desulfo)[association_df[, 2]]
)

# Get matching  labels in the correct order
matching_hosts <- association_df$Host
matching_symbionts <- association_df$Symbiont

# Root trees (if not already)
host_tree <- root(host_tree, outgroup = "FB10N2", resolve.root = TRUE)

# Generate cophylogeny
cophylo_obj <- cophylo(host_tree, symbiont_tree_desulfo, assoc = cbind(matching_hosts, matching_symbionts))

# Save as PNG
png("figures/cophylo_plot_desulfovibrionaceae.png", width = png_width, height = png_hight, res = png_res)
plot(cophylo_obj, fsize = 1.2 ,link.type = "curved", link.lwd = 2, link.col = "blue")
dev.off()