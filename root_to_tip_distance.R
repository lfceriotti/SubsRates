#!/usr/bin/env Rscript

# ---------------------------------------------------------------
# Script: root_to_tip_distance.R
# Description:
#   Computes root-to-tip distances for two phylogenetic trees
#   (one representing dS, the other dN), and merges this information
#   with associated omega, dN, and dS values per taxon.
#   Also outputs annotated tree PDFs and a summary table.
#
# Usage:
#   Rscript root_to_tip_distance.R tree_dS.nwk tree_dN.nwk w_values.txt dN_values.txt dS_values.txt output_table.txt
#
# Arguments:
#   1. tree_dS.nwk   - Rooted Newick tree with branch lengths representing dS
#   2. tree_dN.nwk   - Rooted Newick tree with branch lengths representing dN
#   3. w_values.txt  - Table with omega values per taxon (must have column 'Taxa')
#   4. dN_values.txt - Table with dN values per taxon
#   5. dS_values.txt - Table with dS values per taxon
#   6. output_table.txt - Output table with root-to-tip distances and merged values
# ---------------------------------------------------------------

# Load arguments
args <- commandArgs(trailingOnly = TRUE)

tree_dS_file <- args[1]
tree_dN_file <- args[2]
w_file       <- args[3]
dN_file      <- args[4]
dS_file      <- args[5]
output_file  <- args[6]

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(castor)
  library(ape)
})

# Read rooted Newick trees
tree_dS <- read.tree(file = tree_dS_file)
tree_dN <- read.tree(file = tree_dN_file)

# Read omega, dN, and dS tables
w_table  <- read.table(w_file,  sep = "\t", header = TRUE)
dN_table <- read.table(dN_file, sep = "\t", header = TRUE)
dS_table <- read.table(dS_file, sep = "\t", header = TRUE)

# Initialize result table with taxon labels
Root2tip <- data.frame(Taxa = tree_dS$tip.label)

# Compute root-to-tip distances for terminal nodes (tips) only
Root2tip$dN <- get_all_distances_to_root(tree_dN) %>%
  head(length(tree_dN$tip.label)) - tree_dN$edge.length[1]

Root2tip$dS <- get_all_distances_to_root(tree_dS) %>%
  head(length(tree_dS$tip.label)) - tree_dS$edge.length[1]

# Merge additional values by taxon name
Root2tip <- Root2tip %>%
  full_join(w_table,  by = "Taxa") %>%
  full_join(dN_table, by = "Taxa") %>%
  full_join(dS_table, by = "Taxa")

# Generate annotated tree PDFs
pdf("tree_dS.pdf")
plot.phylo(tree_dS, label.offset = 0.0005, no.margin = TRUE)
add.scale.bar(max(Root2tip$dS, na.rm = TRUE) * 0.85, 1, cex = 0.75)
edgelabels(round(tree_dS$edge.length, 3), frame = "none", cex = 0.5, adj = c(0.5, -0.25))
dev.off()

pdf("tree_dN.pdf")
plot.phylo(tree_dN, label.offset = 0.0005, no.margin = TRUE)
add.scale.bar(max(Root2tip$dN, na.rm = TRUE) * 0.85, 1, cex = 0.75)
edgelabels(round(tree_dN$edge.length, 3), frame = "none", cex = 0.5, adj = c(0.5, -0.25))
dev.off()

# Write final result table
write.table(Root2tip, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
