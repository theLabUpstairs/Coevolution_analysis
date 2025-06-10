# PACo Coevolution Analysis

## Description

This directory contains an R-based workflow for assessing potential
coevolutionary patterns between breviates and associated prokaryotes using the
Procrustes Approach to Cophylogeny (PACo).  Analyses include both
Arcobacteraceae and Desulfovibrionaceae subsets, using phylogenetic distance
matrices and host–symbiont association data.

## Contents

- `01_load.R`: Loads all required data (trees, distance matrices, associations).
- `02_wrangle.R`: Processes and filters data for PACo analysis.
- `03_analysis.R`: Performs PACo analysis and generates tanglegrams.
- `Co-evolution.Rproj`: RStudio project file.
- `data/`: Contains input trees, distance matrices, and association matrices.
  - `_raw/`: Raw phylogenetic trees and distance files from external tools.
- `figures/`: Output folder for tanglegrams (Arcobacteraceae and Desulfovibrionaceae).

## Required packages
- `tidyverse`
- `ape`
- `phytools`
- `paco`

## How to run

1. Open the project in RStudio via `Co-evolution.Rproj`.
2. Run the scripts in order:
   - `01_load.R`
   - `02_wrangle.R`
   - `03_analysis.R`

Figures and results will be saved in the `figures/` directory.

