# Clear workplace
#------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
#------------------------------------------------------------------------------
library(tidyverse)

# Load data
#------------------------------------------------------------------------------

df <- read_csv("data/01_loaded.csv")

# Wrangle and clean data
#------------------------------------------------------------------------------

# Table to long format, to make taxa one column.
df <- df %>%
  pivot_longer(cols = starts_with("d__"),
               names_to = "taxa",
               values_to = "abundance")

# Add columns for each taxonomic level
df <- df %>%
  mutate(domain = str_extract(taxa, "(?<=d__)[^\\.]+"),
         phylum = str_extract(taxa, "(?<=p__)[^\\.]+"),
         class = str_extract(taxa, "(?<=c__)[^\\.]+"),
         order = str_extract(taxa, "(?<=o__)[^\\.]+"),
         family = str_extract(taxa, "(?<=f__)[^\\.]+"),
         genus = str_extract(taxa, "(?<=g__)[^\\.]+"))

# Make sample ID easiser to read
#df <- df %>% 
#  mutate(index = str_replace(index, "NG-32783-", "")) %>% 
#  mutate(index = str_replace(index, "NG-34656-", ""))

# Remove some columns
df <- df %>% select(-ef.barcode, -barcode.sequence, -date.dna.isolation.2,
                    -day.dna.isolation.2, -date.pcr, -date.pcr.purification,
                    -incubation.days, -vol.pcr.product.purified,
                    -vol.eluted.pcr.purification, -date.send.sequencing,
                    -vol.send.sequencing, -taxa)

# Add Family name to missing genus
df <- df %>%
  mutate(genus = ifelse(is.na(genus), paste(family, "(family)"), genus)) %>%
  mutate(genus = ifelse(genus == "P", paste(class, "(class)"), genus))

# Correct excel auto complete mistake
df <- df %>%
  mutate(strain = str_replace(strain, "^LRM2N\\d+", "LRM2N6"))

# Removing the n2o samples
df <- df %>% filter(culture.kno3 != "n2o")

# Write data
#------------------------------------------------------------------------------

df %>% write_csv("data/02_wrangle.csv")
