############################################################
# Script: 01_data_merge.R
# Purpose: Merge marine species datasets (NCBI + GOAT list +
#          WoRMS marine vertebrates), then clean IDs/names and
#          pick one representative row per taxon using release
#          date/assembly presence.
# Inputs  (relative to repo root):
#   - data/raw/marine/supp-table-1.csv           # NCBI assemblies
#   - data/raw/marine/OG_species_goat.csv        # GOAT species list (has '#' header lines)
#   - data/raw/marine/marine-vert-species.csv    # WoRMS marine vertebrates
# Outputs:
#   - data/processed/merged_final.csv
#   - data/processed/species_clean.csv
# Author: Ziheng Wang
# Date  : 2025-08-20
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(lubridate)
})

# paths
in_dir  <- "data/raw"
out_dir <- "data/merged"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# read sources
supp1 <- readr::read_csv(file.path(in_dir, "supp-table-1.csv"))
goat  <- readr::read_csv(file.path(in_dir, "OG_species_goat.csv"), comment = "#")
vert  <- readr::read_csv(file.path(in_dir, "marine-vert-species.csv"))

# rename to normalized keys
# From NCBI assemblies
supp1 <- supp1 %>%
  rename(
    scientific_name = `Organism Name`,
    ncbi_taxon_id   = `Organism Taxonomic ID`
  )

# From GOAT list
goat <- goat %>%
  rename(
    scientific_name = species,
    ncbi_taxon_id   = ncbi_taxon_id
  )

# From WoRMS marine vertebrates
vert <- vert %>%
  rename(
    scientific_name = scientificName
  )

# Ensure IDs are character
supp1 <- supp1 %>% mutate(ncbi_taxon_id = as.character(ncbi_taxon_id))
goat  <- goat  %>% mutate(ncbi_taxon_id = as.character(ncbi_taxon_id))

# Step 1: join NCBI (supp1) + GOAT by ncbi_taxon_id
merged1 <- full_join(
  supp1, goat, by = "ncbi_taxon_id",
  suffix = c("_supp", "_goat") # to coalesce later
  # relationship = "many-to-many"  # uncomment for dplyr >=1.1 to silence warning
)

# Coalesce duplicated columns (scientific_name, family)
merged1 <- merged1 %>%
  mutate(
    scientific_name = coalesce(scientific_name_supp, scientific_name_goat),
    family_name     = coalesce(.data$family_supp, .data$family_goat)
  ) %>%
  select(-scientific_name_supp, -scientific_name_goat, -family_supp, -family_goat)

# Step 2: join WoRMS vertebrates by scientific_name
merged_final <- full_join(
  merged1, vert,
  by = "scientific_name",
  relationship = "many-to-many" 
)

# Coalesce taxonomy dup columns created by the join (class/order/genus)
merge_coalesce <- function(df, col){
  x <- paste0(col, ".x"); y <- paste0(col, ".y")
  if (all(c(x,y) %in% names(df))) {
    df <- df %>% mutate(!!col := coalesce(.data[[x]], .data[[y]])) %>%
      select(-all_of(c(x,y)))
  }
  df
}
for (cname in c("class","order","genus")) {
  merged_final <- merge_coalesce(merged_final, cname)
}

# Step 3: clean IDs and standardize scientific names
clean_id <- function(x){
  x <- as.character(x)
  x[x %in% c("-", "", "NA", "na")] <- NA_character_
  x
}
merged_final <- merged_final %>%
  mutate(
    ncbi_taxon_id = clean_id(ncbi_taxon_id),
    AphiaID       = clean_id(AphiaID),
    sci_name_std  = str_trim(str_replace_all(str_to_lower(scientific_name), "\\s+", " "))
  )

# Build a fallback key: prefer NCBI, else AphiaID, else name
merged_final <- merged_final %>%
  mutate(
    taxon_key  = coalesce(ncbi_taxon_id, AphiaID, sci_name_std),
    key_source = case_when(
      !is.na(ncbi_taxon_id)              ~ "NCBI",
      is.na(ncbi_taxon_id) & !is.na(AphiaID) ~ "AphiaID",
      TRUE                               ~ "sci_name"
    )
  )

# Step 4: pick one representative row per taxon_key
# sort by (has assembly first) then by Assembly Release Date (latest),
# then take the first row per key.
species_clean <- merged_final %>%
  mutate(
    # parse release date safely
    rel_date = suppressWarnings(lubridate::ymd(`Assembly Release Date`)),
    has_assembly = !is.na(`Assembly Accession`)
  ) %>%
  mutate(rel_date = coalesce(rel_date, as.Date("1900-01-01"))) %>%
  arrange(desc(has_assembly), desc(rel_date)) %>% # assemblies first, newest first
  group_by(taxon_key) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-rel_date)

# write outputs
readr::write_csv(merged_final,  file.path(out_dir, "merged_final.csv"))
readr::write_csv(species_clean, file.path(out_dir, "species_clean.csv"))

message("Wrote: ",
        file.path(out_dir, "merged_final.csv"), " (", nrow(merged_final), " rows); ",
        file.path(out_dir, "species_clean.csv"), " (", nrow(species_clean), " rows)")
