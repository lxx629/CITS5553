############################################################
# Script: 03-fishbase-geo-depth-merge.R
# Purpose:
#   Enrich the team's main table with FishBase information:
#   - stocks(): TempMin / TempMax + geographic extents
#   - species(): DepthRangeShallow / DepthRangeDeep
#   Merge into the main table by scientific name and deduplicate
#   per ncbi_taxon_id by keeping the row with the most non-NA fields.
#
# Inputs  (relative to repo root):
#   - data/merged/ncbi_aphialid_merged.csv     # team's main table (adjust name if needed)
#
# Outputs:
#   - data/raw/fishbase_merged_basic.csv
#   - data/merged/main_dedup_by_ncbi.csv
#   - data/raw/fishbase_raw_extended.csv
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(rfishbase)
})

# Paths
in_dir   <- "data/merged"
out_dir  <- "data/merged"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Main input
file_main_candidates <- c(
  file.path(in_dir, "ncbi_aphialid_merged.csv"), 
  file.path(in_dir, "ncbi_aphiaID_merged.csv"),
  file.path(in_dir, "ncbi_aphiaId_merged.csv")
)
file_main <- file_main_candidates[file.exists(file_main_candidates)][1]
if (is.na(file_main)) {
  stop("Main input not found. Put your main table under data/merged/ and update file_main.")
}

# Read main table
main <- read_csv(file_main, show_col_types = FALSE)

# Ensure the scientific name column is present
if (!"scientificname" %in% names(main)) {
  if ("scientific_name" %in% names(main)) {
    main <- main %>% rename(scientificname = scientific_name)
  } else if ("scientificName" %in% names(main)) {
    main <- main %>% rename(scientificname = scientificName)
  } else {
    stop("Column 'scientificname' not found in main table.")
  }
}

# Build species list for FishBase
my_species <- unique(main$scientificname)
valid_names <- validate_names(my_species)

# quick duplicate check by scientificname
dup_scientificname <- main %>% count(scientificname) %>% filter(n > 1)
if (nrow(dup_scientificname)) {
  message("Note: found duplicated scientific names in main (this is OK, we dedup by ncbi later).")
}

# FishBase: stocks()
fields_st <- c(
  "Species","SpecCode",
  "Northernmost","NorthSouthN","Southermost","NorthSouthS",
  "Westernmost","WestEastW","Easternmost","WestEastE",
  "TempMin","TempMax"
)
st <- rfishbase::stocks(species_list = valid_names, fields = fields_st) %>% as_tibble()

# Clean numeric columns that may come as character
num_cols_st <- intersect(c("TempMin","TempMax",
                           "Northernmost","Southermost",
                           "Westernmost","Easternmost"), names(st))
st <- st %>%
  mutate(across(all_of(num_cols_st), ~ suppressWarnings(as.numeric(.x))))

# FishBase: species() basic depth ranges
fields_sp <- c("Species","SpecCode","DepthRangeShallow","DepthRangeDeep")
sp <- rfishbase::species(species_list = valid_names, fields = fields_sp) %>% as_tibble()

num_cols_sp <- intersect(c("DepthRangeShallow","DepthRangeDeep"), names(sp))
sp <- sp %>%
  mutate(across(all_of(num_cols_sp), ~ suppressWarnings(as.numeric(.x))))

# Merge FishBase (basic)
merged_fishbase <- full_join(sp, st, by = c("Species","SpecCode"))
write_csv(merged_fishbase, file.path(out_dir, "fishbase_merged_basic.csv"))

# Join into main by scientific name
merged_depth_range <- main %>%
  left_join(merged_fishbase, by = c("scientificname" = "Species"))

# Deduplicate per ncbi_taxon_id
# Keep the row that has the highest number of non-NA fields (ties: first by original order)
if (!"ncbi_taxon_id" %in% names(merged_depth_range)) {
  merged_depth_range$ncbi_taxon_id <- NA  # create if absent to avoid error
}

main_dedup <- merged_depth_range %>%
  mutate(non_na_count = rowSums(across(everything(), ~ !is.na(.x)))) %>%
  arrange(ncbi_taxon_id, desc(non_na_count)) %>%
  distinct(ncbi_taxon_id, .keep_all = TRUE) %>%
  select(-non_na_count)

write_csv(main_dedup, file.path(out_dir, "main_dedup_by_ncbi.csv"))

# (Optional) FishBase: species() extended fields
fields_sp_ext <- c(
  "Species","SpecCode","Genus","Author","FBname","DemersPelag","Subfamily",
  "DepthRangeShallow","DepthRangeDeep","DepthRangeComShallow","DepthRangeComDeep",
  "Length","CommonLength","Weight"
)
sp1 <- rfishbase::species(species_list = valid_names, fields = fields_sp_ext) %>% as_tibble()

# Numeric cleanup for lengths/weights (they may be character)
num_cols_sp1 <- intersect(c("Length","CommonLength","Weight",
                            "DepthRangeShallow","DepthRangeDeep",
                            "DepthRangeComShallow","DepthRangeComDeep"), names(sp1))
sp1 <- sp1 %>%
  mutate(across(all_of(num_cols_sp1), ~ suppressWarnings(as.numeric(.x))))

# Join extended species table with stocks as a raw extended FishBase dump
raw_fishbase <- full_join(sp1, st, by = c("Species","SpecCode"))
write_csv(raw_fishbase, file.path(out_dir, "fishbase_raw_extended.csv"))

message("Done. Outputs written under data/merged/.")
