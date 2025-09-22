############################################################
# Script: 04_merge_group_fishbase.R
# Purpose:
#   Merge teammate dataset with FishBase stocks (TempMin/TempMax),
#   resolve duplicates (prefer fewer NAs), and append the traits
#   into the base table.
#
# Inputs (CSV under data/merged/):
#   - new_aphia_ncbi_all_18_09_2025.csv   # teammate dataset
#   - non_conflicted_species_final.csv    # teammate dataset
#
# Output:
#   - data/merged/new_final_species.csv   # teammate's + TempMin/TempMax merged
#
# Optional cache:
#   - data/merged/fishbase_stocks_basic.csv
#
# Notes:
#   - Uses rfishbase::validate_names() + stocks(fields = TempMin/TempMax).
#   - Join order:
#       1) Add TempMin/TempMax to teammate table by scientific name
#       2) Fallback join into base table by ncbi OR aphia (key prefers ncbi)
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(rfishbase)
})

# paths
in_dir   <- "data/merged"
out_dir  <- "data/merged"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

file_group <- file.path(in_dir, "new_aphia_ncbi_all_18_09_2025.csv")
file_base  <- file.path(in_dir, "non_conflicted_species_final.csv")
file_cache <- file.path(in_dir, "fishbase_stocks_basic.csv")
file_out   <- file.path(out_dir, "new_final_species.csv")

# safety checks
stopifnot(file.exists(file_group))
stopifnot(file.exists(file_base))

# read inputs
group_df <- read_csv(file_group, show_col_types = FALSE)
base_df  <- read_csv(file_base,  show_col_types = FALSE)

# Standardise column names we rely on
# Expect: scientificname, ncbi_taxon_id, aphia_id
if (!"scientificname" %in% names(group_df) && "scientific_name" %in% names(group_df)) {
  group_df <- group_df %>% rename(scientificname = scientific_name)
}
# Ensure IDs exist (create if absent to avoid errors)
if (!"ncbi_taxon_id" %in% names(group_df)) group_df$ncbi_taxon_id <- NA
if (!"aphia_id"      %in% names(group_df)) group_df$aphia_id      <- NA
if (!"ncbi_taxon_id" %in% names(base_df))  base_df$ncbi_taxon_id  <- NA
if (!"aphia_id"      %in% names(base_df))  base_df$aphia_id       <- NA

# FishBase: validate names + stocks (TempMin/TempMax)
# Pull unique species names from teammate dataset
my_species2   <- group_df %>% distinct(scientificname) %>% pull(scientificname)

# validate to canonical names (FishBase)
valid_names2  <- rfishbase::validate_names(my_species2)

# Try cached stocks first (to avoid repeated API calls)
if (file.exists(file_cache)) {
  st2 <- read_csv(file_cache, show_col_types = FALSE)
} else {
  st2 <- rfishbase::stocks(
    species_list = valid_names2,
    fields = c("Species", "SpecCode", "TempMin", "TempMax")
  ) %>% as_tibble()
  # cache raw stocks for reproducibility
  write_csv(st2, file_cache)
}

# deduplicate FishBase rows by Species
# Rule: prefer rows with fewer NAs in TempMin/TempMax; tie-break by original order
st2_clean <- st2 %>%
  mutate(
    .rowid   = row_number(),
    na_count = rowSums(across(c(TempMin, TempMax), ~ is.na(.x)))
  ) %>%
  group_by(Species) %>%
  arrange(na_count, .rowid, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-na_count, -.rowid)

# bring TempMin/TempMax into teammate dataset
merged_temp <- group_df %>%
  left_join(st2_clean, by = c("scientificname" = "Species"))
# merged_temp now has TempMin, TempMax attached

# helper: robust string key for fallback join
to_chr <- function(x) {
  if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE) else as.character(x)
}
make_key <- function(df) {
  df %>%
    mutate(
      ncbi_chr  = to_chr(ncbi_taxon_id),
      aphia_chr = to_chr(aphia_id),
      join_key  = if_else(!is.na(ncbi_taxon_id),
                          paste0("ncbi:", ncbi_chr),
                          paste0("aphia:", aphia_chr))
    )
}

# fallback join into base_df (prefer ncbi, else aphia)
nf <- make_key(base_df)
mt <- make_key(merged_temp)

# keep only one row per join_key from merged_temp, and only the columns we need
mt_slim <- mt %>%
  distinct(join_key, .keep_all = TRUE) %>%
  select(join_key, TempMin, TempMax)

result <- nf %>%
  left_join(mt_slim, by = "join_key") %>%
  select(-join_key, -ncbi_chr, -aphia_chr)

# sanity checks
n_both_na_before <- base_df %>% filter(is.na(aphia_id) & is.na(ncbi_taxon_id)) %>% nrow()
n_both_na_after  <- result  %>% filter(is.na(aphia_id) & is.na(ncbi_taxon_id)) %>% nrow()
message(sprintf("Rows with BOTH aphia_id & ncbi_taxon_id = NA: before=%d, after=%d",
                n_both_na_before, n_both_na_after))

# write output
write_csv(result, file_out)
message("Wrote: ", file_out)
