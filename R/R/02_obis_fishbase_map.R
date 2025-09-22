############################################################
# Script: 02_obis_map.R
# Purpose:
#   OBIS: pick 1 family, take first N species, fetch occurrences
#             (presence-only, valid coords, coord uncertainty < 10km),
#             drop points on land, and plot global maps (hexbin & points).
#
# Inputs:
#   - data/merged/species_clean.csv   # from 01_data_merge.R
#   - data/merged/merged_final.csv    # for FishBase species list
#
# Outputs:
#   - data/merged/obis_occurrences.csv
#   - data/merged/obis_summary.csv
#   - outputs/maps/obis_hex_<family>_since<date>.png
#   - outputs/maps/obis_points_<family>_since<date>.png
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(lubridate)
  library(robis)
  library(sf)
  library(rnaturalearth); library(rnaturalearthdata)
  library(ggplot2)
  library(rfishbase)
})

# paths
in_merged   <- "data/merged"
out_merged  <- "data/merged"
out_maps    <- "outputs/maps"
dir.create(out_merged, showWarnings = FALSE, recursive = TRUE)
dir.create(out_maps,   showWarnings = FALSE, recursive = TRUE)

# OBIS controls
target_family        <- "Serranidae"      # used Serranidae family
max_species         <- 30                
start_date          <- "2000-01-01"       # OBIS time filter
uncertainty_max_m   <- 10000              # coord uncertainty threshold
hex_bins            <- 60

# read species list
species_clean <- read_csv(file.path(in_merged, "species_clean.csv"), show_col_types = FALSE)

# Take first N species within the family
sp_vec <- species_clean %>%
  filter(!is.na(family), family == !!target_family) %>%
  distinct(scientific_name) %>%
  slice_head(n = max_species) %>%
  pull(scientific_name)

if (length(sp_vec) == 0) stop("No species found under the chosen family.")

# OBIS fetch (per species)
fetch_one <- function(sp) {
  message("Fetching: ", sp)
  out <- tryCatch(robis::occurrence(sp, startdate = start_date), error = function(e) NULL)
  if (is.null(out) || !nrow(out)) return(NULL)
  out$scientificName_query <- sp
  out
}

obis_raw <- purrr::map(sp_vec, fetch_one) %>% list_rbind()
if (is.null(obis_raw) || !nrow(obis_raw))
  stop("No OBIS records returned for the chosen family/species list. Try another family or increase max_species.")

# clean: presence, coords, coord uncertainty
obis_filt <- obis_raw %>%
  filter(is.na(occurrenceStatus) | occurrenceStatus == "present") %>%      # presence-only
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%            # valid coords
  filter(decimalLongitude >= -180, decimalLongitude <= 180,
         decimalLatitude  >=  -90, decimalLatitude  <=  90) %>%
  mutate(coordinateUncertaintyInMeters = suppressWarnings(as.numeric(coordinateUncertaintyInMeters))) %>%
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < uncertainty_max_m)

# drop points on land (sf + rnaturalearth)
world  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |> st_make_valid()
pts_sf <- st_as_sf(obis_filt, coords = c("decimalLongitude","decimalLatitude"), crs = 4326, remove = FALSE)
on_land <- lengths(st_intersects(pts_sf, world)) > 0
obis_sea <- pts_sf[!on_land, ]

# write OBIS CSVs (handy cache for grading)
readr::write_csv(st_drop_geometry(obis_sea), file.path(out_merged, "obis_occurrences.csv"))
obis_summary <- st_drop_geometry(obis_sea) %>%
  group_by(scientificName_query) %>%
  summarise(n_records = n(),
            n_years   = n_distinct(year, na.rm = TRUE),
            .groups = "drop")
readr::write_csv(obis_summary, file.path(out_merged, "obis_summary.csv"))

# maps: hexbin + points
hex_title   <- sprintf("Global OBIS Sightings (%s)", target_family)
hex_sub     <- sprintf("Species: up to %d; Since %s", max_species, start_date)
f_hex_png   <- file.path(out_maps, sprintf("obis_hex_%s_since_%s.png", target_family, start_date))
f_pts_png   <- file.path(out_maps, sprintf("obis_points_%s_since_%s.png", target_family, start_date))

p_hex <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.2) +
  stat_bin_hex(
    data = obis_sea |> st_drop_geometry(),
    aes(x = decimalLongitude, y = decimalLatitude, fill = after_stat(count)),
    bins = hex_bins, color = "white", alpha = 0.85
  ) +
  scale_fill_viridis_c(name = "Sightings (count)", trans = "log10") +
  coord_sf(expand = FALSE) +
  labs(title = hex_title, subtitle = hex_sub, x = NULL, y = NULL) +
  theme_void(base_size = 12) + theme(legend.position = "bottom")

p_points <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.2) +
  geom_point(
    data = obis_sea |> st_drop_geometry(),
    aes(x = decimalLongitude, y = decimalLatitude),
    alpha = 0.3, size = 0.6
  ) +
  coord_sf(expand = FALSE) +
  labs(title = sprintf("Global OBIS Points (%s)", target_family),
       subtitle = hex_sub, x = NULL, y = NULL) +
  theme_void(base_size = 12)

ggsave(f_hex_png, p_hex, width = 11, height = 6.5, dpi = 150)
ggsave(f_pts_png, p_points, width = 11, height = 6.5, dpi = 150)

message("Done. OBIS CSVs and maps written.")
