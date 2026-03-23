############################################################
# Burgwald: Physio-Stratifizierung + Stationskandidaten
############################################################

# ---------------------------------------------------------
# 0) Pakete
# ---------------------------------------------------------
pkgs <- c(
  "here", "fs", "terra", "sf", "dplyr", "tidyr", "tibble",
  "exactextractr"
)

new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)

invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------------------------------------
# 1) Pfade
# ---------------------------------------------------------
aoi_file      <- here::here("data", "processed", "aoi_burgwald.gpkg")
segments_file <- here::here("data", "processed", "burgwald_segments.gpkg")
dem_file      <- here::here("data", "processed", "dem_dgm1_burgwald.tif")
slope_file    <- here::here("data", "processed", "slope_burgwald.tif")
aspect_file   <- here::here("data", "processed", "aspect_burgwald.tif")

out_dir <- here::here("data", "processed", "station_design")
fs::dir_create(out_dir)

attrstack_file      <- file.path(out_dir, "layer0_segments_physio_metrics.gpkg")
station_physio_file <- file.path(out_dir, "station_physio_candidates.gpkg")

# ---------------------------------------------------------
# 2) Daten laden
# ---------------------------------------------------------
aoi_burgwald <- sf::read_sf(aoi_file)
segments_sf  <- sf::read_sf(segments_file)

dem_burgwald <- terra::rast(dem_file)
slope        <- terra::rast(slope_file)
aspect       <- terra::rast(aspect_file)

# ---------------------------------------------------------
# 3) CRS harmonisieren
# ---------------------------------------------------------
crs_ref <- terra::crs(dem_burgwald)

aoi_burgwald <- sf::st_transform(aoi_burgwald, crs_ref)
segments_sf  <- sf::st_transform(segments_sf, crs_ref)

aoi_v <- terra::vect(aoi_burgwald)

dem_burgwald <- terra::mask(terra::crop(dem_burgwald, aoi_v), aoi_v)
slope        <- terra::mask(terra::crop(slope, aoi_v), aoi_v)
aspect       <- terra::mask(terra::crop(aspect, aoi_v), aoi_v)

# ---------------------------------------------------------
# 4) segment_id
# ---------------------------------------------------------
if (!"segment_id" %in% names(segments_sf)) {
  segments_sf$segment_id <- seq_len(nrow(segments_sf))
}
segments_sf$segment_id <- as.integer(segments_sf$segment_id)

# ---------------------------------------------------------
# 5) Physio-Metriken
# ---------------------------------------------------------
southness <- cos((aspect * pi) / 180 - pi)
mean_wind_dir <- 270
windwardness_r <- cos(((aspect - mean_wind_dir) * pi) / 180)

physio_tbl <- tibble(
  segment_id        = segments_sf$segment_id,
  elev_mean         = exactextractr::exact_extract(dem_burgwald, segments_sf, "mean"),
  slope_mean_deg    = exactextractr::exact_extract(slope, segments_sf, "mean"),
  aspect_mean_deg   = exactextractr::exact_extract(aspect, segments_sf, "mean"),
  southness_mean    = exactextractr::exact_extract(southness, segments_sf, "mean"),
  windwardness_mean = exactextractr::exact_extract(windwardness_r, segments_sf, "mean")
)

# ---------------------------------------------------------
# 6) Attrstack
# ---------------------------------------------------------
segments_attr <- segments_sf |>
  left_join(physio_tbl, by = "segment_id")

sf::st_write(segments_attr, attrstack_file, delete_dsn = TRUE, quiet = TRUE)

# ---------------------------------------------------------
# 7) Stratifizierung
# ---------------------------------------------------------
physio_df <- segments_attr |>
  st_drop_geometry() |>
  select(segment_id, elev_mean, slope_mean_deg, southness_mean, windwardness_mean) |>
  tidyr::drop_na() |>
  mutate(
    elev_z      = scale(elev_mean)[,1],
    slope_z     = scale(slope_mean_deg)[,1],
    southness_z = scale(southness_mean)[,1],
    windward_z  = scale(windwardness_mean)[,1]
  )

set.seed(123)
k_physio <- 5

km <- kmeans(
  physio_df[, c("elev_z", "slope_z", "southness_z", "windward_z")],
  centers = k_physio,
  nstart = 50
)

physio_df$physio_stratum <- km$cluster

# Distanz zum Clusterzentrum
centroids <- physio_df |>
  group_by(physio_stratum) |>
  summarise(across(ends_with("_z"), mean), .groups = "drop")

physio_df <- physio_df |>
  left_join(centroids, by = "physio_stratum", suffix = c("", "_c")) |>
  mutate(
    dist_to_center = sqrt(
      (elev_z - elev_z_c)^2 +
        (slope_z - slope_z_c)^2 +
        (southness_z - southness_z_c)^2 +
        (windward_z - windward_z_c)^2
    )
  )

# ---------------------------------------------------------
# 8) Repräsentative Segmente (FIX!)
# ---------------------------------------------------------
target_n_total <- 10
n_per_stratum <- max(1, ceiling(target_n_total / k_physio))

rep_segments <- physio_df |>
  group_by(physio_stratum) |>
  slice_min(dist_to_center, n = n_per_stratum, with_ties = FALSE) |>
  ungroup() |>
  select(segment_id, physio_stratum)

# WICHTIG: Join sauber
rep_segments_sf <- segments_attr |>
  inner_join(rep_segments, by = "segment_id")

# DEBUG
print(names(rep_segments_sf))

# ---------------------------------------------------------
# 9) Punkte
# ---------------------------------------------------------
station_physio_candidates <- rep_segments_sf |>
  st_point_on_surface() |>
  mutate(station_id = paste0("PHYS_", row_number())) |>
  select(
    station_id,
    segment_id,
    physio_stratum,
    elev_mean,
    slope_mean_deg,
    aspect_mean_deg,
    southness_mean,
    windwardness_mean
  )

# ---------------------------------------------------------
# 10) Speichern
# ---------------------------------------------------------
sf::st_write(
  station_physio_candidates,
  station_physio_file,
  delete_dsn = TRUE,
  quiet = TRUE
)

# ---------------------------------------------------------
# 11) Output
# ---------------------------------------------------------
print(
  station_physio_candidates |>
    st_drop_geometry()
)

