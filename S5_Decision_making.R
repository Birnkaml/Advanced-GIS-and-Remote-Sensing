

############################################################
# Burgwald: Physio-basierte Messpunktkandidaten für Luv/Lee
# Angepasst an einfache Projektstruktur ohne outputs.tsv
############################################################

suppressPackageStartupMessages({
  library(here)
  library(fs)
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(exactextractr)
  library(cluster)
})

# ---------------------------------------------------------
# 0) Ordner / Pfade
# ---------------------------------------------------------
message("Project root: ", here::here())

aoi_file      <- here::here("data", "processed", "aoi_burgwald.gpkg")
segments_file <- here::here("data", "processed", "burgwald_segments.gpkg")
dem_file      <- here::here("data", "processed", "dem_dgm1_burgwald.tif")
slope_file    <- here::here("data", "processed", "slope_burgwald.tif")
aspect_file   <- here::here("data", "processed", "aspect_burgwald.tif")

# Optional: falls vorhanden, Kandidaten auf Hull beschränken
hull_file     <- here::here("data", "processed", "gauges_hull.gpkg")

out_dir <- here::here("data", "processed", "station_design")
fs::dir_create(out_dir)

physio_seg_out <- file.path(out_dir, "physio_strata_segments.gpkg")
physio_pts_out <- file.path(out_dir, "physio_candidates_pts.gpkg")

required_files <- c(aoi_file, segments_file, dem_file, slope_file, aspect_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Diese Eingabedateien fehlen:\n", paste(missing_files, collapse = "\n"))
}

# ---------------------------------------------------------
# 1) Parameter
# ---------------------------------------------------------
set.seed(42)

# Zahl der Strata
k_physio <- 5L

# Kandidaten pro Stratum vor Spacing
n_per_stratum <- 4L

# Mindestabstand zwischen finalen Punkten (m)
d_min_m <- 500

# Hauptwindrichtung für Luv/Lee
mean_wind_dir <- 270

# Optionaler Filter auf Segmentebene:
# Beispiele:
#   NULL
#   "elev_mean > 300"
#   "windwardness_mean > 0"
filter_expr <- NULL

# ---------------------------------------------------------
# 2) Daten laden
# ---------------------------------------------------------
aoi_sf      <- sf::read_sf(aoi_file)
segments_sf <- sf::read_sf(segments_file)

dem    <- terra::rast(dem_file)
slope  <- terra::rast(slope_file)
aspect <- terra::rast(aspect_file)

# segment_id sicherstellen
if (!"segment_id" %in% names(segments_sf)) {
  segments_sf$segment_id <- seq_len(nrow(segments_sf))
}
segments_sf$segment_id <- as.integer(segments_sf$segment_id)

# ---------------------------------------------------------
# 3) CRS harmonisieren + AOI anwenden
# ---------------------------------------------------------
crs_ref <- terra::crs(dem)

if (!identical(sf::st_crs(aoi_sf)$wkt, crs_ref)) {
  aoi_sf <- sf::st_transform(aoi_sf, crs_ref)
}
if (!identical(sf::st_crs(segments_sf)$wkt, crs_ref)) {
  segments_sf <- sf::st_transform(segments_sf, crs_ref)
}

aoi_v <- terra::vect(aoi_sf)

dem <- dem |> terra::crop(aoi_v) |> terra::mask(aoi_v)
slope <- slope |> terra::crop(aoi_v) |> terra::mask(aoi_v)
aspect <- aspect |> terra::crop(aoi_v) |> terra::mask(aoi_v)

# ---------------------------------------------------------
# 4) Optionaler Hull-Filter
# ---------------------------------------------------------
if (file.exists(hull_file)) {
  hull <- sf::read_sf(hull_file)
  hull <- sf::st_make_valid(hull)
  hull <- sf::st_union(hull)
  
  if (sf::st_crs(segments_sf) != sf::st_crs(hull)) {
    hull <- sf::st_transform(hull, sf::st_crs(segments_sf))
  }
  
  segments_sf <- segments_sf[sf::st_within(segments_sf, hull, sparse = FALSE), ]
  message("Segmente innerhalb Hull: ", nrow(segments_sf))
} else {
  message("Kein gauges_hull.gpkg gefunden -> arbeite auf allen Segmenten.")
}

stopifnot(nrow(segments_sf) > 0)

# ---------------------------------------------------------
# 5) Physio-Metriken pro Segment
# ---------------------------------------------------------
southness <- cos((aspect * pi) / 180 - pi)
windwardness_r <- cos(((aspect - mean_wind_dir) * pi) / 180)

physio_tbl <- tibble::tibble(
  segment_id        = segments_sf$segment_id,
  elev_mean         = as.numeric(exactextractr::exact_extract(dem, segments_sf, "mean")),
  slope_mean_deg    = as.numeric(exactextractr::exact_extract(slope, segments_sf, "mean")),
  aspect_mean_deg   = as.numeric(exactextractr::exact_extract(aspect, segments_sf, "mean")),
  southness_mean    = as.numeric(exactextractr::exact_extract(southness, segments_sf, "mean")),
  windwardness_mean = as.numeric(exactextractr::exact_extract(windwardness_r, segments_sf, "mean"))
)

segments_attr <- segments_sf |>
  dplyr::left_join(physio_tbl, by = "segment_id")

# ---------------------------------------------------------
# 6) Optionaler Hard Filter (aus Selector-Idee)
# ---------------------------------------------------------
apply_filter_expr <- function(x_sf, filter_expr) {
  if (is.null(filter_expr) || !nzchar(filter_expr)) return(x_sf)
  x_tbl <- x_sf |> sf::st_drop_geometry()
  keep <- eval(parse(text = filter_expr), envir = x_tbl, enclos = baseenv())
  if (!is.logical(keep) || length(keep) != nrow(x_tbl)) {
    stop("filter_expr muss ein logischer Vektor der Länge nrow(data) sein.")
  }
  x_sf[which(keep), ]
}

segments_attr <- apply_filter_expr(segments_attr, filter_expr)
stopifnot(nrow(segments_attr) > 0)

# ---------------------------------------------------------
# 7) Physio-Stratifizierung
# ---------------------------------------------------------
physio_df <- segments_attr |>
  sf::st_drop_geometry() |>
  dplyr::select(
    segment_id,
    elev_mean,
    slope_mean_deg,
    southness_mean,
    windwardness_mean
  ) |>
  tidyr::drop_na()

stopifnot(nrow(physio_df) > 0)

X <- physio_df |>
  dplyr::select(elev_mean, slope_mean_deg, southness_mean, windwardness_mean) |>
  as.matrix()

Xz <- scale(X)

km <- kmeans(Xz, centers = k_physio, nstart = 50)

physio_df$physio_stratum <- as.integer(km$cluster)

# Distanz zum jeweiligen Clusterzentrum
dist_to_center <- numeric(nrow(Xz))
for (i in seq_len(nrow(Xz))) {
  c_id <- km$cluster[i]
  dist_to_center[i] <- sqrt(sum((Xz[i, ] - km$centers[c_id, ])^2))
}
physio_df$physio_dist_to_center <- dist_to_center

segments_attr <- segments_attr |>
  dplyr::left_join(
    physio_df |>
      dplyr::select(segment_id, physio_stratum, physio_dist_to_center),
    by = "segment_id"
  )

# ---------------------------------------------------------
# 8) Kandidaten pro Stratum: repräsentative Segmente
# ---------------------------------------------------------

cand_ids <- physio_df |>
  dplyr::group_by(physio_stratum) |>
  dplyr::slice_min(physio_dist_to_center, n = n_per_stratum, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(segment_id, physio_stratum, physio_dist_to_center)

cand_seg <- segments_attr |>
  dplyr::select(-any_of(c("physio_stratum", "physio_dist_to_center"))) |>
  dplyr::inner_join(cand_ids, by = "segment_id")

cand_pts <- sf::st_point_on_surface(cand_seg) |>
  dplyr::select(
    segment_id,
    physio_stratum,
    physio_dist_to_center,
    elev_mean,
    slope_mean_deg,
    aspect_mean_deg,
    southness_mean,
    windwardness_mean
  )

# ---------------------------------------------------------
# 9) Greedy Spacing (vereinfacht aus Selector-Idee)
# ---------------------------------------------------------
thin_points_greedy_simple <- function(pts_sf, d_min_m = 500) {
  if (nrow(pts_sf) <= 1) return(pts_sf)
  
  pts_proj <- pts_sf
  if (sf::st_is_longlat(pts_proj)) {
    pts_proj <- sf::st_transform(pts_proj, 25832)
  }
  
  keep <- rep(FALSE, nrow(pts_proj))
  chosen <- integer(0)
  
  ord <- order(pts_proj$physio_dist_to_center, decreasing = FALSE)
  
  for (i in ord) {
    if (length(chosen) == 0) {
      keep[i] <- TRUE
      chosen <- c(chosen, i)
    } else {
      d <- as.numeric(sf::st_distance(pts_proj[i, ], pts_proj[chosen, ], by_element = FALSE))
      if (all(d >= d_min_m)) {
        keep[i] <- TRUE
        chosen <- c(chosen, i)
      }
    }
  }
  
  pts_sf[keep, ]
}

cand_pts_thin <- thin_points_greedy_simple(cand_pts, d_min_m = d_min_m)

# Segmente passend zu finalen Punkten
cand_seg_thin <- cand_seg |>
  dplyr::semi_join(cand_pts_thin |> sf::st_drop_geometry(), by = "segment_id")

# ---------------------------------------------------------
# 10) Speichern
# ---------------------------------------------------------
if (file.exists(physio_seg_out)) try(sf::st_delete_dsn(physio_seg_out), silent = TRUE)
if (file.exists(physio_pts_out)) try(sf::st_delete_dsn(physio_pts_out), silent = TRUE)

sf::st_write(cand_seg_thin, physio_seg_out, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(cand_pts_thin, physio_pts_out, delete_dsn = TRUE, quiet = TRUE)

# ---------------------------------------------------------
# 11) Ausgabe
# ---------------------------------------------------------
message("Fertig.")
message("Strata-Segmente: ", physio_seg_out)
message("Kandidaten-Punkte: ", physio_pts_out)

print(
  cand_pts_thin |>
    sf::st_drop_geometry() |>
    dplyr::arrange(physio_stratum, physio_dist_to_center)
)

