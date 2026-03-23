

############################################################

# Burgwald: Luv-Lee & Niederschlag Analyse

############################################################

# ---------------------------------------------------------

# 1) Pakete

# ---------------------------------------------------------

packages <- c("terra", "sf", "here", "fs", "geodata")

new <- packages[!packages %in% installed.packages()[, "Package"]]

if (length(new)) install.packages(new)

invisible(lapply(packages, library, character.only = TRUE))

pkgs <- c(
  "here","terra","sf","stars","dplyr","tidyr",
  "data.table","gstat","automap","pbmcapply",
  "geodata","RANN"
)

# ---------------------------------------------------------

# 2) Projektordner anlegen

# ---------------------------------------------------------

fs::dir_create(
  c(
    here::here("data/raw"),
    here::here("data/processed"),
    here::here("outputs/figures")
  )
)

message("Project root: ", here::here())

# ---------------------------------------------------------

# 3) Helper-Funktionen

# ---------------------------------------------------------

download_if_missing <- function(url, destfile, mode = "wb") {
  
  if (!file.exists(destfile)) {
    
    message("Downloading: ", basename(destfile))
    
    download.file(url, destfile = destfile, mode = mode)

    
  } else {
    
      message("Already exists: ", basename(destfile))
 
    
  }
  
}

# ---------------------------------------------------------

# 4) AOI Burgwald definieren

# ---------------------------------------------------------

burgwald_bbox <- c(
  xmin = 8.7,
  xmax = 9.0,
  ymin = 50.85,
  ymax = 51.05
)

aoi_burgwald_wgs <- sf::st_as_sfc(
  sf::st_bbox(burgwald_bbox, crs = 4326)
)

aoi_file <- here::here(
  "data",
  "processed",
  "aoi_burgwald.gpkg"
)

sf::st_write(
  sf::st_sf(geometry = aoi_burgwald_wgs),
  aoi_file,
  delete_dsn = TRUE,
  quiet = TRUE
)

# ---------------------------------------------------------

# 5) DGM1 herunterladen und mosaikieren

# ---------------------------------------------------------

dem_out_file <- here::here(
  "data",
  "processed",
  "dem_dgm1_burgwald.tif"
)

 if (!file.exists(dem_out_file)) {
  
  options(timeout = 600)
  
  today <- format(Sys.Date(), "%Y%m%d")
  
  base_url_wf <- paste0(
    "https://gds.hessen.de/downloadcenter/",
    today,
    "/3D-Daten/Digitales%20Gel%C3%A4ndemodell%20(DGM1)/Landkreis%20Waldeck-Frankenberg/"
  )
  
  base_url_mr <- paste0(
    "https://gds.hessen.de/downloadcenter/",
    today,
    "/3D-Daten/Digitales%20Gel%C3%A4ndemodell%20(DGM1)/Landkreis%20Marburg-Biedenkopf/"
  )
  
  dgm1_urls <- c(
    
    burgwald = paste0(base_url_wf, "Burgwald%20-%20DGM1.zip"),
    gemuenden = paste0(base_url_wf, "Gem%C3%BCnden%20(Wohra)%20-%20DGM1.zip"),
    rosenthal = paste0(base_url_wf, "Rosenthal%20-%20DGM1.zip"),
    muenchhausen = paste0(base_url_mr, "M%C3%BCnchhausen%20-%20DGM1.zip"),
    rauschenberg = paste0(base_url_mr, "Rauschenberg%20-%20DGM1.zip"),
    coelbe = paste0(base_url_mr, "C%C3%B6lbe%20-%20DGM1.zip"),
    lahntal = paste0(base_url_mr, "Lahntal%20-%20DGM1.zip"),
    wohra = paste0(base_url_mr, "Wohratal%20-%20DGM1.zip"),
    wetter = paste0(base_url_mr, "Wetter%20(Hessen)%20-%20DGM1.zip"),
    haina = paste0(base_url_wf, "Haina%20(Kloster)%20-%20DGM1.zip")
    
    
  )
  
  all_tif_files <- c()
  
  for (nm in names(dgm1_urls)) {
    
  
    url_i <- dgm1_urls[[nm]]
    
    zip_file_i <- here::here(
      "data",
      "raw",
      paste0("dgm1_", nm, ".zip")
    )
    
    unzipped_dir_i <- here::here(
      "data",
      "raw",
      paste0("dgm1_", nm)
    )
    
    download_if_missing(url_i, zip_file_i, mode = "wb")
    
    if (!fs::dir_exists(unzipped_dir_i)) {
      
      fs::dir_create(unzipped_dir_i)
      
      unzip(
        zip_file_i,
        exdir = unzipped_dir_i
      )
      
    }
    
    tif_i <- dir(
      unzipped_dir_i,
      pattern = "\\.tif$",
      full.names = TRUE,
      recursive = TRUE
    )
    
    all_tif_files <- c(all_tif_files, tif_i)
    
    
  }
  
  if (length(all_tif_files) == 0) {
    

    stop("Keine TIF-Dateien in den DGM1-Zips gefunden.")
  
    
  }
  
  message("Gefundene DGM1-TIFs: ", length(all_tif_files))
  
  # ---------------------------------------------------------
  
  # AOI + 3 km Puffer vorbereiten
  
  # ---------------------------------------------------------
  
  aoi_dem_sf <- sf::st_transform(aoi_burgwald_wgs, 25832)
  
  aoi_buffer_sf <- sf::st_buffer(aoi_dem_sf, 3000)
  
  aoi_buffer_v <- terra::vect(aoi_buffer_sf)
  
  # ---------------------------------------------------------
  
  # DEM laden und direkt beschneiden
  
  # ---------------------------------------------------------
  
  dem_list <- lapply(all_tif_files, function(f) {
    

    r <- terra::rast(f)
    
    terra::crs(r) <- "EPSG:25832"
    
    e_int <- terra::intersect(
      terra::ext(r),
      terra::ext(aoi_buffer_v)
    )
    
    if (is.null(e_int)) {
      
      return(NULL)
      
    }
    
    r_crop <- terra::crop(r, aoi_buffer_v)
    
    return(r_crop)
    
    
  })
  
  dem_list <- Filter(Negate(is.null), dem_list)
  
  # ---------------------------------------------------------
  
  # Mosaic
  
  # ---------------------------------------------------------
  
  dem_merged <- do.call(
    terra::mosaic,
    c(dem_list, fun = "min")
  )
  
  # ---------------------------------------------------------
  
  # Finaler Clip
  
  # ---------------------------------------------------------
  
  aoi_dem_v <- terra::vect(aoi_dem_sf)
  
  dem_burgwald <- dem_merged |>
    terra::crop(aoi_dem_v) |>
    terra::mask(aoi_dem_v)
  
  terra::writeRaster(
    dem_burgwald,
    dem_out_file,
    overwrite = TRUE
  )
  
  message("Gespeichertes DGM1: ", dem_out_file)
  
} else {
  
  message("DGM bereits vorhanden, Merge wird übersprungen: ", dem_out_file)
  
  dem_burgwald <- terra::rast(dem_out_file)
  
}

# ---------------------------------------------------------

# 6) Geländeparameter

# ---------------------------------------------------------

slope_file <- here::here("data", "processed", "slope_burgwald.tif")
aspect_file <- here::here("data", "processed", "aspect_burgwald.tif")

if (!file.exists(slope_file)) {
  
  slope <- terra::terrain(dem_burgwald, "slope", unit = "degrees")
  
  terra::writeRaster(slope, slope_file, overwrite = TRUE)
  
} else {
  
  slope <- terra::rast(slope_file)
  
}

if (!file.exists(aspect_file)) {
  
  aspect <- terra::terrain(dem_burgwald, "aspect", unit = "degrees")
  
  terra::writeRaster(aspect, aspect_file, overwrite = TRUE)
  
} else {
  
  aspect <- terra::rast(aspect_file)
  
}

hill <- terra::shade(
  slope,
  aspect,
  angle = 40,
  direction = 315
)

# ---------------------------------------------------------

# 7) Luv / Lee Klassifikation

# ---------------------------------------------------------

wind_dir <- 270

wind_diff <- abs(aspect - wind_dir)

wind_diff <- terra::ifel(
  wind_diff > 180,
  360 - wind_diff,
  wind_diff
)

luv_lee <- terra::classify(
  wind_diff,
  matrix(
    c(
      0, 90, 1,
      90, 180, 2
    ),
    ncol = 3,
    byrow = TRUE
  )
)

names(luv_lee) <- "LuvLee"

luvlee_out_file <- here::here(
  "data",
  "processed",
  "luv_lee_burgwald.tif"
)

terra::writeRaster(
  luv_lee,
  luvlee_out_file,
  overwrite = TRUE
)

# ---------------------------------------------------------

# 8) Niederschlag laden (WorldClim)

# ---------------------------------------------------------

worldclimate_dir <- here::here("data", "raw", "worldclim")

fs::dir_create(worldclimate_dir)

prec <- geodata::worldclim_country(
  country = "DEU",
  var = "prec",
  path = worldclimate_dir
)

prec_dec <- prec[[12]]

aoi_prec_sf <- sf::st_transform(
  aoi_burgwald_wgs,
  terra::crs(prec_dec)
)

aoi_prec_v <- terra::vect(aoi_prec_sf)

prec_crop <- terra::crop(
  prec_dec,
  aoi_prec_v
) |>
  terra::mask(aoi_prec_v)

prec_res <- terra::project(
  prec_crop,
  dem_burgwald,
  method = "bilinear"
)

prec_out_file <- here::here(
  "data",
  "processed",
  "prec_dec_burgwald.tif"
)

terra::writeRaster(
  prec_res,
  prec_out_file,
  overwrite = TRUE
)

# ---------------------------------------------------------

# 9) Luv / Lee Niederschlag trennen

# ---------------------------------------------------------

luv_prec <- terra::mask(
  prec_res,
  luv_lee == 1
)

lee_prec <- terra::mask(
  prec_res,
  luv_lee == 2
)

# ---------------------------------------------------------

# 10) Statistik

# ---------------------------------------------------------

mean_luv <- as.numeric(
  terra::global(
    luv_prec,
    mean,
    na.rm = TRUE
  )
)

mean_lee <- as.numeric(
  terra::global(
    lee_prec,
    mean,
    na.rm = TRUE
  )
)

cat("Mittlerer Niederschlag (Dezember):\n")
cat("Luv:", round(mean_luv, 2), "\n")
cat("Lee:", round(mean_lee, 2), "\n")

# ---------------------------------------------------------

# 11) Karten

# ---------------------------------------------------------

png(
  here::here("outputs", "figures", "relief_burgwald.png"),
  width = 1400,
  height = 1000
)

plot(
  hill,
  col = gray(0:100 / 100),
  legend = FALSE,
  main = "Relief Burgwald"
)

plot(
  dem_burgwald,
  col = terrain.colors(50),
  alpha = 0.4,
  add = TRUE
)

dev.off()

png(
  here::here("outputs", "figures", "luv_lee_burgwald.png"),
  width = 1400,
  height = 1000
)

plot(
  luv_lee,
  col = c("dodgerblue", "tomato"),
  main = "Luv- und Lee-Seiten (Westwind)"
)

dev.off()

png(
  here::here("outputs", "figures", "prec_burgwald.png"),
  width = 1400,
  height = 1000
)

plot(
  prec_res,
  main = "Niederschlag Burgwald (Dezember)",
  col = hcl.colors(40, "YlGnBu")
)

dev.off()

png(
  here::here("outputs", "figures", "prec_luv_burgwald.png"),
  width = 1400,
  height = 1000
)

plot(
  luv_prec,
  main = "Niederschlag Luv",
  col = hcl.colors(40, "YlGnBu")
)

dev.off()

png(
  here::here("outputs", "figures", "prec_lee_burgwald.png"),
  width = 1400,
  height = 1000
)

plot(
  lee_prec,
  main = "Niederschlag Lee",
  col = hcl.colors(40, "YlGnBu")
)

dev.off()

message("Analyse abgeschlossen.")

############################################################

# 12) DWD Tagesniederschlag (RSK) – Interpolation über Hessen

############################################################

pkgs <- c(
  "here",
  "terra",
  "sf",
  "stars",
  "dplyr",
  "tidyr",
  "data.table",
  "gstat",
  "automap",
  "pbmcapply",
  "geodata"
)

new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]

if (length(new_pkgs)) install.packages(new_pkgs)

invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------------------------------------

# fun_clim_data.R laden

# ---------------------------------------------------------

root_folder <- here::here()

candidate_fun_files <- c(
  here::here("src", "r-libs", "fun_clim_data.R")
)

fun_file <- candidate_fun_files[file.exists(candidate_fun_files)][1]

if (is.na(fun_file) || length(fun_file) == 0) {
  
  stop(
    "fun_clim_data.R wurde nicht gefunden.\n",
    "Geprüft wurden:\n",
    paste(candidate_fun_files, collapse = "\n")
  )
  
}

source(fun_file)

message("Funktionsscript geladen: ", fun_file)

needed_funs <- c(
  "kl_daily_core_prepare",
  "kl_daily_extract_sf",
  "sanitize_climate_param"
)

missing_funs <- needed_funs[
  !vapply(needed_funs, exists, logical(1), mode = "function")
]

if (length(missing_funs) > 0) {
  
  stop(
    "Diese Funktionen fehlen nach source(): ",
    paste(missing_funs, collapse = ", ")
  )
  
}

# ---------------------------------------------------------

# Parameter

# ---------------------------------------------------------

epsg <- 3035
res <- 500
start_date <- as.Date("2022-01-01")
end_date <- as.Date("2024-12-31")
minStations <- 20
param <- "RSK"

# ---------------------------------------------------------

# Arbeitsordner

# ---------------------------------------------------------

clim_root <- file.path(
  root_folder,
  "data",
  "processed",
  "dwd_hessen"
)

if (!dir.exists(dirname(clim_root))) {
  
  clim_root <- file.path(
    tempdir(),
    "dwd_hessen"
  )
  
}

root_raw <- file.path(clim_root, "kl_daily_raw")
out_dir <- file.path(clim_root, "kl_daily_extracted")
target_dir <- file.path(clim_root, "RSK_daily_hessen")
clip_dir <- file.path(clim_root, "RSK_daily_burgwald")

dir.create(root_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clip_dir, recursive = TRUE, showWarnings = FALSE)

meta_file <- file.path(
  root_raw,
  "KL_Tageswerte_Beschreibung_Stationen.txt"
)

message("Arbeitsordner: ", clim_root)

# ---------------------------------------------------------

# Hessen AOI

# ---------------------------------------------------------

hessen <- geodata::gadm(
  country = "DEU",
  level = 1,
  path = tempdir()
)

hessen_sf <- sf::st_as_sf(hessen)

hessen_sf <- hessen_sf[hessen_sf$NAME_1 == "Hessen", ]

if (nrow(hessen_sf) == 0) {
  
  stop("Hessen konnte nicht gefunden werden.")
  
}

hessen_sf <- sf::st_transform(hessen_sf, epsg)

# ---------------------------------------------------------

# DEM Hessen

# ---------------------------------------------------------

dem_file <- file.path(
  clim_root,
  "dem_hessen_500m.tif"
)

if (!file.exists(dem_file)) {
  
  dem_de <- geodata::elevation_30s(
    country = "DEU",
    path = tempdir()
  )
  
  dem_de <- terra::project(
    dem_de,
    paste0("EPSG:", epsg)
  )
  
  hessen_v <- terra::vect(hessen_sf)
  
  dem_hessen <- dem_de |>
    terra::crop(hessen_v) |>
    terra::mask(hessen_v)
  
  dem_hessen <- terra::aggregate(
    dem_hessen,
    fact = 2,
    fun = mean,
    na.rm = TRUE
  )
  
  names(dem_hessen) <- "Stationshoehe"
  
  terra::writeRaster(
    dem_hessen,
    dem_file,
    overwrite = TRUE
  )
  
} else {
  
  dem_hessen <- terra::rast(dem_file)
  
}

dem <- stars::st_as_stars(dem_hessen)

gc()

dem_burgwald <- terra::rast(dem_out_file)
slope <- terra::rast(slope_file)
aspect <- terra::rast(aspect_file)
prec_res <- terra::rast(prec_out_file)

############################################################

# 13) Rasterauflösung reduzieren (für Segmentierung)

############################################################

agg_factor <- 40

dem_seg <- terra::aggregate(
  dem_burgwald,
  fact = agg_factor,
  fun = mean
)

slope_seg <- terra::aggregate(
  slope,
  fact = agg_factor,
  fun = mean
)

aspect_seg <- terra::aggregate(
  aspect,
  fact = agg_factor,
  fun = mean
)

prec_seg <- terra::aggregate(
  prec_res,
  fact = agg_factor,
  fun = mean
)

############################################################

# 14) Predictor Stack erstellen

############################################################

stack_file <- here::here(
  "data",
  "processed",
  "predictor_stack_burgwald.tif"
)

if(!file.exists(stack_file)){
  
  predictor_stack <- c(
    dem_seg,
    slope_seg,
    aspect_seg,
    prec_seg
  )
  
  names(predictor_stack) <- c(
    "elevation",
    "slope",
    "aspect",
    "precipitation"
  )
  
  terra::writeRaster(
    predictor_stack,
    stack_file,
    overwrite = TRUE
  )
  
}else{
  
  predictor_stack <- terra::rast(stack_file)
  
}

############################################################

# 15) Z-Score Standardisierung

############################################################

zscore_file <- here::here(
  "data",
  "processed",
  "predictor_stack_zscore.tif"
)

if(!file.exists(zscore_file)){
  
  vals <- terra::spatSample(
    predictor_stack,
    size = 50000,
    method = "random",
    as.df = TRUE,
    na.rm = TRUE
  )
  
  vals <- vals[complete.cases(vals), ]
  
  means <- colMeans(vals)
  
  sds <- apply(vals, 2, sd)
  
  z_stack <- predictor_stack
  
  for(i in 1:nlyr(predictor_stack)){
    
    z_stack[[i]] <- (predictor_stack[[i]] - means[i]) / sds[i]
    
  }
  
  terra::writeRaster(
    z_stack,
    zscore_file,
    overwrite = TRUE
  )
  
}else{
  
  z_stack <- terra::rast(zscore_file)
  
}

############################################################

# 16) PCA Feature Space

############################################################

pca_file <- here::here(
  "data",
  "processed",
  "pca_features_burgwald.tif"
)

if(!file.exists(pca_file)){
  
  vals <- terra::spatSample(
    z_stack,
    size = 30000,
    method = "random",
    as.df = TRUE,
    na.rm = TRUE
  )
  
  vals <- vals[complete.cases(vals), ]
  
  pca <- prcomp(vals)
  
  pc1 <- (z_stack[[1]] * pca$rotation[1,1]) +
    (z_stack[[2]] * pca$rotation[2,1]) +
    (z_stack[[3]] * pca$rotation[3,1]) +
    (z_stack[[4]] * pca$rotation[4,1])
  
  pc2 <- (z_stack[[1]] * pca$rotation[1,2]) +
    (z_stack[[2]] * pca$rotation[2,2]) +
    (z_stack[[3]] * pca$rotation[3,2]) +
    (z_stack[[4]] * pca$rotation[4,2])
  
  pca_raster <- c(pc1, pc2)
  
  names(pca_raster) <- c("PC1","PC2")
  
  terra::writeRaster(
    pca_raster,
    pca_file,
    overwrite = TRUE
  )
  
}else{
  
  pca_raster <- terra::rast(pca_file)
  
}

############################################################
# 17) Segmentierung (k-means + nearest neighbour)
############################################################

segments_file <- here::here(
  "data",
  "processed",
  "burgwald_segments.gpkg"
)

if(!file.exists(segments_file)){
  
  # Stichprobe
  sample_vals <- terra::spatSample(
    pca_raster,
    size = 15000,
    method = "random",
    as.df = TRUE,
    na.rm = TRUE
  )
  
  # k-means
  km <- kmeans(sample_vals, centers = 20, iter.max = 100)
  
  # alle Rasterwerte
  vals <- terra::values(pca_raster)
  
  good <- complete.cases(vals)
  
  clusters <- rep(NA, nrow(vals))
  
  nn <- RANN::nn2(
    data = km$centers,
    query = vals[good,],
    k = 1
  )
  
  clusters[good] <- nn$nn.idx
  
  # Raster erzeugen
  seg_raster <- pca_raster[[1]]
  
  terra::values(seg_raster) <- clusters
  
  # Raster verkleinern
  seg_raster_small <- terra::aggregate(
    seg_raster,
    fact = 5,
    fun = function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
  )
  
  # Polygonisierung
  segments_poly <- terra::as.polygons(
    seg_raster_small,
    dissolve = TRUE
  )
  
  segments_sf <- sf::st_as_sf(segments_poly)
  
  sf::st_write(
    segments_sf,
    segments_file,
    delete_dsn = TRUE
  )
  
}else{
  
  segments_sf <- sf::st_read(segments_file)
  
}
############################################################

# 18) Niederschlag pro Segment

############################################################

seg_raster <- terra::rasterize(
  terra::vect(segments_sf),
  prec_res,
  field = 1
)

segment_prec <- terra::zonal(
  prec_res,
  seg_raster,
  fun = "mean",
  na.rm = TRUE
)

segments_sf$prec_mean <- segment_prec$mean

############################################################

# 19) Segmente speichern

############################################################

segments_prec_file <- here::here(
  "data",
  "processed",
  "burgwald_segments_prec.gpkg"
)

sf::st_write(
  segments_sf,
  segments_prec_file,
  delete_dsn = TRUE
)

message("Segmentanalyse abgeschlossen.")

