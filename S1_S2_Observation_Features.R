
############################################################
# Burgwald: Luv-Lee & Niederschlag Analyse
# Verbessert mit echtem DGM1-Mosaik
############################################################


# ---------------------------------------------------------
# 1) Pakete
# ---------------------------------------------------------
packages <- c("terra", "sf", "here", "fs", "geodata")

new <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

invisible(lapply(packages, library, character.only = TRUE))

# ---------------------------------------------------------
# 2) Projektordner anlegen
# ---------------------------------------------------------
fs::dir_create(c(
  here::here("data/raw"),
  here::here("data/processed"),
  here::here("outputs/figures")
))

message("Project root: ", here::here())

# ---------------------------------------------------------
# 3) Helper-Funktionen
# ---------------------------------------------------------
download_if_missing <- function(url, destfile, mode = "wb", tries = 5) {
  if (!file.exists(destfile)) {
    
    for (i in 1:tries) {
      message("Download Versuch ", i, ": ", basename(destfile))
      
      tryCatch({
        download.file(url, destfile = destfile, mode = mode, quiet = FALSE)
        
        # Check: Datei existiert und ist nicht leer
        if (file.exists(destfile) && file.size(destfile) > 1000000) {
          message("Download erfolgreich!")
          return(TRUE)
        }
        
      }, error = function(e) {
        message("Fehler: ", e$message)
      })
      
      Sys.sleep(5) # kurze Pause
    }
    
    stop("Download fehlgeschlagen nach ", tries, " Versuchen: ", destfile)
    
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

aoi_file <- here::here("data", "processed", "aoi_burgwald.gpkg")
sf::st_write(sf::st_sf(geometry = aoi_burgwald_wgs), aoi_file, delete_dsn = TRUE, quiet = TRUE)

# ---------------------------------------------------------
# 5) DGM1 herunterladen und mosaikieren mit 3km Radius um aoi_burgwald => (nur 757 von 828 )
# ---------------------------------------------------------
dem_out_file <- here::here("data", "processed", "dem_dgm1_burgwald.tif")

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
  burgwald      = paste0(base_url_wf, "Burgwald%20-%20DGM1.zip"),
  gemuenden     = paste0(base_url_wf, "Gem%C3%BCnden%20(Wohra)%20-%20DGM1.zip"),
  rosenthal     = paste0(base_url_wf, "Rosenthal%20-%20DGM1.zip"),
  muenchhausen  = paste0(base_url_mr, "M%C3%BCnchhausen%20-%20DGM1.zip"),
  rauschenberg  = paste0(base_url_mr, "Rauschenberg%20-%20DGM1.zip"),
  coelbe        = paste0(base_url_mr, "C%C3%B6lbe%20-%20DGM1.zip"),
  lahntal       = paste0(base_url_mr, "Lahntal%20-%20DGM1.zip"),
  wohra         = paste0(base_url_mr, "Wohratal%20-%20DGM1.zip"),
  wetter        = paste0(base_url_mr, "Wetter%20(Hessen)%20-%20DGM1.zip"),
  haina         = paste0(base_url_wf, "Haina%20(Kloster)%20-%20DGM1.zip")
)

all_tif_files <- c()

for (nm in names(dgm1_urls)) {
  url_i <- dgm1_urls[[nm]]
  zip_file_i <- here::here("data", "raw", paste0("dgm1_", nm, ".zip"))
  unzipped_dir_i <- here::here("data", "raw", paste0("dgm1_", nm))
  
  download_if_missing(url_i, zip_file_i, mode = "wb")
  
  if (!fs::dir_exists(unzipped_dir_i)) {
    fs::dir_create(unzipped_dir_i)
    unzip(zip_file_i, exdir = unzipped_dir_i)
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
# DEM laden und direkt auf AOI + Puffer beschneiden
# ---------------------------------------------------------

dem_list <- lapply(all_tif_files, function(f) {
  
  r <- terra::rast(f)
  terra::crs(r) <- "EPSG:25832"
  
  # erst prüfen, ob sich die Extents überhaupt schneiden
  e_int <- terra::intersect(terra::ext(r), terra::ext(aoi_buffer_v))
  
  if (is.null(e_int)) {
    return(NULL)
  }
  
  r_crop <- terra::crop(r, aoi_buffer_v)
  
  return(r_crop)
  
})

# NULL-Einträge entfernen
dem_list <- Filter(Negate(is.null), dem_list)

# ---------------------------------------------------------
# Mosaic nur der relevanten Raster
# ---------------------------------------------------------

dem_merged <- do.call(terra::mosaic, c(dem_list, fun = "min"))

# ---------------------------------------------------------
# Finaler Clip auf ursprüngliche AOI
# ---------------------------------------------------------

aoi_dem_v <- terra::vect(aoi_dem_sf)

dem_burgwald <- dem_merged |>
  terra::crop(aoi_dem_v) |>
  terra::mask(aoi_dem_v)

dem_out_file <- here::here("data", "processed", "dem_dgm1_burgwald.tif")
terra::writeRaster(dem_burgwald, dem_out_file, overwrite = TRUE)

message("Gespeichertes DGM1: ", dem_out_file)

} else {
  message("DGM bereits vorhanden, Merge wird übersprungen: ", dem_out_file)
  dem_burgwald <- terra::rast(dem_out_file)
}

# ---------------------------------------------------------
# 6) Geländeparameter
# ---------------------------------------------------------
slope_file  <- here::here("data", "processed", "slope_burgwald.tif")
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
hill   <- terra::shade(slope, aspect, angle = 40, direction = 315)

# ---------------------------------------------------------
# 7) Luv / Lee Klassifikation
# ---------------------------------------------------------
wind_dir <- 270  # Westwind

wind_diff <- abs(aspect - wind_dir)
wind_diff <- terra::ifel(wind_diff > 180, 360 - wind_diff, wind_diff)

luv_lee <- terra::classify(
  wind_diff,
  matrix(
    c(0, 90, 1,
      90, 180, 2),
    ncol = 3,
    byrow = TRUE
  )
)

names(luv_lee) <- "LuvLee"

luvlee_out_file <- here::here("data", "processed", "luv_lee_burgwald.tif")
terra::writeRaster(luv_lee, luvlee_out_file, overwrite = TRUE)

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

# Dezember
prec_dec <- prec[[12]]

# AOI in Niederschlags-CRS transformieren
aoi_prec_sf <- sf::st_transform(aoi_burgwald_wgs, terra::crs(prec_dec))
aoi_prec_v  <- terra::vect(aoi_prec_sf)

prec_crop <- terra::crop(prec_dec, aoi_prec_v) |>
  terra::mask(aoi_prec_v)

# Auf DGM-Raster resamplen
prec_res <- terra::project(prec_crop, dem_burgwald, method = "bilinear")

prec_out_file <- here::here("data", "processed", "prec_dec_burgwald.tif")
terra::writeRaster(prec_res, prec_out_file, overwrite = TRUE)

# ---------------------------------------------------------
# 9) Luv / Lee Niederschlag trennen
# ---------------------------------------------------------
luv_prec <- terra::mask(prec_res, luv_lee == 1)
lee_prec <- terra::mask(prec_res, luv_lee == 2)

# ---------------------------------------------------------
# 10) Statistik
# ---------------------------------------------------------
mean_luv <- as.numeric(terra::global(luv_prec, mean, na.rm = TRUE))
mean_lee <- as.numeric(terra::global(lee_prec, mean, na.rm = TRUE))

cat("Mittlerer Niederschlag (Dezember):\n")
cat("Luv:", round(mean_luv, 2), "\n")
cat("Lee:", round(mean_lee, 2), "\n")


# ---------------------------------------------------------
# 11) Karten
# ---------------------------------------------------------
png(here::here("outputs", "figures", "relief_burgwald.png"), width = 1400, height = 1000)
plot(hill, col = gray(0:100 / 100), legend = FALSE, main = "Relief Burgwald")
plot(dem_burgwald, col = terrain.colors(50), alpha = 0.4, add = TRUE)
dev.off()

png(here::here("outputs", "figures", "luv_lee_burgwald.png"), width = 1400, height = 1000)
plot(luv_lee, col = c("dodgerblue", "tomato"), main = "Luv- und Lee-Seiten (Westwind)")
dev.off()

png(here::here("outputs", "figures", "prec_burgwald.png"), width = 1400, height = 1000)
plot(prec_res, main = "Niederschlag Burgwald (Dezember)", col = hcl.colors(40, "YlGnBu"))
dev.off()

png(here::here("outputs", "figures", "prec_luv_burgwald.png"), width = 1400, height = 1000)
plot(luv_prec, main = "Niederschlag Luv", col = hcl.colors(40, "YlGnBu"))
dev.off()

png(here::here("outputs", "figures", "prec_lee_burgwald.png"), width = 1400, height = 1000)
plot(lee_prec, main = "Niederschlag Lee", col = hcl.colors(40, "YlGnBu"))
dev.off()

message("Analyse abgeschlossen.")

############################################################
# 12) DWD Tagesniederschlag (RSK) – Interpolation über Hessen
# optional: Burgwald-Clip am Ende
############################################################

# ---------------------------------------------------------
# 1) Pakete
# ---------------------------------------------------------
pkgs <- c(
  "here", "terra", "sf", "stars", "dplyr", "tidyr", "data.table",
  "gstat", "automap", "pbmcapply", "geodata"
)

new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)

invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------------------------------------
# 2) fun_clim_data.R robust finden und laden
# ---------------------------------------------------------
root_folder <- here::here()

candidate_fun_files <- c(file.path(root_folder, "src", "r-libs", "fun_clim_data.R"),
                         "fun_clim_data.R"
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

# Prüfen, ob die benötigten Funktionen da sind
needed_funs <- c("kl_daily_core_prepare", "kl_daily_extract_sf", "sanitize_climate_param")
missing_funs <- needed_funs[!vapply(needed_funs, exists, logical(1), mode = "function")]

if (length(missing_funs) > 0) {
  stop("Diese Funktionen fehlen nach source(): ", paste(missing_funs, collapse = ", "))
}

# ---------------------------------------------------------
# 3) Parameter
# ---------------------------------------------------------
epsg <- 3035
res  <- 500

start_date <- as.Date("2022-01-01")
end_date   <- as.Date("2024-12-31")

minStations <- 20
param <- "RSK"

# ---------------------------------------------------------
# 4) Einfache lokale Arbeitsordner
#    -> unabhängig von envrmt / alter Projektlogik
# ---------------------------------------------------------
clim_root  <- file.path(root_folder, "data", "processed", "dwd_hessen")
if (!dir.exists(dirname(clim_root))) {
  clim_root <- file.path(tempdir(), "dwd_hessen")
}

root_raw   <- file.path(clim_root, "kl_daily_raw")
out_dir    <- file.path(clim_root, "kl_daily_extracted")
target_dir <- file.path(clim_root, "RSK_daily_hessen")
clip_dir   <- file.path(clim_root, "RSK_daily_burgwald")

dir.create(root_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clip_dir, recursive = TRUE, showWarnings = FALSE)

meta_file <- file.path(root_raw, "KL_Tageswerte_Beschreibung_Stationen.txt")

message("Arbeitsordner: ", clim_root)

# ---------------------------------------------------------
# 5) Hessen als AOI für die Interpolation
# ---------------------------------------------------------
hessen <- geodata::gadm(country = "DEU", level = 1, path = tempdir())
hessen_sf <- sf::st_as_sf(hessen)
hessen_sf <- hessen_sf[hessen_sf$NAME_1 == "Hessen", ]

if (nrow(hessen_sf) == 0) {
  stop("Hessen konnte in den GADM-Daten nicht gefunden werden.")
}

hessen_sf <- sf::st_transform(hessen_sf, epsg)

# ---------------------------------------------------------
# 6) DEM für Hessen vorbereiten
# ---------------------------------------------------------
dem_file <- file.path(clim_root, "dem_hessen_500m.tif")

if (!file.exists(dem_file)) {
  message("Lade DEM für Deutschland ...")
  dem_de <- geodata::elevation_30s(country = "DEU", path = tempdir())
  dem_de <- terra::project(dem_de, paste0("EPSG:", epsg))
  
  hessen_v <- terra::vect(hessen_sf)
  
  dem_hessen <- dem_de |>
    terra::crop(hessen_v) |>
    terra::mask(hessen_v)
  
  # groberes Raster für Klimainterpolation
  dem_hessen <- terra::aggregate(dem_hessen, fact = 2, fun = mean, na.rm = TRUE)
  
  names(dem_hessen) <- "Stationshoehe"
  terra::writeRaster(dem_hessen, dem_file, overwrite = TRUE)
} else {
  dem_hessen <- terra::rast(dem_file)
}

dem <- stars::st_as_stars(dem_hessen)

# ---------------------------------------------------------
# 7) DWD-Kerndaten vorbereiten
# ---------------------------------------------------------


# fehlende Metadatei gezielt herunterladen
meta_url <- paste0(
  "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/",
  "KL_Tageswerte_Beschreibung_Stationen.txt"
)

if (!file.exists(meta_file)) {
  download.file(meta_url, destfile = meta_file, mode = "wb")
}

if (!file.exists(meta_file)) {
  stop("Metadatei konnte nicht geladen werden: ", meta_file)
}



# Metadatei in UTF-8 bereinigen
raw_lines <- readLines(meta_file, warn = FALSE, encoding = "unknown")
clean_lines <- iconv(raw_lines, from = "", to = "UTF-8", sub = "")
clean_lines[is.na(clean_lines)] <- ""
writeLines(clean_lines, meta_file, useBytes = TRUE)

# Dann erst dt_core erzeugen
dt_core <- kl_daily_core_prepare(
  start_date = start_date,
  end_date   = end_date,
  max_missing_frac = 0.25,
  delta_t = 7L,
  keep_cols = c(
    "STATIONS_ID", "MESS_DATUM", "QN_4",
    "RSK", "RSKF", "SDK", "SHK_TAG", "NM", "VPM",
    "PM", "TMK", "UPM", "TXK", "TNK", "TGK", "eor"
  ),
  root_raw   = root_raw,
  out_dir    = out_dir,
  meta_file  = meta_file,
  do_download = TRUE,
  force_unzip = FALSE,
  force_merge = TRUE,
  auto_resolve_duplicates = TRUE,
  write_merged_csv = TRUE
)

dt_core_rds <- file.path(clim_root, "dt_core.rds")
saveRDS(dt_core, dt_core_rds)

# ---------------------------------------------------------
# 8) RSK als sf extrahieren
# ---------------------------------------------------------
cVar.sf <- kl_daily_extract_sf(
  dt = dt_core,
  param = "RSK"
)

# auf Hessen beschränken
cVar.sf <- sf::st_transform(cVar.sf, epsg)
cVar.sf <- sf::st_intersection(cVar.sf, hessen_sf)

message("DWD-Punkte in Hessen: ", nrow(cVar.sf))
message("Eindeutige Stationen: ", length(unique(cVar.sf$STATIONS_ID)))

if (nrow(cVar.sf) == 0) {
  stop("Nach dem Hessen-Filter sind keine DWD-Daten übrig.")
}

# ---------------------------------------------------------
# 9) Alle Tage bestimmen
# ---------------------------------------------------------
dat_list_all <- sort(unique(as.Date(cVar.sf$MESS_DATUM)))

tifs_existing <- list.files(target_dir, pattern = "\\.tif$", full.names = FALSE)
dates_done <- character(0)

if (length(tifs_existing) > 0) {
  dates_done <- sub(paste0("_", param, "\\.tif$"), "", basename(tifs_existing))
}

dat_list_window <- dat_list_all[
  dat_list_all >= start_date & dat_list_all <= end_date
]

dat_list_todo <- setdiff(as.character(dat_list_window), dates_done)

message(sprintf(
  "[KED] %s: bereits vorhanden: %d | zu berechnen: %d",
  param, length(dates_done), length(dat_list_todo)
))

# ---------------------------------------------------------
# 10) Tagesweise Interpolation über Hessen
# ---------------------------------------------------------
if (length(dat_list_todo) > 0) {
  
  invisible(
    pbmcapply::pbmclapply(seq_along(dat_list_todo), function(n) {
      currentDate <- as.Date(dat_list_todo[n])
      cd <- format(currentDate, "%Y-%m-%d")
      target_file <- file.path(target_dir, paste0(cd, "_", param, ".tif"))
      
      if (file.exists(target_file)) {
        message("Skipping existing file: ", target_file)
        return(NULL)
      }
      
      cVar.sf.day <- cVar.sf[cVar.sf$MESS_DATUM == currentDate, ]
      
      dat <- sanitize_climate_param(cVar.sf.day, param, date = currentDate)
      
      geom_col <- intersect(c("geometry", "geom"), names(dat))
      if (length(geom_col) == 1 && geom_col != "geometry") {
        names(dat)[names(dat) == geom_col] <- "geometry"
      }
      
      sf::st_geometry(dat) <- "geometry"
      dat <- dat[, c("Stationshoehe", "tmp", "geometry")]
      
      dat$tmp <- as.numeric(dat$tmp)
      names(dat)[names(dat) == "tmp"] <- param
      
      data <- tidyr::drop_na(dat)
      
      if (nrow(data) == 0) {
        message("Keine gültigen Daten für: ", cd)
        return(NULL)
      }
      
      if (sf::st_crs(data) != sf::st_crs(dem)) {
        data <- sf::st_transform(data, crs = epsg)
      }
      
      if (sum(!is.na(data[[param]])) > minStations) {
        data <- dplyr::distinct(data, geometry, .keep_all = TRUE)
        data <- sf::st_transform(data, sf::st_crs(dem))
        
        vm.auto <- automap::autofitVariogram(
          as.formula(paste(param, "~1")),
          input_data = data
        )
        
        pred <- gstat::krige(
          formula     = as.formula(paste(param, "~Stationshoehe")),
          locations   = data,
          newdata     = dem,
          model       = vm.auto$var_model,
          debug.level = -1
        )
        
        stars::write_stars(
          pred,
          target_file,
          NA_value = -9999,
          overwrite = TRUE
        )
        
        rm(pred)
        gc()
        
      } else {
        message("Zu wenige Stationen für ", cd, " (", sum(!is.na(data[[param]])), ")")
        stars::write_stars(
          dem * 0 - 9999,
          target_file,
          overwrite = TRUE
        )
      }
      
      message("Written: ", target_file)
      return(target_file)
      
    }, mc.cores = 4, mc.allow.recursive = TRUE)
  )
  
} else {
  message(sprintf(
    "[KED] %s: Keine offenen Tage im Zeitraum %s bis %s.",
    param, start_date, end_date
  ))
}

# ---------------------------------------------------------
# 11) Optional: auf Burgwald clippen
# ---------------------------------------------------------
if (exists("aoi_burgwald_wgs")) {
  message("Clippe Hessen-Raster auf Burgwald ...")
  
  aoi_bw <- sf::st_transform(
    sf::st_as_sf(data.frame(id = 1), geometry = aoi_burgwald_wgs),
    epsg
  )
  aoi_bw_v <- terra::vect(aoi_bw)
  
  tif_files <- list.files(target_dir, pattern = "\\.tif$", full.names = TRUE)
  
  for (f in tif_files) {
    out_bw <- file.path(clip_dir, basename(f))
    
    r <- terra::rast(f)
    r_bw <- r |>
      terra::crop(aoi_bw_v) |>
      terra::mask(aoi_bw_v)
    
    terra::writeRaster(r_bw, out_bw, overwrite = TRUE)
  }
  
  message("Burgwald-Clip gespeichert in: ", clip_dir)
} else {
  message("aoi_burgwald_wgs nicht gefunden -> Hessen-Raster werden nicht auf Burgwald geclippt.")
}


# ---------------------------------------------------------
# 12) Ergebnis anzeigen
# ---------------------------------------------------------
cat("\nFertig.\n")
cat("Hessen-Raster liegen in:\n", target_dir, "\n")
cat("Burgwald-Raster liegen in:\n", clip_dir, "\n")
cat("Kerndaten RDS:\n", dt_core_rds, "\n")


