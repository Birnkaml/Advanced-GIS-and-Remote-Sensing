############################################################
# README.R
# Burgwald – Luv/Lee-Effekte und Niederschlag
############################################################

# ==========================================================
# PROJEKTTITEL
# ==========================================================
# Einfluss von Luv-/Lee-Effekten auf den Niederschlag
# in einem kleinen Mittelgebirge:
# Fallstudie Burgwald

# ==========================================================
# PROJEKTÜBERSICHT
# ==========================================================
# Dieses Projekt untersucht die Forschungsfrage, ob
# Luv-/Lee-Effekte die Niederschlagsverteilung in kleineren
# Mittelgebirgen beeinflussen. Als Untersuchungsgebiet dient
# der Burgwald in Hessen, Deutschland.
#
# Der Workflow kombiniert:
# - Geländeanalyse auf Basis eines hochauflösenden DEM
# - Niederschlagsdaten aus WorldClim und DWD
# - räumliche Segmentierung und Stratifizierung
# - datengetriebene Auswahl repräsentativer Messpunkte
#
# Ziel ist ein transparentes, reproduzierbares Repository,
# das den Weg von den Rohdaten bis zum finalen
# Stationsvorschlag nachvollziehbar dokumentiert.

# ==========================================================
# FORSCHUNGSFRAGE
# ==========================================================
# Zentrale Frage:
# Haben Luv-/Lee-Effekte einen Einfluss auf die
# Niederschlagsverteilung im Burgwald?
#
# Arbeitshypothese:
# Luv-exponierte Hänge erhalten im Mittel mehr Niederschlag
# als Lee-exponierte Hänge, da auf der windzugewandten Seite
# orographische Hebung auftritt.

# ==========================================================
# INHALTSVERZEICHNIS
# ==========================================================
# 1. Erste Schritte
# 2. Installation
# 3. Projektstruktur
# 4. Workflow
# 5. Verwendung
# 6. Konfiguration
# 7. Abhängigkeiten
# 8. Ergebnisse
# 9. Limitationen
# 10. Weiterentwicklung

# ==========================================================
# 1. ERSTE SCHRITTE
# ==========================================================
# Voraussetzungen:
# - R (empfohlen: Version 4.x oder neuer)
# - RStudio oder VS Code
# - Internetverbindung für Paket- und Datendownloads
#
# Das Projekt ist vollständig skriptbasiert. Die Analyse
# kann nachvollzogen werden, indem die Skripte in der unten
# beschriebenen Reihenfolge ausgeführt werden.

# ==========================================================
# 2. INSTALLATION
# ==========================================================
# 1) Repository klonen:
#
#    git clone <DEIN-REPOSITORY-LINK>
#    cd <DEIN-REPOSITORY-ORDNER>
#
# 2) Setup-Skript ausführen:
#
#    source("src/01-setup-burgwald.R")
#
# Dieses Skript:
# - installiert fehlende Pakete
# - erstellt die Projektstruktur
# - definiert das Untersuchungsgebiet (AOI)
# - setzt Pfade und temporäre Arbeitsverzeichnisse
#
# Zugehöriges Skript:
# - 01-setup-burgwald.R

# ==========================================================
# 3. PROJEKTSTRUKTUR
# ==========================================================
# Empfohlene Repository-Struktur:
#
# .
# ├── data/
# │   ├── raw/                  # Rohdaten (DEM, DWD, WorldClim)
# │   ├── processed/            # Abgeleitete Raster und Vektordaten
# │   └── productive/           # Weiterverarbeitete Produkte
# ├── outputs/
# │   └── figures/              # Karten und Abbildungen
# ├── src/
# │   ├── 01-setup-burgwald.R
# │   ├── S1_S2_Observation_Features.R
# │   ├── S3_Structure.R
# │   ├── S4_Station_candidates.R
# │   └── S5_Decision_making.R
# └── README.R

# Zusätzlich verwendete Kurs-Bibliothek (Vorlagen & Funktionen):
#
# src/r-libs/
# ├── 01-fun-data-retrieval.R   # Funktionen zum Datenabruf
# ├── 01-setup-burgwald.R       # Setup (teilweise angepasst)
# ├── fun_clim_data.R           # Klimadaten-Funktionen (DWD etc.)
# ├── metrics-fun.R             # Hilfsfunktionen für Metriken
#
# Diese Skripte stammen aus dem Kursmaterial (Prof. Dr. Reudenbach)
# und werden als Funktionsbasis im Projekt verwendet (z.B. für
# Klimadatenverarbeitung und Preprocessing).

# Hinweis:
# Das README ist hier als R-Skript formuliert, damit es
# direkt im Projekt abgelegt und mit abgegeben werden kann.

# ==========================================================
# 4. WORKFLOW
# ==========================================================
# Der Workflow gliedert sich in fünf Hauptschritte:
#
# ----------------------------------------------------------
# Schritt 1: Setup und Umgebung
# ----------------------------------------------------------
# Skript:
# - 01-setup-burgwald.R
#
# Inhalte:
# - Laden und Installieren der benötigten Pakete
# - Aufbau der Ordnerstruktur
# - Definition der AOI für den Burgwald
# - Vorbereitung der Arbeitsumgebung
#
# ----------------------------------------------------------
# Schritt 2: Beobachtungsmerkmale / Features (S1-S2)
# ----------------------------------------------------------
# Skript:
# - S1_S2_Observation_Features.R
#
# Inhalte:
# - Download und Mosaikierung des DGM1
# - Berechnung von Höhe, Hangneigung und Exposition
# - Klassifikation von Luv- und Lee-Seiten
# - Einbindung von WorldClim-Niederschlag
# - Vergleich des mittleren Niederschlags zwischen Luv
#   und Lee
# - optionale DWD-Interpolation für Tagesniederschlag
#
# Kerngedanke:
# In diesem Schritt wird geprüft, ob sich erste
# Niederschlagsunterschiede zwischen windzugewandten und
# windabgewandten Hängen erkennen lassen.
#
# ----------------------------------------------------------
# Schritt 3: Landschaftsstruktur (S3)
# ----------------------------------------------------------
# Skript:
# - S3_Structure.R
#
# Inhalte:
# - Reduktion der Rasterauflösung für die Segmentierung
# - Aufbau eines Predictor-Stacks aus:
#   elevation, slope, aspect, precipitation
# - Z-Score-Standardisierung
# - PCA zur Dimensionsreduktion
# - k-means-Segmentierung des Untersuchungsraums
#
# Ziel:
# Bildung räumlich homogener Landschaftseinheiten als
# strukturelle Grundlage für die Stationsplanung.
#
# ----------------------------------------------------------
# Schritt 4: Stationskandidaten (S4)
# ----------------------------------------------------------
# Skript:
# - S4_Station_candidates.R
#
# Inhalte:
# - Berechnung physiographischer Metriken pro Segment:
#   elev_mean, slope_mean_deg, aspect_mean_deg,
#   southness_mean, windwardness_mean
# - Clustering der Segmente in physiographische Strata
# - Auswahl repräsentativer Segmente pro Stratum
# - Erzeugung erster Stationskandidaten
#
# Ziel:
# Sicherstellen, dass unterschiedliche Umweltbedingungen
# im Burgwald durch repräsentative Standorte abgedeckt
# werden.
#
# ----------------------------------------------------------
# Schritt 5: Entscheidungsfindung (S5)
# ----------------------------------------------------------
# Skript:
# - S5_Decision_making.R
#
# Inhalte:
# - Verfeinerung der Stationskandidaten
# - optionaler Hull-Filter
# - Auswahl mehrerer Kandidaten pro Stratum
# - Greedy-Spatial-Thinning mit Mindestabstand
# - Export finaler Segment- und Punktlayer
#
# Ziel:
# Entwicklung eines räumlich sinnvollen und möglichst
# redundanzarmen Messnetzes.

# ==========================================================
# 5. VERWENDUNG
# ==========================================================
# Gesamten Workflow ausführen:
#
# source("src/01-setup-burgwald.R")
# source("src/S1_S2_Observation_Features.R")
# source("src/S3_Structure.R")
# source("src/S4_Station_candidates.R")
# source("src/S5_Decision_making.R")
#
# Erwartete zentrale Outputs:
# - data/processed/aoi_burgwald.gpkg
# - data/processed/dem_dgm1_burgwald.tif
# - data/processed/slope_burgwald.tif
# - data/processed/aspect_burgwald.tif
# - data/processed/luv_lee_burgwald.tif
# - data/processed/prec_dec_burgwald.tif
# - data/processed/burgwald_segments.gpkg
# - data/processed/station_design/station_physio_candidates.gpkg
# - data/processed/station_design/physio_candidates_pts.gpkg
# - data/processed/station_design/physio_strata_segments.gpkg
#
# Erwartete Abbildungen:
# - outputs/figures/relief_burgwald.png
# - outputs/figures/luv_lee_burgwald.png
# - outputs/figures/prec_burgwald.png
# - outputs/figures/prec_luv_burgwald.png
# - outputs/figures/prec_lee_burgwald.png

# ==========================================================
# 6. KONFIGURATION
# ==========================================================
# Wichtige Parameter, die in den Skripten angepasst werden
# können:
#
# wind_dir / mean_wind_dir
# - Standardwert: 270
# - interpretiert die dominierende Windrichtung als Westwind
#
# k_physio
# - Anzahl der physiographischen Strata
#
# target_n_total
# - Gesamtzahl repräsentativer Stationskandidaten
#
# n_per_stratum
# - Anzahl der Kandidaten pro Stratum
#
# d_min_m
# - Mindestabstand zwischen finalen Stationspunkten
#
# agg_factor
# - Aggregationsfaktor der Rasterdaten für Segmentierung
#
# start_date / end_date
# - Zeitraum für die optionale DWD-Tagesniederschlagsanalyse
#
# Empfehlung:
# Für eine spätere Weiterentwicklung wäre eine zentrale
# Konfigurationsdatei sinnvoll.

# ==========================================================
# 7. ABHÄNGIGKEITEN
# ==========================================================
# Zentrale R-Pakete:
#
# Geodatenverarbeitung:
# - terra
# - sf
# - stars
#
# Datenverarbeitung:
# - dplyr
# - tidyr
# - data.table
# - tibble
#
# Klimadaten und Interpolation:
# - geodata
# - gstat
# - automap
# - pbmcapply
#
# Zonale Statistik und Hilfspakete:
# - exactextractr
# - here
# - fs
# - RANN
# - cluster
#
# Das Setup-Skript installiert viele dieser Pakete
# automatisch.

# ==========================================================
# 8. ERGEBNISSE
# ==========================================================
# Das Projekt erzeugt:
#
# - Reliefkarten des Burgwalds
# - Karten von Hangneigung und Exposition
# - eine Luv-/Lee-Klassifikation
# - Niederschlagskarten für den Untersuchungsraum
# - segmentierte Landschaftseinheiten
# - physiographische Strata
# - repräsentative Stationskandidaten
# - ein finales, ausgedünntes Stationsnetz
#
# Inhaltliche Kernaussage:
# Der Workflow prüft, ob im Burgwald ein systematischer
# Unterschied zwischen Luv- und Lee-Seiten hinsichtlich des
# Niederschlags erkennbar ist, und nutzt diese Information
# für die Entwicklung eines geeigneten Messnetzes.

# ==========================================================
# 9. LIMITATIONEN
# ==========================================================
# Wichtige Einschränkungen:
#
# - Die Windrichtung wird vereinfacht als konstant
#   angenommen.
# - WorldClim ist für kleinräumige Relief- und
#   Mittelgebirgseffekte relativ grob aufgelöst.
# - k-means benötigt eine vorab definierte Clusterzahl.
# - Eine vollständige Unsicherheitsanalyse ist bislang
#   nicht integriert.
# - Die Aussagekraft hängt auch von der Qualität und
#   Dichte der Niederschlagsdaten ab.

# ==========================================================
# 10. WEITERENTWICKLUNG
# ==========================================================
# Sinnvolle nächste Schritte:
#
# - Einbindung längerer DWD-Zeitreihen zur Validierung
# - Nutzung dynamischer Windfelder anstelle einer fixen
#   Windrichtung
# - Vergleich unterschiedlicher Interpolationsmethoden
# - systematische Unsicherheitsanalyse
# - Multi-Kriterien-Optimierung des Stationsnetzes
# - Ergänzung einer Lizenzdatei und weiterer Metadaten

# ==========================================================
# ZUSAMMENFASSUNG
# ==========================================================
# Dieses Projekt stellt einen transparenten und
# reproduzierbaren Workflow bereit, um den möglichen
# Einfluss von Luv-/Lee-Effekten auf den Niederschlag im
# Burgwald zu analysieren und daraus ein repräsentatives
# Stationsnetz abzuleiten.
############################################################
