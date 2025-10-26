library(dplyr)
library(stringr)
library(sf)
library(leaflet)

# Include Altitude in tibble format
dat_raw <- tibble::tribble(
  ~Locality, ~Altitude, ~Exposition, ~dms,
  1, 600, "W – SW",  "50°32'30.05\"N, 13°26'28.39\"E",
  2, 520, "E – SE",  "50°30'42.72\"N, 13°25'27.43\"E",
  3, 520, "SE",      "50°30'22.64\"N, 13°24'30.10\"E",
  4, 790, "Plain",   "50°31'36.44\"N, 13°19'57.64\"E",
  5, 790, "Plain",   "50°31'40.13\"N, 13°20'24.96\"E",
  6, 810, "Plain",   "50°32'44.70\"N, 13°17'03.77\"E",
  7, 820, "NE",      "50°33'33.68\"N, 13°17'19.52\"E",
  8, 800, "NW",      "50°34'21.58\"N, 13°21'38.15\"E",
  9, 700, "W",       "50°35'16.41\"N, 13°19'36.85\"E",
  10, 760, "W",       "50°35'24.64\"N, 13°21'31.69\"E",
  11, 800, "SW",      "50°35'57.02\"N, 13°22'32.76\"E",
  12, 580, "NE",      "50°37'28.17\"N, 13°23'48.48\"E",
  13, 520, "NE",   "50°37'53.25\"N, 13°24'02.11\"E",
  14, 680, "W – SW",  "50°38'13.17\"N, 13°39'55.43\"E",
  15, 700, "W",       "50°38'08.31\"N, 13°40'05.60\"E",
  16, 800, "N",       "50°40'26.54\"N, 13°32'30.66\"E",
  17, 740, "SW – S",  "50°41'39.80\"N, 13°34'34.51\"E",
  18, 770, "NW – W",  "50°41'32.08\"N, 13°38'06.10\"E",
  19, 390, "S",       "50°46'39\"N,   14°05'27.01\"E",
  20, 480, "S",       "50°46'31\"N,   14°06'11.05\"E"
)

# Convert DMS to decimal degrees
dms_to_dd <- function(x){
  m <- stringr::str_match(x, "(\\d+)°(\\d+)'([0-9.]+)\"?\\s*([NSEW])")
  stopifnot(!any(is.na(m)))
  deg  <- as.numeric(m[,2])
  min  <- as.numeric(m[,3])
  sec  <- as.numeric(m[,4])
  hemi <- m[,5]
  dd <- deg + min/60 + sec/3600
  dd[hemi %in% c("S","W")] <- -dd[hemi %in% c("S","W")]
  dd
}

# Split latitude and longitude, and convert
dat <- dat_raw %>%
  tidyr::separate_wider_delim(
    dms, delim = ",", names = c("lat_dms", "lon_dms"), too_few = "error"
  ) %>%
  mutate(
    lat = dms_to_dd(str_trim(lat_dms)),
    lon = dms_to_dd(str_trim(lon_dms))
  ) %>%
  arrange(Locality)

pts <- st_as_sf(dat, coords = c("lon","lat"), crs = 4326, remove = FALSE)
bb <- sf::st_bbox(pts)

leaflet(pts) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addMiniMap(toggleDisplay = TRUE, minimized = TRUE) |>
  addScaleBar(position = "bottomleft") |>
  addCircleMarkers(
    radius = 6, stroke = TRUE, weight = 1,
    popup = ~sprintf(
      "<b>Locality %s.</b><br>Altitude: %dm<br>Exposition: %s<br>Lat: %.5f, Lon: %.5f",
      Locality, Altitude, Exposition, lat, lon
    )
  ) |>
  addPolylines(
    data = sf::st_cast(sf::st_combine(pts), "LINESTRING"),
    weight = 2, opacity = 0.8
  ) |>
  fitBounds(lng1 = as.numeric(bb["xmin"]), lat1 = as.numeric(bb["ymin"]),
            lng2 = as.numeric(bb["xmax"]), lat2 = as.numeric(bb["ymax"]))
# --- --- --- 
# Map designs with ggplot2
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

# A panel with overall Europe with Czech republic highlighted
#####
# Load Europe and Czech Republic polygons
eu <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
cz <- eu %>% dplyr::filter(admin == "Czechia")

# Project to equal-area CRS for better map proportions
crs_plot <- 3035  # ETRS89 / LAEA Europe
eu_p <- st_transform(eu, crs_plot)
cz_p <- st_transform(cz, crs_plot)

# Theme
theme_map <- theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0))

# Plot: Europe with Czech Republic highlighted in grey
p_eu_cz <- ggplot() +
  geom_sf(data = eu_p, fill = "grey90", color = "grey70", linewidth = 0.2) +
  geom_sf(data = cz_p, fill = "grey40", color = "grey20", linewidth = 0.3) +
  coord_sf(expand = FALSE) +
  theme_map

p_eu_cz

tiff('p_eu_cz.tiff', units = "in", width = 8, height = 10, res = 600)
p_eu_cz
dev.off()

pdf("p_eu_cz.pdf", width = 8, height = 10)
p_eu_cz
dev.off()
#####
# A panel with Central Europe
#####
# Load map data
eu <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
cz <- eu %>% filter(admin == "Czechia")

# Project to equal-area projection for Europe
crs_plot <- 3035  # ETRS89 / LAEA Europe
eu_p <- st_transform(eu, crs_plot)
cz_p <- st_transform(cz, crs_plot)

# Compute centroids for labeling
centroids <- st_centroid(eu_p)

# Extract coordinates and country codes for labeling
label_df <- centroids %>%
  mutate(X = st_coordinates(geometry)[,1],
         Y = st_coordinates(geometry)[,2]) %>%
  st_drop_geometry() %>%
  select(admin, iso_a2, X, Y)

label_df <- label_df %>%
  filter(iso_a2 %in% c("CZ", "DE", "PL", "AT", "SK", "HU", "CH", "SI"))

# Center map around Czech Republic
cz_center <- st_centroid(cz_p)
cz_center_coords <- st_coordinates(cz_center)

# Define zoom window (in meters)
x_range <- 700000   # horizontal distance
y_range <- 550000   # vertical distance
xlim <- c(cz_center_coords[1] - x_range, cz_center_coords[1] + x_range)
ylim <- c(cz_center_coords[2] - y_range, cz_center_coords[2] + y_range)

# Theme
theme_map <- theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0))

# Plot
p_central_labels <- ggplot() +
  geom_sf(data = eu_p, fill = "grey90", color = "grey70", linewidth = 0.25) +
  geom_sf(data = cz_p, fill = "grey40", color = "grey20", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = X, y = Y, label = iso_a2),
    size = 3,
    fontface = "bold",
    color = "black",
    segment.size = 0.2,
    min.segment.length = 0
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_map

p_central_labels

tiff('p_central_labels.tiff', units = "in", width = 8, height = 10, res = 600)
p_central_labels
dev.off()

pdf("p_central_labels.pdf", width = 8, height = 10)
p_central_labels
dev.off()
#####
# B panel with square in Czech republic
#####
# Download basemaps
eu <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
cz <- eu %>% filter(admin == "Czechia")
crs_plot <- 3035
eu_p  <- st_transform(eu, crs_plot)
cz_p  <- st_transform(cz, crs_plot)
pts_p <- st_transform(pts, crs_plot)

eu
# Use the centroid of your sites as before
center_p <- st_centroid(st_union(pts_p))

# Define half-size of the square (in km)
half_size_km <- 35     # adjust as needed
half_size_m  <- half_size_km * 1000

# Convert center to numeric coordinates
cx <- st_coordinates(center_p)[1,1]
cy <- st_coordinates(center_p)[1,2]

# Create a square polygon centered on the site cluster
square_mat <- matrix(c(
  cx - half_size_m, cy - half_size_m,
  cx + half_size_m, cy - half_size_m,
  cx + half_size_m, cy + half_size_m,
  cx - half_size_m, cy + half_size_m,
  cx - half_size_m, cy - half_size_m
), ncol = 2, byrow = TRUE)

square_p <- st_polygon(list(square_mat)) |> 
  st_sfc(crs = st_crs(pts_p))

# Czech overview map with square highlight
p_cz_overview <- ggplot() +
  geom_sf(data = cz_p, fill = "grey95", color = "grey50", linewidth = 0.25) +
  geom_sf(data = square_p, fill = "grey70", color = "grey30", alpha = 0.3, linewidth = 0.4) +
  theme_void()
tiff('p_cz_overview.tiff', units = "in", width = 8, height = 10, res = 600)
p_cz_overview
dev.off()

pdf("p_cz_overview.pdf", width = 8, height = 10)
p_cz_overview
dev.off()
#####
# B panel with describtions
#####
pad <- 0.10
bb  <- st_bbox(circle_p)
xr <- as.numeric(bb["xmax"] - bb["xmin"])
yr <- as.numeric(bb["ymax"] - bb["ymin"])
bb_exp <- st_bbox(
  c(xmin = as.numeric(bb["xmin"]) - pad*xr,
    ymin = as.numeric(bb["ymin"]) - pad*yr,
    xmax = as.numeric(bb["xmax"]) + pad*xr,
    ymax = as.numeric(bb["ymax"]) + pad*yr),
  crs = st_crs(circle_p)
)
win <- st_as_sfc(bb_exp)

line_p <- st_cast(st_combine(pts_p), "LINESTRING")
cz_win   <- st_crop(cz_p,  win)
line_win <- suppressWarnings(st_intersection(line_p, win))
pts_win  <- pts_p[ st_within(pts_p, win, sparse = FALSE), ]

p_cz_detail <- ggplot() +
  geom_sf(data = cz_win,   fill = "grey98", color = "grey70", linewidth = 0.2) +
  geom_sf(data = line_win, color = "grey40", linewidth = 0.5) +
  geom_sf(data = pts_win,  size = 2.2, color = "black") +
  ggrepel::geom_label_repel(
    data = cbind(sf::st_drop_geometry(pts_win), sf::st_coordinates(pts_win)),
    aes(
      X, Y,
      label = sprintf("%s.\n%dm\n%s", Locality, Altitude, Exposition)
    ),
    size = 2.6, label.size = 0.15, label.padding = unit(0.08, "lines")
  ) +
  coord_sf(xlim = c(bb_exp["xmin"], bb_exp["xmax"]),
           ylim = c(bb_exp["ymin"], bb_exp["ymax"]),
           expand = FALSE) +
  theme_void()

p_cz_overview
p_cz_detail

tiff('p_cz_overview.tiff', units = "in", width = 8, height = 10, res = 600)
p_cz_overview
dev.off()

pdf("p_cz_detail.pdf", width = 8, height = 10)
p_cz_detail
dev.off()
