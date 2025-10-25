library(dplyr)
library(stringr)
library(sf)
library(leaflet)

# Locality + DMS strings
dat_raw <- tibble::tribble(
  ~Locality, ~dms,
  1,  "50°32'30.05\"N, 13°26'28.39\"E",
  2,  "50°30'42.72\"N, 13°25'27.43\"E",
  3,  "50°30'22.64\"N, 13°24'30.10\"E",
  4,  "50°31'36.44\"N, 13°19'57.64\"E",
  5,  "50°31'40.13\"N, 13°20'24.96\"E",
  6,  "50°32'44.70\"N, 13°17'03.77\"E",
  7,  "50°33'33.68\"N, 13°17'19.52\"E",
  8,  "50°34'21.58\"N, 13°21'38.15\"E",
  9,  "50°35'16.41\"N, 13°19'36.85\"E",
  10, "50°35'24.64\"N, 13°21'31.69\"E",
  11, "50°35'57.02\"N, 13°22'32.76\"E",
  12, "50°37'28.17\"N, 13°23'48.48\"E",
  13, "50°37'53.25\"N, 13°24'02.11\"E",
  14, "50°38'13.17\"N, 13°39'55.43\"E",
  15, "50°38'08.31\"N, 13°40'05.60\"E",
  16, "50°40'26.54\"N, 13°32'30.66\"E",
  17, "50°41'39.80\"N, 13°34'34.51\"E",
  18, "50°41'32.08\"N, 13°38'06.10\"E",
  19, "50°46'39\"N,   14°05'27.01\"E",
  20, "50°46'31\"N,   14°06'11.05\"E"
)

# Convert DMS to decimal degrees
dms_to_dd <- function(x){
  m <- str_match(x, "(\\d+)°(\\d+)'([0-9.]+)\"?\\s*([NSEW])")
  stopifnot(!any(is.na(m)))
  deg <- as.numeric(m[,2]); min <- as.numeric(m[,3]); sec <- as.numeric(m[,4]); hemi <- m[,5]
  dd <- deg + min/60 + sec/3600
  dd[hemi %in% c("S","W")] <- -dd[hemi %in% c("S","W")]
  dd
}

# Split "lat, lon" and convert
dat <- dat_raw %>%
  tidyr::separate_wider_delim(dms, delim = ",", names = c("lat_dms", "lon_dms"), too_few = "error") %>%
  mutate(lat = dms_to_dd(str_trim(lat_dms)),
         lon = dms_to_dd(str_trim(lon_dms))) %>%
  arrange(Locality)

pts <- st_as_sf(dat, coords = c("lon","lat"), crs = 4326, remove = FALSE)
bb <- sf::st_bbox(pts)

leaflet(pts) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addMiniMap(toggleDisplay = TRUE, minimized = TRUE) |>
  addScaleBar(position = "bottomleft") |>
  addCircleMarkers(
    radius = 6, stroke = TRUE, weight = 1,
    popup = ~sprintf("<b>Locality %s</b><br>Lat: %.5f, Lon: %.5f", Locality, lat, lon)
  ) |>
  addPolylines(
    data = sf::st_cast(sf::st_combine(pts), "LINESTRING"),
    weight = 2, opacity = 0.8
  ) |>
  # use $ instead of ["..."]
  fitBounds(lng1 = as.numeric(bb["xmin"]), lat1 = as.numeric(bb["ymin"]),
            lng2 = as.numeric(bb["xmax"]), lat2 = as.numeric(bb["ymax"]))

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

# Reuse 'pts' from above (sf points, EPSG:4326)
# Download basemaps
eu <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
cz <- eu %>% dplyr::filter(admin == "Czechia")

# Project to an equal-area CRS for nice plotting
crs_plot <- 3035  # ETRS89 / LAEA Europe
eu_p  <- st_transform(eu, crs_plot)
cz_p  <- st_transform(cz, crs_plot)
pts_p <- st_transform(pts, crs_plot)

# Panel A: Europe with CZ highlighted
p_eu <- ggplot() +
  geom_sf(data = eu_p, fill = "grey92", color = "grey60", linewidth = 0.2) +
  geom_sf(data = cz_p, fill = "grey55", color = "grey25", linewidth = 0.3) +
  coord_sf(expand = FALSE) +
  theme_map
p_eu

# Assumes you already have:
#   pts      (sf, EPSG:4326) with your sites
#   cz       (sf, EPSG:4326) Czechia polygon
#   crs_plot <- 3035          # LAEA Europe used above
# and their projected versions:
eu_p  <- st_transform(ne_countries(scale = "medium", continent = "Europe", returnclass = "sf"), crs_plot)
cz_p  <- st_transform(cz, crs_plot)
pts_p <- st_transform(pts, crs_plot)

# Circle over Krušné hory 
# 1) Center = centroid of your sites (good for this dataset)
center_p <- st_centroid(st_union(pts_p))

# If you want to set an explicit center, uncomment (approx Krušné hory):
# center_p <- st_sfc(st_point(c(13.4, 50.55)), crs = 4326) |> st_transform(crs_plot)

# 2) Radius in kilometers (adjust as you like)
radius_km <- 35
circle_p <- st_buffer(center_p, dist = radius_km * 1000)  # buffer in meters

# Panel B: CZ with grey circle (35 km2)
theme_map <- theme_void() +
  theme(plot.title = element_text(hjust = 0, face = "bold"))

p_cz_overview <- ggplot() +
  geom_sf(data = cz_p, fill = "grey95", color = "grey50", linewidth = 0.25) +
  geom_sf(data = circle_p, fill = "grey70", color = "grey30", alpha = 0.25, linewidth = 0.4) +
  coord_sf(expand = FALSE) +
  theme_map

# Panel C: Zoom to circle with points
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
    data = cbind(pts_win, st_coordinates(pts_win)),
    aes(X, Y, label = Locality),
    size = 2.6, label.size = 0.15, label.padding = unit(0.08, "lines")
  ) +
  coord_sf(xlim = c(bb_exp["xmin"], bb_exp["xmax"]),
           ylim = c(bb_exp["ymin"], bb_exp["ymax"]),
           expand = FALSE) +
  theme_void()

p_cz_overview
p_cz_detail

tiff('p_cz_detail.tiff', units = "in", width = 8, height = 10, res = 600)
p_cz_detail
dev.off()
