install.packages(c("sf", "ggplot2", "rnaturalearth", "rnaturalearthdata", "dplyr", "ggspatial"))

library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(lwgeom)

world <- ne_countries(scale = "medium", returnclass = "sf")

### Central Europe area distribution ###
base_bw <- ggplot(world) +
  geom_sf(fill = "white", color = "black", linewidth = 0.2) +
  theme_void()
base_bw

central_europe_names <- c(
  "Austria","Czechia","Germany","Poland","Slovakia","Hungary",
  "Switzerland","Liechtenstein","Slovenia"
)

central_europe <- world |> filter(name %in% central_europe_names)

p_ce <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black", linewidth = 0.2) +
  geom_sf(data = central_europe, fill = "grey80", color = "black", linewidth = 0.3) +
  coord_sf(xlim = c(5, 25), ylim = c(45, 56), expand = FALSE) +
  theme_void()

p_ce

sort(unique(world$name))


# or vector (great for journals)
ggsave("map_central_europe_bw.pdf", p_ce, width = 6, height = 5, bg = "white")

### Circumboreal area distribution ###
# Use GEOS (more forgiving than s2 for messy polygons)
sf::sf_use_s2(FALSE)

# 1) Data ----------------------------------------------------
world  <- ne_countries(scale = "medium", returnclass = "sf") %>% st_transform(4326)
states <- ne_states(returnclass = "sf") %>% st_transform(4326)

# Make valid early
world  <- sf::st_make_valid(world)
states <- sf::st_make_valid(states)

# Helper: polygon north of a latitude
north_of <- function(lat){
  st_as_sfc(st_bbox(c(xmin = -180, ymin = lat, xmax = 180, ymax = 90), crs = 4326))
}

# 2) NORTH AMERICA -------------------------------------------
canada <- world %>% filter(name == "Canada")

# If name_en doesn't exist, replace with name:
# alaska <- states %>% filter(iso_a2 == "US", name == "Alaska")
alaska <- states %>% filter(iso_a2 == "US", name_en == "Alaska")

north_us_states <- c(
  "Washington","Oregon","Idaho","Montana","Wyoming","North Dakota","South Dakota",
  "Minnesota","Wisconsin","Michigan","Vermont","New Hampshire","Maine","New York",
  "Massachusetts","Connecticut","Rhode Island","Pennsylvania","New Jersey"
)

north_us <- states %>% filter(iso_a2 == "US", name_en %in% north_us_states)

circ_na <- bind_rows(
  canada   %>% transmute(region = "Canada", geometry),
  alaska   %>% transmute(region = "Alaska", geometry),
  north_us %>% transmute(region = "N_US",   geometry)
)

# 3) EURASIA -------------------------------------------------
eu_clip_lat <- 44
europe_countries <- world %>% filter(continent == "Europe")
europe_north <- st_intersection(europe_countries, north_of(eu_clip_lat))

russia <- world %>% filter(name == "Russia")

china <- world %>% filter(name == "China")
china_north <- st_intersection(china, north_of(40))

japan <- world %>% filter(name == "Japan")
japan_north <- st_intersection(japan, north_of(38))

circ_eu <- bind_rows(
  europe_north %>% transmute(region = "Europe_N", geometry),
  russia       %>% transmute(region = "Russia",   geometry),
  china_north  %>% transmute(region = "China_N",  geometry),
  japan_north  %>% transmute(region = "Japan_N",  geometry)
)

# 4) MERGE + CLEAN + UNION -----------------------------------
circ_all <- bind_rows(circ_na, circ_eu)

# Robust geometry cleaning BEFORE union
circ_all <- sf::st_make_valid(circ_all)

# Keep only polygons (drop lines/points created by intersections)
circ_all <- st_collection_extract(circ_all, "POLYGON")

# Union into a single circumboreal geometry
circumboreal <- st_union(circ_all) %>%
  st_sf(geometry = .)

# 5) PLOT ----------------------------------------------------
p <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black", linewidth = 0.15) +
  geom_sf(data = circumboreal, fill = "grey80", color = "black", linewidth = 0.25) +
  coord_sf(xlim = c(-170, 170), ylim = c(30, 85), expand = FALSE) +
  theme_void()

print(p)

# 6) SAVE ----------------------------------------------------
ggsave("circumboreal_bw.png", p, width = 12, height = 5, dpi = 600, bg = "white")
ggsave("circumboreal_bw.pdf", p, width = 12, height = 5, bg = "white")
