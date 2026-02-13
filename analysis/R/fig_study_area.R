library(terra)
library(ggplot2)
library(sf)
library(dplyr)
library(ggnewscale)
library(stringr) #str_wrap

#library(USA.state.boundaries)
#data(state_boundaries_wgs84)
# alaska shape in this package has the aleutian islands wrapped around to other hemispehre

# AK CAN AEA crs
ak_can_aea_raw <- readLines("data/ancillary/AK_Canada_aea.prj")
ak_can_aea_noquotes <- paste(ak_can_aea_raw, collapse = "") # Collapse lines into a single string
ak_can_aea_crs <- st_crs(ak_can_aea_noquotes)

# Load shapefiles
alaska <- read_sf("data/raw/ecoregions/alaska/alaska.shp") %>%
    st_transform(ak_can_aea_crs)
canada <- read_sf("data/raw/ecoregions/CANprovinces/world-administrative-boundaries.shp") %>%
    st_transform(ak_can_aea_crs)
ak_can <- st_union(alaska, canada)
ecoregions <- read_sf("data/processed/ecoregions/ecos.shp") %>%
    st_transform(ak_can_aea_crs) 
historical_fires <- read_sf("data/processed/fire/shapefiles/AK_Canada_large_fire_history.shp") %>%
    st_transform(ak_can_aea_crs)  %>%
    mutate(decade=floor(FIREYR/10)*10) %>%  
    st_intersection(ecoregions)
eco_50607 <- ecoregions %>%
    filter(ECO_ID == 50607) 
eco_50609 <- ecoregions %>%
    filter(ECO_ID == 50609) 

# Load raster
tree_cover_raw <- rast("data/processed/veg/mod44b/MOD44B.006_Percent_Tree_Cover_doy2020065_aid0001.tif") 
tree_cover_filt <- clamp(tree_cover_raw, upper=100, value=FALSE) #%>% #ocean >100 for some reason
    #terra::aggregate( fact=100) #TEMP
tree_cover_ecos_crop <- crop(tree_cover_filt, ecoregions, mask = FALSE)
tree_cover_ecos_mask <- mask(tree_cover_ecos_crop, ecoregions, touches = TRUE)
tree_cover_df <- as.data.frame(tree_cover_ecos_mask, xy = TRUE) %>%
    dplyr::rename(perc_tree=MOD44B.006_Percent_Tree_Cover_doy2020065_aid0001) 
    
## Generate label geometries:
# Create 'centroids' of municipalities (you can safely ignore the warning)
sf_centroids <- ecoregions %>% 
    st_centroid() %>%
    mutate(lon = case_when(ECO_ID==50602 ~ st_coordinates(.)[,1] + 300000,
                           ECO_ID==50605 ~ st_coordinates(.)[,1] - 200000,
                           ECO_ID==50606 ~ st_coordinates(.)[,1] - 100000,
                           ECO_ID==50607 ~ st_coordinates(.)[,1] - 100000,
                           ECO_ID==50608 ~ st_coordinates(.)[,1] - 0,
                           ECO_ID==50609 ~ st_coordinates(.)[,1] + 300000,
                           ECO_ID==50610 ~ st_coordinates(.)[,1] + 50000,
                           ECO_ID==50612 ~ st_coordinates(.)[,1] - 0,
                           ECO_ID==50613 ~ st_coordinates(.)[,1] - 0,
                           ECO_ID==50614 ~ st_coordinates(.)[,1] - 0,
                           ECO_ID==50616 ~ st_coordinates(.)[,1] - 0,
                           ECO_ID==51111 ~ st_coordinates(.)[,1] - 0),
           lat = case_when(ECO_ID==50602 ~ st_coordinates(.)[,2] - 0,
                           ECO_ID==50605 ~ st_coordinates(.)[,2] - 200000,
                           ECO_ID==50606 ~ st_coordinates(.)[,2] + 100000,
                           ECO_ID==50607 ~ st_coordinates(.)[,2] + 200000,
                           ECO_ID==50608 ~ st_coordinates(.)[,2] - 50000,
                           ECO_ID==50609 ~ st_coordinates(.)[,2] - 200000,
                           ECO_ID==50610 ~ st_coordinates(.)[,2] - 100000,
                           ECO_ID==50612 ~ st_coordinates(.)[,2] + 0,
                           ECO_ID==50613 ~ st_coordinates(.)[,2] - 0,
                           ECO_ID==50614 ~ st_coordinates(.)[,2] + 200000,
                           ECO_ID==50616 ~ st_coordinates(.)[,2] + 0,
                           ECO_ID==51111 ~ st_coordinates(.)[,2] - 0)) %>%
    st_drop_geometry() %>%
    st_as_sf(coords = c("lon", "lat")) %>%
    st_set_crs(st_crs(ecoregions))
    
# Create points for label locations
sf_labels <- sf_centroids %>%
    mutate(lon = case_when(ECO_ID==50602 ~ st_coordinates(.)[,1] + 500000,
                           ECO_ID==50605 ~ st_coordinates(.)[,1] + 400000,
                           ECO_ID==50606 ~ st_coordinates(.)[,1] - 400000,
                           ECO_ID==50607 ~ st_coordinates(.)[,1],# - 500000,
                           ECO_ID==50608 ~ st_coordinates(.)[,1] - 600000,
                           ECO_ID==50609 ~ st_coordinates(.)[,1],# - 500000,
                           ECO_ID==50610 ~ st_coordinates(.)[,1] - 550000,
                           ECO_ID==50612 ~ st_coordinates(.)[,1] + 100000,
                           ECO_ID==50613 ~ st_coordinates(.)[,1] - 800000,
                           ECO_ID==50614 ~ st_coordinates(.)[,1] + 100000,
                           ECO_ID==50616 ~ st_coordinates(.)[,1] - 300000,
                           ECO_ID==51111 ~ st_coordinates(.)[,1] - 800000),
           lat = case_when(ECO_ID==50602 ~ st_coordinates(.)[,2] - 400000,
                           ECO_ID==50605 ~ st_coordinates(.)[,2] - 500000,
                           ECO_ID==50606 ~ st_coordinates(.)[,2] + 600000,
                           ECO_ID==50607 ~ st_coordinates(.)[,2] + 500000,
                           ECO_ID==50608 ~ st_coordinates(.)[,2] - 500000,
                           ECO_ID==50609 ~ st_coordinates(.)[,2] - 650000,
                           ECO_ID==50610 ~ st_coordinates(.)[,2] - 550000,
                           ECO_ID==50612 ~ st_coordinates(.)[,2] + 800000,
                           ECO_ID==50613 ~ st_coordinates(.)[,2] - 300000,
                           ECO_ID==50614 ~ st_coordinates(.)[,2] + 550000,
                           ECO_ID==50616 ~ st_coordinates(.)[,2] + 900000,
                           ECO_ID==51111 ~ st_coordinates(.)[,2] - 500000)) %>%
    st_drop_geometry() %>%
    st_as_sf(coords = c("lon", "lat")) %>%
    st_set_crs(st_crs(sf_centroids))
# Create lines between label locations and centroids
sf_lines <- rbind(sf_labels[,"ECO_NAME"], 
                  sf_centroids[,"ECO_NAME"]) %>%
    group_by(ECO_NAME) %>%
    summarise(geometry = st_union(geometry)) %>%
    st_cast("LINESTRING")
# Create labels
labels1 <- paste0(ecoregions$ECO_NAME,"\n",ecoregions$ECO_ID)
# labels2 <- ecoregions$ECO_ID


# Tree cover
study_area_tree <- ggplot() +
    geom_sf(data = ak_can, fill = "white", size=0.5) +
    geom_tile(data = tree_cover_df, aes(x, y, fill = perc_tree),
              show.legend=T) +
    scale_fill_gradient(name="Tree Cover (%)", low="white", high="black") +
    geom_sf(data = ecoregions, fill=NA, color="black", size=0.5) + #fill=NA,
    geom_sf(data = sf_lines, linewidth = 0.5) +
    geom_sf(data = sf_labels, shape = 21, fill = "white", size = 0.5) +
    geom_sf(data = sf_centroids, size = 0.75) +
    geom_sf_label(data = sf_labels,
                  aes(label = labels1),
                  size = 3, fontface="italic",#vjust = -0.5,
                  fill="white", alpha=0.8,
                  fun.geometry = st_centroid,
                  colour = "black") +
    coord_sf(clip = "off") +
    theme_void() +
    theme(legend.position = "right", 
          legend.direction = "vertical")
ggsave(file.path("figures","study_area_treecover.png"),
       plot = study_area_tree, width = 12, height = 8, dpi = 300)

# Fire by historical year
study_area_fire <- ggplot() +
    geom_sf(data = ak_can, fill = "white", size=0.5) +
    geom_sf(data = historical_fires,# %>% mutate(decade=ifelse(FIREYR==2020,2010,decade)), 
            aes(fill=decade),
            col=NA, show.legend=T) +
    scale_fill_binned(name = "Historical fires (decade)", 
                      palette = "OrRd", #n.breaks=6, 
                      breaks=c(1940,1960,1980,2000,2020)) +
    guides(fill = guide_coloursteps(show.limits = TRUE)) +
    geom_sf(data = ecoregions, fill=NA, color="black", size=0.5) +
    # geom_sf(data = eco_50607, fill=NA, color="#1c4587", size=0.6) +
    # geom_sf(data = eco_50609, fill=NA, color="#783f04", size=0.6) +
    geom_sf(data = sf_lines, linewidth = 0.5) +
    geom_sf(data = sf_labels, shape = 21, fill = "white", size = 0.5) +
    geom_sf(data = sf_centroids, size = 0.75) +
    geom_sf_label(data = sf_labels,
                  aes(label = labels1),
                  size = 3, fontface="italic",#vjust = -0.5,
                  fill="white", alpha=0.8,
                  fun.geometry = st_centroid,
                  colour = "black") +
    coord_sf(clip = "off") +
    theme_void() +
    theme(legend.title = element_text(face = "bold"),
          legend.position = "right",
          legend.direction = "vertical")
ggsave(file.path("figures","study_area_fire.png"),
       plot = study_area_fire, width = 12, height = 8, dpi = 300)

# Fire by count of overlapping historical fires 
# use the tree‑cover raster as template (same extent, resolution, CRS)
cnt_rast <- rast(tree_cover_ecos_mask)          # copy template
values(cnt_rast) <- 0                           # initialise with 0

# Rasterise each fire polygon, adding 1 for every overlap
fire_vect <- vect(historical_fires)      # fire polygons
cnt_rast <- rasterize(
    fire_vect,
    cnt_rast,
    field = 1,                     # each polygon contributes 1
    fun   = "sum",                 # sum → number of overlapping polygons
    background = 0
)

cnt_rast <- mask(cnt_rast, vect(ecoregions))

# --------------------------------------------------------------
# 2. Convert raster → data.frame (numeric fire count)
# --------------------------------------------------------------
cnt_df <- as.data.frame(cnt_rast, xy = TRUE) %>% 
    rename(fire_cnt = 3)               # columns: x, y, fire_cnt (numeric)

pal_cont <- colorRampPalette(
    colors = c("#FFFFFF", "#FFCCCC", "#FF6666", "#FF0000", "#990000", "#660000"),
    space = "Lab"
)(100)   # 100 steps – enough for a smooth gradient

x_limits <- c(-170, -100)  # Longitude limits (West to East)
y_limits <- c(30, 50) 
fire_count_hex <- ggplot(cnt_df) +
    geom_sf(data = ak_can,
            fill = "white",
            colour = "grey60",
            size = 0.3) +
#   * z = fire_cnt → numeric, needed for max()
#   * after_stat(y) returns the *max* fire count for each hex
    stat_summary_hex( aes(x = x, y = y, z = fire_cnt),
                      fun   = max, 
                      bins  = 70,   # larger → fewer hexes
                      colour = NA ) + 
    scale_fill_gradientn(name = "# of overlapping historical fires",
                         colours = pal_cont,
                         #limits = c(0, 8), #set upper limit to match the cap; remove if no cap
                         oob = scales::squish, #values outside limits are squished to the nearest colour
                         na.value = "transparent" ) +
    geom_sf(data = ecoregions,
            fill = NA,
            colour = "black",
            size = 0.3) +
    geom_sf(data = sf_lines, linewidth = 0.5) +
    geom_sf(data = sf_labels, shape = 21, fill = "white", size = 0.5) +
    geom_sf(data = sf_centroids, size = 0.75) +
    geom_sf_label(data = sf_labels,
                  aes(label = labels1),
                  size = 3, fontface="italic",#vjust = -0.5,
                  fill="white", alpha=0.8,
                  fun.geometry = st_centroid,
                  colour = "black") +
    coord_sf(crs = ak_can_aea_crs) +
    guides(fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 12,
        barheight = 0.4 )) +
    theme_void() +
    theme(
        legend.position = "bottom",
        legend.title    = element_text(size = 10),
        legend.text     = element_text(size = 9),
        plot.margin = margin(5, 5, 5, 5) )
print(fire_count_hex)
ggsave(file.path("figures","study_area_fire_count.png"),
       plot = fire_count_hex, width = 12, height = 8, dpi = 300)


# Plot for CASC Staff meeting - mark 20607 and 50609
plot_07_09 <- ggplot() +
    geom_sf(data = ak_can, fill = "white", linewidth=0.5) +
    geom_sf(data = ecoregions, fill = NA, color = "black", linewidth=0.75) + # Base ecoregions
    geom_sf(data = ecoregions %>% filter(ECO_ID == 50607), fill = "#674ea7", color = "#674ea7") + # Eco 40507
    geom_sf(data = ecoregions %>% filter(ECO_ID == 50609), fill = "#38761d", color = "#38761d") + # Eco 50609
    theme_void() 
ggsave(file.path("figures","ecos_50607_50609_forpres.png"),
       plot = plot_07_09, width = 12, height = 8, dpi = 300)
plot_07_09

### Plot
# source: https://stackoverflow.com/questions/78378959/adding-labels-to-maps-with-lines
# study_area <- ggplot() +
#     geom_sf(data = ak_can, fill = "white", size=0.5) +
#     geom_tile(data = tree_cover_df, aes(x, y, fill = perc_tree),
#                 show.legend=T) +
#     scale_fill_gradient(name="Tree Cover (%)", low="white", high="black") +
#     new_scale_fill() +
#     geom_sf(data = historical_fires, aes(fill=decade),
#             col=NA, show.legend=T) +
#     scale_fill_gradient(name="Historical fires", low="goldenrod1", high="firebrick") +
#     geom_sf(data = ecoregions, fill=NA, color="black", size=0.5) + #fill=NA,
#     geom_sf(data = sf_lines, linewidth = 0.5) +
#     geom_sf(data = sf_labels, shape = 21, fill = "white", size = 0.5) +
#     geom_sf(data = sf_centroids, size = 0.75) +
#     geom_sf_label(data = sf_labels,
#                  aes(label = labels1),
#                  size = 3, fontface="italic",#vjust = -0.5,
#                  fill="white", alpha=0.8,
#                  fun.geometry = st_centroid,
#                  colour = "black") +
#     # geom_sf_label(data = sf_labels,
#     #              aes(label = labels2),
#     #              size = 3,
#     #              vjust = -1,
#     #              fun.geometry = st_centroid,
#     #              colour = "black") +
#     # geom_sf_label(data = st_point_on_surface(ecoregions), 
#     #               aes(label=str_wrap(ECO_NAME, width=21)),
#     #              fill="white", alpha=0.8,
#     #               #check_overlap=T,
#     #              size = 3, color="black", fontface="italic") +
#     coord_sf(clip = "off") +
#     theme_void() +
#     theme(legend.position = "right", 
#           legend.direction = "vertical")
# ggsave(file.path("figures","study_area.png"),
#        plot = study_area, width = 12, height = 8, dpi = 300)
