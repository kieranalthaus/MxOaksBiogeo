# Load packages
library(sf)
library(rmapshaper)
library(ggplot2)
sf_use_s2(F)

# Load data
countries = st_read(dsn = "MxOaks_Prep/data/countries.shp") %>% 
  st_transform(crs = 5070)
bkg = maps::map(database = "world", 
                regions = c("Canada", "Venezuela", "Ecuador", "Peru", "Brazil"),
                plot = FALSE,
                fill = TRUE) %>% 
  st_as_sf(crs = 5070)
cuba = maps::map(database = "world", regions = "Cuba", plot = FALSE, fill = TRUE) %>% 
  st_as_sf() %>% 
  st_transform(crs = 5070) %>% 
  rename(geometry = geom)
countries = rbind(countries, cuba)
regions = st_read(dsn = "MxOaks_DEC/data/DEC_7_bioregions.shp") %>% 
  st_transform(crs = 5070)

# Edit region data
states = regions %>%
  group_by(id) %>% 
  summarize(geometry = st_union(geometry))

states = st_intersection(x = st_make_valid(states), 
                y = st_make_valid(countries))

# Create lat/long data points
points <- data.frame(name = c("40°", "20°", "0°"),
                     x = -68,
                     y = c(40,20,0)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) 

# Get bounding box
bbox = st_bbox(countries)

# Plot map
map = ggplot()+
  # Plot shapes
  geom_sf(data = states,
          aes(fill = as.factor(id)),
          color = "transparent",
          alpha = 0.6)+
  geom_sf(data = countries, 
          fill = "transparent",
          color = "black",
          linewidth = 0.75)+
  geom_sf(data = bkg,
          fill = "transparent",
          color = "black",
          linewidth = 0.5)+
  
  # Plot text
  geom_sf_label(data = points, 
               mapping = aes(label = name),
               size = 10, 
               fill = "transparent",
               label.size = 0,
               angle = c(18,15,10))+
  
  # Change colors
  scale_fill_manual(values =   c("#0D0887FF",
                                 "#E16462FF", 
                                 "#F0F921FF", 
                                 "#6A00A8FF", 
                                 "#B12A90FF", 
                                 "#FCA636FF"),
                    labels = c("Tropics",
                               "Sierra Madre Oriental",
                               "Sierra Madre Occidental",
                               "California Floristic Provicne",
                               "Trans-Mexican Volcanic Belt",
                               "Eastern North America"
                               ))+

  
  # Bounding
  xlim(bbox[1], bbox[3])+
  ylim(bbox[2], bbox[4])+
  
  # Theme 
  theme_bw(base_rect_size = 1, 
           base_size = 3,
           base_line_size = 50)+
  theme(
    # Remove axis labels and ticks 
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    
    # Edit panel
    panel.grid = element_line(colour = "grey",
                              linewidth = 1,
                              linetype = 2),
    
    # Remove legend
    legend.position = "none",
    # legend.position = c(0.2,0.2),
    legend.key.size = unit(1,"cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)

  )+
  coord_sf(); map

  # Save
# ggsave(plot = map,
#        filename = "MX_oaks_DEC/OUT/20240922_studyarea_map.pdf",
#        height = 4.84*2,
#        width = 5.15*2)








