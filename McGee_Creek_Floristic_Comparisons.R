
### import data
# a list of species present in the study area including synonyms
all_names <- read_excel("Names.xlsx")

# a .csv file downloaded from the idigbio portal
dat <- arrow::open_csv_dataset("north_america_idig_bio.csv") 

### initial data cleaning
dat_trimmed <- dat |> 
  rename(canon_name = 6, sci_name = 70, class = 8, 
         phylum = 66, county = 19, state = 73, 
         tax_status = 75, lat_log = 37, uncertainty = 16) |> # renaming columns
  # keeping only relevant species
  filter(sci_name %in% all_names$SpeciesName) |> 
  # keeping only relevant columns
  select(canon_name, sci_name, class, phylum, 
         county, state, tax_status, lat_log, uncertainty) |> 
  collect() |> # loading the data
  # reformatting geographic coordinates
  separate_wider_delim(cols = lat_log, delim = ",", 
                       names = c("latitude", "longitude")) |> 
  separate_wider_delim(cols = latitude, delim = " ", 
                       names = c("leave", "latitude")) |>
  separate_wider_delim(cols = longitude, delim = ":", 
                       names = c("out", "longitude")) |>
  select (-c(leave, out))
# reformatting geographic coordinates
dat_trimmed$longitude <- str_remove(dat_trimmed$longitude, "\\}")
dat_trimmed$longitude <- str_remove(dat_trimmed$longitude, " ")

### turning the data into and sf object and removing points east of the 100th meridian
dat_sf <- st_as_sf(dat_trimmed, 
                   coords = c("longitude", "latitude"), 
                   crs = 4326) %>% 
  st_crop(test, xmin = -130, ymin = 24, xmax=-100, ymax = 50)


### resolving taxonomy
# using tnrs to find an accepted name for each name in the dataset
tnrs_input <- dat_sf$sci_name
unique_input <- unique(tnrs_input)
tnrs_output <- TNRS(taxonomic_names = unique_input)

# replacing each name in the dataset with an accepted name
name_submitted <- as.vector(tnrs_output$Name_submitted)
name_accepted <- as.vector(tnrs_output$Accepted_name)
dat_synonymized <- dat_sf
library(plyr)
dat_synonymized$sci_name <- mapvalues(dat_sf$sci_name, 
                                      name_submitted, name_accepted)
detach(package:plyr) # plyr will interfere with dply functions if not detached

# a list of species present in the study area, excluding synonyms
species_list <- read_excel("McGee_Species.xlsx") 

# removing names remaining in the dataset that do not occur in the study area
species_list_tnrs <- TNRS(taxonomic_names = species_list)
dat_synonymized_trimmed <- dat_synonymized[dat_synonymized$sci_name 
                                           %in% species_list$Accepted_name, ]

### removing occurrences whose coordinates disagree with the county or state on the label
counties_sf <- us_counties()
county_info <- st_join(dat_synonymized_trimmed, counties_sf)

occurrences_geo_clean <- county_info %>% 
  mutate(duplicates = if_else(stateProvince == state_name & county == name,
                              TRUE,
                              FALSE)) %>%
  filter(duplicates == TRUE)

### separating dataset by group of species and updating crs
# NAD83/Conus Albers is used because units are in meters rather than degrees
dat_all <- st_transform(occurrences_geo_clean, crs = 5070)

dat_vascular <- st_transform(dat_synonymized_trimmed, crs = 5070) |>
  filter(phylum == "tracheophyta") 

dat_bryophyte <- st_transform(dat_synonymized_trimmed, crs = 5070) |>
  filter(phylum == "bryophyta" | phylum == "marchantiophyta") 

# importing a list of introduced species present in the study area
introduced <- read_excel("Introduced.xlsx", col_names = FALSE) 

dat_introduced <- dat_all[dat_all$sci_name %in% introduced$...1, ]


### creating a grid with cells 50 km wide
# importing polygons from a folder containing a .shp file of US states
States <- st_read("States") 
State_sf <- st_as_sf(States, 
                     coords = c("long", "lat"), 
                     crs = 4326) 
State_83 <- st_transform(State_sf, crs = 5070)
tempgrid <- st_make_grid(State_83, cellsize = 50000, 
                         square = FALSE, crs = 5070) 
index <- which(lengths(st_intersects(tempgrid, State_83)) > 0)
tomap <- tempgrid[index]             
grid <- st_as_sf(tomap) 
grid_id <- cbind(id = 1:nrow(grid), grid)

### analysis
# join occurrence layer and grid 
occurrences_grid_joined_all <- st_join(dat_all, grid_id, join = st_within)

# count species per grid cell
hex_panel_counts_all <- occurrences_grid_joined_all |> 
  filter(!is.na(id)) |> 
  group_by(id) |> 
  summarise(unique_species_count = n_distinct(sci_name), .groups = 'drop') |> 
  ungroup() |>
  sf::st_drop_geometry()

# join occurrence layer and grid and count species per grid cell: vascular 
occurrences_grid_joined_vascular <- st_join(dat_vascular, 
                                            grid_id, join = st_within)
hex_panel_counts_vascular <- occurrences_grid_joined_vascular |> 
  filter(!is.na(id)) |> 
  group_by(id) |> 
  summarise(unique_species_count = n_distinct(sci_name), .groups = 'drop') |> 
  ungroup() |>
  sf::st_drop_geometry()

# join occurrence layer and grid and count species per grid cell: bryophyte
occurrences_grid_joined_bryophyte <- st_join(dat_bryophyte, 
                                             grid_id, join = st_within)
hex_panel_counts_bryophyte <- occurrences_grid_joined_bryophyte |> 
  filter(!is.na(id)) |> 
  group_by(id) |> 
  summarise(unique_species_count = n_distinct(sci_name), .groups = 'drop') |> 
  ungroup() |>
  sf::st_drop_geometry()

# join occurrence layer and grid and count species per grid cell: introduced
occurrences_grid_joined_introduced <- st_join(dat_introduced, 
                                              grid_id, join = st_within)
hex_panel_counts_introduced <- occurrences_grid_joined_introduced |> 
  filter(!is.na(id)) |> 
  group_by(id) |> 
  summarise(unique_species_count = n_distinct(sci_name), .groups = 'drop') |> 
  ungroup() |>
  sf::st_drop_geometry()

# join species counts back onto hexagon layer
grid_id_all <- grid_id |> 
  left_join(hex_panel_counts_all, by = "id")
grid_id_vascular <- grid_id |> 
  left_join(hex_panel_counts_vascular, by = "id")
grid_id_bryophyte <- grid_id |> 
  left_join(hex_panel_counts_bryophyte, by = "id")
grid_id_introduced <- grid_id |> 
  left_join(hex_panel_counts_introduced, by = "id")

# count number of specimen records in each grid cell

# creating a trimmed dataset with all taxa included
dat_trimmed_all <- dat |> 
  rename(canon_name = 6, sci_name = 70, class = 8, 
         phylum = 66, county = 19, state = 73, 
         tax_status = 75, lat_log = 37, uncertainty = 16) |> # renaming columns
  # keeping only relevant columns
  select(canon_name, sci_name, class, phylum, 
         county, state, tax_status, lat_log, uncertainty) |> 
  collect() |> # loading the data
  # reformatting geographic coordinates
  separate_wider_delim(cols = lat_log, delim = ",", 
                       names = c("latitude", "longitude")) |> 
  separate_wider_delim(cols = latitude, delim = " ", 
                       names = c("leave", "latitude")) |>
  separate_wider_delim(cols = longitude, delim = ":", 
                       names = c("out", "longitude")) |>
  select (-c(leave, out))

# formatting lat/long data for conversion to sf object
dat_trimmed_all$longitude <- str_remove(dat_trimmed$longitude, "\\}")
dat_trimmed_all$longitude <- str_remove(dat_trimmed$longitude, " ")

# converting to sf object 
dat_sf_removed <- st_as_sf(dat_trimmed_all, 
                   coords = c("longitude", "latitude"), 
                   crs = 4326) 
dat_83_removed <- st_transform(dat_sf_removed, crs = 5070) # updating datum to match

# pairing specimen records with grid
occurrences_grid_joined_all <- st_join(dat_83_removed, grid_id, join = st_within)

# counting number of records per grid cell
hex_panel_counts_all <- occurrences_grid_joined_all |> 
  filter(!is.na(id)) |> 
  group_by(id) |> 
  summarise(count_in_cell = n(), .groups = 'drop') |> 
  ungroup() |>
  sf::st_drop_geometry()

