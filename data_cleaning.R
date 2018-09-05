library(tidyverse)
library(sp)
library(raster)
library(rgdal)
library(stringi)
library(stringdist)

# Functions: Convert to unaccented and lowercase characters
ascii_lower <- function (v) {
  stri_trans_general(v, "latin-ascii") %>% tolower()
}

## Data carpentries
dat <- read_csv("data/addresses_google.csv")
dat$ID <- 1:nrow(dat)

dat_complete <- dat %>%
  filter(LAT_HITS != "ZERO_RESULTS")
dat_incomplete <- dat %>%
  filter(LAT_HITS == "ZERO_RESULTS")

loc_urb_lcc <- shapefile("data/shp/loc_urb.shp")
loc_urb_wgs <- loc_urb_lcc  %>%
  spTransform(CRS("+init=epsg:4326"))
loc_rur_lcc <- shapefile("data/shp/loc_rur.shp")
loc_rur_wgs <- loc_rur_lcc %>%
  spTransform(CRS("+init=epsg:4326"))

## Build a table of muncipal & locality
tmp1 <- loc_urb_lcc@data %>% 
  dplyr::select(CVEGEO, MUNICIPIO = NOM_MUN, LOCALIDAD = NOMBRE) %>%
  mutate(URB = "U")

tmp2 <- loc_rur_lcc@data %>% 
  dplyr::select(CVEGEO, MUNICIPIO = NOM_MUN, LOCALIDAD = NOMBRE) %>%
  mutate(URB = "R")

dat_locnombre <- bind_rows(tmp1, tmp2)
dat_locnombre$MUNICIPIO <- stri_trans_general(dat_locnombre$MUNICIPIO, "latin-ascii") %>% tolower()
dat_locnombre$LOCALIDAD <- stri_trans_general(dat_locnombre$LOCALIDAD, "latin-ascii") %>% tolower()

## Merge with school dataset with municipal & locality
dat_loc <- dat_incomplete %>%
  dplyr::select(ID, LOCALIDAD, MUNICIPIO) %>%
  mutate(LOCALIDAD = tolower(LOCALIDAD),
         MUNICIPIO = tolower(MUNICIPIO)) %>%
  left_join(dat_locnombre)
dat_loc

## Data Cleaning
# 46 Unmatched localities
dat_loc_unmatch <- dat_loc %>%
  filter(is.na(CVEGEO))
dat_loc_unmatch

# see if there are fuzzy matches with names within same municipal
nearest_string <- function (loc, mun) {
  nom <- dat_locnombre$LOCALIDAD[dat_locnombre$MUNICIPIO == mun]
  str_dist_min <- stringdist(loc, nom)
  nearest <- nom[str_dist_min == min(str_dist_min)]
  return(nearest)
}

nearest <- apply(dat_loc_unmatch, 1, function (x) nearest_string(x[2], x[3]))
nearest <- nearest %>% 
  map(~ paste(.x, collapse = "; ")) %>%
  unlist()
cbind(dat_loc_unmatch$LOCALIDAD, dat_loc_unmatch$MUNICIPIO, nearest)

# # Match purely by locality names
# nearest_string2 <- function (loc) {
#   nom <- dat_locnombre$LOCALIDAD
#   str_dist_min <- stringdist(loc, nom)
#   ind <- str_dist_min == min(str_dist_min)
#   nearest <- dat_locnombre[ind,] %>%
#     dplyr::select(LOCALIDAD, MUNICIPIO) %>%
#     distinct()
#   return(nearest)
# }
# 
# nearest <- apply(dat_loc_unmatch, 1, function (x) nearest_string2(x[2]))
# nearest <- nearest %>% 
#   map(~ mutate(.x, NAME = paste(LOCALIDAD, MUNICIPIO, sep = ", "))) %>%
#   map(~ .x$NAME) %>%
#   map(~ paste(.x, collapse = "; ")) %>%
#   unlist()
# nearest_df <- cbind(dat_loc_unmatch$LOCALIDAD, dat_loc_unmatch$MUNICIPIO, nearest)

# Fixing names for cases with obvious spelling/naming mismatch
dat_incomplete_v2 <- dat_incomplete %>%
  dplyr::select(ID, LOCALIDAD, MUNICIPIO) %>%
  mutate(LOCALIDAD = tolower(LOCALIDAD),
         MUNICIPIO = tolower(MUNICIPIO))

cond <- (dat_incomplete_v2$LOCALIDAD =="temozon abala") & (dat_incomplete_v2$MUNICIPIO == "abala")
dat_incomplete_v2$LOCALIDAD[cond] <- "temozon"
cond <- (dat_incomplete_v2$LOCALIDAD =="unidades de riego") & (dat_incomplete_v2$MUNICIPIO == "ticul")
dat_incomplete_v2$LOCALIDAD[cond] <- "ticul [unidad de riego]"
cond <- (dat_incomplete_v2$LOCALIDAD =="el alamo") & (dat_incomplete_v2$MUNICIPIO == "tizimin")
dat_incomplete_v2$LOCALIDAD[cond] <- "alamo"
cond <- (dat_incomplete_v2$LOCALIDAD =="x 'lapak") & (dat_incomplete_v2$MUNICIPIO == "yaxcaba")
dat_incomplete_v2$LOCALIDAD[cond] <- "xlapak"

# Rerun previous codes
dat_loc <- dat_incomplete_v2 %>%
  left_join(dat_locnombre)
dat_loc

# Down to 38 unmatched localities
dat_loc_unmatch <- dat_loc %>%
  filter(is.na(CVEGEO))
dat_loc_unmatch

dat_loc_match <- dat_loc %>%
  filter(!is.na(CVEGEO))

# There are localities with same names but multiple points in the shapefile
# Further complication of some localities matching both Urban and Rural polygons
shape_per_sch <- dat_loc %>%
  filter(!is.na(CVEGEO)) %>%
  group_by(ID, URB) %>%
  tally() %>%
  spread(URB, n, fill = 0) %>%
  mutate(TOT = U + R) %>%
  arrange(desc(TOT))
shape_per_sch %>% print(n=40)
sum(shape_per_sch$TOT > 1) # 53 schools with more than one shape
sum(shape_per_sch$U > 0 & shape_per_sch$R > 0) # 21 schools with both urban and rural shapes
