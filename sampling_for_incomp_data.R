library(rgeos)
library(mgcv)

## Augment the nighttime light dataframe
ntl <- read.csv("data/yuc_ntl.csv")
colnames(ntl)[1:2] <- c("long", "lat")
ntl$LOCALIDAD <- NA
ntl$MUNICIPIO <- NA
ntl$URB <- NA
ntl$CVEGEO <- NA

mun_shp <- shapefile("data/shp/municipal.shp") %>%
  spTransform(CRS("+init=epsg:4326"))

ntl_mat <- ntl[,c("long", "lat")] %>% as.matrix()

# For each point in ntl, assign municipality
for (i in 1:nrow(mun_shp)) {
  coords <- mun_shp@polygons[[i]]@Polygons %>%
    map(~ .x@coords)
  cond <- coords %>%
    map(~ in.out(.x, ntl_mat)) %>%
    unlist() %>%
    matrix(nrow = nrow(ntl_mat))
  cond <- rowSums(cond) > 0
  ntl$MUNICIPIO[cond] <- mun_shp$NOMBRE[i]
}

# For each point in ntl, assign urbanity and urban localities
for (i in 1:nrow(loc_urb_wgs)) {
  coords <- loc_urb_wgs@polygons[[i]]@Polygons[[1]]@coords
  cond <- in.out(coords, ntl_mat)
  ntl$LOCALIDAD[cond] <- loc_urb_wgs$NOMBRE[i]
  ntl$MUNICIPIO[cond] <- loc_urb_wgs$NOM_MUN[i]
  ntl$CVEGEO[cond] <- loc_urb_wgs$CVEGEO[i]
  ntl$URB[cond] <- "U"
}

ntl$URB[is.na(ntl$URB)] <- "R"

# Convert accented characters
ntl$LOCALIDAD <- ascii_lower(ntl$LOCALIDAD)
ntl$MUNICIPIO <- ascii_lower(ntl$MUNICIPIO)
# write_rds(ntl, "data/ntl_yucatan_augment.rds")

## Create spatial points objects so CRS is captured
ntl_mat_wgs <- SpatialPoints(ntl_mat, CRS("+init=epsg:4326"))
ntl_mat_lcc <- spTransform(ntl_mat_wgs, crs(loc_rur_lcc))
samppoint_mat <- matrix(NA, nrow = nrow(dat_inc_1sp), ncol = 3,
                        dimnames = list(NULL, c("long", "lat", "ntl")))

## Sampling matched schools
# Sampling criteria - for urban, points must be in urban polygons
# For rural, (1) Within municipal, (2) Not in urban, (3) 5km from ref point
get_valid_point <- function (cvegeo, urb) {
  if(urb == "U") {
    cond <- ntl$CVEGEO == cvegeo
    cond <- ifelse(is.na(cond), F, cond)
  } else {
    mun <- loc_rur_wgs$NOM_MUN[loc_rur_wgs$CVEGEO == cvegeo] %>%
      ascii_lower()
    cond1 <- ntl$MUNICIPIO == mun & ntl$URB == "R"
    cond1 <- ifelse(is.na(cond1), F, cond1)
    ind1 <- which(cond1)
    
    cond <- rep(F, length(cond1))
    pol <- loc_rur_lcc[loc_rur_lcc$CVEGEO == cvegeo,]
    cond2 <- (gDistance(ntl_mat_lcc[ind1], pol, byid = T) %>% as.vector()) < 5000
    cond[ind1] <- cond2
  }
  
  return(cond)
}

# sample a point from all valid points
sample_pt <- function (pts_df) {
  if(nrow(pts_df) == 0) {
    return(data.frame(long=NA, lat=NA, ntl=NA))
  } else {
    samp <- sample(1:nrow(pts_df), size = 1, prob = pts_df$prob)
    return(pts_df[samp,1:3])
  }
}

set.seed(4326)
ids <- dat_loc %>%
  dplyr::filter(!is.na(CVEGEO)) %>%
  .$ID %>% unique()
n <- length(ids)
samppt_df <- data.frame(ID = ids,
                        long = numeric(n),
                        lat = numeric(n),
                        ntl = numeric(n))
for (i in 1:n) {
  id <- ids[i]
  cvegeo <- dat_loc$CVEGEO[dat_loc$ID == id]
  urb <- dat_loc$URB[dat_loc$ID == id]
  
  cond <- map2(cvegeo, urb, ~ get_valid_point(.x, .y)) %>%
    unlist() %>%
    matrix(nrow = nrow(ntl)) %>%
    rowSums() > 0
  
  ntl_sub <- ntl[cond,] %>%
    mutate(ntl2 = ntl+0.0001) %>%
    mutate(prob = ntl2/sum(ntl2)) %>%
    as.data.frame()
  
  samppt_df[i, 2:4] <- sample_pt(ntl_sub)
  print(samppt_df[i,])
}

# 11 unmatched ones, so it'll be total of 49 unmatched ones
samppt_df_comp <- samppt_df[!is.na(samppt_df$long),]

## Sampling unmatched data
# Build unmatched data
ids1 <- samppt_df[is.na(samppt_df$long),"ID"]
ids2 <- dat_loc_unmatch$ID
ids <- c(ids1, ids2) %>% sort

dat_unmatch <- dat_loc[dat_loc$ID %in% ids,] %>%
  arrange(ID)

# Assume rural for all unmatched ones 
get_valid_pt_unmatched <- function (mun) {
  cond <- ntl$MUNICIPIO == mun & ntl$URB == "R"
  return(cond)
}

n <- length(ids)
samppt_unmatch_df <- data.frame(ID = ids,
                                long = numeric(n),
                                lat = numeric(n),
                                ntl = numeric(n))

samppt_unmatch_df <- dat_unmatch$MUNICIPIO %>%
  map(~ get_valid_pt_unmatched(.x)) %>%
  map(~ ntl[.x,] %>%
        mutate(ntl2 = ntl+0.0001) %>%
        mutate(prob = ntl2/sum(ntl2)) %>%
        as.data.frame()) %>%
  map_df(~ sample_pt(.x))
samppt_unmatch_df$ID = ids

samppt_full <- bind_rows(samppt_df_comp, samppt_unmatch_df)
colSums(is.na(samppt_full))

# Jitter the point by +- 7.5 arsecond.
resolution <- 15/3600
x <- runif(nrow(samppt_full), -resolution/2, resolution/2)
y <- runif(nrow(samppt_full), -resolution/2, resolution/2)
samppt_full$long <- samppt_full$long + x
samppt_full$lat <- samppt_full$lat + y

# Merge sampled coordinates with geocoded coordinates
dat_incomplete_samp <- samppt_full %>%
  left_join(dat_incomplete_v2)

dat_final <- dat_complete %>%
  dplyr::select(ID, long = LNG_HITS, lat = LAT_HITS, LOCALIDAD, MUNICIPIO) %>%
  mutate(long = as.numeric(long), lat = as.numeric(lat),
         LOCALIDAD = tolower(LOCALIDAD), MUNICIPIO = tolower(MUNICIPIO)) %>%
  bind_rows(dat_incomplete_samp) %>%
  arrange(ID)
write_csv(dat_final, "data/school_coordinates_full.csv")

# Some graphics
sch_coords <- cbind(dat_final$long, dat_final$lat) %>%
  SpatialPoints(proj4string = CRS("+init=epsg:4326"))

png("google_results.png", width=1200, height=1200)
plot(mun_shp)
plot(sch_coords[is.na(dat_final$ntl),], add=T, col="blue")
dev.off()

png("rand_samp.png", width=1200, height=1200)
plot(mun_shp)
plot(sch_coords[!is.na(dat_final$ntl),], add=T, col="red")
dev.off()
