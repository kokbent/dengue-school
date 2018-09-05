library(rgeos)
library(mgcv)

ntl <- read_tsv("data/all_yucatan_pixels.out", 
                col_names = c("long", "lat", "ntl", "mun_id", "mun"))
ntl$LOCALIDAD <- NA
ntl$MUNICIPIO <- NA
ntl$URB <- NA

mun_shp <- shapefile("data/shp/municipal.shp") %>%
  spTransform(CRS("+init=epsg:4326"))

# For each point in ntl, assign municipality
ntl_mat <- ntl[,c("long", "lat")] %>% as.matrix()

lapply(mun_shp@polygons, function (x) length(x@Polygons)) %>% unlist()

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
  ntl$URB[cond] <- "U"
}

ntl$URB[is.na(ntl$URB)] <- "R"

# Convert accented characters
ntl$LOCALIDAD <- stri_trans_general(ntl$LOCALIDAD, "latin-ascii") %>% tolower()
ntl$MUNICIPIO <- stri_trans_general(ntl$MUNICIPIO, "latin-ascii") %>% tolower()
# write_rds(ntl, "data/ntl_yucatan_augment.rds")

# Schools with exactly one shape
ids <-  shape_per_sch %>%
  filter(TOT == 1) %>%
  .$ID
dat_inc_1sp <- dat_loc[dat_loc$ID %in% ids,]
dat_inc_1sp

ntl_mat_wgs <- SpatialPoints(ntl_mat, CRS("+init=epsg:4326"))
ntl_mat_lcc <- spTransform(ntl_mat_wgs, crs(loc_rur_lcc))
samppoint_mat <- matrix(NA, nrow = nrow(dat_inc_1sp), ncol = 3)

# Sampling criteria
# For urban, points must be in urban polygons

for (i in 1:nrow(dat_inc_1sp)) {
  cvegeo <- dat_inc_1sp$CVEGEO[i]
  if(dat_inc_1sp$URB[i] == "U") {
    cond <- ntl$LOCALIDAD == dat_inc_1sp$LOCALIDAD[i] & 
      ntl$MUNICIPIO == dat_inc_1sp$MUNICIPIO[i]
    cond <- ifelse(is.na(cond), F, cond)
    ntl_sub <- ntl[cond,] %>%
      mutate(prob = ntl/sum(ntl)) %>%
      as.data.frame()
    samp <- sample(1:nrow(ntl_sub), size = 1, prob = ntl_sub$prob)
    samppoint_mat[i,] <- ntl_sub[samp,1:3] %>% unlist
    print(c(i, samppoint_mat[i,]))
  } else {
    cond1 <- ntl$MUNICIPIO == dat_inc_1sp$MUNICIPIO[i] &
      ntl$URB == "R"
    cond1 <- ifelse(is.na(cond1), F, cond1)
    ind1 <- which(cond1)
    ntl_sub1 <- ntl[ind1,]
    
    pol <- loc_rur_lcc[loc_rur_lcc$CVEGEO == cvegeo,]
    cond2 <- (gDistance(ntl_mat_lcc[ind1], pol, byid = T) %>% as.vector()) < 5000
    ntl_sub2 <- ntl_sub1[cond2,] %>%
      mutate(prob = ntl/sum(ntl)) %>%
      as.data.frame()
    if (nrow(ntl_sub2) == 0) {
      samppoint_mat[i,] <- NA
    } else {
      samp <- sample(1:nrow(ntl_sub2), size = 1, prob = ntl_sub2$prob)
      samppoint_mat[i,] <- ntl_sub2[samp,1:3] %>% unlist
    }
    print(c(i, samppoint_mat[i,]))
   }
}

dat_inc_1sp[is.na(samppoint_mat[,2]),] %>% print(n=44)
