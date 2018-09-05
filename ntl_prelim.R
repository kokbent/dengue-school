ntl_ras <- ntl %>%
  dplyr::select(long, lat, ntl) %>%
  rasterFromXYZ(crs = CRS("+init=epsg:4326"))
plot(ntl_ras)

mex <- shapefile("data/shp/gadm36_MEX_1.shp")
yuc <- mex[31,]
plot(yuc)
plot(ntl_ras, add=T)

# Prelim testing of ntl raster and shapes
library(mgcv)
test <- yuc@polygons[[1]]@Polygons[[16]] %>%
  list() %>%
  Polygons(ID=1) %>%
  list() %>%
  SpatialPolygons(proj4string = CRS("+init=epsg:4326"))
plot(ntl_ras)
plot(test, add=T)

test2 <- list()
for (i in 1:16) {
  test2[[i]] <- yuc@polygons[[1]]@Polygons[[i]] %>%
    list() %>%
    Polygons(ID=i)
}
test3 <- test2 %>% SpatialPolygons(proj4string = CRS("+init=epsg:4326"))
plot(test3)
text(coordinates(test3)[c(11, 16),], labels=c(11, 16), col="red")