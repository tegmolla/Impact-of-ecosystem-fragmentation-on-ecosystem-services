
#** Title: Modelling the Impact of Ecosystem Fragmentation on Ecosystem Services in the Degraded Ethiopian Highlands *#
#** Sitotaw, T.M., Willemen, L., Meshesha, D.T., Weldemichael, M.,  & Nelson, A., 2024.  *
#** GAMs model - Nature-based tourism model based on social media geotagged data  *
#** predict NATURE visitation on the basis of ‘reasonable’ predictor variables*
##****************************************************************************
# clear workspace
rm(list=ls())

## Libraries
library(terra)
library(sf)
library(ggplot2)
library(dplyr)

#****************************************************************************
#**                      PART I                                             *
#*
#** Predictors extraction using using the center location of 1x1 km grid    *
#** Predictor variables: natural and cultural sites,                        *
#** road proximity, Proximity, line-of-sight  *
#***************************************************************************

# Set working directory
setwd("D:/ES_ECOINF/nature_tourism_data")

# Import geotagged photo data
visit_data <- read.csv("pud_geotagged_data_2020.csv")

## convert the plot table to points to extract data from raster layers
prj4str <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
visit_data_shp <- st_as_sf(visit_data, coords = c("X", "Y"), crs = prj4str, agr = "constant")

# Export Spatial Points
output_file <- "visit_data.shp"
if (file.exists(output_file)) {file.remove(output_file)}
st_write(visit_data_shp, output_file)

#***************************************************************
#** Load raster predictors for the year 2020
area <- rast("area_pud_2020.tif")
peri_area <- rast("peri_area_pud_2020.tif")
prox <- rast("prox_pud_2020.tif")
cult_sites <- rast("nat_cult_counts_pud_2020.tif")
facility <- rast("facilities_counts_pud_2020.tif")
road_prox <- rast("road_prox_pud_2020.tif")
line_sight <- rast("line_of_sight_pud_2020.tif")


# Combine raster layers
raslist <- list(area, peri_area, prox, cult_sites, facility, road_prox, line_sight)
aligned_rasters <- lapply(raslist, function(r) {
  terra::resample(r, area, method = "bilinear")})

# Stack the rasters
raster_stack <- do.call(c, aligned_rasters)

## Load the geotagged photos point shapefile and check geometry
points <- st_read("visit_data.shp")
if (is.null(st_geometry(points))) {
  stop("Geometry column is not correct.")}

## Convert and and point geometry to X and Y coordinates
coords <- st_coordinates(points)
points <- cbind(points, coords)

## Buffer distance (walking distance of tourists by foot from roads)
buffer_distance <- 500

# Convert to terra object and create buffer
points_terra <- vect(points)
points_buffer <- buffer(points_terra, width = buffer_distance)

# Extract mean raster values
extracted_values <- extract(raster_stack, points_buffer, fun = mean, na.rm = TRUE)

extracted_df <- as.data.frame(extracted_values)
points_df <- as.data.frame(points)

# Combine
combined_df <- cbind(points_df[, c("fid", "X", "Y", "pud")], extracted_df[2:8])  # Append
# Rename columns
colnames(combined_df) <- c("fid", "X", "Y", "pud", "area", "peri_area", "prox", "cult_sites", "facility", "road_prox", "line_sight")

# Save
write.csv(combined_df, "visit_predictors_data_2020.csv", row.names = FALSE)


#****************************************************************************
#**                      PART II                                            *
#*
#**  GAMs Analysis of nature-based tourism serivce using mgcv package       *
#****************************************************************************
# Load the data
visit_data.2020 <- read.csv("visit_predictors_data_2020.csv", stringsAsFactors=FALSE)
# Replace NA using replace() & is.na()
visit_data.2020 <- replace(visit_data.2020, is.na(visit_data.2020), 0)
head(visit_data.2020)

#*******************************************************************
#** Normalize Data columns with Min-Max Scaling*
min_max_norm <- function(x) {(x - min(x)) / (max(x) - min(x))}
visit_data.2020 <- as.data.frame(lapply(visit_data.2020[5:11], min_max_norm))
head(visit_data.2020)

# Combine
visit_data.2020 <- cbind(visit_data[, c("fid", "X", "Y","pud")], visit_data.2020)
names(visit_data.2020) <- c("fid", "X", "Y", "pud", "norm_area", "norm_peri_area", "norm_prox", "norm_cult_sites", "norm_facility", "norm_road_prox", "norm_line_sight")
head(visit_data.2020)


#*******************************************************************
#** GAMs analysis*
library(mgcv)

# GAMs model
visits_mod2020 <- gam(pud ~ s(norm_area, k = 5, bs = 'cr')+
                     s(norm_peri_area, k = 5, bs = 'cr')+
                     s(norm_prox, k=5, bs='cr')+
                     s(norm_cult_sites, k=5, bs='cr')+
                     s(norm_facility, k=5, bs='cr')+
                     s(norm_road_prox, k=5, bs='cr')+
                     s(norm_line_sight, k=5, bs='cr')+
                     ti(X, Y, k=5, bs='cr'),
                   data = visit_data.2020, family = quasipoisson(link = "log"),
                   method = "REML", correlation = corAR1())

## Summary gam model
summary(visits_mod2020)


#**************************************************************************
visit_mod2020_2 <- gam(pud ~ s(norm_area, k = 5, bs = 'cr')+
                     s(norm_peri_area, k = 5, bs = 'cr')+
                     s(norm_prox, k=5, bs='cr')+
                     s(norm_cult_sites, k=5, bs='cr')+
                     s(norm_facility, k=5, bs='cr')+
                     s(norm_road_prox, k=5, bs='cr')+
                     s(norm_line_sight, k=5, bs='cr'),
                   data = visit_data.2020, family = quasipoisson(link = "log"),
                   method = "REML", correlation = corAR1())

## Summary from the gam model
summary(visit_mod2020_2)

######################################################################

#**              Figure 1. area                                      *
#*********************************************************************
#** Nature tourist visitation rate decay due to natural and cultural landscape area variation*
#** GAMs plots using new ecosystem area data generated *
tourism_Data1 <- data.frame(norm_area = seq(from = 0, to = 1, by = 0.02),
                            norm_peri_area = mean(visit_data.2020$norm_peri_area),
                            norm_prox = mean(visit_data.2020$norm_prox),
                            norm_cult_sites = mean(visit_data.2020$norm_cult_sites),
                            norm_facility = mean(visit_data.2020$norm_facility),
                            norm_road_prox = mean(visit_data.2020$norm_road_prox),
                            norm_line_sight = mean(visit_data.2020$norm_line_sight),
                            fid = "1")

## Predict dependent variable results based on GAMs
G1 <- predict(visit_mod2020_2, newdata=tourism_Data1, type="link", se=TRUE)
F1 <- (-G1$fit)
FSEUP1 <- (-G1$fit + 1.96 * G1$se.fit)
FSELOW1 <- (-G1$fit - 1.96 * G1$se.fit)
#plot model results
png(file="Output/Pud_plot_peri_area.png",
   width=10, height=8, units="in", res=300)
plot(tourism_Data1$norm_area, F1, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim = range(c(FSELOW1, FSEUP1)), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#order values for confidence intervals
i.for <- order(tourism_Data1$norm_area)
i.back <- order(tourism_Data1$norm_area, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(tourism_Data1$norm_area[i.for], tourism_Data1$norm_area[i.back])
y.polygon <- c(FSEUP1[i.for], FSELOW1[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(visit_data.2020$norm_area), visit_data.2020$pud, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(tourism_Data1$norm_area, F1, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Area (normalised)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("PUD")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1.0), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1.0), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

#####################################################################

#**              Figure 2. perimeter-area ratio                    *
#********************************************************************
tourism_Data2 <- data.frame(norm_area = mean(visit_data.2020$norm_area),
                            norm_peri_area = seq(from = 0, to = 1, by = 0.02),
                            norm_prox = mean(visit_data.2020$norm_prox),
                            norm_cult_sites = mean(visit_data.2020$norm_cult_sites),
                            norm_facility = mean(visit_data.2020$norm_facility),
                            norm_road_prox = mean(visit_data.2020$norm_road_prox),
                            norm_line_sight = mean(visit_data.2020$norm_line_sight),
                            fid = "1")

## Predict dependent variable results based on GAMs
G2 <- predict(visit_mod2020_2, newdata=tourism_Data2, type="link", se=TRUE)
F2 <- 1/(1+exp(-G2$fit))
FSEUP2 <- 1/(1+exp(-G2$fit + 1.96 * G2$se.fit))
FSELOW2 <- 1/(1+exp(-G2$fit - 1.96 * G2$se.fit))
#plot model results
png(file="Output/Pud_plot_peri_area.png",
    width=10, height=8, units="in", res=300)
plot(tourism_Data2$norm_peri_area, F2, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim = range(c(FSELOW2, FSEUP2)), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#order values for confidence intervals
i.for <- order(tourism_Data2$norm_peri_area)
i.back <- order(tourism_Data2$norm_peri_area, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(tourism_Data2$norm_peri_area[i.for], tourism_Data2$norm_peri_area[i.back])
y.polygon <- c(FSEUP2[i.for], FSELOW2[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(visit_data.2020$norm_peri_area), visit_data.2020$pud, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(tourism_Data2$norm_peri_area, F2, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Perimeter area ratio (normalised)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("PUD")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

#####################################################################

#**        Figure 3. natural and cultural landscape proximity       *
#********************************************************************
tourism_Data3 <- data.frame(norm_area = mean(visit_data.2020$norm_area),
                            norm_peri_area = mean(visit_data.2020$norm_peri_area),
                            norm_prox = seq(from = 0, to = 1, by = 0.02),
                            norm_cult_sites = mean(visit_data.2020$norm_cult_sites),
                            norm_facility = mean(visit_data.2020$norm_facility),
                            norm_road_prox = mean(visit_data.2020$norm_road_prox),
                            norm_line_sight = mean(visit_data.2020$norm_line_sight),
                            fid = "1")

## Predict dependent variable results based on GAMs
G3 <- predict(visit_mod2020_2, newdata=tourism_Data3, type="link", se=TRUE)
F3 <- 1/(1+exp(-G3$fit))
FSEUP3 <- 1/(1+exp(-G3$fit + 1.96 * G3$se.fit))
FSELOW3 <- 1/(1+exp(-G3$fit - 1.96 * G3$se.fit))
# plot model results
png(file="Output/Pud_plot_prox.png",
    width=10, height=8, units="in", res=300)
plot(tourism_Data3$norm_prox, F3, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim = range(c(FSELOW3, FSEUP3)), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
# Confidence intervals
i.for <- order(tourism_Data3$norm_prox)
i.back <- order(tourism_Data3$norm_prox, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(tourism_Data3$norm_prox[i.for], tourism_Data3$norm_prox[i.back])
y.polygon <- c(FSEUP3[i.for], FSELOW3[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(visit_data.2020$norm_prox), visit_data.2020$pud, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(tourism_Data3$norm_prox, F3, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Proximity (normalised)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("PUD")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

#****************************************************************************
#**                      PART III                                           *
#*
#** Nature visitation rate prediction in space and time                     *                                                     * 

#****************************************************************************
#** Load raster predictors for the year 2020
area20 <- rast("area_pud_2020.tif")
peri_area20 <- rast("peri_area_pud_2020.tif")
prox20 <- rast("prox_pud_2020.tif")
cult_sites20 <- rast("nat_cult_counts_pud_2020.tif")
facility20 <- rast("facilities_counts_pud_2020.tif")
road_prox20 <- rast("road_prox_pud_2020.tif")
line_sight20 <- rast("line_of_sight_pud_2020.tif")

# Combine raster layers
raslist <- list(area20, peri_area20, prox20, cult_sites20, facility20, road_prox20, line_sight20)

# Normalising the raster datasets to 0 and 1
normalize_raster <- function(ras) {
  ras_min <- minmax(ras)[1]
  ras_max <- minmax(ras)[2]
  (ras - ras_min) / (ras_max - ras_min)
}

# Resample and normalise predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resampled_ras <- resample(crop_ras, area20, method = "bilinear")
  normalize_raster(resampled_ras)
})

# Stack raster layers and rename them to make it similar to the GAMs model names
pud_preds_stack20 <- rast(pred_rasters)
names(pud_preds_stack20) <- c("area", "peri_area", "prox", "cult_sites", "facility", "road_prox", "line_sight")

#** Spatial prediction (scaling-up) for the whole Lake Tana basin for 2020*
pud_pred2020 <- terra::predict(pud_preds_stack20, model = visit_mod2020_2, type = "response")
# Plot
plot(pud_pred2020, col = terrain.colors(100), main = "PUD Prediction Map 2020")
# Export predicted map
filename20 <- file.path("Output", "PUD map 2020.tif")
writeRaster(pud_pred2020, filename, overwrite = TRUE)


##***********************************************************
##** Load spatial predictor rasters of 2000*
area00 <- rast("area_pud_2020.tif")
peri_area20 <- rast("peri_area_pud_2020.tif")
prox20 <- rast("prox_pud_2020.tif")
cult_sites20 <- rast("nat_cult_counts_pud_2020.tif")
facility20 <- rast("facilities_counts_pud_2020.tif")
road_prox20 <- rast("road_prox_pud_2020.tif")
line_sight20 <- rast("line_of_sight_pud_2020.tif")

# Combine raster layers
raslist <- list(area00, peri_area00, prox00, cult_sites00, facility00, road_prox00, line_sight00)

# Normalising the raster datasets to 0 and 1
normalize_raster <- function(ras) {
  ras_min <- minmax(ras)[1]
  ras_max <- minmax(ras)[2]
  (ras - ras_min) / (ras_max - ras_min)
}

# Resample and normalise predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resampled_ras <- resample(crop_ras, area00, method = "bilinear")
  normalize_raster(resampled_ras)
})

# Stack raster layers and rename them to make it similar to the GAMs model names
pud_preds_stack00 <- rast(pred_rasters)
names(pud_preds_stack00) <- c("area", "peri_area", "prox", "cult_sites", "facility", "road_prox", "line_sight")

#** Temporal prediction (scaling-up) for the whole Lake Tana basin for 2000*
pud_pred2000 <- terra::predict(pud_preds_stack00, model = visit_mod2020_2, type = "response")
# Plot
plot(pud_pred2000, col = terrain.colors(100), main = "PUD Prediction Map 2000")
# Export predicted map
filename20 <- file.path("Output", "PUD map 2000.tif")
writeRaster(pud_pred2000, filename, overwrite = TRUE)


##****************************************************************************
##**                      PART IV                                            *
##*
##**  Nature visitation contour plot using vis.gam function                                *
##*
##***************************************************************************

#jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

myvis.gam <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                       col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                       plot.type = "persp", zlim = NULL, nCol = 50, ...) 
{
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(paste(c("view variables must be one of", v.names), 
                 collapse = ", "))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(paste("View variables must contain more than one value. view = c(", 
               view[1], ",", view[2], ").", sep = ""))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 3
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    ### customized here
    else if (color == 'jet') {
      pal <- jet.colors(nCol)
      con.col = 1
    }
    ####
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit, n.grid, n.grid)
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
        eval(parse(text = txt))
      }
      else {
        txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        op <- par(bg = "white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
        eval(parse(text = txt))
        par(op)
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      z.max <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      z.min <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
    }
    zlim <- c(z.min, z.max)
    z <- fv$fit - fv$se.fit * se
    z <- matrix(z, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit + se * fv$se.fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
    eval(parse(text = txt))
  }
}


#**********************************************************************
#** Contour plot for nature visitation  *
#*
########################################################################
# GAMs model
visit_mod2020_2 <- gam(pud ~ s(norm_area, k = 5, bs = 'cr')+
                     s(norm_peri_area, k = 5, bs = 'cr')+
                     s(norm_prox, k=5, bs='cr')+
                     s(norm_cult_sites, k=5, bs='cr')+
                     s(norm_facility, k=5, bs='cr')+
                     s(norm_road_prox, k=5, bs='cr')+
                     s(norm_line_sight, k=5, bs='cr')+
                     ti(X, Y, k=5, bs='cr'),
                   data = visit_data.2020, family = gaussian(link = "identity"),
                   method = "REML", correlation = corAR1())

## Summary gam model
summary(visit_mod2020_2)

# colors palette definition
jet.colors <- colorRampPalette(c("#fe9929","#ffffb2", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"))
## Saving the plot output
#png(file="Output/Pud_plot_2023_2.png",
#    width=10, height=8, units="in", res=300)
myvis.gam(visit_mod2020_2, view=c("X", "Y"), color = "jet",
          type = "response", xlab="", ylab="",
          main = "", cex.main=1.5, font.main=1,
          plot.type="contour", n.grid=600, too.far=0.5,
          ylim = c(min(visit_data.2020$Y),max(visit_data.2020$Y)), xlim = c(min(visit_data.2020$X),max(visit_data.2020$X)),
          cex.lab=1.5, cex.axis=1.5, las=1, yaxt="n", xaxt="n",
          contour.col="black", labcex=1.3, lwd=1.5, nlevels = 30)
## Create box bound
box(which="plot", lwd = 2)
## Plot y-axis for latitude
axis(side=2, las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
## Plot x-axis for longitude
axis(1, cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)
## Axis and figure labels
mtext("UTM Easting", side=1, line=3, cex=1.5)
mtext("UTM Northing", side=2, line=2.5, cex=1.5)

## Finally, Run dev.off() to create the file!
#dev.off()




