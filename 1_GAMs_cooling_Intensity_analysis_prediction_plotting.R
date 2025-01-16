
#** Title: Modelling the Impact of Ecosystem Fragmentation on Ecosystem Services in the Degraded Ethiopian Highlands *#
#** Sitotaw, T. M., Willemen, L., Meshesha, D. T., & Nelson, A., 2024.      *
#** GAMs model for 2020 for cooling intensity model based on field air temperature data             *
#
#**************************************************************************
# clear workspace
rm(list=ls())

# Load libraries
library(terra)
library(sf)
library(ggplot2)
library(dplyr)

#****************************************************************************
#**                      PART I                                             *
#**  Predictors extraction using the locations of air temperature records   *
#**  Predictor variables: fvc, impermeable_surface, perimeter-area ratio, proximity  *
#*
#***************************************************************************

# Set working directory
setwd("D:/ES_ECOINF/heat_stress_regulation_data")

# Cooling intensity data
temp_data_plot <- read.csv("airtemp_data_raw_2020.csv")
head(temp_data_plot)

## convert the plot table to points to extract data from raster layers
prj4str <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
temp_data_shp <- st_as_sf(temp_data_plot, coords = c("X", "Y"), crs = prj4str, agr = "constant")

#** Load predictors*
fvc20 <- rast("fvc_cooling_2020.tif")
imper_surf20 <- rast("Impersurf_cooling_2020.tif")
area20 <- rast("area_cooling_2020.tif")
peri_area20 <- rast("peri_area_cooling_2020.tif")
prox20 <- rast("prox_cooling_2020.tif")

# Read boundary as spatial data
bound <- st_read("D:/ES_Paper_3/Basin_Boundary.json")

# Combine raster layers
raslist <- list(fvc20, imper_surf20, area20, peri_area20, prox20)

## Resample predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resample(crop_ras, fvc20, method = "bilinear")
})

## Stack raster layers
raster_stack <- rast(pred_rasters)

# Define the extraction function within a 15 me buffer around each point
extract.fun <- function(rast, temp_data_shp) {
  result <- extract(rast, temp_data_shp, buffer = 1000, fun = mean, na.rm = TRUE)
  return(result[, -1])  # Extracted values
}

# Apply the extraction function to each raster in the stack
rasextract <- extract.fun(raster_stack, temp_data_shp)

# Assign column names based on raster layer names
layer_names <- names(raster_stack)
colnames(rasextract) <- layer_names

head(temp_data_plot)
## Combine extracted predictor values with plot data
temp_preds_data <- cbind(temp_data_plot[, c("transect_id", "fid", "veg_type", "X", "Y", "distance_m", "ref_temp_max", "temp_oC", "cooling_inten_oC")], rasextract)
head(temp_preds_data)

## Rename
names(temp_preds_data) <- c("transect_id", "fid", "veg_type", "X", "Y", "distance_m", "ref_temp_max", "temp_oC", "cooling_inten_oC", "fvc",
                           "imper_surf", "area", "peri_area", "prox")
head(temp_preds_data)
max(temp_preds_data$imper_surf)

#Export
#write.table(temp_preds_data, file="Output/temp_predictors_data_2020.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


##***************************************************************************
##** GAMs analysis per transect*
library(mgcv)

temp_modeldata <- read.csv("airtemp_data_raw_2020_final.csv", stringsAsFactors=FALSE)
head(temp_modeldata)

# Cleaning No Data values. Change NA values to 0 values.
temp_modeldata[is.na(temp_modeldata[,"fvc"]), "fvc"] <- 0
temp_modeldata[is.na(temp_modeldata[,"imper_surf"]), "imper_surf"] <- 0
temp_modeldata[is.na(temp_modeldata[,"area"]), "area"] <- 0
temp_modeldata[is.na(temp_modeldata[,"peri_area"]), "peri_area"] <- 0
temp_modeldata[is.na(temp_modeldata[,"prox"]), "prox"] <- 0

# Define a function for min-max normalization
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

max(temp_modeldata$prox)

# Apply normalization to the specified columns
temp_modeldata$norm_fvc <- normalize(temp_modeldata$fvc)
temp_modeldata$norm_imper_surf <- normalize(temp_modeldata$imper_surf)
temp_modeldata$norm_area <- normalize(temp_modeldata$area)
temp_modeldata$norm_peri_area <- normalize(temp_modeldata$peri_area)
temp_modeldata$norm_prox <- normalize(temp_modeldata$prox)

## Factorization
temp_modeldata$veg_type <- factor(temp_modeldata$veg_type)
temp_modeldata$transect_id <- factor(temp_modeldata$transect_id)

# Fit the GAM model with normalized data
cooling.gam2020 <- gam(cooling_inten_oC ~ s(distance_m, k = 6, bs = 'cr') +
                         s(norm_imper_surf, k = 5, bs = 'cr') +
                         s(norm_fvc, k = 5, bs = 'cr') +
                         s(norm_area, k = 5, bs = 'cr') +
                         s(norm_peri_area, k = 5, bs = 'cr') +
                         s(norm_prox, k = 5, bs = 'cr'),
                       data = temp_modeldata, random=list(transect_id = ~ 1),
                       family = gaussian(link = "log"), method = "REML", 
                       correlation = corGaus())

# Summary of GAMs result
summary(cooling.gam2020)

# Check the variance of the random effects
gam.vcomp(cooling.gam2020)


#***********************************************************************
#*
# GAM model
cooling.gam2020 <- gam(cooling_inten_oC ~ s(distance_m, k = 6, bs = 'cr') +
                         s(norm_imper_surf, k = 5, bs = 'cr') +
                         s(norm_fvc, k = 5, bs = 'cr') +
                         s(norm_area, k = 5, bs = 'cr') +
                         s(norm_peri_area, k = 5, bs = 'cr') +
                         s(norm_prox, k = 5, bs = 'cr'),
                       data = temp_modeldata, random=list(transect_id = ~ 1),
                       family = gaussian(link = "log"), method = "REML", 
                       correlation = corGaus())

########################################################################
##**              Figure 1. area                                      *

##*********************************************************************
#** Cooling intensity decay due to green area variation along distance*
cooling_data1 <- data.frame(
      distance_m = mean(temp_modeldata$distance_m),
      norm_fvc = mean(temp_modeldata$norm_fvc),
      norm_imper_surf = mean(temp_modeldata$norm_imper_surf),
      norm_area = seq(from = 0, to = 1, by = 0.1),
      norm_peri_area = mean(temp_modeldata$norm_peri_area),
      norm_prox = mean(temp_modeldata$norm_prox),
      transect_id = "1")

## Predict dependent variable results based on GAMs
G1 <- predict(cooling.gam2020, newdata=cooling_data1, type="link", se=TRUE)
F1 <- exp(G1$fit)
#calculate confidence intervals around mean model results
FSEUP1 <- exp(G1$fit + 1.96 * G1$se.fit)
FSELOW1 <- exp(G1$fit - 1.96 * G1$se.fit)
#plot model results

png(file="Output/Cooling_intensity_area.png",
    width=10, height=8, units="in", res=300)
plot(cooling_data1$norm_area, F1, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,7), xlim = c(min(temp_modeldata$norm_area), max(temp_modeldata$norm_area)), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#order values for confidence intervals
i.for <- order(cooling_data1$norm_area)
i.back <- order(cooling_data1$norm_area, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(cooling_data1$norm_area[i.for], cooling_data1$norm_area[i.back])
y.polygon <- c(FSEUP1[i.for], FSELOW1[i.back])
polygon(x.polygon, y.polygon, col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(temp_modeldata$norm_area), temp_modeldata$cooling_inten_oC, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(cooling_data1$norm_area, F1, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Green patch area (norm)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("Cooling intensity")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, at=c(0,1,2,3,4,5,6,7), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0.0,0.2,0.4,0.6,0.8,1.0), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

########################################################################
##**              Figure 2. perimeter-area ratio                       *

##*********************************************************************
#** Cooling intensity decay due to perimeter-area ratio variation along distance*
cooling_data2 <- data.frame(
  distance_m = mean(temp_modeldata$distance_m),
  norm_fvc = mean(temp_modeldata$norm_fvc),
  norm_imper_surf = mean(temp_modeldata$norm_imper_surf),
  norm_area = mean(temp_modeldata$norm_area),
  norm_peri_area = seq(from = 0, to = 1, by = 0.01),
  norm_prox = mean(temp_modeldata$norm_prox),
  transect_id = "1")

## Predict dependent variable results based on GAMs
G2 <- predict(cooling.gam2020, newdata=cooling_data2, type="link", se=TRUE)
F2 <- exp(G2$fit)
#calculate confidence intervals around mean model results
FSEUP2 <- exp(G2$fit + 1.96 * G2$se.fit)
FSELOW2 <- exp(G2$fit - 1.96 * G2$se.fit)
#plot model results
png(file="Output/Cooling_intensity_peri_area.png",
    width=10, height=8, units="in", res=300)
plot(cooling_data2$norm_peri_area, F2, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,7), xlim = c(min(temp_modeldata$norm_peri_area), max(temp_modeldata$norm_peri_area)), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#order values for confidence intervals
i.for <- order(cooling_data2$norm_peri_area)
i.back <- order(cooling_data2$norm_peri_area, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(cooling_data2$norm_peri_area[i.for], cooling_data2$norm_peri_area[i.back])
y.polygon <- c(FSEUP2[i.for], FSELOW2[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(temp_modeldata$norm_peri_area), temp_modeldata$cooling_inten_oC, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(cooling_data2$norm_peri_area, F2, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Green patch perimeter-area ratio (norm)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("Cooling intensity")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, at=c(0,1,2,3,4,5,6,7), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0.0,0.2,0.4,0.6,0.8,1.0), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()
########################################################################
##**              Figure 3. proximity                                 *

##*********************************************************************
#** Cooling intensity decay due to mean proximity to surrounding green patches variation along distance*
cooling_data3 <- data.frame(
  distance_m = mean(temp_modeldata$distance_m),
  norm_fvc = mean(temp_modeldata$norm_fvc),
  norm_imper_surf = mean(temp_modeldata$norm_imper_surf),
  norm_area = mean(temp_modeldata$norm_area),
  norm_peri_area = mean(temp_modeldata$norm_peri_area),
  norm_prox = seq(from = 0, to = 1, by = 0.01),
  transect_id = "1")

## Predict dependent variable results based on GAMs
G3 <- predict(cooling.gam2020, newdata=cooling_data3, type="link", se=TRUE)
F3 <- exp(G3$fit)
#calculate confidence intervals around mean model results
FSEUP3 <- exp(G3$fit + 1.96 * G3$se.fit)
FSELOW3 <- exp(G3$fit - 1.96 * G3$se.fit)
#plot model results
png(file="Output/Cooling_intensity_prox.png",
    width=10, height=8, units="in", res=300)
plot(cooling_data3$norm_prox, F3, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,7), xlim = c(min(temp_modeldata$norm_prox), max(temp_modeldata$norm_prox)), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#order values for confidence intervals
i.for <- order(cooling_data3$norm_prox)
i.back <- order(cooling_data3$norm_prox, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(cooling_data3$norm_prox[i.for], cooling_data3$norm_prox[i.back])
y.polygon <- c(FSEUP3[i.for], FSELOW3[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(temp_modeldata$norm_prox), temp_modeldata$cooling_inten_oC, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(cooling_data3$norm_prox, F3, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Green patch proximity (norm)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("Cooling intensity")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, at=c(0,1,2,3,4,5,6,7), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0.0,0.2,0.4,0.6,0.8,1.0), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

##****************************************************************************
##**                      PART III                                           *
##*
##** GAMs cooling intensity prediction in space and time                     *
##** Predict across the benefiting areas in the study area                                                      * 
##*
##****************************************************************************
# Fit the GAM model with normalized data
cooling.gam2020 <- gam(cooling_inten_oC ~ s(norm_fvc, k = 5, bs = 'cr') +
                         s(norm_imper_surf, k = 5, bs = 'cr') +
                         s(norm_area, k = 5, bs = 'cr') +
                         s(norm_peri_area, k = 5, bs = 'cr') +
                         s(norm_prox, k = 5, bs = 'cr'),
                       data = temp_modeldata, 
                       random = list(transect_id = ~1,  distance_m=~0),
                       family = gaussian(link = "log"), method = "REML", 
                       correlation = corGaus())



#** Load raster predictors 0f 2020*
fvc20 <- rast("fvc_cooling_2020.tif")
imper_surf20 <- rast("Impersurf_cooling_2020.tif")
area20 <- rast("area_cooling_2020.tif")
peri_area20 <- rast("peri_area_cooling_2020.tif")
prox20 <- rast("prox_cooling_2020.tif")

# Combine raster layers
raslist <- list(fvc20, imper_surf20, area20, peri_area20, prox20)

# Normalising the raster datasets to 0 and 1
normalize_raster <- function(ras) {
  ras_min <- minmax(ras)[1]
  ras_max <- minmax(ras)[2]
  (ras - ras_min) / (ras_max - ras_min)
}

# Resample and normalise predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resampled_ras <- resample(crop_ras, fvc20, method = "bilinear")
  normalize_raster(resampled_ras)
})

# Stack raster layers and rename them to make it similar to the GAMs model names
preds_stack20 <- rast(pred_rasters)
names(preds_stack20) <- c("norm_fvc", "norm_imper_surf", "norm_area", "norm_peri_area", "norm_prox")

#** Spatial prediction for 2020*
cooling_pred2020 <- terra::predict(preds_stack20, model = cooling.gam2020, type = "response")
# Plot
plot(cooling_pred2020, col = terrain.colors(100), main = "Cooling Prediction Map 2020")
# Export
cooling20 <- file.path("Output", "Cooling intensity map 2020.tif")
writeRaster(cooling_pred2020, cooling20, overwrite = TRUE)


##************************************************
##** Temporal prediction for 2000               *
#** Load raster predictors of 2020*
fvc00 <- rast("fvc_cooling_2000.tif")
imper_surf00 <- rast("Impersurf_cooling_2000.tif")
area20 <- rast("area_cooling_2000.tif")
peri_area00 <- rast("peri_area_cooling_2000.tif")
prox00 <- rast("prox_cooling_2000.tif")

# Combine raster layers
raslist <- list(fvc00, imper_surf00, area00, peri_area00, prox00)

# Normalising the raster datasets to 0 and 1
normalize_raster <- function(ras) {
  ras_min <- minmax(ras)[1]
  ras_max <- minmax(ras)[2]
  (ras - ras_min) / (ras_max - ras_min)
}

# Resample and normalise predictor rasters
pred_rasters00 <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resampled_ras <- resample(crop_ras, fvc00, method = "bilinear")
  normalize_raster(resampled_ras)
})

# Stack raster layers and rename them to make it similar to the GAMs model names
preds_stack00 <- rast(pred_rasters00)
names(preds_stack00) <- c("norm_fvc", "norm_imper_surf", "norm_area", "norm_peri_area", "norm_prox")


#** Temporal prediction for 2000*
cooling_pred2000 <- terra::predict(preds_stack00, model = cooling.gam2020, type = "response")
# Plot
plot(cooling_pred2000, col = terrain.colors(100), main = "Cooling Prediction Map 2000")
# Export
cooling20 <- file.path("Output", "Cooling intensity map 2020.tif")
writeRaster(cooling_pred2000, cooling00, overwrite = TRUE)


##****************************************************************************
##**                      PART IV   vis.gam plot                             *
##**  GAMs plotting of grass biomass using vis.gam function                  *
##*
##***************************************************************************

## Custom color for plotting
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

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


##**********************************************************


##** Model result visualization using vis.gam              *
##*
##**********************************************************
temp_modeldata$veg_type <- factor(temp_modeldata$veg_type)

# Fit the GAM model with normalized data
cooling.gam2020 <- gam(cooling_inten_oC ~ s(distance_m, k = 6, bs = 'cr') +
                         s(norm_imper_surf, k = 5, bs = 'cr') +
                         s(fvc, k = 5, bs = 'cr') +
                         s(norm_area, k = 5, bs = 'cr') +
                         s(norm_peri_area, k = 5, bs = 'cr') +
                         s(norm_prox, k = 5, bs = 'cr')+
                         ti(X, Y, k = 20),
                       data = temp_modeldata, random=list(transect_id = ~1),
                       family = quasipoisson(link = "identity"), method = "REML", 
                       correlation = corGaus())

# Summary of GAMs result
summary(cooling.gam2020)

## colors palette definition
jet.colors <- colorRampPalette(c("#fe9929","#ffffb2", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"))

## Saving the plot output - 
png(file="Output/Cooling_intensity_plot_2020.png",
    width=10, height=8, units="in", res=300)
myvis.gam(cooling.gam2020, view=c("X", "Y"), color = "jet",
        type = "response", xlab="", ylab="",
        main = "", cex.main=1.5, font.main=1,
        plot.type="contour", n.grid=600, too.far=0.5,
        ylim = c(min(temp_modeldata$Y),max(temp_modeldata$Y)), xlim = c(min(temp_modeldata$X),max(temp_modeldata$X)),
        cex.lab=1.5, cex.axis=1.5, las=1, yaxt="n", xaxt="n",
        contour.col="black", labcex=1.3, lwd=1.5, nlevels = 30)
## Create box bound
box(which="plot", lwd = 2)
## Plot y-axis for latitude
axis(side=2, las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)
## Axis and figure labels
mtext("UTM Easting", side=1, line=3, cex=1.5)
mtext("UTM Northing", side=2, line=2.5, cex=1.5)

## Finally, Run dev.off() to create the file!
dev.off()

##** OR *
##*******

## Saving the plot output
#png(file="Output/Cooling_intensity.png",
#    width=10, height=8, units="in", res=300)
myvis.gam(cooling.gam2020, view=c("X", "Y"), color = "heat",
          type = "response", xlab="", ylab="",
          main = "", cex.main=1.5, font.main=1,
          plot.type="contour", n.grid=600, too.far=0.0,
          ylim = c(1278000,1305000), xlim = c(300000,340500),
          cex.lab=1.5, cex.axis=1.5, las=1, yaxt="n", xaxt="n",
          contour.col="black", labcex=1.3, lwd=1.5, nlevels = 40)
## Create box bound
box(which="plot", lwd = 2)
## Plot y-axis for latitude
axis(side=2, las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
## Plot x-axis for longitude
axis(1, cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)
## Axis and figure labels
mtext("Longitude", side=1, line=3, cex=1.5)
mtext("Latitude", side=2, line=2.5, cex=1.5)
## Overlay polygon data
plot.sf(bound, col = NA, bg =  NULL, lwd = 2, add = TRUE)

## Finally, Run dev.off() to create the file!
#dev.off()


