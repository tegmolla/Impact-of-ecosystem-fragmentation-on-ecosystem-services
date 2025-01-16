#####################################################################################
##** Modelling the impact of ecosystem fragmentation on ecosystem services in the degraded Ethiopian Highlands    *
##** Tegegne Molla, Louise Willemen, Derege Tsegaye Meshesha, Martha Weldemichael, and Andrew Nelson  * 
##** ES Type: Crop pollination services based on flower visitation rate records  *
##** Journal of Ecological Informatics                                           *
##** GAMMs Analysis                                            *

## clear workspace
rm(list=ls())

## Libraries
library(terra)
library(sf)
library(ggplot2)
library(dplyr)


#****************************************************************************
#**                      PART I                                             *
#**  Raster predictors extraction using the locations of pollinators visitation rate data*
#**  Predictor variables: fvc, area, perimeter-area ratio, proximity  *
#***************************************************************************

# Set working directory
setwd("D:/ES_ECOINF/crop_pollination_data")

# Crop pollination data
visit_data <- read.csv("pollinator_visitation_data_raw.csv")
head(visit_data)

## convert the plot table to points to extract data from raster layers
prj4str <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
visit_data_shp <- st_as_sf(visit_data, coords = c("X", "Y"), crs = prj4str, agr = "constant")


#** Load predictors *
fvc20 <- rast("fvc_pollin_2020.tif")
area20 <- rast("area_pollin_2020.tif")
peri_area20 <- rast("peri-area_pollin_2020.tif")
prox20 <- rast("prox_pollin_2020.tif")

# Read boundary as spatial data
bound <- st_read("Basin_Boundary.json")

# Combine raster layers
raslist <- list(fvc20, area20, peri_area20, prox20)

## Resample predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resample(crop_ras, fvc20, method = "bilinear")
})

## Stack raster layers
raster_stack_poll <- rast(pred_rasters)

# Define the extraction function
extract.fun <- function(rast, visit_data_shp) {
  result <- extract(rast, visit_data_shp, buffer = 1500, method = "max", na.rm = TRUE)
  return(result[, -1])
}

# Apply the extraction function to each raster in the stack
rasextract <- extract.fun(raster_stack_poll, visit_data_shp)

# Assign column names based on raster layer names
layer_names <- names(raster_stack_poll)
colnames(rasextract) <- layer_names

head(visit_data)
## Combine extracted predictor values with plot data
pollin_preds_data <- cbind(visit_data[, c("transect", "plot_name", "X", "Y", "crop_type", "distance_m", "visitation_rate")], rasextract)
head(pollin_preds_data)

## Rename
names(pollin_preds_data) <- c("transect", "plot_name", "X", "Y", "crop_type", "distance_m", "visitation_rate", "fvc",
                            "area", "peri_area", "prox")
head(pollin_preds_data)
pollin_preds_data20 <- pollin_preds_data
#write.table(pollin_preds_data, file="pollin_predictors_data_2020.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


#****************************************************************************
#**                      PART II                                            *
#
#**  GAMs Analysis of crop pollination services using mgcv package          *
#
#***************************************************************************
library(mgcv)
pollin_preds_data20 <- read.csv("pollin_predictors_data_2020.csv", stringsAsFactors = FALSE)

## Replace NA using replace() & is.na()
pollin_preds_data20 <- replace(pollin_preds_data20, is.na(pollin_preds_data20), 0)

# Normalise predictor variables between 0 and 1
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# Apply normalization to the selected columns
pollin_preds_data20$norm_fvc <- normalize(pollin_preds_data20$fvc)
pollin_preds_data20$norm_area <- normalize(pollin_preds_data20$area)
pollin_preds_data20$norm_prox <- normalize(pollin_preds_data20$prox)
pollin_preds_data20$norm_peri_area <- normalize(pollin_preds_data20$peri_area)


# Replace NA values with zeros across the entire dataframe
pollin_preds_data20[is.na(pollin_preds_data20)] <- 0

# Convert Crop_type to factor
pollin_preds_data20$crop_type2 <- factor(pollin_preds_data20$crop_type)

#************************************************************************
library(mgcv)
library(car)  # check multicollinearity

# Calculate Variance Inflation Factor (VIF) to assess collinearity (VIF < 10)
vif_model <- lm(visitation_rate ~ distance_m + norm_area + norm_fvc + 
                  norm_peri_area + norm_prox, data = pollin_preds_data20)
vif(vif_model)

# GAM model
pollin_mod20 <- gam(visitation_rate ~ s(distance_m, k = 4, bs = 'cr') + crop_type2 +
                    s(norm_area, k = 6, bs = 'cr') +
                    s(norm_fvc, k = 6, bs = 'cr') +
                    s(norm_peri_area, k = 6, bs = 'cr') +
                    s(norm_prox, k = 6, bs = 'cr')+
                    ti(X, Y, k = 10),
                  data = pollin_preds_data20,
                  random = list(transect = ~ 1),
                  family = quasipoisson(link = "log"),
                  method = "REML", correlation = corGaus(), select = TRUE)

# Model summary for p-values
summary(pollin_mod20)


#*****************************************************************************
## Fit the GAM model
pollin_mod20 <- gam(visitation_rate ~ s(distance_m, k = 4, bs = 'cr') +
                      s(norm_area, k = 6, bs = 'cr') +
                      s(norm_fvc, k = 6, bs = 'cr') +
                      s(norm_peri_area, k = 6, bs = 'cr') +
                      s(norm_prox, k = 6, bs = 'cr'),
                    data = pollin_preds_data20,
                    random = list(transect = ~ 1),
                    family = quasipoisson(link = "log"),
                    method = "REML", correlation = corGaus(), select = TRUE)


#######################################################################

#**              Figure 1. area                                      *
#*********************************************************************
#** Crop flower visitation rate decay due to habitat area variation along distance*

##** GAMs plots using new data generated                            *
## Change in crop flower visitation with habitat area variation
polliData1 <- data.frame(distance_m = mean(pollin_preds_data20$distance_m),
                         norm_fvc = mean(pollin_preds_data20$fvc),
                         norm_area = seq(from = 0, to = 1, by = 0.01),
                         norm_peri_area = mean(pollin_preds_data20$norm_peri_area),
                         norm_prox = mean(pollin_preds_data20$norm_prox),
                         transect = "2020-01")

## Predict dependent variable results based on GAMs
G1 <- predict(pollin_mod20, newdata=polliData1, type="link", se=TRUE)
F1 <- 1/(1+exp(-G1$fit))
FSEUP1 <- 1/(1+exp(-G1$fit + 1.96 * G1$se.fit))
FSELOW1 <- 1/(1+exp(-G1$fit - 1.96 * G1$se.fit))

#plot predicted results
#png(file="Output/Pollination_plot_area.png",
#   width=10, height=8, units="in", res=300)
plot(polliData1$norm_area, F1, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,1), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#Confidence intervals
i.for <- order(polliData1$norm_area)
i.back <- order(polliData1$norm_area, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(polliData1$norm_area[i.for], polliData1$norm_area[i.back])
y.polygon <- c(FSEUP1[i.for], FSELOW1[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(pollin_preds_data20$norm_area), pollin_preds_data20$visitation_rate, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(polliData1$norm_area, F1, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Habitat area (normalised)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("Crop visitation rate")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

#####################################################################

#**              Figure 2. perimeter-area ratio                    *
#********************************************************************
#** Crop flower visitation rate decay due to habitat perimeter-area ratio variation along distance*

#** New data generated                            *
# Change in crop flower visitation along distance gradient from church forest habitats
polliData2 <- data.frame(distance_m = mean(pollin_preds_data20$distance_m),
                         norm_fvc = mean(pollin_preds_data20$fvc),
                         norm_area = mean(pollin_preds_data20$norm_area),
                                          norm_peri_area = seq(from = 0, to = 1, by = 0.01),
                         norm_prox = mean(pollin_preds_data20$norm_prox),
                         transect = "2020-01")

## Predict dependent variable results based on GAMs
G2 <- predict(pollin_mod20, newdata=polliData2, type="link", se=TRUE)
F2 <- 1/(1+exp(-G2$fit))
FSEUP2 <- 1/(1+exp(-G2$fit + 1.96 * G2$se.fit))
FSELOW2 <- 1/(1+exp(-G2$fit - 1.96 * G2$se.fit))
#plot model results
#png(file="Output/Pollination_plot_peri_area.png",
#    width=10, height=8, units="in", res=300)
plot(polliData2$norm_peri_area, F2, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,1), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#order values for confidence intervals
i.for <- order(polliData2$norm_peri_area)
i.back <- order(polliData2$norm_peri_area, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(polliData2$norm_peri_area[i.for], polliData2$norm_peri_area[i.back])
y.polygon <- c(FSEUP2[i.for], FSELOW2[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(pollin_preds_data20$norm_peri_area), pollin_preds_data20$visitation_rate, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(polliData2$norm_peri_area, F2, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Habitat perimeter-area ratio (normalised)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("Crop visitation rate")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

#####################################################################

#**              Figure 3. Pollinator habitat proximity          *
#********************************************************************
#** Crop flower visitation rate decay due to habitat proximity variation along distance*
##** New data generated                            *
polliData3 <- data.frame(distance_m = mean(pollin_preds_data20$distance_m),
                         norm_fvc = mean(pollin_preds_data20$fvc),
                         norm_area = mean(pollin_preds_data20$norm_area),
                         norm_peri_area = mean(pollin_preds_data20$norm_peri_area),
                         norm_prox = seq(from = 0, to = 1, by = 0.01),
                         transect = "2020-01")

## Predict dependent variable results based on GAMs
G3 <- predict(pollin_mod20, newdata=polliData3, type="link", se=TRUE)
F3 <- 1/(1+exp(-G3$fit))
FSEUP3 <- 1/(1+exp(-G3$fit + 1.96 * G3$se.fit))
FSELOW3 <- 1/(1+exp(-G3$fit - 1.96 * G3$se.fit))
#plot model results
#png(file="Pollination_plot_prox.png",
#    width=10, height=8, units="in", res=300)
plot(polliData3$norm_prox, F3, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,1), cex.lab=2.5, cex.axis=2.0, las=1, tck=0.02)
#Confidence intervals
i.for <- order(polliData3$norm_prox)
i.back <- order(polliData3$norm_prox, decreasing=TRUE)
#plot confidence intervals as shaded polygons
x.polygon <- c(polliData3$norm_prox[i.for], polliData3$norm_prox[i.back])
y.polygon <- c(FSEUP3[i.for], FSELOW3[i.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
#plot data points and jitter them to reduce overlap
points(jitter(pollin_preds_data20$norm_prox), pollin_preds_data20$visitation_rate, pch=19, col="black", cex=1.25)
#plot mean model results as line
lines(polliData3$norm_prox, F3, lty=1, lwd=4, col="black")
#axis labels
mtext(text=expression(paste("Habitat proximity (normalised)")), side=1, line=3, cex=1.5)
mtext(text=expression(paste("Crop visitation rate")), side=2, line=2.5, cex=1.5)
#add thicker line around plot
box(which="plot", lwd=4)
#add custom axis labels
axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

#***************************************************************************
#**                      PART III                                          *

#** Crop pollination service prediction in space and time                 * 
#***************************************************************************
# Fitting GAM model
pollin_mod20 <- gam(visitation_rate ~ s(norm_area, k = 4, bs = 'cr') +
                      s(norm_fvc, k = 4, bs = 'cr') +
                      s(norm_peri_area, k = 4, bs = 'cr') +
                      s(norm_prox, k = 4, bs = 'cr') ,
                    #te(X, Y, k = 4, bs = 'cr')
                    data = pollin_preds_data20,
                    random=list(transect = ~1), 
                    family = gaussian(link = "log"), 
                    method = "REML", correlation = corGaus())

#** Load raster predictors for 2020 and stack together*
fvc20 <- rast("fvc_pollin_2020.tif")
area20 <- rast("area_pollin_2020.tif")
peri_area20 <- rast("peri-area_pollin_2020.tif")
prox20 <- rast("prox_pollin_2020.tif")

# Combine raster layers
raslist <- list(fvc20, area20, peri_area20, prox20)

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
pollin_preds_stack20 <- rast(pred_rasters)
names(pollin_preds_stack20) <- c("norm_fvc", "norm_area", "norm_peri_area", "norm_prox")

#** Spatial prediction (scaling-up) for the whole Lake Tana basin for 2020*
pollination_pred2020 <- terra::predict(pollin_preds_stack20, model = pollin_mod20, type = "response")
# Plot
plot(pollination_pred2020, col = terrain.colors(100), main = "Pollination Prediction Map 2020")
# Export predicted map
filename20 <- file.path("Output", "Crop pollination map 2020.tif")
writeRaster(pollination_pred2020, filename, overwrite = TRUE)


##***********************************************************
##** Load spatial predictor rasters of 2000*
fvc00 <- rast("fvc_pollin_2000.tif")
area00 <- rast("area_pollin_2000.tif")
peri_area00 <- rast("peri-area_pollin_2000.tif")
prox00 <- rast("prox_pollin_2000.tif")

# Combine raster layers
raslist <- list(fvc00, area00, peri_area00, prox00)

# Normalising the raster datasets to 0 and 1
normalize_raster <- function(ras) {
  ras_min <- minmax(ras)[1]
  ras_max <- minmax(ras)[2]
  (ras - ras_min) / (ras_max - ras_min)
}

# Resample and normalise predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resampled_ras <- resample(crop_ras, fvc00, method = "bilinear")
  normalize_raster(resampled_ras)
})

# Stack raster layers and rename them to make it similar to the GAMs model names
pollin_preds_stack00 <- rast(pred_rasters)
names(pollin_preds_stack00) <- c("norm_fvc", "norm_area", "norm_peri_area", "norm_prox")

##****************************************************************************
#** Spatial prediction for 2020*
pollination_pred2000 <- terra::predict(pollin_preds_stack00, model = pollin_mod20, type = "response")
# Plot
plot(pollination_pred2000, col = terrain.colors(100), main = "Pollination Prediction Map 2020")

# Export predicted map
filename00 <- file.path("Crop pollination map 2000.tif")
writeRaster(pollination_pred2000, filename00, overwrite = TRUE)


##****************************************************************************
##**                      PART IV                                            *
##*
##**  GAMs plotting of crop pollination using vis.gam function               *
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



##***********************************************************************
library(mgcv)

#** Pollination service contour map Using vis.gam to show interactions*
# Fit the GAM model
poll.mod2020 <- gam(visitation_rate ~ s(distance_m, k = 4, bs = 'cr') + crop_type2 +
                    s(norm_area, k = 6, bs = 'cr') +
                    s(norm_fvc, k = 6, bs = 'cr') +
                    s(norm_peri_area, k = 6, bs = 'cr') +
                    s(norm_prox, k = 6, bs = 'cr')+
                    ti(X, Y, k = 10),
                  data = pollin_preds_data20,
                  random = list(transect = ~ 1),
                  family = quasipoisson(link = "log"), method = "REML",
                  correlation = corGaus(), select = TRUE)

# Model summary for p-values
summary(pollin_mod20)


# colors palette
jet.colors <- colorRampPalette(c("#fe9929","#ffffb2", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"))

## Saving the plot output
png(file="Output/Pollination_plot_2020.png",
    width=10, height=8, units="in", res=300)
myvis.gam(poll.mod2020, view=c("X", "Y"), color = "jet",
        type = "response", xlab="", ylab="",
        main = "", cex.main=1.5, font.main=1,
        plot.type="contour", n.grid=400, too.far=0.5,
        #ylim = c(),
        cex.lab=1.5, cex.axis=1.5, las=1, yaxt="n", xaxt="n",
        contour.col="black", labcex=1.3, lwd=1.5, nlevels = 40)
# Create box bound
box(which="plot", lwd = 2)
# Plot x- and y-axis
axis(side=2, las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

# Axis and figure labels
mtext("UTM Easting", side=1, line=3, cex=1.5)
mtext("UTM Northing", side=2, line=2.5, cex=1.5)

## Finally, Run dev.off() to create the file!
dev.off()
