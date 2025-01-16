
#** Title: Modelling the Impact of Ecosystem Fragmentation on Ecosystem Services in the Degraded Ethiopian Highlands *#
#** Sitotaw, T. M., Willemen, L., Meshesha, D. T., & Nelson, A., 2024.            *
#** GAMs Analysis for 2020 for grass biomass model based on field biomass records *

#**************************************************************************
# clear workspace
rm(list=ls())

# Libraries
library(terra)
library(sf)
library(ggplot2)
library(dplyr)

##****************************************************************************
##**                      PART I                                             *
##**  Predictors extraction using the plot location from rasters             *
##**  Predictor variables: ndvi, fvc, area, perimeter-area ratio, proximity  *
##****************************************************************************
# Set working directory
setwd("D:/ES_ECOINF/grass_biomass_data")

# Load grass biomass data and geographic locations
grassBiomass <- read.csv("Plot_GrassAGB_Raw.csv", stringsAsFactors=FALSE)

# Cleaning No Data values and change to 0 values
grassBiomass[is.na(grassBiomass[,"freshAGB_gm_0.25m2"]), "freshAGB_gm_0.25m2"] <- 0
grassBiomass[is.na(grassBiomass[,"dryAGB_gm_0.25m2"]), "dryAGB_gm_0.25m2"] <- 0

## Convert the grass biomass (gm/0.25m^2) to biomass (ton per 15m radius = 705.5 m^2)
grassBiomass$freshAGB_ton_15mRad <- ((grassBiomass$freshAGB_gm_0.25m2*705.5)/0.25)*0.000001
grassBiomass$dryAGB_ton_15mRad <- ((grassBiomass$dryAGB_gm_0.25m2*705.5)/0.25)*0.000001
head(grassBiomass)

## convert the plot table to points to extract data from raster layers
prj4str <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
plot_data <- st_as_sf(grassBiomass, coords = c("X", "Y"), crs = prj4str, agr = "constant")

# Export Spatial Points
output_file <- "biomass_plot.shp"
if (file.exists(output_file)) {file.remove(output_file)}
st_write(plot_data, output_file)


#**  Predictors extraction using the plot location from rasters*
##***************************************************************
# Read grass biomass plot data
plot_biomass <- st_read("biomass_plot.shp")
plot_biomass <- plot_data
# Read styudy area boundary as spatial data
bound <- st_read("D:/ES_ECOINF/crop_pollination_data/Basin_Boundary.json")

#** Load predictors - 2020*
ndvi20 <- rast("ndvi_biomass_2020.tif")
area20 <- rast("area_biomass_2020.tif")
peri_area20 <- rast("peri-area_biomass_2020.tif")
prox20 <- rast("prox_biomass_2020.tif")

# Combine raster layers
raslist <- list(ndvi20, area20, peri_area20, prox20)

# Resample predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resample(crop_ras, ndvi20, method = "bilinear")
})

# Stack raster layers
raster_stack <- rast(pred_rasters)

# Define the extraction function within a 15 me buffer around each point
extract.fun <- function(rast, plot_biomass) {
  result <- extract(rast, plot_biomass, buffer = 15, fun = mean, na.rm = TRUE)
  return(result[, -1])  # Extracted values
}

# Apply the extraction function to each raster in the stack
rasextract <- extract.fun(raster_stack, plot_biomass)

# Assign column names based on raster layer names
layer_names <- names(raster_stack)
colnames(rasextract) <- layer_names

head(plot_biomass)
## Combine extracted predictor values with plot data
modeldata <- cbind(plot_biomass[, c("site", "plot_id", "freshAGB_ton_15mRad", "dryAGB_ton_15mRad")], rasextract)
head(modeldata)

## Append X, Y coordinates
biomass_pred_data <- cbind(grassBiomass[, c("X", "Y")], modeldata)
biomass_pred_data <- as.data.frame(biomass_pred_data)
## Re-order column data
biomass_pred_data <- biomass_pred_data[, c("site", "plot_id", "X", "Y", "freshAGB_ton_15mRad", "dryAGB_ton_15mRad",
                      "ndvi_biomass_2000", "area_ha", "perimeter_area_ratio", "prox_biomass_2020")]
## Rename
names(biomass_pred_data) <- c("site", "plot_id", "X", "Y", "fAGB_ton", "dAGB_ton", "ndvi", 
                              "area", "Peri_area", "prox")
head(biomass_pred_data)
# Export combined result for further GAMs analysis
write.table(biomass_pred_data, file="Output/biomass_predictors_data_2020.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


#****************************************************************************

#**                      PART II                                            *
#**  GAMs Analysis of grass biomass supply using mgcv package               *
#***************************************************************************
biomass_frag_2020 <- read.csv("Output/biomass_predictors_data_2020.csv", stringsAsFactors=FALSE)
biomass_frag_2020 <- biomass_pred_data
## Replace NA using replace() & is.na()
biomass_frag_2020 <- replace(biomass_frag_2020, is.na(biomass_frag_2020), 0)

#** Normalise Data columns with Min-Max Scaling*
cols_to_normalize <- (ncol(biomass_frag_2020) - 3):ncol(biomass_frag_2020)
biomass_frag_2020[cols_to_normalize] <- lapply(biomass_frag_2020[cols_to_normalize], function(x) {
  (x - min(x)) / (max(x) - min(x))
})
names(biomass_frag_2020) <- c("site", "plot_id", "X", "Y", "fAGB_ton", "dAGB_ton", "ndvi", "norm_area", "norm_peri_area", "norm_prox")
head(biomass_frag_2020)


##** Fit GAMs. Both fixed effects and random effects.*
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(mgcv)

biomass_gam2020 <- gam(dAGB_ton ~ s(ndvi, k = 6, bs = 'cr') +
                         s(norm_area, k = 6, bs = 'cr') + 
                         s(norm_peri_area, k = 6, bs = 'cr') +
                         s(norm_prox, k = 6, bs = 'cr')+
                         ti(X, Y, k = 10),
                       data = biomass_frag_2020,
                       random=list(site = ~1),
                       family = quasipoisson(link = "log"),
                       method = "REML", correlation = corAR1())

summary(biomass_gam2020)

#**********************************************************************
biomass_gam2020 <- gam(dAGB_ton ~ s(ndvi, k = 6, bs = 'cr') +
                         s(norm_area, k = 6, bs = 'cr') + 
                         s(norm_peri_area, k = 6, bs = 'cr') +
                         s(norm_prox, k = 6, bs = 'cr'),
                       data = biomass_frag_2020,
                       random=list(site = ~1),
                       family = quasipoisson(link = "log"),
                       method = "REML", correlation = corAR1())


########################################################################
##**              Figure 1. Area                                        *
##*
########################################################################
# Variation in grass biomass
predict_data1 <- data.frame(ndvi = mean(biomass_frag_2020$ndvi),
                            norm_area = seq(from = 0, to = 1, by = 0.01),
                            norm_peri_area = mean(biomass_frag_2020$norm_peri_area),
                            norm_prox = mean(biomass_frag_2020$norm_prox),
                            Site = "1")

G2 <- predict(biomass_gam2020, newdata=predict_data1, type="link", se=TRUE)
F2 <- 1/(1+exp(-G2$fit))
FSEUP2 <- 1/(1+exp(-G2$fit - 1.96 * G2$se.fit))
FSELOW2 <- 1/(1+exp(-G2$fit + 1.96 * G2$se.fit))

## Saving the plot output
#png(file="Output/Biomass_area_plot.png",
#    width=8, height=6, units="in", res=300)
plot(predict_data1$norm_area, F2, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,1), cex.lab=2.5, cex.axis=1.7, las=1, tck=0.02)
i.for <- order(predict_data1$norm_area)
i.back <- order(predict_data1$norm_area, decreasing=TRUE)
x.polygon <- c(predict_data1$norm_area[i.for], predict_data1$norm_area[i.back])
y.polygon <- c(FSEUP2[i.for], FSELOW2[i.back]) 

polygon(x.polygon,y.polygon,col="grey", border=NA)
points(jitter(biomass_frag_2020$norm_area), biomass_frag_2020$dAGB_ton, pch=19, col="black", cex=1.25)
lines(predict_data1$norm_area, F2, lty=1, lwd=3, col="black")

mtext(text="Grass biomass", side=2, line=2.5, cex=1.5)
mtext(text=expression(paste("Wetland fragment area")), side=1, line=3, cex=1.5)
box(which="plot", lwd=4)
axis(side=2, at=c(0,0.4,0.6,0.8,1.2), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()
##************

#######################################################################
#**              Figure 2. Perimeter:area ratio                      *
########################################################################
## Create new dataset
predict_data2 <- data.frame(ndvi = mean(biomass_frag_2020$ndvi),
                            norm_area = mean(biomass_frag_2020$norm_peri_area),
                            norm_peri_area = seq(from=0, to=1, by=0.01),
                            norm_prox = mean(biomass_frag_2020$norm_prox),
                            Site="1")

## predict predictor variable results
pred.1 <- predict(biomass_gam2020, newdata=predict_data2, type="link", se=TRUE)
fit.1 <- 1/(1+exp(-pred.1$fit))
conUpper.1 <- 1/(1+exp(-pred.1$fit - 1.96 * pred.1$se.fit))
conLower.1 <- 1/(1+exp(-pred.1$fit + 1.96 * pred.1$se.fit))

## Saving the plot output
png(file="Output/Grass_biomass_perimeter-area_2020.png",
   width=8, height=6, units="in", res=300)
## plot predicted model results
plot(predict_data2$norm_peri_area, fit.1, ylab="", xlab="", yaxt="n", xaxt="n",
     cex.main=1.5, font.main=1,
     type="n", ylim=c(0,1), cex.lab=1.5, cex.axis=1.5, las=1, tck=0.02)
# Confidence intervals
order.for <- order(predict_data2$norm_peri_area)
order.back <- order(predict_data2$norm_peri_area, decreasing=TRUE)
## plot confidence intervals as shaded polygons
x.polygon <- c(predict_data2$norm_peri_area[order.for], predict_data2$norm_peri_area[order.back])
y.polygon <- c(conUpper.1[order.for], conLower.1[order.back])
polygon(x.polygon, y.polygon,col="grey", border=NA)
# plot data points and jitter them to reduce overlap
points(jitter(biomass_frag_2020$norm_peri_area), biomass_frag_2020$dAGB_ton, pch=19, col="grey20", cex=1.25)
## plot mean model results as line
lines(predict_data2$norm_peri_area, fit.1, lty=1, lwd=3, col="black")
# add thicker line around plot
box(which="plot", lwd=2)
# add axis labels
axis(side=2, at=c(0,0.4,0.8,1.2), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)
# add axis labels
mtext(text=expression(paste("Grass biomass")), side=2, line=2.5, cex=1.5)
mtext(text=expression(paste("Wetland perimeter-area ratio")), side=1, line=3, cex=1.5)

## Finally, Run dev.off() to create the file!
dev.off()

########################################################################
#**              Figure 3. Proximity                                  *
########################################################################
## Change in grass biomass
predict_data3 <- data.frame(ndvi = mean(biomass_frag_2020$ndvi),
                            norm_area = mean(biomass_frag_2020$norm_area),
                            norm_peri_area = mean(biomass_frag_2020$norm_peri_area),
                            norm_prox = seq(from = 0, to = 1, by = 0.01),
                            Site = "1")

G2 <- predict(biomass_gam2020, newdata=predict_data3, type="link", se=TRUE)
F2 <- 1/(1+exp(-G2$fit))
FSEUP2 <- 1/(1+exp(-G2$fit - 1.96 * G2$se.fit))
FSELOW2 <- 1/(1+exp(-G2$fit + 1.96 * G2$se.fit))

png(file="Output/Grass_biomass_prox_2020.png",
    width=8, height=6, units="in", res=300)
plot(predict_data3$norm_prox, F2, ylab="", xlab="", yaxt="n", xaxt="n",
     main="", cex.main=2.5, font.main=1,
     type="n", ylim=c(0,1), cex.lab=2.5, cex.axis=1.7, las=1, tck=0.02)
i.for <- order(predict_data3$norm_prox)
i.back <- order(predict_data3$norm_prox, decreasing=TRUE)
x.polygon <- c(predict_data3$norm_prox[i.for], predict_data3$norm_prox[i.back])
y.polygon <- c(FSEUP2[i.for], FSELOW2[i.back])
polygon(x.polygon,y.polygon,col="grey", border=NA)
points(jitter(biomass_frag_2020$norm_prox), biomass_frag_2020$dAGB_ton, pch=19, col="black", cex=1.25)
lines(predict_data3$norm_prox, F2, lty=1, lwd=3, col="black")

mtext(text="Grass biomass", side=2, line=2.5, cex=1.5)
mtext(text=expression(paste("Wetland fragment proximity")), side=1, line=3, cex=1.5)
box(which="plot", lwd=4)
axis(side=2, at=c(0,0.4,0.6,0.8,1.2), las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)

dev.off()

##************

#****************************************************************************
#**                      PART III                                           *
#** GAMs grass biomass prediction in space and time                         *                                                            * 
#*
#****************************************************************************
biomass_gam2020 <- gam(dAGB_ton ~ s(ndvi, k = 6, bs = 'cr') +
                         s(norm_area, k = 6, bs = 'cr') + 
                         s(norm_peri_area, k = 5, bs = 'cr') +
                         s(norm_prox, k = 5, bs = 'cr'),
                       data = biomass_frag_2020,
                       random=list(site = ~1),
                       family = quasipoisson(link = "log"),
                       method = "REML", correlation = corAR1())

# Read boundary as spatial data
bound <- st_read("D:/ES_ECOINF/crop_pollination_data/Basin_Boundary.json")

#** Load predictors - 2020*
ndvi20 <- rast("ndvi_biomass_2020.tif")
area20 <- rast("area_biomass_2020.tif")
peri_area20 <- rast("peri-area_biomass_2020.tif")
prox20 <- rast("prox_biomass_2020.tif")

# Combine raster layers
raslist <- list(ndvi20, area20, peri_area20, prox20)

# Normalising the raster datasets to 0 and 1
normalize_raster <- function(ras) {
  ras_min <- minmax(ras)[1]
  ras_max <- minmax(ras)[2]
  (ras - ras_min) / (ras_max - ras_min)
}

# Resample and normalise predictor rasters
pred_rasters20 <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resampled_ras <- resample(crop_ras, ndvi20, method = "bilinear")
  normalize_raster(resampled_ras)
})

# Stack raster layers and rename them to make it similar to the GAMs model names
preds_stack20 <- rast(pred_rasters20)
names(preds_stack20) <- c("ndvi", "norm_area", "norm_peri_area", "norm_prox")

#** Spatial prediction for 2020*
biomass_pred20 <- terra::predict(preds_stack20, model = biomass_gam2020, type = "response")
# Plot
plot(biomass_pred2020, col = terrain.colors(100), main = "Grass biomass map 2020")
# Export
biomass20 <- file.path("Output", "Grass biomass map 2020.tif")
writeRaster(biomass_pred20, biomass00, overwrite = TRUE)


#***********************************************************
#** Grass biomass prediction for 2000                      *
#** Load predictor rasters of 2000                         *
#****************************************************************************

##** Load predictors - 2000*
ndvi00 <- rast("ndvi_biomass_2000.tif")
area00 <- rast("area_biomass_2000.tif")
peri_area00 <- rast("peri-area_biomass_2000.tif")
prox00 <- rast("prox_biomass_2000.tif")

# Combine raster layers
raslist <- list(ndvi00, area00, peri_area00, prox00)

# Normalising the raster datasets to 0 and 1
normalize_raster <- function(ras) {
  ras_min <- minmax(ras)[1]
  ras_max <- minmax(ras)[2]
  (ras - ras_min) / (ras_max - ras_min)
}

# Resample and normalise predictor rasters
pred_rasters <- lapply(raslist, function(ras) {
  crop_ras <- crop(ras, ext(bound))
  resampled_ras <- resample(crop_ras, ndvi00, method = "bilinear")
  normalize_raster(resampled_ras)
})

# Stack raster layers and rename them to make it similar to the GAMs model names
preds_stack00 <- rast(pred_rasters)
names(preds_stack00) <- c("ndvi", "norm_area", "norm_peri_area", "norm_prox")
## GAMs model

#** Temporal prediction for 2000*
biomass_pred00 <- terra::predict(predstack00, model = biomass_gam2020, type = "response")
# Plot
plot(biomass_pred2000, col = terrain.colors(100), main = "Grass biomass map 2000")
# Export predicted map
biomass00 <- file.path("Output", "Grass biomass map 2000.tif")
writeRaster(biomass_pred00, biomass00, overwrite = TRUE)


#**************************************************************************
#**                      PART IV                                          *
#**  Grass biomass contour maps using vis.gam function                    *
#*
#**************************************************************************

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

##*******************************************************************

library(mgcv)

biomass_gam2020 <- gam(dAGB_ton ~ s(ndvi, k = 6, bs = 'cr') +
                         s(norm_area, k = 6, bs = 'cr') + 
                         s(norm_peri_area, k = 6, bs = 'cr') +
                         s(norm_prox, k = 6, bs = 'cr')+
                         ti(X, Y, k = 10),
                       data = biomass_frag_2020,
                       random=list(site = ~1),
                       family = quasipoisson(link = "log"),
                       method = "REML", correlation = corAR1())

summary(biomass_gam2020)


summary(biomass_gam2020)

## colors palette definition
#jet.colors <- colorRampPalette(c("#fe9929","#ffffb2", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"))

## Saving the plot output
png(file="Output/Grass_bomass.png",
    width=8, height=6, units="in", res=300)
myvis.gam(biomass_gam2020, view=c("X", "Y"), color = "jet",
          type = "response", xlab="", ylab="",
          main = "", cex.main=1.5, font.main=1,
          plot.type="contour", n.grid=600, too.far=0,
          cex.lab=1.5, cex.axis=1.5, las=1, yaxt="n", xaxt="n",
          contour.col="black", labcex=1.3, lwd=1.5, nlevels = 60)
## Create box bound
box(which="plot", lwd = 2)
## Plot y-axis for latitude
axis(side=2, las=1, cex.axis=1.5, tck=0.02, hadj=0.8, lwd=2.5, las=3)
## Plot x-axis for longitude
axis(1, cex.axis=1.5, tck=0.02, padj=0.2, lwd=2.5)
## Axis and figure labels
mtext("Easting", side=1, line=3, cex=1.5)
mtext("Northing", side=2, line=2.5, cex=1.5)

## Finally, Run dev.off() to create the file!
dev.off()

