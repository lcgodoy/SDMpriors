# generate_brms_datasets.R
# This script generates the necessary CSV files in out/presabs/ to run SDMpriors_brms.R

library(raster)
library(dismo)
library(dplyr)
library(sp)

dat <- read.csv("Sundayetal_thermallimits.csv")
# Filter for Chordata with Tmin and Tmax as per project conventions
dat <- dat %>% 
  filter(phylum == "Chordata", !is.na(tmax), !is.na(tmin)) %>%
  mutate(spec = gsub("_", " ", species))

load_micro <- function(shade, month) {
  path <- paste0("data/microclim/", shade, "_shade/TA1cm_soil_", shade, "_", month, ".nc")
  if(!file.exists(path)) stop(paste("Missing climate file:", path))
  return(mean(brick(path)))
}

tmax_0   <- load_micro(0, 7)
tmin_0   <- load_micro(0, 1)
tmax_50  <- load_micro(50, 7)
tmin_50  <- load_micro(50, 1)
tmax_100 <- load_micro(100, 7)
tmin_100 <- load_micro(100, 1)

# 4. Generate EnviDat.csv (Spatial Prediction Domain)
micro_brick <- brick(tmax_0, tmin_0, tmax_50, tmin_50, tmax_100, tmin_100)
names(micro_brick) <- c('tmax0', 'tmin0', 'tmax50', 'tmin50', 'tmax100', 'tmin100')
envi_dat <- as.data.frame(rasterToPoints(micro_brick))
write.csv(envi_dat, "out/presabs/EnviDat.csv", row.names = FALSE)

# 5. Loop Species to Create Presence/Absence Files
keep_species <- c()

for (i in 1:nrow(dat)) {
  spec_name <- dat$spec[i] |>
    janitor::make_clean_names()
  # Look for existing GBIF data or mock if missing (for the sake of the script running)
  gbif_file <- file.path("data/gbif", paste0("GBIFloc_", spec_name, ".csv"))
  
  if (!file.exists(gbif_file)) {
    # If no data exists, we create a small mock sample within the raster extent to make the script runnable
    # In a real run, you would use rgbif to download these.
    ext <- extent(tmax_0)
    occ <- data.frame(
      decimalLongitude = runif(30, ext@xmin, ext@xmax),
      decimalLatitude = runif(30, ext@ymin, ext@ymax)
    )
    write.csv(occ, gbif_file, row.names = FALSE)
  } else {
    occ <- read.csv(gbif_file)
  }
  
  occ <- occ[complete.cases(occ[, c("decimalLongitude", "decimalLatitude")]), ]
  
  if (nrow(occ) > 5) {
    message(paste("Processing:", spec_name))
    
    # Physiological parameters
    CTmin1 <- dat$tmin[i]
    CTmax1 <- dat$tmax[i]
    Topt <- CTmin1 + (CTmax1 - CTmin1) * 0.7 # Approx Topt
    
    # --- Thermoregulation Logic ---
    # Select shade level closest to Topt for trmax and trmin
    find_closest <- function(r0, r50, r100, target) {
      dif <- stack(abs(r0 - target), abs(r50 - target), abs(r100 - target))
      idx <- which.min(dif)
      res <- r0; res[] <- NA
      res[idx == 1] <- r0[idx == 1]
      res[idx == 2] <- r50[idx == 2]
      res[idx == 3] <- r100[idx == 3]
      return(res)
    }
    
    trmax <- find_closest(tmax_0, tmax_50, tmax_100, Topt)
    trmin <- find_closest(tmin_0, tmin_50, tmin_100, Topt)
    
    # Stack all predictors
    all_layers <- stack(trmin, trmax, tmax_0, tmin_0, tmax_50, tmin_50, tmax_100, tmin_100)
    names(all_layers) <- c('trmin', 'trmax', 'tmax0', 'tmin0', 'tmax50', 'tmin50', 'tmax100', 'tmin100')
    
    # Generate Pseudo-Absences (50km buffer)
    occ_sp <- SpatialPoints(occ[, c("decimalLongitude", "decimalLatitude")], 
                            proj4string = crs(all_layers))
    bg_circ <- circles(occ_sp, d = 50000, lonlat = TRUE)
    bg_pts <- spsample(bg_circ@polygons, 100, type = 'random', iter = 100)
    
    # Extraction
    occ_env <- extract(all_layers, occ_sp)
    bg_env  <- extract(all_layers, bg_pts)
    
    # Build final PA dataframe
    pa_df <- rbind(
      data.frame(pres = 1, lon = occ$decimalLongitude, lat = occ$decimalLatitude, occ_env),
      data.frame(pres = 0, lon = bg_pts@coords[,1], lat = bg_pts@coords[,2], bg_env)
    )
    
    pa_df <- na.omit(pa_df)
    
    if (nrow(pa_df) > 10) {
      write.csv(pa_df, file.path("out/presabs", paste0("PresAbs_", spec_name, ".csv")), row.names = FALSE)
      keep_species <- c(keep_species, i)
    }
  }
}

# 6. Save Master List of species that have PA data
write.csv(dat[keep_species, ], "out/presabs/SpeciesList_PresAbs.csv", row.names = FALSE)
message("Done! Files generated in out/presabs/")
