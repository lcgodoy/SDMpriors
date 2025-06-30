
#Try out lepidoptera models using IBIS iSDM

desktop<- "y"

library(ggplot2)
library(reshape)
library(viridis)
library(patchwork)
library(raster)
library(splitTools) #divide into test and training
library(dismo) #for pseudo absences

#Data access: https://github.com/RebeccaLovell/OrangeTipAnalyses/blob/main/DataPreparation.R

#=================
#Assess data for European species

#toggle between desktop (y) and laptop (n)
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/Proposals/2024_NSF_BoCP_Sep4/data/EuropeSpecies/")
if(desktop=="n") setwd("/Users/lbuckley/My Drive/Buckley/Work/Proposals/2024_NSF_BoCP_Sep4/data/EuropeSpecies/")

#read species list
ukb<- read.csv("UKchecklist.csv")
ukb$GenSpec<- paste(ukb$Genus, ukb$Species, sep=" ")

finb<- read.csv("Hallfors_Data_Shifts_NicheMetrics_Traits.csv")
finb<- finb[which(finb$Taxonomic.group=="Butterfly"),]

#total species count: 140
spc<- unique(c(ukb$GenSpec, finb$Species))

#combine
ukb$location<- "uk"
finb$location<-"finland"
finb$GenSpec<- finb$Species
spc<- rbind(ukb[,c("GenSpec","location")], finb[,c("GenSpec","location")] )

#---------------
#check match with development data
dev.dat<-read.csv("AppendixS3_SeasonalityDatabase.csv")

spc$devmatch<- match(spc$GenSpec, dev.dat$Species)

#lep match
#read trait data
lept<- read.csv("LepTraits/consensus/consensus.csv")

spc$lepmatch <- match(spc$GenSpec, lept$verbatimSpecies)

#thermal database
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/SDMpriors/data/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/My Drive/Buckley/Work/SDMpriors/data/")

ct.dat<-read.csv("PapilioTherm/PapilioThermMultiContinentalCT_2.csv")
spc$ctmatch<- match(spc$GenSpec, ct.dat$Latin_binomial)

kd.dat<-read.csv("PapilioTherm/PapilioThermMultiContinentalHKDT.csv")
spc$kdmatch<- match(spc$GenSpec, kd.dat$Latin_binomial)

#Pieris rapae, Vanessa atalanta 


#Buckley et al LDT: Aglais urticae, Aricia agestis, Inachis io, Pieris brassicae, Pieris rapae, Polygonia c-album

#================
#try out moth data
#https://repository.rothamsted.ac.uk/item/988z5/yearly-occurrence-of-544-species-of-moths-uk-1990-2019-with-trait-values-and-putative-environmental-drivers
#https://gitlab.com/Yo-B/ann_trait_msdm

envi.dat<-read.csv("UKmoths/gridXData.csv")

#traits
tr.dat<-read.csv("UKmoths/moths/processed/TrTaxoData.csv")
#envi data
x.dat<-read.csv("UKmoths/moths/processed/XData.csv")
#presence absence
y.dat<-read.csv("UKmoths/moths/processed/YData.csv")

library(ncdf4)
#================
#Update UK data
# HadUK grid daily temperature and precipitation data is available here https://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb 
#https://github.com/RebeccaLovell/OrangeTipAnalyses/blob/main/DataPreparation.R

#use bioclim? https://www.worldclim.org/data/bioclim.html

#https://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb/
#tmax<- nc_open("./UKclimate/tasmax_hadukgrid_uk_60km_seas-20y_198101-200012.nc")
#tmax<- terra::rast("./UKclimate/tasmax_hadukgrid_uk_60km_seas-20y_198101-200012.nc")

#HadUK provisional data
#https://www.metoffice.gov.uk/hadobs/hadukgrid/data/download_2025-05.html

tmax<- nc_open("./UKclimate/tasmax_hadukgrid_uk_5km_day_20250501-20250531.nc")

time <- ncvar_get(tmax,"time");# extracts time, days 
lats<- ncvar_get(tmax,"projection_x_coordinate")
lons<- ncvar_get(tmax,"projection_y_coordinate")
tasmax<- ncvar_get(tmax,"tasmax")

tmax<- terra::rast("./UKclimate/tasmax_hadukgrid_uk_5km_day_20250501-20250531.nc")
head(names(tmax))
tasmax<- tmax[["tasmax_1"]]

plot(tasmax)

#provisional data
#https://www.metoffice.gov.uk/hadobs/hadukgrid/data/download_2025-05.html

#================
#UK butterfly analysis from Buckley et al 2011 Ecology

if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/SDMpriors/data/")
if(desktop=="n") setwd("/Users/lbuckley/My Drive/Buckley/Work/SDMpriors/data/")

#load butterfly physiolgy data
Spec<-read.csv("./UKBMS/ButterflyDegreeDays/LDTs_0504.csv")

yr.range<-c("1970-82","1995-9","2000-4")
years<-c("19741978", "19951999", "20002004")

# climate data: mean temp of coldest month and mean annual precip (Luoto 2006) 
temp<-raster("./UKBMS/ButterflyDegreeDays/climate/Temp_coldestmonth.asc")
pre<-raster("./UKBMS/ButterflyDegreeDays/climate/AnnualPrecip.asc")

#load degree day estimates
dd.gen<-function(dd){dd/ddgen}

yeark<-1
speciesk<-5

  #load distribution data
  file<- paste("./UKBMS/ButterflyDegreeDays/distribution/DistData_MaxentFormat_UKgrid/",Spec$CommonName[speciesk]," ",yr.range[yeark]," .csv", sep='',collapse=NULL)
  xy.dat<-read.csv(file)
  xy.dat<-(unique(xy.dat, MARGIN=1))
  file<- paste("./UKBMS/ButterflyDegreeDays/distribution/DistData_MaxentFormat_UKgrid/",Spec$CommonName[speciesk]," ",yr.range[2]," .csv", sep='',collapse=NULL)
  xy.dat1995<-read.csv(file)
  xy.dat1995<-(unique(xy.dat1995, MARGIN=1))
  file<- paste("./UKBMS/ButterflyDegreeDays/distribution/DistData_MaxentFormat_UKgrid/",Spec$CommonName[speciesk]," ",yr.range[3]," .csv", sep='',collapse=NULL)
  xy.dat2000<-read.csv(file)
  xy.dat2000<-(unique(xy.dat2000, MARGIN=1)) 
  
  #load degree day data
  file<-paste("./UKBMS/ButterflyDegreeDays/climate/degreeday/UKDegreeDays_Apr1Oct1_5yrMean_7082/",Spec$Species[speciesk],".asc",sep="",collapse=NULL)
  dd<-raster(file) 
  #mask to get rid of zeros
  dd<-mask(dd, pre) 
  
  file<-paste("./UKBMS/ButterflyDegreeDays/climate/degreeday/UKDegreeDays_Apr1Oct1_5yrMean_9599/",Spec$Species[speciesk],".asc",sep="",collapse=NULL)
  dd1995<-raster(file) 
  #mask to get rid of zeros
  dd1995<-mask(dd1995, pre) 
  
  file<-paste("./UKBMS/ButterflyDegreeDays/climate/degreeday/UKDegreeDays_Apr1Oct1_5yrMean_0004/",Spec$Species[speciesk],".asc",sep="",collapse=NULL)
  dd2000<-raster(file) 
  #mask to get rid of zeros
  dd2000<-mask(dd2000, pre) 
  
  #divide by dd per generation
  ddgen<-Spec$DDlarval[speciesk]
  dd<-calc(dd, dd.gen, ddgen=ddgen)
  dd1995<-calc(dd1995, dd.gen, ddgen=ddgen)
  dd2000<-calc(dd2000, dd.gen, ddgen=ddgen)
  
#plot
  plot(dd,axes=FALSE, nlevel=50, zlim=range(2,8), main="1970-1982",col=rev(heat.colors(22)), horizontal=TRUE, xlim=range(000,680000), ylim=range(0,1250000))
  points(xy.dat[,-1],col="#0000ff22") #color makes transparent points

#==============================  
 #IBIS iSDM model
  
  # Load the package
  library(ibis.iSDM)
  library(inlabru)
  library(xgboost)
  library(terra)
  library(uuid)
  library(assertthat)
  library(sf)
  
  # Don't print out as many messages
  options("ibis.setupmessages" = FALSE)
  
  #convert to terra spatial raster
  temp<- terra::rast(temp)
  pre<- terra::rast(pre)
  dd<- terra::rast(dd)
  
  #make predictor stack
  #predictors<- rast(temp, pre, dd)
  predictors<- c(temp, pre, dd) #makes list
  
  # First we define a distribution object using the background layer
  background<- terra::mask(temp, temp, inverse=TRUE, maskvalue=NA, updatevalue=1)
  
  #set projection
  crs(background) <- "+proj=utm +zone=30"
  
  mod <- ibis.iSDM::distribution(background)
  
  # Load species points
  xy.dat$Observed=1
  
  pts=st_as_sf(xy.dat[,-1], coords=c("Longitude","Latitude"), crs="+proj=utm +zone=30")
  
  # This data needs to be in sf format and key information is that
  # the model knows where occurrence data is stored (e.g. how many observations per entry) as
  # indicated by the field_occurrence field.
  mod <- add_biodiversity_poipo(mod, pts,
                                name = "Virtual test species",
                                field_occurrence = "Observed")
  
  # Then lets add predictor information
  # Here we are interested in basic transformations (scaling), but derivates (like quadratic)
  # for now, but check options
  mod <- add_predictors(mod, 
                        env = predictors,
                        transform = "scale", derivates = "none")
  
  # Finally define the engine for the model
  # This uses the default data currently backed in the model,
  # !Note that any other data might require an adaptation of the default mesh parameters used by the engine!
  mod <- engine_inlabru(mod)
  
  # Print out the object to see the information that is now stored within
  print(mod)
  
  #Create model
  mod <- distribution(background) |> 
    add_biodiversity_poipo(pts,
                           name = "Virtual test species",
                           field_occurrence = "Observed") |>  
    add_predictors(env = predictors, transform = "scale", derivates = "none") |>
    engine_inlabru() 
  
  # Make visualization of the contained biodiversity data
  plot(mod$biodiversity)
  
  # Other options to explore
  names(mod)
  
  #--------------------------
  # Define prior
  # In this case and for the INLA engine we define normal prior on the mean and precision
  # https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.html
  # Required parameters are a mean and a precision estimate provided to "hyper". Note that precision is not equivalent (rather the inverse) to typical standard deviation specified in Gaussian priors. 
  p <- INLAPrior(variable = "layer",type = "normal",hyper = c(2, 10))
  # This is then wrapped in a PriorList
  pp <- priors(p)
  print( pp )
  
  ##Create model
  #mod <- distribution(background) |> 
  #  add_biodiversity_poipo(pts,
  #                         name = "Virtual test species",
  #                         field_occurrence = "Observed") |>  
  #  add_predictors(env = predictors, transform = "scale", derivates = "none") |>
  #  add_priors(priors = pp) |>
  #  engine_inlabru() 
  #---------------------------
  
  # Finally train
  fit <- train(mod,
               runname =  "Test INLA run",
               aggregate_observations = FALSE, # Don't aggregate point counts per grid cell
               verbose = TRUE # Don't be chatty
  )
  
  # Plot the mean of the posterior predictions
  plot(fit, "mean")
  
  # Print out some summary statistics
  summary(fit)
  
  # Show the default effect plot from inlabru
  effects(fit)
  
  # To calculate a partial effect for a given variable
  o <- partial(fit, x.var = "layer", plot = TRUE)
  
  # The object o contains the data underlying this figure
  
  # Similarly the partial effect can be visualized spatially as 'spartial'
  s <- spartial(fit, x.var = "layer")
  plot(s[[1]], col = rainbow(10), main = "Marginal effect of dd on the relative reporting rate")
  
  # Calculate a threshold based on a 50% percentile criterion
  fit <- threshold(fit, method = "percentile", value = 0.5)
  
  # Notice that this is now indicated in the fit object
  print(fit)
  
  # There is also a convenient plotting function
  fit$plot_threshold()
  
  # It is also possible to use truncated thresholds, which removes non-suitable areas
  # while retaining those that are suitable. These are then normalized to a range of [0-1]
  fit <- threshold(fit, method = "percentile", value = 0.5, format = "normalize")
  fit$plot_threshold()
  
 #------------
#Compare with and without priors

  mod1 <- train(mod,
               runname =  "Test INLA run",
               aggregate_observations = FALSE, # Don't aggregate point counts per grid cell
               verbose = TRUE # Don't be chatty
  )
  
      # Add Priors
  mod2 <- train(mod |> add_priors(pp), only_linear = TRUE)
  
  # Compare the difference in effects
  p1 <- partial(mod1, pp$varnames(), plot = TRUE)
  p2 <- partial(mod2, pp$varnames(), plot = TRUE)
  
  #summarize and evaluate
  plot(mod1)
  #summarize coefficients
  summary(mod1)
  
  #Continuous validations use error metrics (e.g. RMSE) to infer prediction precision (Jung, 2022), while discrete validations can be calculated on a-priori mapped thresholded distributions with a range of different options from binary to normalized estimation
  #area under the curve (AUC) and true skill statistic (TSS)
  validate(mod1)
  validate(mod2)
  #validate the model with independent data
  #validate(mod1, point=XXXX)
  
  coef(mod1)
  limiting(mod1)
  
  
  #Jung 2023
  #https://doi.org/10.1016/j.ecoinf.2023.102127
  #Another idea is to enable support for dedicated equations, for example for population growth or microclimatic thresholds (Schouten et al., 2020), and integrate them into inference and projections (Talluto et al., 2016).

  #-------
  #make pseudoabsent for test, train
  
  #----------------------------
  #generate pseudo absence
 
  # define circles with a radius of 50 km around the subsampled points
  x = circles(xy.dat[,c("Longitude","Latitude")], d=50000, lonlat=T)
  # draw random points that must fall within the circles in object x
  bg = spsample(x@polygons, 1000, type='random', iter=100)
  
  # #extract climate variables
  # # pulling bioclim values
  # occ_bc = extract(BClim, occ[,c("decimalLongitude","decimalLatitude")] ) # for the subsampled presence points
  # bg_bc = extract(BClim, bg) # for the pseudo-absence points
  # occ_bc = data.frame(lon=occ$decimalLongitude, lat=occ$decimalLatitude, occ_bc)
  # bgpoints = bg@coords
  # colnames(bgpoints) = c("lon","lat")
  # bg_bc = data.frame(cbind(bgpoints,bg_bc))
  # 
  # # Create dataframe from bioclim and presense/absance.
  # pres<-rep(1,dim(occ_bc)[1])
  # temp1<-data.frame(pres,occ_bc[,3:21])
  # pres<-rep(0,dim(bg_bc)[1])
  # temp2<-data.frame(pres,bg_bc[,3:21])
  # df<-rbind(temp1,temp2)
  # head(df,5)
  # 
  # locs= occ_bc
  # 
  # #------------
  # #split into test and train
  # set.seed(3451)
  # inds <- partition(xy.dat$Observed, p = c(train = 0.8, test = 0.2))
  # 
  # #or use folds: folds <- create_folds(train$Sepal.Length, k = 5, seed = 2734)
  # 
  # str(inds)
  # #> List of 3
  # #>  $ train: int [1:81] 2 3 6 7 8 10 11 18 19 20 ...
  # #>  $ valid: int [1:34] 1 12 14 15 27 34 36 38 42 48 ...
  # #>  $ test : int [1:35] 4 5 9 13 16 17 25 39 41 45 ...
  # 
  # train <- iris[inds$train, ]
  # valid <- iris[inds$valid, ]
  # test <- iris[inds$test, ]
  # 
  