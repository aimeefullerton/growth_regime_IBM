#*************************************************************************************************
# Growth Regime IBM, VERSION 1.0.1
#     An individual based model to evaluate the contribution of seasonally warm habitat to  
#     Oncorhynchus mykiss growth, and used to generate an example presented in:
#
#     Armstrong, J.B., A.H. Fullerton, C.E. Jordan, J.L. Ebersole, J.R. Bellmore, I. Arismendi, 
#         B. Penaluna, and G.H. Reeves. The significance of warm habitat to the growth regime of 
#         coldwater fishes.
#
# Adapted from:
#     Fullerton, A.H., B.J. Burke, J.J. Lawler, C.E. Torgersen, J.L. Ebersole, and S.G. Leibowitz. 
#         2017. Simulated juvenile salmon growth and phenology respond to altered thermal regimes 
#         and stream network shape. Ecosphere 8(12):e02052. [Original model.]
#     Hawkins, B.L., A.H. Fulleton. B.L. Sanderson, and E.A. Steel. 2020. Individual-based simulations 
#         suggest mixed impacts of warmer temperatures and a non-native predator on Chinook salmon.
#         Ecosphere 11(8):e03218. [Updated movement rules.]
#
# Last functional updates 13 June 2019; editorial edits 15 Oct 2020
# See development notes at the end of this script
#*************************************************************************************************

# === SETUP =====================================================================================
rm(list = ls()) #clear workspace
start.time <- proc.time() #get initial time stamp for calculating processing time
options(scipen = 999)

# Load required libraries
library(SSN)
library(rgdal)
library(plyr)
library(RColorBrewer)
library(gam)

# Load Functions
source("code/growth_regime_functions_v1.0.1.R")

# List of simulation replicates (iterations) to run (if longer than one, will need to un-comment the brackets for looping below)
iter.list <- 10 #1:100  
iter <- iter.list[1] 
plot.iter <- 99 #select iteration(s) to plot maps for fish and water temperature at each time step and some diagnostic plots for certain fish
run <- fncGetRun()

# DIRECTORY STRUCTURE
  loadDir<- "data.in"
    ifelse(!dir.exists(file.path(loadDir)), dir.create(file.path(loadDir)), FALSE)
  dataDir<- "data.out"
    ifelse(!dir.exists(file.path(dataDir)), dir.create(file.path(dataDir)), FALSE)
  plotDir<- "plots"
    ifelse(!dir.exists(file.path(plotDir)), dir.create(file.path(plotDir)), FALSE)
  aniDir<- "plots.ani"
    ifelse(!dir.exists(file.path(aniDir)), dir.create(file.path(aniDir)), FALSE)

# SSN-specific fields (different for GIS-based and generated streams)
length.field <- "Length" #reach length
wt.field <- "WT" #water temperature
so.field <- "shreve" #stream order
xlat <- "NEAR_X" #x-latitude
ylon <- "NEAR_Y" #y-longitude

# Color palette for plotting water temperatures
cb <- fncColBrewPlus(n = 14, paint = F)


# === DEFINE PARAMETER VALUES ===================================================================

# FISH PARAMETERS
  spp <- "steelhead"
  JulianDate.Begin <- 1 #1 January
  numFish <- 500 #number of fish in a network per scenario
  corr.factor <- 5000 #expansion from numFish ("superorganism") to estimate real fish count

# Movement parameters used in fish movement functions
  mvmt.multiplier <- 1 #for sensitivity analysis, to expand or shrink movement distances
  mvdst.sdlog <- 4 #standard deviation for lognormal distribution of target movement distances
  move.cost <- 0.2 #cost to growth of moving (movement uses energy)

# Bioenergetics 
  Constants <- fncReadConstants.steelhead_ration()
  st.wt <- 100 #starting weight of individuals
  MaxDensity4Growth <- 1 #growth is not depressed further above this density
  ration.lo <- 0.02 #lowest ration in the network (headwaters)
  ration.hi <- ration.lo + 0.05 # highest ration in the network (river outlet)
  pred.en.dens <- 5900 #predator energy density for bioenergetics
  oxy <- 13560 #oxygen consumed for bioenergetics
  pff <- 0.1 #percent indigestible prey for bioenergetics
  prey.en.dens <- 3500 #prey energy density for bioenergetics
  #this value from Beauchamp et al. 2007, which was for marine prey
  #Thompson and Beauchap 2016 have similarly high energy densities for aquatic prey, mean 3276, 3422 in 2 Skagit River streams
  #2500 J/g in Railsback and Rose 1999 

# THERMAL REGIMES
  # Based on 'cool slow warm stable' thermal regime in Fullerton et al. 2017, which is based on observations from the Snoqualmie River, western WA
  # increased to represent a warmer system for this application, and smoothed
  tr730 <- read.csv("data.in/thermal.regime.730ts.csv") # values for each 12-h time step (at outlet reach)
  ifelse(JulianDate.Begin > 1, timeStep <- c(seq(JulianDate.Begin + 0.1, 365.6, 0.5), seq(1.1, JulianDate.Begin, 0.5)), timeStep <- seq(1.1, 365.6, 0.5))

# NETWORK CONFIGURATIONS
  netnm <- "network-swh" 
  missed <- NULL #list of segments with no water temperature predictions (may happen for GIS-generated streams)

  # #Network(s) already exist and will be loaded below
  #  #Optional: re-create network
  #     numSegs <- 151 #number of segments in a network
  #     s <- 123456792; set.seed(s); netnm <- "swh.net"
  #     ssn1 <- fncGenerateNetwork(numSegs, lmin = 0.5, lmax = 1, sampleprob = 1, ps = NA, ordermax = 3, wt.spc = 0.1, fish.spc = 0, numFish = 10000, plotit = "file", path = paste0(loadDir, "/", netnm, ".ssn"), seed = s)
  #     #this takes a long time with many fails, but produces a good network in the end!
  #   #To recreate lists that give, for each segment, (1) all downstream segments, (2) all upstream segments, and (3) list of all segments at each junction
  #      ssn1 <- importSSN(paste0(loadDir, "/", netnm, ".ssn")); obs.df <- getSSNdata.frame(ssn1, "Obs")
  #      dnsegs <- fncDnSegs(ssn1, path = paste0(loadDir, "/", netnm, ".ssn"))
  #      upsegs <- fncUpSegs(ssn1, path = paste0(loadDir, "/", netnm, ".ssn"))
  #      jct.list <- fncJctLst(ssn1, path = paste0(loadDir, "/", netnm, ".ssn"))
  
# === SCENARIOS & STARTUP ================================================

  # FOOD AVAILABILITY SCENARIOS
  food.scenarios <- c("VariFood", "ConstFood")
  fs <- food.scenarios[1]
  
  # MANAGEMENT SCENARIOS
  mgmt.scenarios <- c("Base", "DivestSWH", "EnhancePCH", "DivSWH.EnhPCH")
  ms <- mgmt.scenarios[1]
  divest.prop <- 0.1 #proportion to reduce ration in seasonally warm habitat reaches (if this scenario is used; has no effect on 'Base')
  enhance.prop <- 0.1 #proportion to increase ration in perennially cool habitat reaches (if this scenario is used; has no effect on 'Base')

  (scenario.nm <- paste0(ms, ".", fs))
  
  # Import the appropriate existing SSN network to get structure info
  ssn1 <- importSSN(paste0(loadDir, "/", netnm, ".ssn"))
  nloc <- nrow(getSSNdata.frame(ssn1, "Obs"))
  rm(ssn1)
  
  # Load topology lists for the network (previously created)
  load(paste0(loadDir, "/", netnm, ".ssn/dnsegs.RData"))
  load(paste0(loadDir, "/", netnm, ".ssn/upsegs.RData"))
  load(paste0(loadDir, "/", netnm, ".ssn/jct.list.RData"))
  
  # Load pre-calculated growth
  load(paste0(loadDir, "/wt.growth.array.RData")) 
  # dimensions are 500 WT (0-25 C), 151 rations (0.02-0.17 g/g/d), 700 fish weights (1-700 g)
  # created with 'pre-calculate_growth.R'
  
  # Create arrays that will store data:
  fa <- array(NA, dim = c(numFish, 13, 730, length(iter.list))) #250 fish, no. variables, 730 time steps, no. of iterations 
  WT <- array(NA, dim = c(730, nloc, length(iter.list))) #730 time steps, N locations, no. of iterations
  
  # Loop over iterations in simulation 
  #  for(iter in iter.list){ #(un-comment this line and its closing bracket if running more than one scenario at a time)

  
# === INITIALIZE HABITAT ========================================================================

# Seed for reproducible results
  set.seed(1) 
    
# 1. Create spatial structure of water temperature observations on the network (stochastic)
  ssn1 <- importSSN(paste0(loadDir, "/", netnm, ".ssn")) #load ssn
  ssn1 <- importPredpts(ssn1, "preds", "ssn") #import the water temperature points
  ssn1 <- fncInitializeWT(ssn = ssn1, plotit = "screen") #initialize water temperature for first time step
  WT.init <- obs.df <- ssn1@obspoints@SSNPoints[[1]]@point.data
  # change factor/character fields to numeric (format depends on version of RStudio)
  if(is.factor(WT.init$rid)) WT.init$rid <- as.numeric(levels(WT.init$rid)[as.numeric(WT.init$rid)])  # factor
  if(is.character(WT.init$rid)) WT.init$rid <- as.numeric(WT.init$rid) # character
  
  # Add some especially cold reaches
  CWR <- c(2, 6, 11, 31, 32, 39, 56, 79, 93)
  WT.init$scaledWT[WT.init$rid %in% CWR] <- WT.init$scaledWT[WT.init$rid %in% CWR] * 0.2
  obs.df <- WT.init
  ssn1 <- putSSNdata.frame(WT.init, ssn1, "Obs")
 
  # Get reach widths, adjusted to represent the proportion useable by fish, 
  # which decreases as width increases
  ssn1@data$WIDTH_M <- ssn1@data$addfunccol * 500
  useable.widths <- fncUseableWidths(dat = ssn1@data[, c("rid", "WIDTH_M")])
  ssn1@data$UseableWidth <- useable.widths[useable.widths[,1] == ssn1@data$rid, 3]
  rm(useable.widths)
  
  # Add column to denote which reaches are accessible to fish (for blocking movement into ephemeral headwater reaches)
  ssn1@data$accessible <- 0
  ssn1@data$accessible[ssn1@data$shreve >= 16] <- 1 #Accessible based on shreve order
  ssn1@data$accessible[ssn1@data$rid %in% CWR] <- 1 #Additional accessible cold reaches
  
  
# 2. Define ration based on the food and management scenarios
  
  # Food scenarios
  if(fs == "ConstFood"){
    # Constant ration (g/g/d) available across stream network
    obs.df$ration <- ssn1@obspoints@SSNPoints[[1]]@point.data$ration <- rep(0.055, nrow(obs.df)) 
  }
  
  if(fs == "VariFood"){
    # Calculate ration (g/g/d) available to fish which varies across the network
    # Assumption: more productive food webs in lower elevation habitats
    ration_1 <- fncRescale(log(obs.df$shreve + (1 - obs.df$ratio)), c(ration.lo, ration.hi)) #linearly relate to stream order (shreve) and position within reach (ratio)
    ration_2 <- fncRescale(obs.df$upDist, c(ration.hi, ration.lo)) #inverse-linearly relate to distance upstream from river outlet (upDist)
    obs.df$ration <- ssn1@obspoints@SSNPoints[[1]]@point.data$ration <- (ration_1 + ration_2) / 2 #mean of these two processes
    #plot(ssn1, "ration", breaktype = "even", nclasses = 6, xlab = "", ylab = "")
  }
  
  # Management scenarios
  if(! ms %in% c("DivestSWH", "EnhancePCH", "DivSWH.EnhPCH")){ divest.prop <- enhance.prop <- 0}
  
  # "Divest" in seasonally warm habitats by cranking down food availability by 'divest.prop' in those reaches to represent habitat quality degradation
  if(ms %in% c("DivestSWH", "DivSWH.EnhPCH")){
    swh.segs <- c(1, 3, 10, 7, 22, 33, 15, 42, 49, 14, 57, 21, 66, 78, 30)
    obs.df$ration[obs.df$rid %in% swh.segs] <- obs.df$ration[obs.df$rid %in% swh.segs] * (1 - divest.prop)
    ssn1@obspoints@SSNPoints[[1]]@point.data$ration[ssn1@obspoints@SSNPoints[[1]]@point.data$rid %in% swh.segs] <- 
      ssn1@obspoints@SSNPoints[[1]]@point.data$ration[ssn1@obspoints@SSNPoints[[1]]@point.data$rid %in% swh.segs] * (1 - divest.prop)
  }
  
  # "Enhance" coldwater habitats by cranking up food availability by 'enhance.prop' in those reaches to represent habitat quality improvement
    # (e.g., with nutrient additions, reduced sportfishing, improved physical habitat, etc.)
  if(ms %in% c("EnhancePCH", "DivSWH.EnhPCH")){
    pch.segs <- CWR
    obs.df$ration[obs.df$rid %in% pch.segs] <- obs.df$ration[obs.df$rid %in% pch.segs] * (1 + enhance.prop)
    ssn1@obspoints@SSNPoints[[1]]@point.data$ration[ssn1@obspoints@SSNPoints[[1]]@point.data$rid %in% pch.segs] <- 
      ssn1@obspoints@SSNPoints[[1]]@point.data$ration[ssn1@obspoints@SSNPoints[[1]]@point.data$rid %in% pch.segs] * (1 + enhance.prop)
  }
  
  # Save this 'base case' ration because we will be updating ration later depending on conspecific fish density in a reach
  obs.df$ration_base <- ssn1@obspoints@SSNPoints[[1]]@point.data$ration_base <- obs.df$ration 
  
  # Initialize fish density
  obs.df$density <- ssn1@obspoints@SSNPoints[[1]]@point.data$density <- 0
    
  # Store water temperatures through time for each site
  WT.ts <- matrix(NA, length(timeStep), length(obs.df[,1]))
  WT.ts[1,] <- t(WT.init$WT)
  


# === INITIALIZE FISH ===========================================================================

# 1. Select initial locations and import the locations into the SSN
  
  # Use pred points as basis for selecting locations
  fish.pts <- readOGR(paste0(loadDir, "/", netnm, ".ssn"), "preds") # load preds shapefile
  if(is.factor(fish.pts@data$rid)) fish.pts@data$rid <- as.numeric(levels(fish.pts@data$rid)[as.numeric(fish.pts@data$rid)])  #changing factor fields to numeric, if necessary
  if(is.character(fish.pts@data$rid)) fish.pts@data$rid <- as.numeric(fish.pts@data$rid)  #changing character fields to numeric, if necessary
  fish.pts$seg <- fish.pts$rid
  # limit to accessible reaches
  fish.pts <- fish.pts[fish.pts$rid %in% ssn1@data$rid[ssn1@data$accessible == 1],]
  # assign closest water temp values
  fish.pts$WT <- as.numeric(fncGetNearestWT(WT.init[,c("pid", "rid", "ratio", wt.field)], fish.pts@data[,c("pid", "seg", "ratio")], ssn1)[,"WT"]) #assign nearest water temperature to each fish
  fish.pts$scaledWT <- as.numeric(fncGetNearestWT(WT.init[,c("pid", "rid", "ratio", "scaledWT")], fish.pts@data[,c("pid", "seg", "ratio")], ssn1)[,"WT"]) #scaled WT (before diel variation)
  
  # Generate a distribution of locations where fish will start (based on scaledWT)
  ss <- rnorm(numFish, 0.5, 0.2) 
  
  # Find closest location that matches the value from the generated distribution:
  fswt <- fish.pts$scaledWT; names(fswt) <- fish.pts$pid
  sss <- NULL; idx2 <- NULL; i <- 1
  for(s in ss){
    idx <- which.min(abs(fswt - s)) 
    sss[i] <- fswt[idx]
    idx2[i] <- as.numeric(names(fswt[idx]))
    fswt <- fswt[-idx] #remove, so can't be sampled twice
    i <- i + 1
  }
  fish.pts <- fish.pts[fish.pts$pid %in% idx2,]
  rownames(fish.pts@data) <- fish.pts@data$pid 
  rm(fswt, idx, idx2, s, ss, sss)

  # Save/overwrite "fish" shapefile in SSN folder & remove temporary shapefile
  writeOGR(fish.pts, paste0(loadDir, "/", netnm, ".ssn"), driver = "ESRI Shapefile", "fish", overwrite_layer = TRUE)
  rm(fish.pts)

  # Load into SSN
  ssn1 <- importPredpts(ssn1, "fish", "ssn")
  
# 2. Create fish data frame for manipulating during simulation
  #  pid: unique identification for each fish (prediction site id)
  #  seg: which network segment is the fish in (same as rid)
  #  ratio: relative distance along the segment (from bottom to top [0-1])
  #  xloc: x coordinate of fish; yloc: y coordinate of fish
  #  length2segBase is the distance from the current location to the base of that segment
  #  upDist is the distance from the base of the whole network
  
  DFfish<-getSSNdata.frame(ssn1, "fish") #get data frame from SSN
  # Changing factor fields to numeric, if necessary
  if(is.factor(DFfish$rid)) DFfish$rid <- as.numeric(levels(DFfish$rid)[as.numeric(DFfish$rid)])
  if(is.factor(DFfish$WT)) DFfish$WT <- as.numeric(levels(DFfish$WT)[as.numeric(DFfish$WT)])
  if(is.character(DFfish$rid)) DFfish$rid <- as.numeric(DFfish$rid)
  if(is.character(DFfish$WT)) DFfish$WT <- as.numeric(DFfish$WT)
  
  fish <- data.frame(pid = 1:nrow(DFfish), seg = 0, ratio = 0, xloc = 0, yloc = 0, length2segBase = 0, upDist = 0, WT = 0, TU = 0, emrg = 1, weight = st.wt, 
          growth = 0, consCum = 0, direction = rep("up", nrow(DFfish)), dateEm = 0, dateOm = 0, density = 0, movedist = 0, consInst = 0, pvals = 0, survive = 1, 
          shreve = 0, ration = 0, stringsAsFactors = FALSE) #create new dataframe that will hold results
  # Fill in what we can from the ssn
  fish$seg <- DFfish$rid
  fish$ratio <- DFfish$ratio
  fish$xloc <- DFfish[, xlat]; fish$yloc <- DFfish[, ylon]
  fish$length2segBase <- ddply(fish, .(pid), function(x) data.frame(dist=segRatio2length(x$seg, x$ratio)))$dist
  fish$upDist <- DFfish$upDist
  fish$WT <- DFfish[, wt.field]
  fish$shreve <- DFfish[, so.field]
  fish$dateSp <- JulianDate.Begin
  fish$growth2 <- fish$growth
  fish$ration <- fncGetNearestWT(obs.df[, c("pid", "rid", "ratio", "ration")], fish[, c("pid", "seg", "ratio")], ssn1)[, "WT"] #assign nearest ration to each fish
  
  plot(ssn1, PredPointsID = "fish", addWithLegend = T, cex = 0.7)
  rm(DFfish)
  
# === RUN 1 YEAR SIMULATION ======================================================================

  # Make temporary array to hold fish data for this time step (number of fish x number of time steps)
  fa.tmp <- array(NA, dim =c (numFish, dim(fa)[2], 730))

  # For each time step, where 3am (t = JulianDay.1) and 3pm (t = JulianDay.6)
  for(tt in 1:length(timeStep)){

  # 1. UPDATE WATER TEMPERATURE.
    # Alter stream temperature in each reach at time step tt
    WT.init <- ssn1@obspoints@SSNPoints[[1]]@point.data 
    ssn1 <- fncUpdateWaterTemps(ts = timeStep[tt], WT.init, tr = tr730, ssn = ssn1, type = "Obs", diel = T, plotit = "none")
    obs.df <- ssn1@obspoints@SSNPoints[[1]]@point.data 
    if(is.factor(obs.df$rid)) obs.df$rid <- as.numeric(levels(obs.df$rid)[as.numeric(obs.df$rid)])
    if(is.character(obs.df$rid)) obs.df$rid <- as.numeric(obs.df$rid) 
    WT.ts[tt,] <- t(obs.df$WT)

    # Update water temperature in 'fish'
    fish$WT <- fncGetNearestWT(obs.df[, c("pid", "rid", "ratio", wt.field)], fish[, c("pid", "seg", "ratio")], ssn = ssn1)[,"WT"]
    
    # Plot updated water temperatures and fish locations (only for select iteration due to storage constraints)
    if(iter == plot.iter){
      
      aniDir1 <- paste0("ani.", iter, ".", spp, ".", scenario.nm)
      ifelse(!dir.exists(file.path(paste0(aniDir, "/", aniDir1))), dir.create(file.path(paste0(aniDir, "/", aniDir1))), FALSE)

      breaks <- c(0, seq(4, 26, length.out = 12), 40)
      cb <- fncColBrewPlus(n = 14, paint = F)
      left <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
      rght <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 27)
      nd <- "d"; if(timeStep[tt] - floor(timeStep[tt]) < 0.2) {nd <- "n"}
      
      if(nd == "d"){
        anifile <- paste0(aniDir, "/", aniDir1, "/", nd, sprintf("%03d", tt), ".png")
        png(anifile, width = 6, height = 6, units = "in", res = 150)
        
        plot(ssn1, "WT", breaktype = "user", brks = breaks, lwdLineCol  =  "addfunccol", lwdLineEx = 8, lineCol = "darkgray", xlab = "", ylab = "", color.palette = cb, main = timeStep[tt])
        alive <- which(fish$survive == 1)
        plotSSN.mod(ssn1, PredPointsID = "fish", cex = 0.8, addWithLegend = TRUE, myvar = alive, col = 1, pch = 20)
        emerged <- which(fish$emrg == 1 & fish$survive == 1)
        if(sum(emerged) > 0){
          plotSSN.mod(ssn1, PredPointsID  =  "fish", cex = 1, addWithLegend = TRUE, myvar = emerged, col = 1, pch = 21)
        }
        # #track a few individuals:
        if(1 %in% emerged) {plotSSN.mod(ssn1, PredPointsID = "fish", cex = 1.8, addWithLegend = TRUE, myvar = 1, col = "purple", pch = 19); plotSSN.mod(ssn1, PredPointsID = "fish", cex = 1.8, addWithLegend = TRUE, myvar = 1, pch = 21)}
        if(10 %in% emerged) {plotSSN.mod(ssn1, PredPointsID = "fish", cex=1.8, addWithLegend = TRUE, myvar = 10, col = 4, pch = 19); plotSSN.mod(ssn1, PredPointsID = "fish", cex = 1.8, addWithLegend = TRUE, myvar = 10, pch = 21)}
        if(20 %in% emerged) {plotSSN.mod(ssn1, PredPointsID = "fish", cex=1.8, addWithLegend = TRUE, myvar = 20, col = 5, pch = 19); plotSSN.mod(ssn1, PredPointsID = "fish", cex = 1.8, addWithLegend = TRUE, myvar = 20, pch = 21)} 
        if(30 %in% emerged) {plotSSN.mod(ssn1, PredPointsID = "fish", cex=1.8, addWithLegend = TRUE, myvar = 30, col = 6, pch = 19); plotSSN.mod(ssn1, PredPointsID = "fish", cex = 1.8, addWithLegend = TRUE, myvar = 30, pch = 21)} 
        if(40 %in% emerged) {plotSSN.mod(ssn1, PredPointsID = "fish", cex=1.8, addWithLegend = TRUE, myvar = 40, col = 2, pch = 19); plotSSN.mod(ssn1, PredPointsID = "fish", cex = 1.8, addWithLegend = TRUE, myvar = 40, pch = 21)}
        dev.off()
      }
    }

  # 2. MOVE FISH.
    # Assess probability of movement and move fish if appropriate
    step <- fncMoveOneStep(ssn1)
    fish <- step[[1]]
    ssn1 <- step[[2]]

    
  # 3. UPDATE HABITAT QUALITY.
    # Re-assess reach quality and fish density after all fish have moved
    
    # Update water temperature in 'fish' at fish's new location
    fish$WT <- fncGetNearestWT(obs.df[, c("pid", "rid", "ratio", wt.field)], fish[, c("pid", "seg", "ratio")], ssn = ssn1)[,"WT"]

    # Update fish density in 'fish'
    fish$density <- fncFishDensity(fish1.seg = fish[,"seg"], fish2.seg = NULL, corr.factor,
             fish1.alive.emerged.segs = seq(1:nrow(fish)),fish2.alive.emerged.segs = NULL, ssn = ssn1)[,1]

    # Get density effect on ration (ration linearly decreases with fish density at the new location)
    fdens <- fish$density; fdens[fdens > MaxDensity4Growth] <- MaxDensity4Growth # high densities cap out at this parameter value
    density.effect <- fncRescale((1 - c(fdens, 0.01, MaxDensity4Growth)), c((0.5),1))
      density.effect <- density.effect[-c(length(density.effect), (length(density.effect) - 1))] #remove the last 2 temporary values
      density.effect <- cbind("rid" = fish$seg, "pid" = fish$pid, density.effect)
      rm(fdens)
    
    # Update ration in reaches occupied by fish after fish have moved and after accounting for density
    for(r in unique(fish$seg)){
      ration1 <- obs.df$ration_base[obs.df$rid ==r]
      ssn1@obspoints@SSNPoints[[1]]@point.data$ration[ssn1@obspoints@SSNPoints[[1]]@point.data$rid ==r] <- 
        obs.df$ration[obs.df$rid ==r] <- ration1 * density.effect[,"density.effect"][density.effect[,"rid"] ==r][1]
    }
    
    # Update ration in 'fish'
    fish$ration <- fncGetNearestWT(obs.df[,c("pid", "rid", "ratio", "ration")], fish[,c("pid", "seg", "ratio")], ssn1)[,"WT"] #this gets ration even though variable says "WT" because 'ration' is passed to function
    
    
  # 4. GROW FISH. (Wisconsin Bioenergetics model)

    # Get parameters and inputs for Bioenergetics model
      # give it nearest water temperature, fish weight (grams), and ration
    Input <- fncGetBioEParms("steelhead", pred.en.dens, prey.en.dens, oxy, pff,
                            wt.nearest = fish[,c("pid", wt.field)],
                            startweights = fish$weight, pvals = rep(1, nrow(fish)), ration = fish$ration) 
    
    # Run bioenergetics for all fish at this 12-hour time step
    Results <- BioE(Input, Constants)
    
    # Add a cost to moving by reducing growth
    gr.ins <- t(Results$Growth / 2) / fish$weight #instantaneous growth (g/g/d)
    gr.mov <- gr.ins - fncMoveCost(md = fish$movedist, gr = gr.ins)

    # Put results back into fish dataframe
      # note: dividing results in half because bioenergetics model operates on 
      # a daily time step and we are using a half-day time step
    fish$growth <- gr.mov #instantaneous growth, discounted for movement
    fish$growth2 <- gr.ins #instantaneous growth, no effect of movement
    
    fish$consInst <- t(Results$Consumption / 2) #instantaneous consumption for that time step only (g/g/d)
    fish$weight <- fish$weight + gr.mov * fish$weight #cumulative growth (g/d)
    fish$consCum <- fish$consCum + fish$consInst #cumulative amount consumed

    # Save fish in a multidimensional array through time
    cols2save <- c("pid", "seg", "WT", "density", "movedist", "shreve", "ration", "consInst", "consCum", "growth", "growth2", "weight", "survive") 
    fa.tmp[,,tt]<- as.matrix(fish[,cols2save])
    
  } #end time steps (tt)

  # Store results
    iii <- iter - iter.list[1] + 1
    fa[,,,iii] <- fa.tmp
    colnames(fa) <- colnames(fish[,cols2save])
    WT[,,iii] <- WT.ts
    
  # Make some diagnostic plots
    filename <- paste0(plotDir, "/", iter, ".", spp, ".", scenario.nm, ".png")
    fncPlotDiagnostics(fa[,,,iii], filename, plotit = "file")
  
  #} #end iteration (un-comment if running more than one iteration)
    
  #Save data
  save("fa", file = paste0(dataDir, "/fa.", iter, ".", spp, ".", scenario.nm, ".RData"))
  save("WT", file = paste0(dataDir, "/WT.", iter, ".", spp, ".", scenario.nm, ".RData"))

  end.time <- proc.time()
  runtime <- (end.time[3] - start.time[3]) / 60

# Store info on no. replicates, scenarios,  runtime, parameters, etc.
  textDir <- paste0(dataDir, "/run.info.", iter, ".", spp, ".", scenario.nm, ".txt")
  file.create(textDir)

  run.info <- c(paste0("netnwork: ", netnm), 
            paste0("species: ", spp), 
            paste0("replicates: ", length(iter.list)), 
            paste0("management scenario: ", ms), 
            paste0("SWH ration reduced by: ", divest.prop),
            paste0("PCH ration enhanced by: ", enhance.prop),
            paste0("food scenario: ", fs), 
            paste0("max density: ", MaxDensity4Growth), 
            paste0("movement multiplier: ", mvmt.multiplier), 
            paste0("move-distance (SDlog): ", mvdst.sdlog), 
            paste0("upper ration: ", ration.hi), 
            paste0("lower ration: ", ration.lo), 
            paste0("first date:", JulianDate.Begin), 
            paste0("runtime (minutes): ", runtime), 
            paste0("iter: ", iter), 
            paste0("run: ", run))
  write(c("RUN INFO:", run.info),file = textDir)


# # Animate IBM results images
# # Mac Version; Requires ffmpeg to be installed already (developer tool)
# ani.path<- paste0(aniDir, "/", aniDir1)
# (dir(ani.path))
# setwd(ani.path)
# system("ffmpeg -framerate 5 -pattern_type glob -i '*.png' -c:v libxvid -b:v 2400k 1Fish.ani.mp4") #other libraries: libx264, libx265, postproc
# setwd("../../")

  
# === END OF FILE ================================================================================


# DEVELOPMENT NOTES:

# Changes in this version that differ from the Fullerton et al. (2017) version are outlined below.
# Updated movement rules are presented in Hawkins et al. (2020).

# SPECIES
# Oncorhynchus mykiss, with species-specific parameters
# All fish begin on same date at 100g
# Starting locations are spread across accessible reaches (shreve order >= 16 plus several smaller really cold reaches)

# NETWORK SHAPE & THERMAL REGIME
# 'Virtual' network shape with short (maxlength = 1) and even-sized stream reaches; N = 151
# Attributed each reach with widths based on drainage area and estimated 'useable widths', which are smaller in mainstems than in tribs
# Thermal regime is based on Snoqualmie River, begins 1 Jan
# Set constant seed to always create same thermal landscape
# Created some extra-cold perennially cold tributary reaches
# Diel temperature variation enabled

# MOVEMENT
# Target movement distance is determined by growth potential (longer if conditions are poor for growth, shorter if closer to optimum) 
# Fish can "sense" optimal habitat and stop to grow there instead of moving their pre-allocated target distance (operates each time step)
# Sensitivity analysis showed that the ability to stop early is more influential than initial target movement distance
# Distance of rare (and all) fish movements increased
# Removed effect of fish size on movement (had assumed bigger fish could swim farther, but doesn't mean they will and larger fish are territorial) 
# Fish always choose the reach with the highest growth potential at tributary junctions (this is the only deterministic movement rule; others are stochastic)
# Fish density affects movement via its effect on growth potential (fish move more if it is crowded in a reach)
# Parameters affecting density tuned to produce reasonable results
# Fish cannot move into inaccessible reaches
# Outmigration functionality disabled; all fish remain in the stream network for the duration of simulations

# GROWTH
#	Maps of food availability simulated as either:
# 'ConstFood': same ration available in all reaches
# 'VariFood': scales with distance from river outlet and stream order
# Fish density affects growth via its effect on ration (reduced prey availability if there are more fish in a reach)
# Growth penalized by movement distance (swim cost)
# Used adult steelhead parameters for bioenergetics because modeled fish are 100g

# SURVIVAL
#	Disabled; all fish survive (if enabled, survival is adjusted to be slightly higher when fish bellies are full and lower when empty based on last time step)

# ACCOUNTING & DISPLAY
#	Track thermal characteristics of the habitats fish use so we could quantify how much of their growth comes from particular habitats
# New diagnostic output plots

