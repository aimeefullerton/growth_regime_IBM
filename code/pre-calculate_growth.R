## Pre-calculate O. mykiss growth at each temperature that can be looked up instead of having to run bioenergetics on the fly each time
# Last updated 15 October 2020

# Load Functions
source("code/growth_regime_functions_v1.0.1.R")

# Parameters
spp <- "steelhead"
numFish <- 500 #number of fish in a network per scenario
Constants <- fncReadConstants.steelhead_ration()
st.wt <- 100 #starting weight
ration.lo <- 0.02
ration.hi <- ration.lo + 0.05
pred.en.dens <- 5900 #predator energy density
oxy <- 13560 #oxygen consumed
pff <- 0.1 #percent indigestible prey
prey.en.dens <- 3500 # see Railsback and Rose 1999, Beauchamp et al. 2007

  
# Set up arrays
waterTemps <- cbind("pid" = seq(1, length(seq(0.05, 25, 0.05))),"WT" = seq(0.05, 25, 0.05))
rations <- seq(0.02, 0.17, 0.001)
weights <- seq(1, 1500, 1)
wt.growth <- array(NA, dim = c(nrow(waterTemps), length(rations), length(weights))) #dim = c(500 WTs, 151 rations, 700 weights)
#respiration <- consumption <- wastes <- wt.growth

# Loop through arrays to record growth for each case
for(w in 1:dim(wt.growth)[3]){
  ww <- weights[w]
  for(r in 1:dim(wt.growth)[2]){
    rr <- rations[r]
    Input<- fncGetBioEParms("steelhead", pred.en.dens, prey.en.dens, oxy, pff,
                            waterTemps, 
                            startweights = rep(ww, nrow(waterTemps)), 
                            pvals = rep(1, nrow(waterTemps)),
                            ration = rep(rr, nrow(waterTemps)))
    
    
    Results <- BioE(Input, Constants) #run BioE code given input and constants
    wt.growth[,r,w] <- c(Results$Gg_WinBioE)
    # consumption[,r,w] <- c(Results$Consumption)
    # wastes[,r,w] <- c(Results$Egestion + Results$Excretion + Results$S.resp)
    # respiration[,r,w] <- c(Results$Respiration)
    
    ## In Joules instead of grams:
    # wt.growth[,r,w] <- c(Results$Growth_j)
    # consumption[,r,w] <- c(Results$Consumption_j)
    # wastes[,r,w] <- c(Results$Egestion_j + Results$Excretion_j + Results$Sj.resp)
    # respiration[,r,w] <- c(Results$Respiration_j)
  }
}

# Save outputs. Only 'wt.growth' is needed for simulations. The others are merely for plotting below.
save("wt.growth", file = "data.in/wt.growth.array.RData")
# save("respiration", file = "data.in/respiration.array.RData")
# save("wastes", file = "data.in/wastes.array.RData")
# save("consumption", file = "data.in/consumption.array.RData")



# Load back in to plot
load("data.in/wt.growth.array.RData")
# load("data.in/consumption.array.RData")
# load("data.in/wastes.array.RData")
# load("data.in/respiration.array.RData")


# Create multi-panel plot of growth curves for different sized fish
png("growth.curves.png", 8, 6, "in", res = 300)
par(mfrow = c(2, 3), las = 1)
ylb <- "Rates (g/g/d)"
for(fish.mass in c(1, 10, 100, 250, 500, 1100)){
  plot(waterTemps[,"WT"], wt.growth[, 151, fish.mass], type='l', ylab = ylb, xlab="Water temperature (C)", 
       col = 4, ylim = c(-0.04, 0.04), main = paste0("fish weight =", fish.mass," g")) #max ration
  lines(waterTemps[,"WT"], wt.growth[, 1, fish.mass], col = 2) #min ration
  for(i in seq(3, 150, length.out = 10)){ lines(waterTemps[,"WT"], wt.growth[,i, fish.mass], col = "gray30", lty = 2)}
  abline(h = 0, lty = 2)
}
dev.off()

# # Create multi-panel 'growth banana' plots for different sized fish
# # need to have calculated consumption, wastes, and respiration arrays (un-comment to run)
# png("growth.banana.png", 8, 6, "in", res = 300)
# par(mfrow = c(2,3), las = 1)
# ylb <- "Rates (g/g/d)"
# for(fish.mass in c(1, 10, 100, 250, 500, 1100)){
#   plot(waterTemps[,"WT"], consumption[,151, fish.mass], type='l', ylab = ylb, xlab="Water temperature (C)",
#        col = 4, ylim = c(-0.02, 0.18), main = paste0("fish weight =", fish.mass," g")) #max ration
#   lines(waterTemps[,"WT"], wastes[,151, fish.mass], col = 3) #min ration
#   lines(waterTemps[,"WT"], respiration[,151, fish.mass], col = 2) #min ration
#   abline(h = 0, lty = 2)
#   if(fish.mass == 1) legend("topleft", legend = c("C-Max", "Wastes", "Respiration"), lty = 1, col = c(4, 3 , 2), bty = 'n')
# }
# dev.off()
# 
# # Create single panel plot of 'growth banana'
# png("growth.banana2.png", 8, 6, "in", res = 300)
# par(mfrow = c(1,1), las = 1)
# ylb <- "Rates (g/g/d)"
# fish.mass <- 500
# plot(waterTemps[,"WT"], consumption[,151, fish.mass], type = 'l', ylab = ylb, xlab = "Water temperature (C)",
#      ylim = c(-0.02, 0.17), main = paste0("fish weight =", fish.mass," g")) #max ration
# polygon(c(waterTemps[1:480, "WT"], rev(waterTemps[1:480, "WT"])), c(wastes[1:480, 151, fish.mass], rev(respiration[1:480, 151, fish.mass])),
#         border = NA, col = "gray70")
# lines(waterTemps[,"WT"], wastes[,151, fish.mass]) #min ration
# lines(waterTemps[,"WT"], respiration[,151, fish.mass]) #min ration
# dev.off()

