# Manuscript Figures for Armstrong et al. Growth Regimes
# This script creates individual panels that were assembled into the final figures using GIMP software
# Last updated 15 October 2020

# === SETUP #####
library(SSN)
library(RColorBrewer)
library(viridis)

source("code/growth_regime_functions_v1.0.1.R") # need 'fncColBrewPlus'
data.dir <- "data.out"
plotDir <- "plots"

# CHOOSE the iteration/scenario in the next 3 lines before proceeding:
iter <- 1
ms <- "Base" # "DivestSWH", "EnhancePCH", "DivSWH.EnhPCH"
fs <- "VariFood" #ConstFood

# Load data
scenario.nm <-  paste0(ms, ".", fs)
load(paste0("data.out/fa.",iter,".steelhead.", scenario.nm, ".RData")); fa=fa[,,,1]
load(paste0("data.out/WT.",iter,".steelhead.", scenario.nm, ".RData")); WT=WT[,,1]

ssn1 <- importSSN(paste0("data.in/network-swh.ssn")) #load ssn
ssn1 <- importPredpts(ssn1,"preds", "ssn") #import the water temperature points
  #add column to denote which reaches are accessible to fish (for blocking movement)
  ssn1@data$accessible <- 0
  ssn1@data$accessible[ssn1@data$shreve >= 16] <- 1 #Accessible based on shreve order
  ssn1@data$accessible[ssn1@data$rid %in% c(6, 11, 2, 32, 31, 39, 56, 79, 93)] <- 1 #Additional accessible cold reaches

nfish <- dim(fa)[1]; nvars <- dim(fa)[2]; ntime <- dim(fa)[3] #dimensions of arrays
d <- seq(2, ntime, 2); n<- seq(1, ntime, 2)
x.labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan")


# Define annual maximum water temperature
annual.max <- apply(t(WT), 1, quantile, probs = 1, na.rm = T)
max.tt <- which.max(WT[, 3])
reaches <- ssn1@data
obs <- ssn1@obspoints@SSNPoints[[1]]@point.data
obs$annual.max <- annual.max

# Next, aggregate temperature points for each reach
if(is.factor(obs$rid)) obs$rid <- as.numeric(levels(obs$rid)[as.numeric(obs$rid)])
annual <- aggregate(obs$annual.max, list(obs[, "rid"]), mean)
colnames(annual) <- c("seg", "annual.max")
annual$seg <- as.numeric(annual$seg)
annual <- annual[order(annual$seg), ]
reaches <- cbind(reaches, annual)

# Define reaches as seasonally warm habitat (SWH) or perennially cool habitat (PCH)
wtAug15 <- WT[454, ] #time step 454 is Julian Day 225 temps (Aug 15); or use max.tt
sw <- obs[obs$pid %in%which(wtAug15 >= 20), ] 
pc <- obs[obs$pid %in%which(wtAug15 < 20), ] 
sw.numpts <- tapply(sw$pid, sw$rid, length); sw.numpts <- cbind.data.frame("seg" = names(sw.numpts), sw.numpts)
pc.numpts <- tapply(pc$pid, pc$rid, length); pc.numpts <- cbind.data.frame("seg" = names(pc.numpts), pc.numpts)
dat <- data.frame("seg" = seq(1:151))
dat <- merge(dat, sw.numpts, by = "seg", all.x = T)
dat <- merge(dat, pc.numpts, by = "seg", all.x = T)
dat <- merge(dat, ssn1@data[, c("rid", "accessible")], by.x = "seg", by.y = "rid", all.x = T)
dat$sw.numpts[is.na(dat$sw.numpts)] <- 0
dat$pc.numpts[is.na(dat$pc.numpts)] <- 0
dat$habtype <- NA
dat$habtype[dat$sw.numpts > dat$pc.numpts] <- "swh"
dat$habtype[dat$sw.numpts < dat$pc.numpts] <- "pch"
swh <- dat$seg[dat$habtype == "swh" & dat$accessible == 1]
pch <- dat$seg[dat$habtype == "pch" & dat$accessible == 1]
sum(ssn1@data$Length[ssn1@data$rid %in%swh])
sum(ssn1@data$Length[ssn1@data$rid %in%pch])
pch; swh

segs.swh <- segs.pch <- fa[, "seg", ]
segs.pch[segs.pch == 1] <- NA #exclude lowest reach
for(i in 1:nfish){
  segs.swh[i, which(segs.swh[i, ] %in% swh)] <-1; segs.swh[i, segs.swh[i, ] > 1] <- NA
  segs.pch[i, which(segs.pch[i, ] %in% pch)] <-1; segs.pch[i, segs.pch[i, ] > 1] <- NA
}


# === FIGURE 1 #####

# Figure 1a: Map of summer max as points
ssn1@obspoints@SSNPoints[[1]]@point.data$annual.max <- obs$annual.max 
ssn1@data$col.class <- NA
breaks<- c(seq(18, 22, length.out = 13), 31)
left <- breaks[1:(length(breaks) - 1)]; rght<- breaks[2:length(breaks)]
cb <- fncColBrewPlus(n = 14, paint = F)

png(paste0(plotDir, "/SummerMaxMap_Points.png"), width = 6, height = 6, units = "in", res = 300) 
  plot(ssn1, "annual.max", breaktype = "user", brks = breaks, lwdLineCol = "addfunccol", lwdLineEx = 8, lineCol = "darkgray", xlab = "", ylab = "", color.palette = cb) 
dev.off()

# Figure 1c: Map of prioritize vs. divest, plots lines
ssn1@data$annual.max <- annual$annual.max
ssn1@data$Prioritize <- NA
ssn1@data$Prioritize[ssn1@data$annual.max >= 20] <- 1
cb <- rep("orange", nrow(ssn1@data))
ssn1@obspoints@SSNPoints[[1]]@point.data$constant <- 1

png(paste0(plotDir,"/Map_prioritize.png"), width = 6, height = 6, units = "in", res = 300) 
plot(ssn1, "xyz", axes = F, xlab = "", ylab = "", linecol = "white", lwdLineEx = 0.1) #this gives an error, just ignore - needs a base plot to work with
# plot gray stream lines (for places where no ST data exist)
for (i in 1:length(ssn1@lines)) for (j in 1:length(ssn1@lines[[i]])) lines(ssn1@lines[[i]]@Lines[[j]]@coords, col = "gray60", lwd = 10 * ssn1@data[i, "addfunccol"] + 0.4)
# plot water temperature as colored stream lines with line thickness based on cumulative drainage area
for (i in 1:length(ssn1@lines)) {
  for (j in 1:length(ssn1@lines[[i]])) {
    lines(ssn1@lines[[i]]@Lines[[j]]@coords, col = cb[ssn1@data[i, "Prioritize"]], lwd = 7 * (ssn1@data[i, "addfunccol"] + 0.4))
  }
}
dev.off()


# === FIGURE S1 #####

# (a) Temperature map
ssn1@obspoints@SSNPoints[[1]]@point.data$annual.max <- obs$annual.max 
ssn1@data$col.class <- NA
breaks<- c(seq(18, 27, length.out = 6), 31)
left <- breaks[1:(length(breaks) - 1)]; rght<- breaks[2:length(breaks)]
cb <- fncColBrewPlus(n = 7, paint = F)

png(paste0(plotDir, "/WT.Map.png"), width = 6, height = 6, units = "in", res = 300) 
  plot(ssn1, "annual.max", breaktype = "user", brks = breaks, lwdLineCol = "addfunccol", lwdLineEx = 8, lineCol = "darkgray", xlab = "", ylab = "", color.palette = cb) 
dev.off()

# (b) Habitat map
dat <- getSSNdata.frame(ssn1, "Obs")
swh.idx <- which(dat$rid %in% swh)
pch.idx <- which(dat$rid %in% pch)  
ssn1@obspoints@SSNPoints[[1]]@point.data$HabType <- NA
ssn1@obspoints@SSNPoints[[1]]@point.data$HabType[pch.idx] <- 1
ssn1@obspoints@SSNPoints[[1]]@point.data$HabType[swh.idx] <- 2
breaks <- c(1, 2, 3)
cb.hab <- c("#2171B5", "#F16913")

png(paste0(plotDir, "/Habitat.Map.png"), width = 6, height = 6, units = "in", res = 300) 
  plot(ssn1, "HabType", breaktype = "user", brks = breaks, lwdLineEx = 2, lineCol = "darkgray", xlab = "", ylab = "", color.palette = cb.hab)
dev.off()


# (c) Ration map
ration.lo <- 0.02
obs.df <- ssn1@obspoints@SSNPoints[[1]]@point.data
ration_1 <- fncRescale(log(obs.df$shreve + (1 - obs.df$ratio)), c(ration.lo, ration.lo + 0.05)) #linearly relate to stream order (shreve) and position within reach (ratio)
ration_2 <- fncRescale(obs.df$upDist, c(ration.lo + 0.05, ration.lo)) #inverse-linearly relate to distance upstream from river outlet (upDist)
obs.df$ration <- ssn1@obspoints@SSNPoints[[1]]@point.data$ration <- (ration_1 + ration_2) / 2 #mean of these two processes

png(paste0(plotDir, "/Ration.Map.png"), width = 6, height = 6, units = "in", res = 300) 
  plot(ssn1, "ration", breaktype = "even", nclasses = 6, xlab = "", ylab = "")
dev.off()


# === FIGURE 5 (and FIGURE S2) ####

# Growth/weight trajectories 

# make color matrix 
cb <- cb.wt <- fncColBrewPlus(n = 14, paint = F)
breaks <- c(0, seq(4, 26, length.out = 12), 27) 
left <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
rght <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 27)
col.mat <- fa[,"WT",]; col.mat <- col.mat[,d]
for(n in 1:length(cb)) {col.mat[col.mat >= left[n] & col.mat <= rght[n]] <- n}
col.mat2 <- col.mat
for(r in 1:nrow(col.mat2)){
  for(c in 1:ncol(col.mat2)){
    col.mat2[r, c]<- cb[col.mat[r, c]]
  }}
col.mat.wt <- col.mat2

# make color matrix based on SWH or PCH occupancy
cb.hab <- c("#9E9AC8", "#2171B5", "#F16913")
col.mat <- fa[,"seg",d]
col.mat[!fa[,"seg",d] %in% swh & !fa[,"seg", d] %in% pch] <- cb.hab[1]
col.mat[fa[,"seg",d] %in% pch] <- cb.hab[2]
col.mat[fa[,"seg",d] %in% swh] <- cb.hab[3]
col.mat.hab <- col.mat
rm(col.mat)

# Figure 5a: Growth trajectory colored by WT
png(paste0(plotDir, "/Trajectories_growth_WT.", iter, ".png"), width = 5.25, height = 4, units = "in", res = 300) 
par(las = 1, mar = c(4, 5, 2, 2), oma = c(0, 0, 0, 0))
  var <- "growth"
  yls <- c(-0.022, 0.01); abval <- 0; ylb <- "Growth (g/g*d)"; yline <- 4 
  xlb <- "Time"; xline <- 2.5
  
  # WT colors
  # plot individual colored trajectories (captures min/max)
  dat <- fa[,var,]; dat <- dat[,d]
  plot(dat[1,], type = 'n', ylim = yls, xaxt = 'n', xlab = "", ylab = "", cex = 0.1)
  axis(1, at = seq(1, 365, length.out = 13), labels = x.labels)
  mtext(ylb, 2, yline, las = 3, cex = 1.1)
  mtext(xlb, 1, xline, cex = 1.1)
  for(i in 1:nrow(dat)){points(dat[i,], col = col.mat.wt[i,], cex = 0.1)}
  leglabs <- paste(seq(4, 24, 4), "to", seq(8, 28, 4))
  mytitle<- expression("Temperature "~(degree~C))
  legend("bottomleft", legend = leglabs, title = mytitle, lwd = 2, col = cb.wt[c(4, 6, 8, 10, 12, 14)], bty = 'n', cex = 0.9)
  legend("top",legend = "Median trajectory", lwd = 2, bty = 'n')
 
  # overlay interquartile range and median
  dat<- fa[,var,]; dat<- apply(dat, 2, quantile, na.rm = T)
  dat<- t(dat)
  dat.d<- dat[d,]
  lines(dat.d[,3], lwd = 2, col = rgb(5, 5, 5, 255, NULL, 255)) #Median

dev.off()


# Figure 5b: Weight trajectory colored by habitat type
png(paste0(plotDir,"/Trajectories_weight_habitat.", iter, ".png"), width = 5.25, height = 4, units = "in", res = 300) 
par(las = 1, mar = c(4, 5, 2, 2), oma = c(0, 0, 0, 0))
  var <- "weight"
  yls <- c(0, 1000); abval = 100; ylb = "Mass (g)"; yline = 3
  xlb <- "Time"; xline <- 2.5

  # SWH or PCH colors
  # plot individual colored trajectories (captures min/max)
  dat <- fa[,var,]; dat <- dat[,d]
  plot(dat[1,], type = 'n', ylim = yls, xaxt = 'n', xlab = "", ylab = "", cex = 0.1)
  axis(1, at = seq(1, 365, length.out = 13), labels = x.labels)
  mtext(ylb, 2, yline, las = 3, cex = 1.1)
  mtext(xlb, 1, xline, cex = 1.1)
  for(i in 1:nrow(dat)){points(dat[i,], col = col.mat.hab[i,], cex = 0.1)}
    legend(1, 900, legend = c("In seasonally warm habitat", "In perennially cold habitat", "Median trajectory", "Interquartile range"),
           lwd = c(3, 3, 3, 12), col = c(cb.hab[3], cb.hab[2], 1, rgb(216, 216, 216, 160, NULL, 255)), bty = 'n')

  # overlay interquartile range and median
  dat <- fa[,var,]; dat <- apply(dat, 2, quantile, na.rm = T)
  dat <- t(dat)
  dat.d <- dat[d,]
  polygon(c(1:365, rev(1:365)), c(dat.d[,2], rev(dat.d[,4])), border = NA, col = rgb(216, 216, 216, 160, NULL, 255)) #Q1/Q3
  lines(dat.d[,3], lwd = 2, col = rgb(5, 5, 5, 255, NULL, 255)) #Median
  
dev.off()


# Figure 5c: Seasonal occupancy and trends

# Reaches over time, with aggregated fish info
fnc <- sum
byreach <- array(NA, dim = c(151, nvars, 730))
colnames(byreach) = colnames(fa)
for(i in 1:730){
  f1 <- aggregate(fa[,,i], list(fa[,"seg",i]), fnc) #all fish in the reach
  f1$seg <- f1$Group.1
  f1 <- f1[order(f1$seg), -1]
  byreach[f1$seg,,i] <- as.matrix(f1)
}

# Aggregate fish occupancy over time (sum)
  # sum of how many fish used each reach during a given season, divided by total number of fish
var <- "survive"; fnc <- sum

wi.days <- c((320 * 2):(365 * 2), 1:(81 * 2)); wi.days <- wi.days[seq(1, length(wi.days), 2)]
sp.days <- (82 * 2):(202 * 2); sp.days <- sp.days[seq(1, length(sp.days), 2)]
su.days <- (203 * 2):(265 * 2); su.days <- su.days[seq(1 ,length(su.days), 2)]
fa.days <- (266 * 2):(319 * 2); fa.days <- fa.days[seq(1 ,length(fa.days), 2)]
winter.occ <- apply(byreach[,var,wi.days], 1,fnc,na.rm=T) / length(wi.days) / nfish # / (ntime * 0.25) / nfish
spring.occ <- apply(byreach[,var,sp.days], 1,fnc,na.rm=T) / length(sp.days) / nfish
summer.occ <- apply(byreach[,var,su.days], 1,fnc,na.rm=T) / length(su.days) / nfish
fall.occ <- apply(byreach[,var,fa.days], 1,fnc,na.rm=T) / length(fa.days) / nfish
occ.dat <- cbind.data.frame(fall.occ, winter.occ, spring.occ, summer.occ)
reaches2 <- cbind(reaches, occ.dat)
reaches2$swh <- reaches2$pch <- 0
reaches2$swh[reaches2$rid %in% swh] <- 1
reaches2$pch[reaches2$rid %in% pch] <- 1
reaches2 <- reaches2[reaches2$accessible == 1,]

lm.win <- lm(reaches2$winter.occ ~ reaches2$annual.max); summary(lm.win)
lm.spr <- lm(reaches2$spring.occ ~ reaches2$annual.max); summary(lm.spr)
lm.sum <- lm(reaches2$summer.occ ~ reaches2$annual.max); summary(lm.sum)
lm.aut <- lm(reaches2$fall.occ ~ reaches2$annual.max); summary(lm.aut)


png(paste0(plotDir,"/pOccupancy_Season_Habitat",iter,".png"), width = 6, height = 4.5, units = "in", res = 300) 
par(las = 1, mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))

n <- 10; cb <- viridis(n)[c(2, 4, 7, 10)]; image(1:n, 1, as.matrix(1:n), col = cb, axes = F, ylab = "", xlab = "")
colrs <- cb

  plot(reaches2$annual.max, reaches2$winter.occ, ylim = c(0, 0.3), type = "n", pch = 19, cex = 0.8, ylab = "", xlab = "", cex.axis = 1.2)
  mtext("Reach occupancy", side = 2, line = 3.4, las = 3, cex = 1.4)
  mtext(expression("Maximum temperature in reach"~(degree ~ C)), side = 1, line = 2.9, cex = 1.4)
  points(reaches2$annual.max, reaches2$winter.occ, pch = 22, cex = 0.7, bg = colrs[2], col = "gray20")
  points(reaches2$annual.max, reaches2$spring.occ, pch = 21, cex = 0.7, bg = colrs[3], col = "gray20")
  points(reaches2$annual.max, reaches2$summer.occ, pch = 24, cex = 0.7, bg = colrs[4], col = "gray20")
  points(reaches2$annual.max, reaches2$fall.occ, pch = 23, cex = 0.7, bg = colrs[1], col = "gray20")
  lines(reaches2$annual.max, lm.win$fitted.values, col = colrs[2], lwd = 4)
  lines(reaches2$annual.max, lm.aut$fitted.values, col = colrs[1], lwd = 4)
  lines(reaches2$annual.max, lm.spr$fitted.values, col = colrs[3], lwd = 4)
  lines(reaches2$annual.max, lm.sum$fitted.values, col = colrs[4], lwd = 4, lty = 2)
  
  legend("top",legend = c("winter", "spring", "summer", "fall"), col = colrs[c(2, 3, 4, 1)], lwd = 6, bty = 'n', cex = 1.1)
dev.off()



# Figure 5d: Scenario boxplots of production by habitat

# Mass (g) accrued in each habitat
dat.swh <- fa[,"growth",] * fa[,"weight",] * segs.swh
dat.pch <- fa[,"growth",] * fa[,"weight",] * segs.pch
dat.all <- fa[,"growth",] * fa[,"weight",]

swh.production <- apply(dat.swh, 1, sum, na.rm = T)
pch.production <- apply(dat.pch, 1, sum, na.rm = T)
all.production <- apply(dat.all, 1, sum, na.rm = T)
  
summary(swh.production)
summary(pch.production)
summary(all.production)


png(paste0(plotDir,"/Production_Box_habitat",iter,".png"), width = 4, height = 6, units = "in", res = 300) 
  par(las = 1, oma = c(0, 1.5, 0, 0), cex = 1.2)
  boxplot(list(all.production, swh.production, pch.production), col = c("gray70", "#F16913", "#2171B5"),
          ylim = c(-90, 1090), cex = 0.2, ylab = "Mass accrued (g)", names = c("Total", "Warm", "Cool"))
  abline(h = 0, lty = 3)
dev.off()

production.df <- cbind.data.frame(all.production, swh.production, pch.production)
write.csv(production.df, paste0("data.out/production_",iter,".csv"))

prod.base <- read.csv("data.out/production_1.csv")[,-1]
prod.divest <- read.csv("data.out/production_2.csv")[,-1]
prod.enhance <- read.csv("data.out/production_3.csv")[,-1]
prod <- cbind.data.frame(prod.base,NA,prod.divest,NA,prod.enhance)

png(paste0(plotDir,"/Production_Box_habitat.All.png"), width = 6, height = 4.5, units = "in", res = 300)
  par(las = 1, mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0), cex = 1.3)
  colrs <- c("gray70", "#F16913", "#2171B5", NA)
  nms <- c("Baseline", "Scenario 1", "Scenario 2")
  boxplot(prod, col = rep(colrs, 3), ylim = c(-90, 1090), cex = 0.3, ylab = "Mass accrued (g)", xaxt = 'n', cex.axis = 0.9)
  axis(side = 1, at = c(2, 6, 10), labels = nms)
  legend("topright", legend = c("Total (all habitats)", "In seasonally warm habitat", "In perennially cold habitat"),
         col = colrs[1:3], bty = 'n', cex = 0.9, lwd = 12)
  abline(h = 0, lty = 3)
dev.off()

#Plot scenario boxplots together including constant ration
prod.constration <- read.csv("data.out/production_8.csv")[,-1]
prod <- cbind.data.frame(prod.base, NA, prod.divest, NA, prod.enhance, NA, prod.constration)

png(paste0(plotDir,"/Production_Box_habitat.All2.png"), width = 8, height = 4.5, units = "in", res = 300)
  par(las = 1, mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0), cex = 1.3)
  colrs <- c("gray70", "#F16913", "#2171B5", NA)
  nms <- c("Baseline", "Scenario 1", "Scenario 2", "ConstRation")
  boxplot(prod, col = rep(colrs, 4), ylim = c(-90, 1090), cex = 0.3, ylab = "Mass accrued (g)", xaxt = 'n', cex.axis = 0.9)
  axis(side = 1, at = c(2, 6, 10, 14), labels = nms)
  legend("topright", legend = c("Total (all habitats)", "In seasonally warm habitat", "In perennially cold habitat"),
         col = colrs[1:3], bty = 'n', cex = 0.9, lwd = 12)
  abline(h = 0, lty = 3)
dev.off()



# === FIGURE S3 #####
# Sensitivity Analysis
sens.anal.dir <- "SensitivityAnalysis"

plotDir <- paste0(sens.anal.dir, "/plots.sa")
dataDir <- paste0(sens.anal.dir, "/data.out.sa")
ifelse(!dir.exists(file.path(plotDir)), dir.create(file.path(plotDir)), FALSE)
ifelse(!dir.exists(file.path(dataDir)), dir.create(file.path(dataDir)), FALSE)

sa.parms <- read.csv(paste0(sens.anal.dir, "/", "sens.anal.parms.csv"), stringsAsFactors = F)

# Define  scenarios
spp <- "steelhead"
scenarios <- c("density.hi", "density.lo", "movement.hi", "movement.lo", "food.hi", "food.lo", "nfish.hi", "nfish.lo", "constration")  
sign <- c(rep(c("hi", "lo"), length(scenarios) / 2))
vars <- gsub(".lo", "", scenarios); vars <- gsub(".hi", "", vars)

# # Un-comment this section below to re-create the summarized data files used below. Original data runs need to be available.
# prod.base <- read.csv("data.out/production_1.csv")[,-1]
# prod.density.hi <- read.csv("data.out/production_10.csv")[,-1]
# prod.density.lo <- read.csv("data.out/production_11.csv")[,-1]
# prod.mvmt.hi <- read.csv("data.out/production_12.csv")[,-1]
# prod.mvmt.lo <- read.csv("data.out/production_13.csv")[,-1]
# prod.ration.hi <- read.csv("data.out/production_14.csv")[,-1]
# prod.ration.lo <- read.csv("data.out/production_15.csv")[,-1]
# prod.nfish.hi <- read.csv("data.out/production_16.csv")[,-1]
# prod.nfish.lo <- read.csv("data.out/production_17.csv")[,-1]
# prod.constration <- read.csv("data.out/production_8.csv")[,-1]
# 
# prod.all <- cbind.data.frame(prod.base[,1], prod.density.hi[,1], prod.density.lo[,1], prod.mvmt.hi[,1], prod.mvmt.lo[,1],
#                              prod.ration.hi[,1], prod.ration.lo[,1], prod.nfish.hi[,1], prod.nfish.lo[,1], prod.constration[,1])
# prod.swh <- cbind.data.frame(prod.base[,2], prod.density.hi[,2], prod.density.lo[,1], prod.mvmt.hi[,2], prod.mvmt.lo[,1],
#                              prod.ration.hi[,2], prod.ration.lo[,2], prod.nfish.hi[,2], prod.nfish.lo[,2], prod.constration[,2])
# colnames(prod.all) <- colnames(prod.swh) <- c("baseline", scenarios)
# 
# write.csv(prod.all, paste0(sens.anal.dir, "/data.out.sa/SWH.IBM.SA.TotalProduction.csv"))
# write.csv(prod.swh, paste0(sens.anal.dir, "/data.out.sa/SWH.IBM.SA.SWHproduction.csv"))


# Reload data and plot
ibm.all.sa <- read.csv(paste0(sens.anal.dir, "/data.out.sa/SWH.IBM.SA.TotalProduction.csv"), header = T, row.names = 1, stringsAsFactors = F)
ibm.swh.sa <- read.csv(paste0(sens.anal.dir, "/data.out.sa/SWH.IBM.SA.SWHproduction.csv"), header = T, row.names = 1, stringsAsFactors = F)

ibm.all.sa.mn <- apply(ibm.all.sa, 2, median)
ibm.swh.sa.mn <- apply(ibm.swh.sa, 2, median)

png(paste0(plotDir,"/SensitivityAnalysis.png"), width = 7, height = 6, units = "in", res = 600)
par(las = 1, mar = c(5, 4, 1, 1), oma = c(0, 0, 0, 0), cex = 1.1)

# two response columns
hi.plot <- c(ibm.all.sa.mn[1] - ibm.all.sa.mn[2], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[2], NA,
             ibm.all.sa.mn[1] - ibm.all.sa.mn[4], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[4], NA,
             ibm.all.sa.mn[1] - ibm.all.sa.mn[6], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[6], NA,
             ibm.all.sa.mn[1] - ibm.all.sa.mn[8], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[8])
names(hi.plot) <- NULL

lo.plot <- c(ibm.all.sa.mn[1] - ibm.all.sa.mn[3], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[3], NA,
             ibm.all.sa.mn[1] - ibm.all.sa.mn[5], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[5], NA,
             ibm.all.sa.mn[1] - ibm.all.sa.mn[7], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[7], NA,
             ibm.all.sa.mn[1] - ibm.all.sa.mn[9], ibm.swh.sa.mn[1] - ibm.swh.sa.mn[9])
names(lo.plot) <- NULL

hi.plot <- -hi.plot; lo.plot <- -lo.plot

barplot(hi.plot, ylim = c(-225, 225), col = "gray70", ylab = "Difference in growth (g)",
        names.arg = rep(c("All", "SWH", NA), 4)[-12])
barplot(hi.plot, add = T, density = c(rep(c(0, 4, NA), 4)[-12]), col = 1)
barplot(lo.plot, add = T, col = "white")
barplot(lo.plot, add = T, density = c(rep(c(0, 4, NA), 4)[-12]), col = 1)
box(); abline(h = 0)
legend("top", legend = c("High parameter value", "Low parameter value"), fill = c("gray70", "white"), bty = 'n')
mtext(c("Density", "Movement", "Food", "# Fish"), side = 1, line = 2, at = c(1.5, 5, 8.5, 12), cex = 1.1)

dev.off()

