## Pre-calculate O. mykiss growth at each temperature that can be Looked Up instead of having to run bioenergetics on the fly each time

# Need to go run 'setup' and 'define parameters' in 'growth_regime_functions_v1.0.1.R' first


Constants = fncReadConstants.steelhead_ration()

waterTemps=cbind("pid"=seq(1,length(seq(0.05,25,0.05))),"WT"=seq(0.05,25,0.05))
rations = seq(0.02,0.17,0.001)
weights = seq(1,1500,1)
wt.growth=array(NA,dim=c(nrow(waterTemps),length(rations),length(weights))) #dim=c(500 WTs, 151 rations, 700 weights)
respiration = consumption = wastes = wt.growth
for(w in 1:dim(wt.growth)[3]){
  ww = weights[w]
  for(r in 1:dim(wt.growth)[2]){
    rr = rations[r]
    Input<- fncGetBioEParms("steelhead", pred.en.dens, prey.en.dens, oxy, pff,
                            waterTemps, 
                            startweights=rep(ww,nrow(waterTemps)), 
                            pvals=rep(1,nrow(waterTemps)),
                            ration=rep(rr,nrow(waterTemps)))
    
    
    Results = BioE(Input, Constants) #run BioE code given input and constants
    wt.growth[,r,w] = c(Results$Gg_WinBioE)
    consumption[,r,w] = c(Results$Consumption)
    wastes[,r,w] = c(Results$Egestion + Results$Excretion + Results$S.resp)
    respiration[,r,w] = c(Results$Respiration)
    
    #in Joules instead:
      #wt.growth[,r,w] = c(Results$Growth_j)
      #consumption[,r,w] = c(Results$Consumption_j)
      #wastes[,r,w] = c(Results$Egestion_j + Results$Excretion_j + Results$Sj.resp)
      #respiration[,r,w] = c(Results$Respiration_j)
  }
}
save("wt.growth",file="data.in/wt.growth.array.RData")
save("respiration",file="data.in/respiration.array.RData")
save("wastes",file="data.in/wastes.array.RData")
save("consumption",file="data.in/consumption.array.RData")


# Load back in and plot
load("data.in/wt.growth.array.RData")
load("data.in/consumption.array.RData")
load("data.in/wastes.array.RData")
load("data.in/respiration.array.RData")

fish.mass = 1 #1-g fish
png("growth.curves.png",8,6,"in",res=300)
par(mfrow=c(2,3))
for(fish.mass in c(1,10,100,250,500,1100)){
  plot(waterTemps[,"WT"],wt.growth[,151,fish.mass],type='l',ylab="Growth (g/g/d)",xlab="Water temperature (C)",col=4,ylim=c(-0.04,0.04),main=paste0("fish weight =",fish.mass," g")) #max ration
  lines(waterTemps[,"WT"],wt.growth[,1,fish.mass],col=2) #min ration
  for(i in seq(3,150,length.out = 10)){ lines(waterTemps[,"WT"],wt.growth[,i,fish.mass],col="gray30",lty=2)}
  abline(h=0,lty=2)
}
dev.off()


png("growth.banana.png",8,6,"in",res=300)
par(mfrow=c(2,3))
for(fish.mass in c(1,10,100,250,500,1100)){
  plot(waterTemps[,"WT"],consumption[,151,fish.mass],type='l',ylab="Rates (g/g/d)",xlab="Water temperature (C)",col=4,ylim=c(-0.02,0.18),main=paste0("fish weight =",fish.mass," g")) #max ration
  lines(waterTemps[,"WT"],wastes[,151,fish.mass],col=3) #min ration
  lines(waterTemps[,"WT"],respiration[,151,fish.mass],col=2) #min ration
  abline(h=0,lty=2)
  if(fish.mass ==1) legend("topleft",legend=c("C-Max","Wastes","Respiration"),lty=1,col=c(4,3,2),bty='n')
}
dev.off()

png("growth.banana2.png",8,6,"in",res=300)
fish.mass = 500
plot(waterTemps[,"WT"],consumption[,151,fish.mass],type='l',ylab="Rates (g/g/d)",xlab="Water temperature (C)",ylim=c(-0.02,0.17),main=paste0("fish weight =",fish.mass," g")) #max ration
polygon(c(waterTemps[1:480,"WT"],rev(waterTemps[1:480,"WT"])),c(wastes[1:480,151,fish.mass],rev(respiration[1:480,151,fish.mass])),border=NA,col="gray70")
lines(waterTemps[,"WT"],wastes[,151,fish.mass]) #min ration
lines(waterTemps[,"WT"],respiration[,151,fish.mass]) #min ration
dev.off()

png("growth.banana_100g.png",8,6,"in",res=300)
par(las=1)
fish.mass = 100
ylm = c(0,250) #100-g
plot(waterTemps[,"WT"],consumption[,151,fish.mass],type='n',ylim=ylm,ylab="Rates (J/g/d)",xlab = expression("Temperature "~degree(C)),main=paste0("fish weight =",fish.mass," g")) #max ration
ii = 447 #100g
polygon(c(waterTemps[1:ii,"WT"],rev(waterTemps[1:ii,"WT"])),c(wastes[1:ii,151,fish.mass],rev(respiration[1:ii,151,fish.mass])),border=NA,col="gray70")
lines(waterTemps[,"WT"],wastes[,151,fish.mass]) #min ration
lines(waterTemps[,"WT"],respiration[,151,fish.mass]) #min ration
dev.off()

