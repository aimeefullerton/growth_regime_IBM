# Plot Sensitivity Analysis
# Last updated 15 June 2019

# Set Working Directory
sens.anal.dir= "SensitivityAnalysis"

plotDir <- paste0(sens.anal.dir,"/plots.sa")
dataDir <- paste0(sens.anal.dir,"/data.out.sa")
ifelse(!dir.exists(file.path(plotDir)), dir.create(file.path(plotDir)), FALSE)
ifelse(!dir.exists(file.path(dataDir)), dir.create(file.path("data.out.SA")), FALSE)
theDirs <- dir(dataDir)

sa.parms <- read.csv(paste0(sens.anal.dir,"/","sens.anal.parms.csv"),stringsAsFactors=F)

# Define  scenarios
spp <- "steelhead"
scenarios <- c("density.hi","density.lo","movement.hi","movement.lo","food.hi","food.lo","nfish.hi","nfish.lo","constration")  
sign <- c(rep(c("hi","lo"),length(scenarios)/2))
vars <- gsub(".lo","",scenarios); vars<- gsub(".hi","",vars)

#Plot scenario boxplots together
prod.base = read.csv("data.out.1-1/production_1.csv")[,2:3]
prod.density.hi = read.csv("data.out.5-5/production_5.csv")[,2:3]
prod.density.lo = read.csv("data.out.6-6/production_6.csv")[,2:3]
prod.mvmt.hi = read.csv("data.out.7-7/production_7.csv")[,2:3]
prod.mvmt.lo = read.csv("data.out.8-8/production_8.csv")[,2:3]
prod.ration.hi = read.csv("data.out.9-9/production_9.csv")[,2:3]
prod.ration.lo = read.csv("data.out.10-10/production_10.csv")[,2:3]
prod.nfish.hi = read.csv("data.out.11-11/production_11.csv")[,2:3]
prod.nfish.lo = read.csv("data.out.12-12/production_12.csv")[,2:3]
prod.constration = read.csv("data.out.4-4/production_4.csv")[,2:3]

prod.all = cbind.data.frame(prod.base[,1],prod.density.hi[,1],prod.density.lo[,1],prod.mvmt.hi[,1],prod.mvmt.lo[,1],
                            prod.ration.hi[,1],prod.ration.lo[,1],prod.nfish.hi[,1],prod.nfish.lo[,1],prod.constration[,1])
prod.swh = cbind.data.frame(prod.base[,2],prod.density.hi[,2],prod.density.lo[,1],prod.mvmt.hi[,2],prod.mvmt.lo[,1],
                            prod.ration.hi[,2],prod.ration.lo[,2],prod.nfish.hi[,2],prod.nfish.lo[,2],prod.constration[,2])
colnames(prod.all) = colnames(prod.swh) = c("baseline",scenarios)

write.csv(prod.all,paste0(dataDir, "SA.TotalProduction.csv"))
write.csv(prod.swh,paste0(dataDir, "SA.SWHproduction.csv"))



#Load data and plot####
ibm.all.sa <- read.csv(dataDir, "SA.TotalProduction.csv",header=T,row.names = 1,stringsAsFactors=F)
ibm.swh.sa <- read.csv(dataDir, "SA.SWHproduction.csv",header=T,row.names = 1,stringsAsFactors=F)

ibm.all.sa.mn = apply(ibm.all.sa,2,median)
ibm.swh.sa.mn = apply(ibm.swh.sa,2,median)

png(paste0(plotDir,"/SensitivityAnalysis.png"),width=7,height=6,units="in",res=600)
par(las=1,mar=c(5,4,1,1),oma=c(0,0,0,0),cex=1.1)

    #two response columns
    hi.plot = c(ibm.all.sa.mn[1]-ibm.all.sa.mn[2],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[2],NA,
                ibm.all.sa.mn[1]-ibm.all.sa.mn[4],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[4],NA,
                ibm.all.sa.mn[1]-ibm.all.sa.mn[6],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[6],NA,
                ibm.all.sa.mn[1]-ibm.all.sa.mn[8],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[8])
    names(hi.plot)=NULL
    
    lo.plot = c(ibm.all.sa.mn[1]-ibm.all.sa.mn[3],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[3],NA,
                ibm.all.sa.mn[1]-ibm.all.sa.mn[5],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[5],NA,
                ibm.all.sa.mn[1]-ibm.all.sa.mn[7],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[7],NA,
                ibm.all.sa.mn[1]-ibm.all.sa.mn[9],ibm.swh.sa.mn[1]-ibm.swh.sa.mn[9])
    names(lo.plot)=NULL
    
    hi.plot = -hi.plot; lo.plot = -lo.plot
    
    barplot(hi.plot,ylim=c(-225,225),col="gray70",ylab="Difference in growth (g)",
            names.arg=rep(c("All","SWH",NA),4)[-12])
    barplot(hi.plot,add=T,density=c(rep(c(0,4,NA),4)[-12]),col=1)
    barplot(lo.plot,add=T,col="white")
    barplot(lo.plot,add=T,density=c(rep(c(0,4,NA),4)[-12]),col=1)
    box();abline(h=0)
    legend("top",legend=c("High parameter value","Low parameter value"),fill=c("gray70","white"),bty='n')
    mtext(c("Density","Movement","Food", "# Fish"),side=1,line=2,at=c(1.5,5,8.5,12),cex=1.1)
      
dev.off()



#Plot scenario boxplots together
prod.base = read.csv("data.out.1-1/production_1.csv")[,-1]
prod.divest = read.csv("data.out.13-13/production_13.csv")[,-1]
prod.enhance = read.csv("data.out.14-14/production_14.csv")[,-1]
prod.constration = read.csv("data.out.4-4/production_4.csv")[,-1]
prod = cbind.data.frame(prod.base,NA,prod.divest,NA,prod.enhance,NA,prod.constration)


png(paste0(plotDir,"/Production_Box_habitat.png"),width=7.5,height=4.5,units="in",res=300)
  par(las=1,mar=c(4,5,1,1),oma=c(0,0,0,0),cex=1.3)
  colrs = c("gray70","#F16913","#2171B5",NA)
  nms = c("Baseline","Scenario 1","Scenario 2","ConstRation")
  boxplot(prod,col=rep(colrs,4),ylim=c(-90,1090),cex=0.3,ylab="Mass accrued (g)",xaxt='n',cex.axis=0.9)
  axis(side=1,at=c(2,6,10,14),labels=nms)
  legend("topright",legend=c("Total (all habitats)","In seasonally warm habitat","In perennially cold habitat"),col=colrs[1:3],bty='n',cex=0.9,lwd=12)
  abline(h=0,lty=3)
dev.off()
