#=================================================================================================
# Functions called by Growth Regime IBM, VERSION 1.0.1
#
#     Armstrong, J.B., A.H. Fullerton, C.E. Jordan, J.L. Ebersole, J.R. Bellmore, I. Arismendi, 
#         B. Penaluna, and G.H. Reeves. The significance of warm habitat to the growth regime of 
#         coldwater fishes.
#
# Last functional updates 13 June 2019; editorial edits 29 July 2020
#=================================================================================================

# === BASIC FUNCTIONS ============================================================================

# Get the four coordinates for a single segment
getSegmentCoords <- function(seg = seg, ssn = ssn1) {
  idx <- which(ssn@data$rid == seg)
  return(ssn@lines[[idx]]@Lines[[1]]@coords)
}

# Distance to network base for a segment (including the length of this segment)
distance2base4seg <- function(seg = seg, ssn = ssn1) {
  return(ssn@data$upDist[ssn@data$rid == seg])
}

# Distance to network base for a fish
distance2base4fish <- function(fishID = 1, ssn = ssn1) {
  return(fish$length2segBase[fish$pid == fishID] - ssn@data[fish$seg[fish$pid == fishID], length.field] + distance2base4seg(fish$seg[fish$pid == fishID]))
}

# Convert segment ratio (proprtional length) to length downstream
segRatio2length <- function(seg = seg, segRatio = 0.5, ssn = ssn1) {
  return (ssn@data[ssn@data[, "rid"] == seg, length.field] * segRatio)
}

# Convert length downstream to segment ratio
length2segRatio <- function(seg = seg, lengthDownstream = 0.5, ssn = ssn1) {
  return (lengthDownstream / ssn@data[ssn@data[, "rid"] == seg, length.field])
}

# Returns the downstream segment
getSegDownstream <- function(seg = seg, ssn = ssn1) {
  if (seg == min(ssn@data$rid)) return(0)
  # loop over each segment to see if the first coordinate set matches this one's second coordinate set
  for (x in min(ssn@data$rid):max(ssn@data$rid)) {
    coords.x <- getSegmentCoords(seg = x)
    coords.seg <- getSegmentCoords(seg = seg)
    nrows.seg <- nrow(coords.seg)
    if(coords.x[1, 1] == coords.seg[nrows.seg, 1] & coords.x[1, 2] == coords.seg[nrows.seg, 2]) return(x)
  }  
}

# Returns the upstream segment(s)
getSegUpstream <- function(seg = seg, ssn = ssn1) {
  mySegs <- c()
  for (x in min(ssn@data$rid):max(ssn@data$rid)) {
    coords.x <- getSegmentCoords(seg = x)
    coords.seg <- getSegmentCoords(seg = seg)
    nrows.x <- nrow(coords.x)
    if  (coords.x[nrows.x, 1] == coords.seg[1, 1] & coords.x[nrows.x, 2] == coords.seg[1, 2]) mySegs <- c(mySegs, x)
  }
  return(mySegs)
}

# Get the x, y coordinate for a particular location along a straight-line segment in generated stream networks
#    For line coordinates, the higher values (upstream) are the first row in the coordinate set, as if going
#    from x1 to x2 is moving downstream, but foor other things, like ratio, it's more like moving upstream
#    Called automatically when 'netnm' begins with "network-".
getXY <- function(seg = seg, ratio = 0.5, ssn = ssn1) {
  # Function is (x1-x2) * proportion + x2
  rowtoget <- nrow(getSegmentCoords(seg))
  idx <- which(ssn@data$rid == seg)
  xloc <- (ssn1@lines[[idx]]@Lines[[1]]@coords[1, 1] - ssn1@lines[[idx]]@Lines[[1]]@coords[rowtoget, 1]) * ratio + ssn1@lines[[idx]]@Lines[[1]]@coords[rowtoget, 1]
  yloc <- (ssn1@lines[[idx]]@Lines[[1]]@coords[1, 2] - ssn1@lines[[idx]]@Lines[[1]]@coords[rowtoget, 2]) * ratio + ssn1@lines[[idx]]@Lines[[1]]@coords[rowtoget, 2]
  return(cbind(xloc, yloc))
}

# This is the curvy line version for segments in stream networks generated in GIS
#   It operates by moving fish to closest existing node within a reach (i.e., coordinates of the points along
#   a line that was created in GIS); This version is called if 'netnm' does not begin with "network-".
getXY.arc <- function(seg = seg, ratio = 0.5, ssn = ssn1) {
  idx <- which(ssn@data$rid == seg)
  ratio <- 1 - ratio #may need to check for your network if it was generated in GIS
  ll <- LineLength(ssn1@lines[[idx]]@Lines[[1]]@coords, sum=F)
  lL <- LineLength(ssn1@lines[[idx]]@Lines[[1]]@coords, sum=T)
  l.idx <- 0; counter <- 1
  while(l.idx < ratio * lL){l.idx <- l.idx + ll[counter]; counter <- counter+1}
  foo <- c(lL - sum(ll[1:(counter-1)]), lL - sum(ll[1:(counter-2)])); l.new <- which.min(foo)
  ifelse(l.new == 1, new.idx <- counter, new.idx <- counter-2)
  xloc <- ssn1@lines[[idx]]@Lines[[1]]@coords[new.idx, 1]
  yloc <- ssn1@lines[[idx]]@Lines[[1]]@coords[new.idx, 2]
  return(cbind(xloc, yloc))
}

# Get stream order for a seg
getSO <- function(seg = seg, ssn = ssn1) {
  #return(ssn@data[, so.field][ssn@data$rid == seg])
  return(obs.df[, so.field][obs.df$rid == seg][1])
}

# Trace upstream
fncTraceUp <- function(seg = seg, ssn = ssn1, plotit = T, col = "turquoise", no2get=50000){
  
  ups <- getSegUpstream(seg = seg, ssn = ssn)
  set <- c(seg, ups)
  new <- ups
  
    while(length(ups)>0 & length(set)<no2get){
    for(u in 1:length(ups)){
        new <- c(new, getSegUpstream(seg=ups[u], ssn = ssn))
      }
      new <- new[!new %in% set]
      ups <- new
      set <- c(set, new)
    }

  if(plotit == T){for(x in set){highlightSegment(x, col = col)}}
  
  return(set)
}    
# To plot:
   #pp <- plot(ssn1, lwdLineCol = "addfunccol", lwdLineEx=15, linecol = "gray50", xlab = "", ylab = "", cex = 0)
   #op <- par(usr=pp$usr)
   #fncTraceUp(seg=5, ssn = ssn1, plotit = T, col = "purple")

# Trace downstream
fncTraceDn <- function(seg = seg, ssn = ssn1, plotit = T, col = "yellow", no2get = 50000){
  
  dn <- getSegDownstream(seg = seg, ssn = ssn)
  set <- c(seg, dn)
  new <- dn
  
  while(length(dn) > 0 & dn > 1 & length(set) < no2get){
    for(u in 1:length(dn)){
      new <- c(new, getSegDownstream(seg = dn[u], ssn = ssn))
    }
    new <- new[!new %in% set]
    dn <- new
    set <- c(set, new)
  }
  
  if(plotit == T){for(x in set){highlightSegment(x, col = col)}}
  
  return(set)
}

# Highlight a particular fish
highlightFish <- function(fishID, col = "turquoise", cex = 2, pch = 20) {
  points(fish$xloc[fish$pid == fishID], fish$yloc[fish$pid == fishID], col = col, cex = cex, pch = pch)
}

# Highlight chosen segment
highlightSegment <- function(seg = seg, col = "turquoise", lwd = 2) {
  lines(getSegmentCoords(seg = seg)[, 1], getSegmentCoords(seg = seg)[, 2], col = col, lwd = lwd)
}

# Label chosen segment
labelSegment <- function(seg = seg, col = 1, lwd = 1) {
  lines(getSegmentCoords(seg = seg)[, 1], getSegmentCoords(seg = seg)[, 2], col = col, lwd = lwd)
  #label seg:
  text(mean(getSegmentCoords(seg = seg)[, 1]), mean(getSegmentCoords(seg = seg)[, 2]), seg, cex = 0.8)
  #label shreve:
  #text(mean(getSegmentCoords(seg = seg)[, 1]), mean(getSegmentCoords(seg = seg)[, 2]), ssn1@data$shreve[ssn1@data$rid == seg], cex = 0.8)
}
# To plot:
  #pp <- plot(ssn1, lwdLineCol = "addfunccol", lwdLineEx = 15, linecol = "gray50", xlab = "", ylab = "", cex = 0)
  #op <- par(usr = pp$usr)
  #for(i in 1:nrow(ssn1@data)){labelSegment(i)}


# Get closest pred points for all branches of a junction at the upper end of a given segment where water temperature values exist
  # only needed for stream networks created in GIS, not generated networks
  # for junctions at confluences, function will return the top point on the segment itself and lowest points on each of the two upstream segments for which predictions exist
  # for a junction of two arcs (no branches), function will return the the top point on the segment and the lowest point on the upstream segment for which predictions exist
  # if a segment is a headwater segment, function will return the top point on the segment itself for which predictions exist
  # for segments without prediction values: 
    # the function will look to see if there are predictions on any segments upstream, and if so, will return the downstream-most prediction point
    # if there are no segments upstream that have predictions, the function will return either the top point from the segment if a prediction point exists, 
    # or will grab the top point of the nearest downstream prediction
fncGetJunctionPreds <- function(dat.df, seg = seg, ssn = ssn1){
  if(seg == 0) seg <- 1
  mySegs <- c(seg, getSegUpstream(seg = seg, ssn = ssn))
  out <- matrix(NA, length(mySegs), 2)
  nn <- 1
  for(r in mySegs){
    if(r %in% missed){ #if r is one of the reaches with no prediction points...
      ups <- upsegs[[r]][1:min(length(upsegs[[r]]), 10)]; ups <- ups[!ups %in% missed]; if(length(ups) == 0) ups <- NA
      dns <- dnsegs[[r]][1:min(length(dnsegs[[r]]), 6)];dns <- dns[!dns %in% missed]; if(length(dns) == 0) dns <- 0
      if(is.na(ups)[1]) { #if no more upstream segs with temperatures, use value from nearest downstream seg
        rr <- dns[which(!dns %in% missed)][1]
        pt <- dat.df[dat.df$rid == rr, ]; pt <- pt[which.max(pt$ratio), "pid"] #if reach has more than one obs, use most upstream one
      } else{ #if any of upstream segs have water temperatures, return the lowest point for each seg
        pt <- NULL
        if(!any(ups %in% mySegs)){ #if none of these matches the next seg in mySegs
        for(rr in ups){
          tmp.df <- dat.df[dat.df$rid == rr, ]
          p1 <- tmp.df[which.min(tmp.df$ratio), "pid"]; pt <- c(pt, p1)
        }
        } else {
          rr <- ups[which(ups %in% mySegs)][1]
          tmp.df <- dat.df[dat.df$rid == rr, ]
          ifelse(rr == mySegs[1], pt <- tmp.df[which.max(tmp.df$ratio), "pid"], pt <- tmp.df[which.min(tmp.df$ratio), "pid"])
        }
      }
    } else{
      # for the trunk, we want the upstream-most point, for branches we want the downstream-most point
      tmp.df <- dat.df[dat.df$rid == r, ]
      ifelse(r == mySegs[1], pt <- tmp.df[which.max(tmp.df$ratio), "pid"], pt <- tmp.df[which.min(tmp.df$ratio), "pid"])
    }
    out[nn, ] <- cbind(r, pt)
    nn <- nn+1
  }
  colnames(out) <- c("r", "pt")
  return(out)
}


# Get list of downstream segs for each seg
fncDnSegs <- function(ssn = ssn1, path){
  seglist <- sort(unique(ssn@data$rid))
  seglist <- seglist[-c(1)] #remove the base one because can't index a list on zero
  
  dnsegs <- list()
  for(seg in seglist){
    sg <- fncTraceDn(seg, plotit = F)
    cat(seg, "\n")
    dnsegs[[seg]] <- sg[-1]
  }
  
  save("dnsegs", file=paste0(path, "/dnsegs.RData"))
  return(dnsegs)
}

# Get list of Upstream segs for each seg
fncUpSegs <- function(ssn = ssn1, path){
  seglist <- sort(unique(ssn@data$rid))
  seglist <- seglist[-c(1)] #remove the base one because can't index a list on zero
  
  upsegs <- list()
  for(seg in seglist){
    sg <- fncTraceUp(seg, plotit = F)
    cat(seg, "\n")
    upsegs[[seg]] <- sg[-1]
  }
  
  save("upsegs", file=paste0(path, "/upsegs.RData"))
  return(upsegs)
}

# Get list of segs at junctions for each seg
fncJctLst <- function(ssn = ssn1, path){
  seglist <- sort(unique(ssn@data$rid))
  if(min(seglist) == 0) seglist <- seglist[-c(1)] #remove the base one because can't index a list on zero
  
  jct.list <- list()
  for(seg in seglist){
    if(seg>0) jct.list[[seg]] <- data.frame(fncGetJunctionPreds(obs.df, seg, ssn = ssn1))
  }
  
  save("jct.list", file=paste0(path, "/jct.list.RData"))
  return(jct.list)
}

# Rescale a vector to a specific
fncRescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
    (x - from[1]) / diff(from) * diff(to) + to[1]
}

#Get Julian month from Julian day
fncJulianMonth <- function(x){
  if(x >= 0){
    if(x <= 31) return(1)
    else if(x >= 32 & x < 60) return (2)
    else if(x >= 60 & x < 91) return (3)
    else if(x >= 91 & x < 121) return (4)
    else if(x >= 121 & x < 152) return (5)
    else if(x >= 152 & x < 182) return (6)
    else if(x >= 182 & x < 213) return (7)
    else if(x >= 213 & x < 244) return (8)
    else if(x >= 244 & x < 274) return (9)
    else if(x >= 274 & x < 305) return (10)
    else if(x >= 305 & x < 335) return (11)
    else if(x >= 335 & x < 366) return (12)
  } else if(x < 0){
    if(-x <= 31) return(12)
    else if(-x >= 32 & -x < 60) return (11)
    else if(-x >= 60 & -x < 91) return (10)
    else if(-x >= 91 & -x < 121) return (9)
    else if(-x >= 121 & -x < 152) return (8)
    else if(-x >= 152 & -x < 182) return (7)
    else if(-x >= 182 & -x < 213) return (6)
    else if(-x >= 213 & -x < 244) return (5)
    else if(-x >= 244 & -x < 274) return (4)
    else if(-x >= 274 & -x < 305) return (3)
    else if(-x >= 305 & -x < 335) return (2)
    else if(-x >= 335 & -x < 366) return (1)
    
  }
}

# Get iteration
fncGetRun <- function(){
  run =  as.character(Sys.time())
  run =  gsub("-", ".", run)
  run =  gsub(":", ".", run)
  run =  gsub(" ", ".", run)
  return(run)
}
#****************************************************************************************

# === CREATE STREAM NETWORKS =====================================================================


# These functions were provided by Nick Som, US Fish & Wildlife Service
# and Pascal Monestiez, INRA France, and adapted for our use
# in generating SSN stream networks with more flexible topology/shape
# changes denoted with "AHF"
# see also: Som, N. A., P. Monestiez, J. M. Ver Hoef, D. L. Zimmerman, and E. E. Peterson. 2014. 
# Spatial sampling on streams: principles for inference on aquatic networks. Environmetrics 25:306-323.


# Simulate tree
simul_tree7toSSNBIDforce2trib <- function(n = 200, ordermax = 3, alph = pi / 3, beta = pi / 8, lmin = 0.6, lmax = 1.1, sampleprob = 1, ps = 1){
  # version 7.0 : angle control added (relatively to the others already existing tributaries)
  # proximity constraint between segments to avoid crossing or edge contact, 
  # branching order limited to 3 (or 4) i.e. two (or 3) maximum tribs
  # introduction of a weight for the sampling of nodes (no longer equiprobable) 
  
  nb <- n ## number of stream segments to simulate
  if( ordermax != 3 | ordermax != 4 ) ordermax = 3 
  # ordermax <- 3 # Two maximum tributaries on the same edge extremity.
  # parameter ordermax has only 2 possible values 3 an 4 (if 2, the tree resumes to one broken line)
  # we do not consider more than three tributaries on the same edge extremity.  
  
  # sampleprob <- 1  # equiprobability to sample any of the previous edge extremities to branch a new one
  
  # sampleprob <- 2  # the sources (extremities) have a prob ps times greater to be individualy selected
  # than inner existing nodes 
  
  # sampleprob <- 3  # the probablity to sample a source is proportional to the distance from the origin
  # and is ps times greater than an inner node (giving elongated networks) 
  
  # sampleprob <- 4  # the probablity to sample a source is inversely-proportional to the distance from the origin 
  # and is ps times greater than an inner node (giving very compact networks if ps < 1)               
  # it can become very long to simulate because the larger probability corresponds to impossible
  # branchings (for example for ps = 0.2)
  
  pex <- 5
  # sampleprob <- 5  # as 4 at power pex if close enough of the extremity (farthest point) 
  
  tree <- matrix(NA, nb, 10)
  tree[, 1] <- 1:nb
  colnames(tree) <- c("indseg", "x", "y", "ancestor", "nbtrib", "length", "theta", "tribind1", "tribind2", "dist2orig")
  tree[1, 1] <- 1   # indseg
  tree[1, 2] <- 0   # x of extremity
  tree[1, 3] <- 1   # y of extremity
  tree[1, 4] <- 0   # ind of ancestor
  tree[1, 5] <- 0   # nb of tributatries
  tree[1, 6] <- 1   # length
  tree[1, 7] <- 0   # theta relatively to the y direction
  tree[1, 8] <- NA
  tree[1, 9] <- NA
  tree[1, 10] <- 1  # total length 
  
  i <- 2
  tree[i, 4] <- 1
  tree[i, 5] <- 0
  tree[tree[i, 4], 5] <- 1
  length <- runif(1, lmin, lmax)
  theta <- runif(1, -alph, +alph)                   
  tree[i, 2] <- tree[tree[i, 4], 2] + length * sin(theta)
  tree[i, 3] <- tree[tree[i, 4], 3] + length * cos(theta)           
  tree[i, 6] <- length
  tree[i, 7] <- theta
  tree[i, 8] <- NA
  tree[i, 9] <- NA
  tree[i, 10] <- tree[tree[i, 4], 10] + length
  tree[tree[i, 4], 8] <- i             
  
  for (i in 3:nb)
  { proba <- rep(1, (i-1)) # default choice
  if(sampleprob == 2) {proba[tree[1:(i-1), 5] == 0] <- ps} ## i is next, so prob of samp for previous nodes
  if(sampleprob == 3) {maxdist <- max(tree[1:(i-1), 10][tree[1:i-1, 5] == 0])
  proba[tree[1:(i-1), 5] == 0] <- ps * 
    (tree[1:(i-1), 10][tree[1:i-1, 5] == 0]/maxdist)} 
  if(sampleprob == 4) {maxdist <- max(tree[1:(i-1), 10][tree[1:i-1, 5] == 0])
  proba[tree[1:(i-1), 5] == 0] <- ps * 
    ((maxdist - tree[1:(i-1), 10][tree[1:i-1, 5] == 0])/maxdist)}    
  if(sampleprob == 5) {maxdist <- max(tree[1:(i-1), 10][tree[1:i-1, 5] == 0])
  if(maxdist < 10) dista <- tree[1:(i-1), 10]
  if(maxdist >= 10) dista <- tree[1:(i-1), 10] -(maxdist-10)
  dista[dista < 0] <- 0  
  proba[tree[1:(i-1), 5] == 0] <- ps * 
    (dista[tree[1:i-1, 5] == 0]/10)**pex}                  
  
  if(i %% 2 != 0){ #i not ending in 0
    last.anc <- tree[(i-1), 4]	
    proba <- rep(0, (i-1))
    proba[last.anc] <- 1
  }
  
  tree[i, 4] <- sample(seq(1, (i-1), 1), 1, prob = proba)  
  while(tree[tree[i, 4], 5] >= (ordermax-1)) {tree[i, 4] <- sample(seq(1, (i-1), 1), 1, prob = proba)} 
  
  tree[i, 5] <- 0
  theta <- runif(1, -alph, +alph)
  length <- runif(1, lmin, lmax)
  tree[i, 2] <- tree[tree[i, 4], 2] + length * sin(theta)
  tree[i, 3] <- tree[tree[i, 4], 3] + length * cos(theta)           
  tree[i, 6] <- length
  tree[i, 7] <- theta
  tree[i, 10] <- tree[tree[i, 4], 10] + length
  
  # computation of minimum distance with all other segments excepted
  # it ancestor and other tributaries on the same node if there exists   
  dsi <- Inf
  for (k in 2:(i-1))
  {if ( k != tree[i, 4] & tree[k, 4] != tree[i, 4] ) 
  {dsi <- min(dsi, ds(tree[i, 2], tree[tree[i, 4], 2], 
                     tree[i, 3], tree[tree[i, 4], 3], 
                     tree[k, 2], tree[tree[k, 4], 2], 
                     tree[k, 3], tree[tree[k, 4], 3]))
  } 
  }
  crit1 <- dsi 
  
  # computation of the angle between the new edge and the tributaries
  # on the same node if there exists
  
  crit2 <- pi
  nbotrib <- tree[tree[i, 4], 5] 
  if (nbotrib  ==  1)  crit2 <- abs(tree[i, 7] - tree[tree[tree[i, 4], 8], 7]) 
  if (nbotrib  ==  2)  {crit2 <- min(abs(tree[i, 7] - tree[tree[tree[i, 4], 8], 7]), 
                                   abs(tree[i, 7] - tree[tree[tree[i, 4], 9], 7]))}  
  
  nbloops <- 1   
  while(crit1 <= 0.3 | crit2 <= beta)
  {tree[i, 4] <- sample(seq(1, (i-1), 1), 1, prob = proba)	
  while(tree[tree[i, 4], 5]>=2) {tree[i, 4] <- sample(seq(1, (i-1), 1), 1, prob = proba)} 
  
  tree[i, 5] <- 0
  theta <- runif(1, -alph, +alph)
  length <- runif(1, lmin, lmax)
  tree[i, 2] <- tree[tree[i, 4], 2] + length * sin(theta)
  tree[i, 3] <- tree[tree[i, 4], 3] + length * cos(theta)           
  tree[i, 6] <- length
  tree[i, 7] <- theta
  tree[i, 10] <- tree[tree[i, 4], 10] + length
  
  dsi <- Inf
  for (k in 2:(i-1))
  {if ( k != tree[i, 4] & tree[k, 4] != tree[i, 4] ) 
  {dsi <- min(dsi, ds(tree[i, 2], tree[tree[i, 4], 2], 
                     tree[i, 3], tree[tree[i, 4], 3], 
                     tree[k, 2], tree[tree[k, 4], 2], 
                     tree[k, 3], tree[tree[k, 4], 3]))
  } 
  }
  crit1 <- dsi
  
  crit2 <- pi
  nbotrib <- tree[tree[i, 4], 5] 
  if (nbotrib  ==  1)  crit2 <- abs(tree[i, 7] - tree[tree[tree[i, 4], 8], 7]) 
  if (nbotrib  ==  2)  {crit2 <- min(abs(tree[i, 7] - tree[tree[tree[i, 4], 8], 7]), 
                                   abs(tree[i, 7] - tree[tree[tree[i, 4], 9], 7]))}    
  
  nbloops <- nbloops + 1
  if(nbloops > 1000) { stop("no more possibilities to add a new branch without excessive time")}
  }  
  tree[tree[i, 4], 5] <- tree[tree[i, 4], 5] + 1
  if(tree[tree[i, 4], 5] == 1){tree[tree[i, 4], 8] <- i}
  if(tree[tree[i, 4], 5] == 2){tree[tree[i, 4], 9] <- i}
  
  }
  
  # the matrix to use is named tree (10 columns, nb rows) 
  
  tree[, 1] <- 1:nrow(tree)
  out.mat <- tree
  
  out.mat <- out.mat[order(out.mat[, 4], out.mat[, 1]), ]
  out.mat2 <- out.mat
  
  while(!identical(as.integer(out.mat2[, 1]), 1:n)){
    
    for(i in 2:nrow(out.mat)){
      sb <- i
      iS <- out.mat[i, 1]
      if(sb == iS){next}
      else{
        which.sb <- which(out.mat[, 1] == sb)
        which.ac.iS <- which(out.mat[, 4] == iS)
        which.ac.sb <- which(out.mat[, 4] == sb)
        which.t1.iS <- which(out.mat[, 8] == iS)
        which.t2.iS <- which(out.mat[, 9] == iS)
        which.t1.sb <- which(out.mat[, 8] == sb)
        which.t2.sb <- which(out.mat[, 9] == sb)
        out.mat[i, 1] <- sb
        out.mat[which.sb, 1] <- iS
        out.mat[which.ac.iS, 4] <- sb
        out.mat[which.ac.sb, 4] <- iS
        
        out.mat[which.t1.iS, 8] <- sb;out.mat[which.t1.sb, 8] <- iS
        out.mat[which.t2.iS, 9] <- sb;out.mat[which.t2.sb, 9] <- iS
      }
    }
    out.mat <- out.mat[order(out.mat[, 4], out.mat[, 1]), ]
    out.mat2 <- out.mat
  }
  
  ### make binID table
  RID <- as.numeric(rep(NA, nrow(out.mat)))
  binaryID <- rep(NA, nrow(out.mat))
  
  RID[1] <- 1
  binaryID[1] <- "1"
  for(i in 1:nrow(out.mat)){
    if(sum(is.na(out.mat[i, 8:9])) == 2){next}
    if(sum(is.na(out.mat[i, 8:9])) == 1){
      Trib1 <- out.mat[i, 8]
      RID[Trib1] <- Trib1
      binaryID[Trib1] <- paste(binaryID[i], "1", sep="")
    }
    if(sum(is.na(out.mat[i, 8:9])) == 0){
      Trib1 <- out.mat[i, 8]
      RID[Trib1] <- Trib1
      binaryID[Trib1] <- paste(binaryID[i], "0", sep="")
      Trib2 <- out.mat[i, 9]
      RID[Trib2] <- Trib2
      binaryID[Trib2] <- paste(binaryID[i], "1", sep="")
    }
  }
  
  bin.table <- data.frame(rid = RID, binaryID = binaryID)
  
  gr <- c(0, 1)
  locs.x <- c(0, out.mat[1, 2])
  locs.y <- c(0, out.mat[1, 3])
  
  for(i in 2:nrow(out.mat)){
    gr <- c(gr, out.mat[i, 4], out.mat[i, 1]) 
    locs.x <- c(locs.x, out.mat[i, 2])
    locs.y <- c(locs.y, out.mat[i, 3])
  }
  
  gr <- gr + 1 #AHF: to avoid value of zero, which is not a valid index name
  gr <- graph(gr, directed=T)
  locs <- matrix(cbind(locs.x, locs.y), ncol = 2)
  DendNet <- vector(mode="list")
  DendNet$graph <- gr
  DendNet$locations <- locs
  DendNet$initialPoint <- 1 #0
  return(list(DendNet = DendNet, bin.table = bin.table))
  
}

# Tree function 1 (called by simul_tree...)
ds <- function(x1, x2, y1, y2, x3, x4, y3, y4){
  xp <- seq(x1, x2, length=10)
  yp <- seq(y1, y2, length=10)
  xq <- seq(x3, x4, length=10)
  yq <- seq(y3, y4, length=10)
  dps <- distpq(xp, yp, xq, yq)
  mds <- min(dps)
  return(mds)
}

# Tree function 2 (called by simul_tree...)
distpq <- function(xp, yp, xq, yq){ 
  np <- length(xp)
  nq <- length(xq)
  dist <- matrix(0, np, nq)  
  dist <- dist + outer(xp, xq, "-")^2
  dist <- dist + outer(yp, yq, "-")^2
  dist <- sqrt(dist)
  return(dist)
}

# Import simulated tree into SSN
createSSNldotsBID <- function (n, obsDesign, predDesign = noPoints, path, importToR = FALSE, treeFunction = igraphKamadaKawai, ...){
  if (!require(igraph))
    stop("simulate requires access to the igraph package")
  if (!require(maptools))  #AHF: added (needed for writeSpatialShape function)
    stop("simulate requires access to the maptools package")
  if (missing(obsDesign)) 
    stop("Input obsDesign cannot be missing")
  if (missing(path)) 
    stop("Path cannot be missing")
  if (missing(n))
    stop("Input n cannot be missing")
  if (length(path) != 1) 
    stop("Please enter a single path")
  info <- file.info(path)
  isdir <- info$isdir
  if (is.na(isdir)) {
    dir.create(path)
  }
  else if (isdir  ==  FALSE) {
    stop("Unable to create directory")
  }
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(path)
  n_networks <- length(n)
  edges <- vector(mode = "list", length = n_networks)
  tree.graphs <- edges
  locations <- edges
  rids <- edges
  initial_points <- n_edges <- vector(mode = "numeric", length = n_networks)
  edge_lengths <- edges
  cumulative_nedges <- 0
  max_x <- vector(mode = "numeric", length = n_networks)
  min_x <- vector(mode = "numeric", length = n_networks)
  for (i in 1:n_networks) {
    graph1 <- treeFunction(n[i], ...)
    graph <- graph1$DendNet
    edges_this_network <- get.edgelist(graph$graph)
    reordering <- order(edges_this_network[, 1])
    edges_this_network <- edges_this_network[reordering, 
                                             ]
    locations_this_network <- graph$locations[reordering]
    tree.graphs[[i]] <- graph.edgelist(edges_this_network)
    initial_points[i] <- graph$initialPoint
    locations_this_network <- graph$locations
    locations[[i]] <- locations_this_network
    edges[[i]] <- edges_this_network
    n_edges[i] = nrow(edges_this_network)
    edge_lengths_function <- function(indicies) {
      return(sqrt(sum((locations_this_network[indicies[1] + 0, ] #AHF: changed +1 to +0 bc of changes made in simul_tree
             - locations_this_network[indicies[2] + 0, ])^2)))
    }
    edge_lengths[[i]] <- apply(edges_this_network, 1, edge_lengths_function)
    rids[[i]] <- (1:n_edges[i]) + cumulative_nedges
    names(edge_lengths[[i]]) <- rids[[i]] 
    cumulative_nedges <- cumulative_nedges + n_edges[i]
    min_x[i] <- min(locations_this_network)
    max_x[i] <- max(locations_this_network)
  }
  cumulative_x <- max_x[1]
  if (n_networks > 1) {
    for (i in 2:n_networks) {
      locations[[i]][, 1] <- locations[[i]][, 1] + cumulative_x - in_x[i] + 0.1
      cumulative_x <- cumulative_x + max_x[i] - min_x[i] + 0.1
    }
  }
  spatial_edges <- vector(mode = "list", length = sum(n_edges))
  cumulative_nedges <- 0
  for (netid in 1:n_networks) {
    locations_this_network <- locations[[netid]]
    edges_this_network <- edges[[netid]]
    for (edge.index in 1:n_edges[netid]) {
      edge <- edges_this_network[edge.index, ]
      first.location = locations_this_network[edge[1] + 0, ] #AHF: changed +1 to +0
      second.location = locations_this_network[edge[2] + 0, ] #AHF: changed +1 to +0
      spatial_edges[[edge.index + cumulative_nedges]] = Lines(list(Line(rbind(second.location, 
           first.location))), ID = as.character(edge.index + cumulative_nedges))
    }
    cumulative_nedges = cumulative_nedges + n_edges[netid]
  }
  sl <- SpatialLines(spatial_edges)
  edge_updist <- vector(mode = "list", length = n_networks)
  line_data <- data.frame()
  shreve <- vector(mode = "list", length = n_networks)
  for (netid in 1:n_networks) {
    rids_this_network <- rids[[netid]]
    edges_this_network <- edges[[netid]]
    edge_updist_this_network = 0
    known_points = initial_points[netid]
    remaining_edges <- edges_this_network
    edge_lengths_this_network <- edge_lengths[[netid]]
    remaining_edge_lengths <- edge_lengths_this_network
    known_rids <- c()
    remaining_rids <- rids_this_network
    while (TRUE) {
      can_calculate <- (remaining_edges[, 1] %in% known_points) & 
        (!(remaining_edges[, 2] %in% known_points))
      upstream_point_indicies <- remaining_edges[can_calculate, 2]
      downstream_point_indicies <- remaining_edges[can_calculate, 1]
      edge_updist_this_network <- c(edge_updist_this_network, 
            edge_updist_this_network[match(downstream_point_indicies, known_points)] + remaining_edge_lengths[can_calculate])
      known_points <- c(known_points, upstream_point_indicies)
      remaining_edges <- remaining_edges[!can_calculate, , drop = FALSE]
      remaining_edge_lengths <- remaining_edge_lengths[!can_calculate]
      known_rids <- c(known_rids, remaining_rids[can_calculate])
      remaining_rids <- remaining_rids[!can_calculate]
      if (length(remaining_edges)  ==  0) 
        break
    }
    edge_updist_this_network <- edge_updist_this_network[-1][order(known_rids)]
    names(edge_updist_this_network) <- sort(known_rids)
    edge_updist[[netid]] <- edge_updist_this_network
    shreve_this_network <- vector(mode = "numeric", length = n_edges[netid])
    is_initial <- !(edges_this_network[, 2] %in% edges_this_network[, 1])
    known_points <- edges_this_network[which(is_initial), 2]
    remaining_points <- edges_this_network[which(!is_initial), 2]
    shreve_values_points <- rep(1, length(known_points))
    shreve_values_rids <- c()
    known_rids <- c()
    remaining_rids <- rids_this_network - min(rids_this_network) + 1
    remaining_edges <- edges_this_network
    while (TRUE) {
      can_calculate <- (remaining_edges[, 2] %in% known_points)
      shreve_values_rids <- c(shreve_values_rids, shreve_values_points[match(remaining_edges[can_calculate, 2], known_points)])
      remaining_edges <- remaining_edges[!can_calculate, , drop = FALSE]
      known_rids <- c(known_rids, remaining_rids[can_calculate])
      remaining_rids <- remaining_rids[!can_calculate]
      if (length(remaining_edges)  ==  0) 
        break
      can_calculate <- !(remaining_points %in% remaining_edges[, 1])
      can_calculate_function <- function(index) {
        relevant_edges <- which(edges_this_network[, 1]  ==  remaining_points[index])
        return(all(relevant_edges %in% known_rids))
      }
      can_calculate[can_calculate] <- sapply(which(can_calculate), can_calculate_function)
      calculate_shreve <- function(index) {
        relevant_edges <- which(edges_this_network[, 1]  ==  remaining_points[index])
        relevant_known_rids <- match(relevant_edges, known_rids)
        return(sum(shreve_values_rids[relevant_known_rids]))
      }
      new_shreve_values_points <- sapply(which(can_calculate), calculate_shreve)
      shreve_values_points <- c(shreve_values_points, new_shreve_values_points)
      known_points <- c(known_points, remaining_points[can_calculate])
      remaining_points <- remaining_points[!can_calculate]
    }
    known_rids_ordering <- order(known_rids)
    shreve_values_rids <- shreve_values_rids[known_rids_ordering]
    additive_function_values <- shreve_values_rids/max(shreve_values_rids)
    additional_line_data <- data.frame(rid = rids_this_network, netID = as.factor(rep(netid, n_edges[netid])), 
           upDist = edge_updist_this_network, shreve = shreve_values_rids, Length = edge_lengths_this_network, addfunccol = additive_function_values)
    rownames(additional_line_data) <- rids_this_network
    line_data <- rbind(additional_line_data, line_data)
  }
  sldf <- SpatialLinesDataFrame(sl, data = line_data, match.ID = TRUE)
  writeSpatialShape(sldf, "edges") #AHF: requires maptools library
  setwd(old_wd)
  
  
  binary_ids_tables <- list()
  for (netid in 1:n_networks) {
    binary_ids_tables[[netid]] <- graph1$bin.table 
    write.table(binary_ids_tables[[netid]], file = file.path(path, paste("netID", netid, ".dat", sep = "")), col.names = T, sep = ", ", row.names = FALSE)
  }
  setwd(path)
  distance_matrices <- list()
  for (netid in 1:n_networks) {
    edge_updist_this_network <- edge_updist[[netid]]
    edges_this_network <- edges[[netid]]
    binary_id_table <- binary_ids_tables[[netid]]
    distance_matrix <- matrix(0, nrow(binary_id_table), nrow(binary_id_table))
    colnames(distance_matrix) <- rownames(distance_matrix) <- binary_id_table$rid
    partial_match_function <- function(binary_id1, binary_id2) {
      min_len <- min(nchar(binary_id1), nchar(binary_id2))
      for (j in 1:min_len) {
        if (substr(binary_id1, j, j) != substr(binary_id2, j, j)) 
          return(j - 1)
      }
      return(min_len)
    }
    character_binary_ids <- as.character(binary_id_table$binaryID)
    for (i in 1:nrow(binary_id_table)) {
      current_binary_id <- binary_id_table$binaryID[i]
      current_rid <- as.character(binary_id_table$rid[i])
      current_updist <- edge_updist_this_network[current_rid]
      matching_characters <- sapply(character_binary_ids, partial_match_function, as.character(current_binary_id))
      matching_substring <- substr(binary_id_table$binaryID, 1, matching_characters)
      indices <- match(matching_substring, binary_id_table$binaryID)
      if (any(is.na(indices))) 
        stop("Internal Error")
      downstream_rids <- binary_id_table$rid[indices]
      downstream_updists <- edge_updist_this_network[as.character(downstream_rids)]
      distance_matrix[as.character(current_rid), ] <- pmax(current_updist - downstream_updists, rep(0, length(downstream_updists)))
    }
    reindex <- match(binary_id_table$rid, rids[[netid]])
    distance_matrix <- distance_matrix[reindex, reindex]
    distance_matrices[[netid]] <- distance_matrix + t(distance_matrix)
  }
  obs_sites <- obsDesign(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
  
  pred_sites <- predDesign(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
  max_observed_locID <- max(unlist(lapply(obs_sites, function(x) max(x$locID))))
  for (i in 1:length(pred_sites)) {
    pred_sites[[i]]$locID <- pred_sites[[i]]$locID + max_observed_locID
  }
  n_obs_sites <- unlist(lapply(obs_sites, function(x) return(dim(x)[1])))
  n_pred_sites <- unlist(lapply(pred_sites, function(x) return(dim(x)[1])))
  sites_data <- data.frame()
  combined_site_location_data <- c()
  pred_data <- data.frame()
  combined_pred_location_data <- c()
  cumulative_pids <- 0
  for (netid in 1:n_networks) {
    edge_lengths_this_network <- edge_lengths[[netid]]
    edges_this_network <- edges[[netid]]
    edge_updist_this_network <- edge_updist[[netid]]
    locations_this_network <- locations[[netid]]
    n_locations_this_network <- n_obs_sites[netid] + n_pred_sites[netid]
    rids_this_network <- rids[[netid]]
    pred_sites_this_network <- pred_sites[[netid]]
    obs_sites_this_network <- obs_sites[[netid]]
    f <- function(row) {
      rid <- as.character(row[1])
      edge_id <- match(rid, rids_this_network)
      proportion <- as.numeric(row[2])
      downstream_point <- edges_this_network[edge_id, 1] + 0 #AHF: replaced 1 with 0 due to igraph pkg version
      upstream_point <- edges_this_network[edge_id, 2] + 0 #AHF: replaced 1 with 0 due to igraph pkg version
      downstream_location <- locations_this_network[downstream_point, ]
      upstream_location <- locations_this_network[upstream_point, ]
      location <- downstream_location + proportion * (upstream_location - downstream_location)
      ret <- c(location, (sqrt(sum((location - downstream_location)^2)) + edge_updist_this_network[rid] - edge_lengths_this_network[rid]))
      names(ret) <- c("NEAR_X", "NEAR_Y", "upDist")
      return(ret)
    }
    obs_location_data <- data.frame(rid = obs_sites_this_network$edge, ratio = obs_sites_this_network$ratio, locID = obs_sites_this_network$locID, stringsAsFactors = FALSE)
    pred_location_data <- data.frame(rid = pred_sites_this_network$edge, ratio = pred_sites_this_network$ratio, locID = pred_sites_this_network$locID, stringsAsFactors = FALSE)
    if (n_locations_this_network > 0) {
      pred_location_data_this_network <- t(apply(pred_location_data, 1, f))
      if (length(pred_location_data_this_network)  ==  0) {
        pred_location_data_this_network <- matrix(0, 0, 3)
        colnames(pred_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
      }
      obs_location_data_this_network <- t(apply(obs_location_data, 1, f))
      obs_pids <- (1:n_obs_sites[netid]) + cumulative_pids
      if (n_pred_sites[netid] > 0) {
        pred_pids <- n_obs_sites[netid] + (1:n_pred_sites[netid]) + cumulative_pids
      }
      else pred_pids <- integer(0)
    }
    else {
      obs_location_data_this_network <- pred_location_data_this_network <- data.frame(NEAR_X = numeric(0), NEAR_Y = numeric(0), upDist = numeric(0))
      obs_pids <- integer(0)
      pred_pids <- integer(0)
    }
    cumulative_pids <- cumulative_pids + n_locations_this_network
    obs_data_this_network <- data.frame(locID = obs_location_data[, "locID"], upDist = obs_location_data_this_network[, "upDist"], 
          pid = obs_pids, netID = rep(netid, length(obs_pids)), rid = obs_location_data[, "rid"], ratio = obs_location_data[,  "ratio"], 
          shreve = line_data[match(obs_location_data[, "rid"], line_data[, "rid"]), "shreve"], addfunccol = line_data[match(obs_location_data[, "rid"], 
          line_data[, "rid"]), "addfunccol"], stringsAsFactors = FALSE)
    if (ncol(obs_sites_this_network) > 3) {
      obs_data_this_network <- cbind(obs_data_this_network, obs_sites_this_network[, -match(c("edge", "ratio", "locID"), colnames(obs_sites_this_network)), drop = FALSE])
    }
    rownames(obs_data_this_network) <- obs_pids
    rownames(obs_location_data_this_network) <- obs_pids
    pred_data_this_network <- data.frame(locID = pred_location_data[, "locID"], upDist = pred_location_data_this_network[, "upDist"], 
          pid = pred_pids, netID = rep(netid, length(pred_pids)), rid = pred_location_data[, "rid"], ratio = pred_location_data[,  "ratio"], 
          shreve = line_data[match(pred_location_data[, "rid"], line_data[, "rid"]), "shreve"], addfunccol = line_data[match(pred_location_data[, "rid"], 
          line_data[, "rid"]), "addfunccol"], stringsAsFactors = FALSE)
    if (ncol(pred_sites_this_network) > 3) {
      pred_data_this_network <- cbind(pred_data_this_network, pred_sites_this_network[, -match(c("edge", "ratio", "locID"), colnames(pred_sites_this_network)), drop = FALSE])
    }
    rownames(pred_data_this_network) <- pred_pids
    rownames(pred_location_data_this_network) <- pred_pids
    if (n_obs_sites[netid] > 0) {
      sites_data <- rbind(obs_data_this_network[1:n_obs_sites[netid], , drop = FALSE], sites_data)
      combined_site_location_data <- rbind(obs_location_data_this_network[1:n_obs_sites[netid], , drop = FALSE], combined_site_location_data)
    }
    if (n_pred_sites[netid] > 0) {
      pred_data <- rbind(pred_data_this_network[(1:n_pred_sites[netid]), , drop = FALSE], pred_data)
      combined_pred_location_data <- rbind(pred_location_data_this_network[(1:n_pred_sites[netid]), , drop = FALSE], combined_pred_location_data)
    }
  }
  if (length(combined_site_location_data)  ==  0) 
    stop("At least one observation site must be present")
  sites <- SpatialPointsDataFrame(combined_site_location_data[, c("NEAR_X", "NEAR_Y"), drop = FALSE], sites_data, match.ID = TRUE)
  writeSpatialShape(sites, "sites") #requires maptools library
  if (length(combined_pred_location_data) > 0) {
    preds <- SpatialPointsDataFrame(combined_pred_location_data[, c("NEAR_X", "NEAR_Y"), drop = FALSE], pred_data, match.ID = TRUE)
    writeSpatialShape(preds, "preds") 
  }
  setwd(old_wd)
  if (importToR) {
    if (sum(n_pred_sites) > 0) 
      return(importSSN(path, predpts = "preds", o.write = TRUE))
    else return(importSSN(path, o.write = TRUE))
  }
  else return(invisible(NULL))
}

# Uses above functions to generate a simulated stream network in the SSN framework
fncGenerateNetwork <- function(numSegs, lmin = 0.5, lmax = 4, sampleprob = 2, ps = 5, ordermax = 3, wt.spc = 0.5, 
  fish.spc = 0, numFish = 200, plotit = "file", path = paste0(getwd(), "/data.in/simulated", netnm, ".ssn"), seed = seed){
  
  #evenly and densely spaced fish; randomly distributed, can define how many fish to create
  ifelse(fish.spc > 0, pd <- systematicDesign(fish.spc), pd <- binomialDesign(numFish))
  
  set.seed(seed)
  
  ssn = list()
  while(is.list(ssn) == TRUE){ #FALSE means it created a network; TRUE means the algorithm got stuck attempting to add segments that don't cross
    try(ssn <- createSSNldotsBID(n = numSegs, 
               obsDesign = systematicDesign(wt.spc), #evenly and densely spaced obs points
               predDesign = pd, 
               importToR = TRUE, path=path, 
               treeFunction = simul_tree7toSSNBIDforce2trib, 
               ordermax=ordermax, alph=pi/3, beta=pi/8, lmin = lmin, lmax=lmax, sampleprob = sampleprob, ps=ps), 
        #see simul_tree7toSSNBIDforce2trib() for how parameters affect network shape
        silent=F)
  }
  
  #Import a second set of prediction points called "sites" for water temperature. We'll use "preds" for fish.
  ssn <- importPredpts(ssn, "sites", "ssn") #here, we are just duplicating the obs sites
  
  # Plot it
  if(plotit != "none"){
    if(plotit == "file") {png(paste0(path, "/locs.", netnm, ".png"), width = 6, height = 6, units = "in", res = 300)}
    pp <- plot(ssn, lwdLineCol = "addfunccol", lwdLineEx=8, linecol = "cyan3", xlab = "X", ylab = "Y", cex = 0) #cex = 0 makes obs points disappear
    #add fish points:
    ifelse(fish.spc > 0, pp2 <- plot(ssn, PredPointsID = "preds", add = TRUE, col = "blue", cex = 0.2), pp2 <- plot(ssn, PredPointsID = "preds", add = TRUE, col = "blue", cex = 0.7))
    #add temperature prediction points:
    pp3 <- plot(ssn, PredPointsID = "sites", add = TRUE, col = "orange", cex = 0.4)
    if(plotit == "file") {dev.off()}
    
    if(plotit == "file") {png(paste0(path, "/net.", netnm, ".png"), width = 6, height = 6, units =  "in", res = 300)}
    pp <- plot(ssn, lwdLineCol = "addfunccol", lwdLineEx = 8, linecol = "cyan3", xlab = "X", ylab = "Y", cex = 0) #cex = 0 makes obs points disappear
    if(plotit == "file") {dev.off()}
  }
  
  return(ssn)
}

# === CREATE THERMAL REGIMES =====================================================================

# Generate altered thermal regimes for analyses
fncGenAlteredTRs <- function(){
  #generate altered thermal regimes based on thermal regimes observed in the Snoqualmie River (Puget Sound, WA state)
  #temperatures are later spread across the network using fncUpdateWaterTemps()
  
  #contemporary: Snoqualmie mean across 4 normal years (2012, 2013, 2014, 2016)
  cscs <- c(13.59, 13.3, 13, 12.7, 12.6, 12.2, 11.84, 11.71, 11.58, 11.46, 11.33, 11.2, 11.08, 10.95, 10.83, 10.7, 10.58, 10.45, 10.33, 10.2, 10.08, 9.96, 9.83, 9.71, 9.59, 9.47, 9.34, 9.22, 9.1, 8.98, 8.87, 8.75, 8.63, 8.51, 8.4, 8.28, 8.17, 8.06, 7.95, 7.84, 7.73, 7.62, 7.52, 7.41, 7.31, 7.21, 7.11, 7.01, 6.91, 6.82, 6.72, 6.63, 6.55, 6.46, 6.37, 6.29, 6.21, 6.13, 6.05, 5.98, 5.91, 5.84, 5.77, 5.7, 5.64, 5.58, 5.52, 5.46, 5.41, 5.36, 5.3, 5.26, 5.21, 5.16, 5.12, 5.08, 5.04, 5.01, 4.97, 4.94, 4.9, 4.87, 4.85, 4.82, 4.79, 4.77, 4.75, 4.73, 4.71, 4.69, 4.67, 4.66, 4.65, 4.63, 4.62, 4.61, 4.61, 4.6, 4.6, 4.59, 4.59, 4.59, 4.59, 4.59, 4.59, 4.6, 4.6, 4.61, 4.62, 4.63, 4.64, 4.65, 4.67, 4.68, 4.7, 4.71, 4.73, 4.75, 4.77, 4.79, 4.82, 4.84, 4.86, 4.89, 4.91, 4.94, 4.97, 5, 5.02, 5.05, 5.08, 5.11, 5.15, 5.18, 5.21, 5.24, 5.28, 5.31, 5.35, 5.38, 5.42, 5.46, 5.49, 5.53, 5.57, 5.61, 5.65, 5.69, 5.73, 5.77, 5.81, 5.85, 5.89, 5.94, 5.98, 6.03, 6.07, 6.12, 6.17, 6.21, 6.26, 6.31, 6.36, 6.41, 6.46, 6.52, 6.57, 6.62, 6.67, 6.73, 6.78, 6.84, 6.9, 6.95, 7.01, 7.07, 7.13, 7.18, 7.24, 7.3, 7.36, 7.42, 7.48, 7.54, 7.6, 7.66, 7.72, 7.78, 7.85, 7.91, 7.97, 8.03, 8.09, 8.15, 8.22, 8.28, 8.34, 8.41, 8.4, 8.53, 8.6, 8.66, 8.72, 8.79, 8.85, 8.92, 8.99, 9.05, 9.12, 9.19, 9.25, 9.32, 9.39, 9.46, 9.53, 9.6, 9.67, 9.74, 9.81, 9.88, 9.95, 10.02, 10.09, 10.17, 10.24, 10.31, 10.39, 10.47, 10.54, 10.62, 10.7, 10.77, 10.85, 10.93, 11.02, 11.1, 11.18, 11.27, 11.36, 11.44, 11.53, 11.62, 11.72, 11.81, 11.91, 12, 12.1, 12.2, 12.3, 12.4, 12.51, 12.61, 12.72, 12.83, 12.94, 13.05, 13.16, 13.27, 13.39, 13.51, 13.63, 13.75, 13.87, 13.99, 14.12, 14.24, 14.37, 14.5, 14.62, 14.75, 14.88, 15.01, 15.15, 15.28, 15.41, 15.54, 15.67, 15.8, 15.93, 16.06, 16.19, 16.32, 16.44, 16.57, 16.69, 16.81, 16.93, 17.05, 17.16, 17.27, 17.38, 17.49, 17.6, 17.7, 17.79, 17.89, 17.98, 18.07, 18.15, 18.23, 18.31, 18.38, 18.45, 18.51, 18.57, 18.62, 18.67, 18.72, 18.75, 18.79, 18.82, 18.84, 18.86, 18.87, 18.88, 18.88, 18.87, 18.86, 18.84, 18.82, 18.79, 18.76, 18.72, 18.68, 18.63, 18.57, 18.51, 18.44, 18.37, 18.29, 18.21, 18.12, 18.03, 17.93, 17.83, 17.72, 17.61, 17.49, 17.37, 17.25, 17.13, 17, 16.87, 16.73, 16.6, 16.46, 16.31, 16.17, 16.02, 15.88, 15.73, 15.58, 15.43, 15.28, 15.13, 14.97, 14.82, 14.67, 14.51, 14.36, 14.21, 14.05, 13.9, 13.74, 13.59)
  
  #a warm dry year in the Snoqualmie with low snowpack (2015)
  wrws <- c(14.17, 14.02, 13.87, 13.72, 13.57, 13.42, 13.26, 13.11, 12.96, 12.81, 12.66, 12.51, 12.36, 12.21, 12.06, 11.91, 11.76, 11.61, 11.46, 11.31, 11.16, 11.01, 10.87, 10.72, 10.57, 10.43, 10.29, 10.14, 10, 9.86, 9.73, 9.59, 9.45, 9.32, 9.19, 9.06, 8.93, 8.8, 8.68, 8.56, 8.44, 8.32, 8.21, 8.09, 7.98, 7.88, 7.77, 7.67, 7.57, 7.48, 7.38, 7.29, 7.21, 7.12, 7.04, 6.96, 6.89, 6.81, 6.74, 6.68, 6.61, 6.55, 6.49, 6.43, 6.38, 6.33, 6.28, 6.23, 6.19, 6.14, 6.1, 6.06, 6.03, 5.99, 5.96, 5.93, 5.9, 5.88, 5.85, 5.83, 5.81, 5.79, 5.77, 5.75, 5.74, 5.73, 5.72, 5.71, 5.7, 5.7, 5.69, 5.69, 5.69, 5.69, 5.69, 5.7, 5.7, 5.71, 5.72, 5.73, 5.74, 5.76, 5.77, 5.79, 5.8, 5.82, 5.84, 5.86, 5.89, 5.91, 5.93, 5.96, 5.98, 6.01, 6.04, 6.07, 6.1, 6.13, 6.16, 6.19, 6.22, 6.25, 6.29, 6.32, 6.35, 6.39, 6.42, 6.46, 6.5, 6.53, 6.57, 6.61, 6.65, 6.68, 6.72, 6.76, 6.8, 6.84, 6.88, 6.92, 6.96, 7, 7.04, 7.08, 7.13, 7.17, 7.21, 7.26, 7.3, 7.35, 7.39, 7.44, 7.49, 7.54, 7.59, 7.64, 7.69, 7.74, 7.79, 7.84, 7.9, 7.95, 8.01, 8.06, 8.12, 8.18, 8.24, 8.3, 8.36, 8.43, 8.49, 8.56, 8.62, 8.69, 8.76, 8.83, 8.91, 8.98, 9.06, 9.13, 9.21, 9.29, 9.38, 9.46, 9.55, 9.64, 9.73, 9.82, 9.91, 10.01, 10.11, 10.21, 10.31, 10.42, 10.52, 10.63, 10.74, 10.86, 10.97, 11.09, 11.21, 11.33, 11.46, 11.58, 11.71, 11.84, 11.97, 12.11, 12.24, 12.38, 12.52, 12.66, 12.81, 12.95, 13.1, 13.25, 13.4, 13.55, 13.71, 13.86, 14.02, 14.18, 14.33, 14.49, 14.66, 14.82, 14.98, 15.14, 15.31, 15.47, 15.64, 15.81, 15.97, 16.14, 16.31, 16.48, 16.65, 16.81, 16.98, 17.15, 17.31, 17.48, 17.64, 17.81, 17.97, 18.13, 18.29, 18.45, 18.61, 18.76, 18.91, 19.06, 19.21, 19.36, 19.5, 19.65, 19.78, 19.92, 20.05, 20.18, 20.31, 20.43, 20.55, 20.67, 20.78, 20.89, 21, 21.1, 21.2, 21.29, 21.38, 21.47, 21.55, 21.62, 21.7, 21.76, 21.83, 21.88, 21.94, 21.99, 22.03, 22.07, 22.11, 22.14, 22.17, 22.19, 22.21, 22.22, 22.23, 22.23, 22.23, 22.23, 22.22, 22.21, 22.19, 22.17, 22.14, 22.11, 22.08, 22.04, 22, 21.95, 21.9, 21.85, 21.79, 21.73, 21.67, 21.6, 21.52, 21.45, 21.37, 21.28, 21.2, 21.1, 21.01, 20.91, 20.81, 20.7, 20.6, 20.48, 20.37, 20.25, 20.13, 20.01, 19.88, 19.75, 19.62, 19.48, 19.34, 19.2, 19.06, 18.92, 18.77, 18.62, 18.47, 18.32, 18.17, 18.01, 17.86, 17.7, 17.54, 17.38, 17.22, 17.06, 16.9, 16.74, 16.57, 16.41, 16.25, 16.08, 15.92, 15.75, 15.59, 15.42, 15.25, 15.09, 14.92, 14.75, 14.59, 14.42, 14.25, 14.09, 13.92, 13.75, 13.59)
  
  #warmer winters
  wscs <- cscs; wscs[30:165] <- wscs[30:165] + 1; wscs[90:120] <- wscs[90:120] + 1
  wscs2 <- gam(wscs ~ s(1:365, spar=0.9))$fitted
  foo <- cbind(round(cscs, 1), round(wscs2, 1))
  idx <- which(foo[, 1] == foo[, 2]); idx1 <- idx[1]; idx2 <- min(idx[idx > 200])
  wscs <- cscs; wscs[idx1:idx2] <- wscs2[idx1:idx2]
  hscs <- cscs; hscs[30:165] <- hscs[30:165] + 2; hscs[90:120] <- hscs[90:120] + 1.5
  hscs2 <- gam(hscs ~ s(1:365, spar=0.9))$fitted
  foo <- cbind(round(cscs, 1), round(hscs2, 1))
  idx <- which(foo[, 1] == foo[, 2]); idx1 <- 15; idx2 <- min(idx[idx > 200])
  hscs <- cscs; hscs[idx1:idx2] <- hscs2[idx1:idx2]
  
  #warmer summers
  csws <- cscs; csws[295:335] <- csws[295:335] * 1.3; csws[c(275:294, 336:355)] <- csws[c(275:294, 336:355)] * 1.15
  csws[355:365] <- cscs[355:365] * 0.6
  csws2 <- gam(csws ~ s(1:365, spar = 0.9))$fitted
  foo <- cbind(round(cscs, 1), round(csws2, 1))
  idx <- which(foo[, 1] == foo[, 2]); idx1 <- min(idx[idx > 200])
  csws <- cscs; csws[idx1:365] <- csws2[idx1:365]
  cshs <- cscs; cshs[295:335] <- cshs[295:335] * 1.5; cshs[c(275:294, 336:355)] <- cshs[c(275:294, 336:355)] * 1.3
  cshs[355:365] <- cscs[355:365] * 0.5
  cshs2 <- gam(cshs ~ s(1:365, spar = 0.9))$fitted
  foo <- cbind(round(cscs, 1), round(cshs2, 1))
  idx <- which(foo[, 1] == foo[, 2]); idx1 <- 218
  cshs <- cscs; cshs[idx1:365] <- cshs2[idx1:365]
  
  #more rapid spring warming
  crcs <- cscs; crcs[165:285] <- crcs[165:285] * 1.18; crcs[145:164] <- crcs[145:164] * 1.05; crcs[which.max(cscs)] <- max(cscs) * 3
  crcs2 <- gam(crcs ~ s(1:365, spar = 0.9))$fitted
  foo <- cbind(round(cscs, 1), round(crcs2, 1))
  idx <- which(foo[, 1] == foo[, 2]); idx1 <- min(idx[idx >= 100]); idx2 <- min(idx[idx > 300])
  crcs <- cscs; crcs[idx1:idx2] <- crcs2[idx1:idx2]
  crcs[which.max(crcs):which.max(cscs)] <- max(crcs)
  
  cics <- crcs; cics[165:275] <- cics[165:275] * 1.18; cics[145:164] <- cics[145:164] * 1.05
  cics2 <- gam(cics ~ s(1:365, spar = 0.9))$fitted
  foo <- cbind(round(cscs, 1), round(cics2, 1))
  idx <- which(foo[, 1] == foo[, 2]); idx1 <- min(idx[idx >= 100]); idx2 <- min(idx[idx > 300])
  cics <- cscs; cics[idx1:idx2] <- cics2[idx1:idx2]
  cics[which.max(cics):which.max(cscs)] <- max(cics)
  
  #variability increased
  cscv <- cscs * runif(365, 0.85, 1.15)
  csci <- cscv * runif(365, 0.85, 1.15)
  
  #combination scenarios
  hshs <- apply(cbind(hscs, cshs), 1, max)
  hihs <- apply(cbind(hshs, cics), 1, max)
  hihs[1:10] <- cscs[1:10] * 1.1; hihs[295:320] <- hihs[295:320] * 1.2; hihs[355:365] <- cscs[355:365] * 0.65
  hihs <- gam(hihs ~ s(1:365, spar = 0.9))$fitted
  
  out1 <- cbind(cscs, wscs, hscs, crcs, cics, csws, cshs, cscv, csci, wrws, hihs)
  out1 <- out1 * 1.3
  
  # Translating days to 12-h time steps
  out2 <- matrix(NA, 730, (dim(out1)[2] + 1)); colnames(out2) <- c("JD", colnames(out1))
  jd <- c("JD"=c(274:365, 1:273))
  jd1 <- jd + 0.1; jd2 <- jd + 0.6
  jd <- c(jd1, jd2)
  ts1 <- 273 * 2; ts2 <- ts1 + 1
  out2[, 1] <- sort(jd)[c(ts2:730, 1:ts1)]
  
  i <- 2
  for(cs in colnames(out1)){
    cc <- out1[, cs]
    out <- cbind("JD" = jd, c(cc, cc))
    out <- out[order(out[, "JD"]), ]
    out <- out[c(ts2:730, 1:ts1), ]
    out2[, i] <- out[, 2]
    i <- i + 1
  }
  
  return(list(out1, out2[, -1], out2[, 1]))
}

# Plot all climate scenarios together
fncPlotCC <- function(dat = therm.regimes365, plotit = "file"){
  
  if(JulianDate.Begin > 1){xvals <- c(JulianDate.Begin:730, 1:(JulianDate.Begin - 1))}
  
  cscens.list <- NULL
  #increase winter
  cscens.list[[1]] <- c("cscs", "wscs", "hscs")
  #increase summer
  cscens.list[[2]] <- c("cscs", "csws", "cshs")
  #increase spring warming
  cscens.list[[3]] <- c("cscs", "crcs", "cics")
  #increase variability
  cscens.list[[4]] <- c("cscs", "cscv", "csci")
  names(cscens.list) <- c("Warmer Winter", "Warmer Summer", "More Rapid Spring Warming", "Increased Variability")
  
  #col1 <- 1; col2 <- "gray50"
  col3 <- rgb(244, 155, 2, 255, NULL, 255) #orange
  col2 <- rgb(206, 193, 8, 255, NULL, 255) #yellow
  col1 <- rgb(3, 38, 178, 255, NULL, 255) #blue
  
  if(plotit == "file") png(paste0(getwd(), local.path, plotDir, "/1ThermalRegimes.png"), width = 4, height = 10, units = "in", res = 600)
  par(mfrow=c(4, 1), las=1, mar=(c(2, 4, 3, 0)+0.5), oma=c(0, 0, 0, 1))
  
  cscens <- cscens.list[[1]]
  plot(dat[, "cscs"], ylim = c(2, 38), type = 'n', ylab = expression("Water temperature "~(degree ~ C)), xlab = "Date", axes=F, main = names(cscens.list[1]))
  axis(side=2, ylim = c(0, 30)); axis(side=1, at=c(0, 74, 147, 220, 293, 365)*2, labels=c("1 Oct", "13 Dec", "24 Feb", "8 May", "20 Jul", "30 Sep")) 
  lwd = 2; lty = 1
  nms <- cscens[cscens!="cscs"]
  lines(dat[, "wscs"], col = col2, lty = lty, lwd = lwd)
  lines(dat[, "hscs"], col = col3, lty = lty, lwd = lwd)
  lines(dat[, "cscs"], lwd = 2, col = col1)
  abcol <- "gray10"; abln <- 3
  abline(h=10, lty = abln, col = abcol); abline(h=20, lty = abln, col = abcol); abline(h=30, lty = abln, col = abcol)
  leg <- c(nms[2], nms[1], "cscs")
  legend("bottomright", legend=leg, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c(col3, col2, col1), box.col = "white", bty = 'o', cex = 0.9)
  box()
  text(10, 35, "(a)", cex = 1.1)
  
  cscens <- cscens.list[[3]]
  plot(dat[, "cscs"], ylim = c(2, 38), type = 'n', ylab = expression("Water temperature "~(degree ~ C)), xlab = "Date", axes=F, main = names(cscens.list[3]))
  axis(side=2, ylim = c(0, 30)); axis(side=1, at=c(0, 74, 147, 220, 293, 365)*2, labels=c("1 Oct", "13 Dec", "24 Feb", "8 May", "20 Jul", "30 Sep")) 
  lwd = 2; lty = 1
  nms <- cscens[cscens!="cscs"]
  lines(dat[, "crcs"], col = col2, lty = lty, lwd = lwd)
  lines(dat[, "cics"], col = col3, lty = lty, lwd = lwd)
  lines(dat[, "cscs"], lwd = 2, col = col1)
  abcol <- "gray10"; abln <- 3
  abline(h=10, lty = abln, col = abcol); abline(h=20, lty = abln, col = abcol); abline(h=30, lty = abln, col = abcol)
  leg <- c(nms[2], nms[1], "cscs")
  legend("bottomright", legend=leg, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c(col3, col2, col1), box.col = "white", bty = 'o', cex = 0.9)
  box()
  text(10, 35, "(b)", cex = 1.1)
  
  cscens <- cscens.list[[2]]
  plot(dat[, "cscs"], ylim = c(2, 38), type = 'n', ylab = expression("Water temperature "~(degree ~ C)), xlab = "Date", axes=F, main = names(cscens.list[2]))
  axis(side=2, ylim = c(0, 30)); axis(side=1, at=c(0, 74, 147, 220, 293, 365)*2, labels=c("1 Oct", "13 Dec", "24 Feb", "8 May", "20 Jul", "30 Sep")) 
  col <- col1; lwd = 2; lty = 1
  nms <- cscens[cscens!="cscs"]
  lines(dat[, "csws"], col = col2, lty = lty, lwd = lwd)
  lines(dat[, "cshs"], col = col3, lty = lty, lwd = lwd)
  lines(dat[, "cscs"], col = col1, lwd = 2)
  abcol <- "gray10"; abln <- 3
  abline(h=10, lty = abln, col = abcol); abline(h=20, lty = abln, col = abcol); abline(h=30, lty = abln, col = abcol)
  leg <- c(nms[2], nms[1], "cscs")
  legend("bottomright", legend=leg, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c(col3, col2, col1), box.col = "white", bty = 'o', cex = 0.9)
  box()
  text(10, 35, "(c)", cex = 1.1)
  
  cscens <- cscens.list[[4]]
  plot(dat[, "cscs"], ylim = c(2, 38), type = 'n', ylab = expression("Water temperature "~(degree ~ C)), xlab = "Date", axes=F, main = names(cscens.list[4]))
  axis(side=2, ylim = c(0, 30)); axis(side=1, at=c(0, 74, 147, 220, 293, 365)*2, labels=c("1 Oct", "13 Dec", "24 Feb", "8 May", "20 Jul", "30 Sep")) 
  col <- col1; lwd = 1; lty = 1
  nms <- cscens[cscens!="cscs"]
  lines(dat[, "csci"], col = col3, lty = lty, lwd = lwd)
  lines(dat[, "cscv"], col = col2, lty = lty, lwd = lwd)
  lines(dat[, "cscs"], col = col1, lwd = 2)
  abcol <- "gray10"; abln <- 3
  abline(h=10, lty = abln, col = abcol); abline(h=20, lty = abln, col = abcol); abline(h=30, lty = abln, col = abcol)
  leg <- c(nms[2], nms[1], "cscs")
  legend("bottomright", legend=leg, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c(col3, col2, col1), box.col = "white", bty = 'o', cex = 0.9)
  box()
  text(10, 35, "(d)", cex = 1.1)
  
  if(plotit == "file") dev.off()
  
}

# === SPATIOTEMPORAL WATER TEMPERATURE PATTERNS ==================================================

# Generate initial stream temperature observations
#   Used when no data are available (i.e., from GIS or models)
fncInitializeWT <- function(cc1 = cc[1], ssn = ssn1, plotit = "file"){
  
  # Get obs points from SSN
  tmpvar <- getSSNdata.frame(ssn, "Obs") 
  
  # Create initial stream temperature data
  #create spatial autocorrelation based on stream order (cooler upstream, warmer downstream)
  #v <- fncRescale(log(tmpvar[, so.field] + (1 - tmpvar[, "ratio"])), c(0, 1)) #accounts for position of point within reach using 'ratio'
  #here, instead using exponential relationship to keep headwaters cool and mainstems warm
  z <- log(tmpvar[, so.field] + (1 - tmpvar[, "ratio"]))
  v <- fncRescale(z^2, c(0, 1))

  tmpvar[, wt.field] <- (runif(nrow(tmpvar), 15, 20) * v) + 8  #runif adds some stochasticity; +8 to make a warm network
  
  # Make some tributary watersheds cooler and some warmer
  cooltribs <- NULL
  idx <- as.numeric(unique(as.character(tmpvar$rid[tmpvar[, so.field] < 15 & tmpvar[, so.field] > 4])))
  (cools <- sample(idx, 3))
  for(i in 1:length(cools)){cooltribs <- fncTraceUp(seg = cools[i], ssn = ssn, plotit = F, col = "blue")}
  tmpvar[tmpvar$rid %in% as.factor(cools), wt.field] <- tmpvar[tmpvar$rid %in% as.factor(cools), wt.field] - 2
  
  warmtribs <- NULL
  idx <- as.numeric(unique(as.character(tmpvar$rid[tmpvar[, so.field] < 15 & tmpvar[, so.field] > 4])))
  (warms <- sample(idx, 3))
  for(i in 1:length(warms)){warmtribs <- fncTraceUp(seg=warms[i], ssn = ssn, plotit = F, col = "red")}
  tmpvar[tmpvar$rid %in% as.factor(warms), wt.field] <- tmpvar[tmpvar$rid %in% as.factor(warms), wt.field] + 2
  
  # To plot:
  #   pp <- plot(ssn1, lwdLineCol = "addfunccol", lwdLineEx = 15, linecol = "dodgerblue", xlab = "X", ylab = "Y", cex = 0)
  #   op <- par(usr=pp$usr)
  #   for(i in 1:length(cools)){cooltribs <- fncTraceUp(seg = cools[i], ssn = ssn, plotit = T, col = "blue")}
  #   for(i in 1:length(warms)){warmtribs <- fncTraceUp(seg = warms[i], ssn = ssn, plotit = T, col = "red")}
  
  # Add some local spatial variation by making some points cooler/warmer than others
  # these represent cold water refuges or thermal barriers that are caused by local
  # drivers and not picked up by spatial autocorrelation or climate covariates
  (cooler <- sample(unique(tmpvar$pid), 35))
  (warmer <- sample(unique(tmpvar$pid), 35))
  tmpvar[cooler, wt.field] <- tmpvar[cooler, wt.field] - 2
  tmpvar[warmer, wt.field] <- tmpvar[warmer, wt.field] + 2
  
  # Ensure none go below zero (Bioenergetics doesn't handle negative temps)
  tmpvar[tmpvar[, wt.field] < 0.1, wt.field] <- 0.1
  
  # Add scaled water temperature column
  # This is the spatial scaffolding of relative water temperatures over which temporal patterns in climate scenario thermal regimes will be draped
  tmpvar$scaledWT <- (tmpvar[, wt.field] - min(tmpvar[, wt.field])) / (max(tmpvar[, wt.field]) - min(tmpvar[, wt.field]))
  #tmpvar$scaledWT <- fncRescale(tmpvar[, wt.field], c(0, 1)) #does exact same thing
  
  # Alter temperature data for each point according to climate scenario for the first time step
  #shift all values up/down by a constant proportion of what they were
  #gives values ranging across the network from min to WTcc.jd
  tmpvar[, wt.field] <- (cc1-min(tmpvar[, wt.field]))*tmpvar$scaledWT+ min(tmpvar[, wt.field])
  
  # Put it back in ssn
  ssn <- putSSNdata.frame(tmpvar, ssn, "Obs") 
  
  # Plot
  if(plotit!="none"){
    cb <- brewer.pal(10, "Spectral"); cb <- cb[c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)]
    if(plotit == "file"){png(paste0(plotDir, "/WTi", scenario, ".png"), width = 6, height = 6, units = "in", res = 300)}
    plot(ssn, "WT", breaktype = "user", brks = seq(5, 27, length.out = 11), lwdLineCol = "addfunccol", lwdLineEx = 8, linecol = "darkgray", xlab = "X", ylab = "Y", color.palette = cb)
    if(plotit == "file") {dev.off()}
  }
  
  return(ssn)
}

# Update each reach with new temperature data based on climate scenario and day
#   Used when predictions are based on fncInitializeWT() and fncWTcc(); not used with empirical or modeled data
fncUpdateWaterTemps <- function(ts, WTu, cc = cc1, ssn = ssn1, type = "Obs", diel = T, plotit = "file"){
  #WTu is a matrix of spatially distributed water temperatures from ssn
  #ts is the timestep; cc is the thermal regime scenario at the warmest location
  
  #ts = timeStep[tt]
  WTcc.jd <- cc[, "cc"][cc[, "JD"] == ts]
  ind <- which(cc[, "JD"] == ts)
  ifelse(ts - floor(ts) < 0.2, night <- TRUE, night <- FALSE)
  daily.chg <- ifelse(night == T, (WTcc.jd - cc[ind - 1, "cc"]) / 2, (WTcc.jd - cc[ind - 2, "cc"]) / 2)
  if(is.na(daily.chg)) daily.chg <- 0
  
  WTu[, wt.field] <- ((WTu[, wt.field] + daily.chg) + (WTu[, wt.field] + daily.chg * WTu[, "scaledWT"])) / 2
  
  # Add diel variability and apply it disproportionally more to mid-order reaches
  #(evidence suggests most diel variability occurs in mid-order reaches)
  if(diel  == T){
    diel.var <- 0.2 * WTcc.jd 
    if(night == T) diel.var <- -diel.var 
    WTu[WTu$addfunccol <= 0.3 & WTu$addfunccol >= 0.1, wt.field] <- WTu[WTu$addfunccol <= 0.3 & WTu$addfunccol >= 0.1, wt.field] + diel.var
    WTu[WTu$addfunccol > 0.3 | WTu$addfunccol < 0.1, wt.field] <- WTu[WTu$addfunccol > 0.3 | WTu$addfunccol < 0.1, wt.field] + diel.var / 2
  }
  
  # Keep water from freezing (bioenergetics calcs can't handle zero temps)
  WTu[WTu[, wt.field] < 0.01, wt.field] <- 0.01
  
  # Put new data back into SSN
  ssn <- putSSNdata.frame(WTu, ssn, "Obs")
  
  if(plotit!="none"){
    if(plotit == "file"){png(paste0(plotDir, "/WTu", scenario, "-", t, ".png"), width = 6, height = 6, units = "in", res = 300)}
    plot(ssn, "WT", breaktype = "user", brks = c(0, seq(5, 29, length.out = 9), 34), lwdLineCol = "addfunccol", lwdLineEx = 8, linecol = "darkgray", xlab = "X", ylab = "Y", color.palette = cb)
    if(plotit == "file"){dev.off()}
  }
  
  return(ssn)
}

# === HABITAT FUNCTIONS ==========================================================================
# Get water temperature nearest to a fish
fncGetNearestWT <- function(dat.df = dat.df[, c("pid", "rid", "ratio", wt.field)], fish = fish[, c("pid", "seg", "ratio")], ssn = ssn1){
  fish <- as.matrix(fish)
  dat <- as.matrix(dat.df)
  wt.out <- matrix(NA, length(unique(fish[, "pid"])), 4)
  missed.tmp <- c(missed, 0) #add seg = 0 to the list of segs without temperature data
  thesegs <- unique(fish[, "seg"]) #each segment that has fish in it
  ifelse(length(missed) > 0, r.missed <- thesegs %in% missed.tmp, r.missed <- FALSE)
  nn <- 1
  for(r in thesegs[r.missed]){ #for segs with fish but no WT data 
    n.fish <- length(fish[fish[, "seg"] == r, "pid"])
    tmp.fish <- matrix(fish[fish[, "seg"] == r], n.fish, 3) #get fish in that seg
    colnames(tmp.fish) <- c("pid", "seg", "ratio")
    for(i in 1:nrow(tmp.fish)){ #for each of these fish
      if(tmp.fish[, "ratio"][i] >= 0.5){ #if the fish is closer to the top of the seg
        #get WTs from seg or tribs flowing into seg (handles internally if these don't have WT values) and record the mean water temp of these segs
        wt <- mean(fncGetJunctionWTs(dat[, c("pid", "rid", wt.field)], r, ssn)[, "WT"])
      }
      else if(tmp.fish[, "ratio"][i] < 0.5){ #if the fish is closer to the bottom of the seg
        dnseg <- getSegDownstream(r, ssn) #get next seg downstream
        if(!dnseg %in% missed.tmp){ #if the downstream seg has a water temp value
          tmp.dat <- dat[dat[, "rid"] == dnseg, ]
          if(is.vector(tmp.dat)) tmp.dat <- matrix(tmp.dat, 1, 4); colnames(tmp.dat) <- c("pid", "rid", "ratio", wt.field)
          wt <- tmp.dat[which.max(tmp.dat[, "ratio"]), wt.field] #take the highest value from that seg
        } else { #if the downstream seg does not have a water temp value
          #get WTs from seg or tribs flowing into seg (handles internally if these don't have WT values) and record the mean water temp of these segs
          wt <- mean(fncGetJunctionWTs(dat[, c("pid", "rid", wt.field)], dnseg, ssn)[, "WT"])
        }
      }
      wt.out[nn, ] <- cbind(tmp.fish[i, "pid"], r, tmp.fish[i, "ratio"], wt)
      nn <- nn + 1
    }
  }
  for(r in thesegs[!r.missed]){ #for segs with fish that do have WT data
    n.fish <- length(fish[fish[, "seg"] == r, "pid"])
    tmp.fish <- matrix(fish[fish[, "seg"] == r], n.fish, 3) #get fish in that seg
    colnames(tmp.fish) <- c("pid", "seg", "ratio")
    tmp.dat <- dat[dat[, "rid"] == r, ] #get WT data for the seg
    if(is.vector(tmp.dat)) tmp.dat <- matrix(tmp.dat, 1, 4); colnames(tmp.dat) <- c("pid", "rid", "ratio", wt.field)
    for(i in 1:nrow(tmp.fish)){ #for each of these fish
      idx <- which.min(abs(tmp.fish[i, "ratio"]-tmp.dat[, "ratio"]))
      wt <- tmp.dat[idx, wt.field] #record the water temp value in the seg that is closest to the fish's position
      wt.out[nn, ] <- cbind(tmp.fish[i, "pid"], r, tmp.fish[i, "ratio"], wt)
      nn <- nn + 1
    }
  }
  
  colnames(wt.out) <- c("pid", "seg", "ratio", wt.field)
  wt.out <- wt.out[order(wt.out[, "pid"]), ]
  return(wt.out)
}

# Get water temperatures for each of the (1 to 3) options at each junction/confluence
fncGetJunctionWTs <- function(dat.df = dat.df[, c("pid", "rid", wt.field)], seg = seg, ssn = ssn){
  dat <- as.matrix(dat.df)
  if(seg == 0) seg <- 1
  out <- as.matrix(jct.list[[seg]]) #get stored list of segments at the upper junction of 'seg'
  if(length(unique(out[, "pt"])) == 1) { #if there is only one water temp point, return this value
    ret <- cbind("r" = out[, "r"], "WT" = dat[dat[, "pid"] %in% out[, "pt"], wt.field])
  } else { #if there are multiple water temp points per segment
    out <- cbind(out, "WT" = rep(NA, nrow(out)))
    for(i in 1:nrow(out)){ #collect them all
      out[i, "WT"] <- dat[dat[, "pid"] == out[i, "pt"], wt.field]
    }
    ret <- tapply(out[, "WT"], out[, "r"], mean, na.rm=T)  #return the mean of all water temp values on a seg, for each seg in junction
    ret <- cbind("r" = as.integer(names(ret)), "WT" = ret); rownames(ret) <- NULL
  }
  return(ret)
}

# Fill water temperatures for segs with no predictions, based on surrounding values
#(sometimes happens in GIS datasets; shouldn't happen in generated virtual networks)
fncFillNullPreds <- function(seglist = missed, ssn = ssn){
  out <- matrix(NA, length(seglist), 2)
  for(i in 1:length(seglist)){
    out[i, ] <- cbind(seglist[i], mean(fncGetJunctionWTs(obs.df[, c("pid", "rid", wt.field)], seglist[i], ssn)[, "WT"]))
  }
  colnames(out) <- c("seg", "WT")
  return(out)
}

# Calculate the density of species1, species2, and total fish in a segment that
# a species1 fish inhabits
# Input:
#   fish1.seg = vector of seg, length is number of species1. Required.
#   fish2.seg = vector of seg, length is number of species2. Can set to NULL.
#   fish1.alive.emerged.segs = vector of seg, for species1 that are alive and emerged. Required.
#   fish2.alive.emerged.segs = same as above for species2. Can set to NULL.
#   rid = vector of segment rid from ssn object, length is number of segments
#   length = vector of segment length (km) from ssn object, length is number of segments
# Output:
#   matrix where nrow = number of fish for species1, ncol = 3
#   first column is density of species 1 fish in that fish's segment
#   second column is density of species 2 fish in that fish's segment
#   third column is density of total fish in that fish's segment
fncFishDensity <- function(fish1.seg = fish1[, "seg"], fish2.seg = fish2[, "seg"], corr.factor=corr.factor, 
                          fish1.alive.emerged.segs= fish1[, "seg"][fish1[, "emrg"] == 1 & fish1[, "survive"] == 1], 
                          fish2.alive.emerged.segs= fish2[, "seg"][fish2[, "emrg"] == 1 & fish2[, "survive"] == 1], 
                          ssn = ssn){
  rid <- ssn@data$rid
  length <- ssn@data[, length.field]
  width <- ssn@data$UseableWidth
  num.fish1 <- length(fish1.seg)
  density.species1 <- vector("numeric", num.fish1)
  density.species2 <- vector("numeric", num.fish1)
  density.total <- vector("numeric", num.fish1)
  # multiply by 1000 because length is in km, width is in m
  # area is in m^2
  area <- length * 1000 * width
  for(i in 1:num.fish1){
    num.fish1.in.seg <- sum(fish1.alive.emerged.segs == fish1.seg[i]) #Here, we only include fish that are alive and emerged
    num.fish2.in.seg <- sum(fish2.alive.emerged.segs == fish1.seg[i]) #here too
    seg.area = area[which(rid  ==  fish1.seg[i])]
    density.species1[i] <- num.fish1.in.seg / seg.area
    density.species2[i] <- num.fish2.in.seg / seg.area
  }
  density.total <- density.species1 + density.species2
  density <- matrix(c(density.species1, density.species2, density.total), 
                    nrow = num.fish1, ncol = 3)
  #apply a correction factor:
  density <- density * corr.factor 
  #this treats each individual like a "superindividual" to reduce the resources required
  #to run (doing so for this many individual fish could be computationally expensive). 
  #This adjustment makes densities more representative of observed fish densities.
  
  return(density)
}

# Estimate the width of different reaches that is useable to fish
fncUseableWidths <- function(dat=ssn@data[, c("rid", "WIDTH_M")]){
  widths <- dat[, 2]
  rid <- dat[, 1]
  prop.useable <- 1 - fncRescale(log(widths), to = c(0, 0.9)) #assuming 10% of really wide rivers is still useable, and 100% of tribs is useable habitat
  useable.widths <- widths * prop.useable
  return(cbind(rid, widths, useable.widths))
}

# === MOVEMENT FUNCTIONS =========================================================================

# Calculates movement distance potentials based on fish densities and fish sizes
# returns data frame of max movement distances for each fish and whether they will move or not (zi)
fncMoveDistance <- function(ssn){
    
  # get global variables
  mvmt.multiplier <- get("mvmt.multiplier")
  mvdst.sdlog <- get("mvdst.sdlog")
  wt.field <- get("wt.field")
  fpids <- fish[, "pid"]
  fdens <- fish[, "density"]
  nFish <- length(fdens)
  
  # look up growth at fish’s location based on its WT, ration, and weight
  growths <- matrix(NA, nFish, 3, dimnames = list(fpids, c("growth", "growth.pot", "growth.min")))
  growths[, "growth"] <- fncGrowthFish(fpids)[, "growth"]
  
  # look up max growth possible for each fish across all accessible reaches during this time step
  for(fpid in fpids){
    growths[, "growth.pot"][fpid] <- max(fncGrowthPossible(fish[, "weight"][fpid], ssn)[, "growth"])
    growths[, "growth.min"][fpid] <- min(fncGrowthPossible(fish[, "weight"][fpid], ssn)[, "growth"])
    }
  
  #movement probability and distance is influenced by how close fish's current growth is to max growth potential
    #bigger values mean higher probability of longer movements
  gr.diff <- apply(growths[, c("growth", "growth.pot")], 1, diff) / apply(growths[, c("growth.min", "growth.pot")], 1, diff)

  #zero-inflated part: will a fish move at all?
  zi <- rep(1, length(fdens))
  #if growth is within 10% of max possible growth, fish does not move
  growth.move.idx <- which(abs(gr.diff) <= 0.1)
  zi[growth.move.idx] <- 0
  
  #lognormal move distance (if fish is moving)
  moveDist <- rlnorm(nFish, meanlog = (gr.diff * mvmt.multiplier), sdlog = mvdst.sdlog)
  moveDist <- cbind(fpids, moveDist, zi)
  colnames(moveDist) <- c("pid", "moveDist", "zi")
  
  return(moveDist)
}

# Determines if a fish can stop before its assigned movement distance if it encounters good habitat; returns list(T/F, pStop, dist2move)
fncStopEarly <- function(fpid, seg, remainingDist, moveLength, length2segBase, ssn){
  # the probability of stopping increases as growth potential increases
  
  length.field <- get("length.field")
  
  #look up growth possible for the fish at all locations in this reach
  result <- fncGrowthPossible(fish[, "weight"][fpid], ssn, seg)
  growth.here <- max(result[, "growth"])
  best.obs <- result[, "pid"][which.max(result[, "growth"])]
  
  #look up growth possible for each fish across accessible reaches during this time step
  result2 <- fncGrowthPossible(fish[, "weight"][fpid], ssn)

  #growth in reach relative to max growth potential in stream network
    #bigger values mean higher probability of stopping
  probStop <- 1 - diff(c(growth.here, max(result2[, "growth"]))) / diff(range(result2[, "growth"]))
  
  #determine position of fish relative to obs to stop at
  obs.seg.ratio <- obs.df$ratio[obs.df$rid == seg & obs.df$pid == best.obs] #seg ratio at best.obs
  obs.dist2Base <- ssn@data[ssn@data[, "rid"] == seg, length.field] * obs.seg.ratio
  ifelse(length2segBase > obs.dist2Base, moveDirection <- 'down', moveDirection <- 'up') #if equal, keep original moveDirection
  
  dist2move <- abs(length2segBase - obs.dist2Base)
  
  pStop <- c(probStop, (1-probStop)) #prob. of stopping, prob. of continuing
  StopEarly <- sample(c(T, F), size = 1, prob = pStop) #T = stop, F = continue
  
  return(list(StopEarly, pStop, dist2move, moveDirection))
}

# Determines if there's room within the segment to move (segLength is required if moving upstream); returns T/F
fncRoom2Move <- function(seg, moveDirection, remainingDist, length2segBase, ssn) {
  
  length.field = get("length.field")
  if (moveDirection == 'down' & remainingDist <= length2segBase) {
    return(TRUE)
  } else if (moveDirection == 'up' & remainingDist <= (ssn@data[ssn@data$rid == seg, length.field] - length2segBase))  {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Determines which segment a fish moves into as it encounters junctions; returns 'seg'
fncJunctionDecisions <- function(fpid, seg, currMoveDirection, stochastic = TRUE, ssn) {
  # As a fish gets to a junction, whether heading up or down
  # it looks at all possible options for movement and chooses the best one
  
  seg.orig <- seg # store the original segment
  
  # We will always look at the junction upstream of a particular segment, 
  # So, if a fish is heading downstream
  # it's easiest to go to the downstream segment and look at all options from there
  if(currMoveDirection  ==  "down"){
    if(seg > min(ssn@data$rid)){ seg <- getSegDownstream(seg)}
  }
  
  # Get conditions at downstream segment and at each of the upstream segments
  
  # first get which reaches participate in the junction
  rids = jct.list[[seg]][, "r"]
  #(inaccessible reaches will be removed in the next step)
  
  #Look up previously calculated growth values
  growth.confluence <- fncGrowthPossible(fish[fpid, "weight"], ssn, rids)
  
  if(length(growth.confluence[, "growth"]) == 1 & is.na(growth.confluence[, "growth"][1])){
    bestSeg <- seg
  } else {
    gr.segs <- tapply(growth.confluence[, "growth"], growth.confluence[, "rid"], mean)
    if(any(gr.segs < 0)) gr.segs <- gr.segs + 2 * abs(min(gr.segs)) #ensure no negative probabilities

  #choose best seg
    if(stochastic  ==  TRUE){
      # base probabilities on relative growth potential in each reach
      ifelse(sum(gr.segs) > 0, probs <- gr.segs / sum(gr.segs), probs <- rep(0.5, length(gr.segs)))
      probs <- gr.segs/sum(gr.segs)
      ifelse(length(gr.segs) > 1, bestSeg <- sample(as.integer(names(gr.segs)), size = 1, prob = probs), bestSeg <- seg)
    } else {
      # choose seg deterministically
      bestSeg <- as.integer(names(which.max(gr.segs)))
    }
  }
  
  return(bestSeg)
}

# Determines a fish's position based on it's current circumstances; returns list(dist2move, remainingDist, length2segBase, segRatio)
fncUpdatePosition <- function(newSeg, seg, moveDirection, length2segBase, dist2move, remainingDist, StopEarly, Rm2Move, ssn){
  #where seg was original seg and newSeg is the result of fncJunctionDecisions
  
  if(StopEarly == T | Rm2Move == T){ #if fish stopped in the reach, update its position
    
    if (moveDirection == 'down') {
      length2segBase <- length2segBase - dist2move 
    } else if (moveDirection == 'up') {
      length2segBase <- length2segBase + dist2move 
    }
    # Then update it's segRatio
    segRatio <- length2segRatio(seg, length2segBase)
    remainingDist <- 0
    #seg and moveDirection stay the same
    
  } else { #move fish to next upstream or downstream reach and update its position
    
    if(moveDirection == 'up'){ # if heading upstream
      
      dist2move <- ssn@data[ssn@data$rid == seg, length.field] - length2segBase #travel the rest of the way up the original reach
      
      # if no accessible upstream reaches or if best choice is to stay in same reach, 
      # then locate fish at top of existing reach facing downstream
      if(sum(ssn@data$accessible[ssn@data$rid %in% getSegUpstream(seg)]) == 0 | newSeg == seg){
        segRatio = 1
        length2segBase = ssn@data[ssn@data$rid == seg, length.field] #length of current reach
        moveDirection = 'down'
        
      } else{ #otherwise, 
        
        #if newSeg is upstream of seg or if both seg and newSeg are upper branches
        ifelse (getSegDownstream(seg) > 0, test.seg <- getSegUpstream(getSegDownstream(seg)), test.seg <- getSegUpstream(seg)) #(no seg zero)
        if(newSeg %in% getSegUpstream(seg) | newSeg %in% test.seg){ 
          # then locate fish at the base of the new reach and maintain original direction
          segRatio <- 0 
          length2segBase <- 0
        } 
      }
      
    } else if(moveDirection == 'down'){ # if heading downstream
      
      dist2move <- length2segBase #travel to base of the reach
      
      # if no accessible downstream reaches or if best choice is to stay in same reach, 
      # then locate fish at bottom of existing reach facing upstream
      if(newSeg <= min(ssn@data$rid) | seg == newSeg){
        segRatio <- 0 
        length2segBase <- 0
        moveDirection <- 'up'
        if(newSeg < min(ssn@data$rid)) newSeg <- min(ssn@data$rid)
        
      } else{ #otherwise, 
        
        if(newSeg %in% getSegDownstream(seg)){ #if newSeg is downstream of seg
          # then locate fish at the top of the new reach and maintain original direction
          segRatio <- 1
          length2segBase <- ssn@data[ssn@data$rid == newSeg, length.field]
          
        } else { #i.e., newSeg and seg are both upper branches
          # then locate fish at the bottom of the new reach facing upstream
          segRatio <-0
          length2segBase <- 0
          moveDirection <- 'up'
        }
      }
    } #end move up or down
    
    #update distance remaining (i.e., that a fish can still move this time step)
    remainingDist = remainingDist - dist2move
    
  } #end stop early or continue
  
  return(list(dist2move, remainingDist, length2segBase, segRatio, moveDirection)) #removed seg
}

# Run one full step of movement (all fish, one time step)
fncMoveOneStep <- function(ssn, plotit = "none") {
  
  #1. Get starting info
  corr.factor <- get("corr.factor")
  fish.all <- fish #store whole original fish df

  #Get fish density in each reach
  fish$density <- fncFishDensity(fish1.seg = fish[, "seg"], fish2.seg = NULL, corr.factor, 
                  fish1.alive.emerged.segs = seq(1:nrow(fish)), 
                  fish2.alive.emerged.segs = NULL, ssn = ssn1)[, 1]

  #Get movement distance (uses fish density and growth potential)
  moveDist <- fncMoveDistance(ssn)

  
  #2. Move each fish
  for (fpid in fish$pid) {
    
    set.seed(fpid) 

    #Get starting reach and location info
    seg <- fish$seg[fish$pid == fpid]
    segLength <- ssn@data[ssn@data$rid == seg, length.field]
    length2segBase <- fish$length2segBase[fish$pid == fpid]
    segRatio <- length2segRatio(seg, length2segBase)
    initDirection <- moveDirection <- fish$direction[fish$pid == fpid]

    # Only move if zero-inflated part of likelihood to move is 1:
    if(moveDist[fpid, "zi"] == 1){
      
      remainingDist <- moveLength <- moveDist[, "moveDist"][moveDist[, "pid"] == fpid]
      dist2move <- dist.moved.total <- 0
      
      while(remainingDist>0){ #Continue as long as the fish has not yet moved its allocated distance
        
        seg.orig <- seg #for testing
        
        #Determine if there is room within the reach to move
        Rm2Move <- fncRoom2Move(seg, moveDirection, remainingDist, length2segBase, ssn)
        if(Rm2Move == T){
          dist2move <- remainingDist
        }
        
        #Determine if fish will stop early due to good thermal habitat somewhere in the reach
        #takes precedent over room to move
        StopResult <- fncStopEarly(fpid, seg, remainingDist, moveLength, length2segBase, ssn)
        StopEarly <- StopResult[[1]] #Did the fish stop early in good habitat? T = stop, F = continue
        if(StopEarly == T) {
          dist2move <- StopResult[[3]]
          moveDirection <- StopResult[[4]]
        }
        
        #If fish is not stopping in the reach, determine which will be the next reach it enters
        #(only allowed to enter accessible reaches; allows fish to change direction at any segment edge, not just confluences)
        if (StopEarly == F & Rm2Move == F) {
          if(seg <= min(ssn@data$rid)) {seg <- min(ssn@data$rid)}
          newSeg <- fncJunctionDecisions(fpid, seg, moveDirection, stochastic = FALSE, ssn)
        } else {newSeg <- seg}
        
        #First move the fish within the current reach and update its position
        #either somewhere within the reach, or the rest of the way up/down the reach to the junction   
        Position <- fncUpdatePosition(newSeg, seg, moveDirection, length2segBase, dist2move, remainingDist, StopEarly, Rm2Move, ssn)
        dist2move <- Position[[1]]
        remainingDist <- Position[[2]]
        length2segBase <- Position[[3]]
        segRatio <- Position[[4]]
        seg <- newSeg
        moveDir.orig <- moveDirection
        moveDirection <- Position[[5]]
        dist.moved.total <- dist.moved.total + dist2move
        
        #stop if errors
        if(remainingDist<0) browser()
        if(segRatio > 1 | segRatio < 0) browser()
        
        #temp.df=rbind(temp.df, cbind(seg, moveDir.orig, Rm2Move, StopEarly, dist2move, dist.moved.total, remainingDist, newSeg, moveDirection, length2segBase, segRatio)) #testing only
        
      } #end while remaining dist loop
  
      #Put info back in fish table
      fish$seg[fish$pid == fpid] <- seg
      fish$length2segBase[fish$pid == fpid] <- length2segBase
      fish$ratio[fish$pid == fpid] <- segRatio
      fish$direction[fish$pid == fpid] <- moveDirection
      fish$upDist[fish$pid == fpid] <- distance2base4fish(fpid)
      fish$movedist[fish$pid == fpid] <- dist.moved.total
      fish[, so.field][fish$pid == fpid] <- getSO(seg)
    
    } #end zero-inflated part of likelihood to move
  } #end for each fish
  
  
  #3. 
  #Update x and y coordinates
  #if this was a generated network, use straight-line distance calcs
  #if this network was created in GIS, use different approach for moving along a potentially curvy reach
  ifelse(length(grep("network-", netnm)) > 0, locs <- ddply(fish, .(pid), function(x) data.frame(getXY(x$seg, x$ratio))), locs <- ddply(fish, .(pid), function(x) data.frame(getXY.arc(x$seg, x$ratio))))
  fish$xloc <- locs$xloc
  fish$yloc <- locs$yloc
  
  #put fish that moved back into whole original dataframe
  fish.all[fish.all$pid %in% fish$pid, ] <- fish
  
  #Put updated data into ssn
  DFpred <- getSSNdata.frame(ssn, "fish")
  DFpred$upDist <- fish.all$upDist
  DFpred$ratio <- fish.all$ratio
  DFpred[, xlat] <- fish.all$xloc
  DFpred[, ylon] <- fish.all$yloc
  DFpred$rid <- fish.all$seg
  DFpred[, so.field] <- fish.all[, so.field]
  ssn <- putSSNdata.frame(DFpred, ssn, "fish")
  ssn@predpoints@SSNPoints[[2]]@point.coords[, 1] <- DFpred[, xlat]
  ssn@predpoints@SSNPoints[[2]]@point.coords[, 2] <- DFpred[, ylon]
  
  
  if(plotit != "none"){
    if(plotit != "screen"){png(paste0(plotDir, "/", plotit, ".png"), width = 6, height = 6, units = "in", res = 300)}
    plot(ssn, lwdLineCol = "addfunccol", lwdLineEx = 8, linecol = "blue", xlab = "X", ylab = "Y", cex = 0)
    plot(ssn, PredPointsID = "fish", add = TRUE, cex = 0.7)
    if(plotit != "screen"){dev.off()}
  }
  
  return(list(fish.all, ssn))
}



# === BIOENERGETICS FUNCTIONS ==================================================================== 

#look up growth based on actual water temperature, ration, and weight of fish
  #PIDs is vector of fish or a single fish
fncGrowthFish <- function(PIDs = NA){
  
  #these need to match what was used when pre-calculating growth 
  #(they define the dimensions of the array that the precalculated data are stored in)
  wt.seq <- seq(0.05, 25, 0.05) #water temperature
  ra.seq <- seq(0.02, 0.17, 0.001) #ration
  ma.seq <- seq(1, 1500, 1) #fish mass
  
  #get data
  if(any(is.na(PIDs))){
    td <- fish[, c("pid", "weight", "ration", wt.field)] #all fish
  } else {
    td <- fish[fish$pid %in% PIDs, c("pid", "weight", "ration", wt.field)] #specific fish
  }
  td[, "weight"][td[, "weight"] < 1] <- 1 #minimum value that can be looked up is 1 g
  
  if(nrow(td) > 0){ #ensure there are data before proceeding
    growth <- watemp <- vector(length = length(PIDs)) #create empty vectors
    #for each fish, lookup pre-calculated growth
    for(x in 1:length(PIDs)){
      wt.idx <- which.min(abs(wt.seq - td[, wt.field][x])) #which water temperature is closest to fish's?
      ra.idx <- which.min(abs(ra.seq - td[, "ration"][x])) #which ration is closest to fish's?
      ma.idx <- which.min(abs(ma.seq - td[, "weight"][x])) #which weight is closest to fishs?
      growth[x] <- wt.growth[wt.idx, ra.idx, ma.idx] #use these indices to look up pre-calculated growth
      watemp[x] <- wt.seq[wt.idx]
    }
    if(any(is.na(PIDs))) PIDs <- fish[, "pid"]
    lookup = cbind("pid" = PIDs, "weight" = td[, "weight"], "ration" = td[, "ration"], "WT.actual" = td[, wt.field], "WT" = watemp, growth)
  
    return(lookup)
  } else{
    message("No fish IDs provided")
    return(NA)
  }
  
}
#result <- fncGrowthFish(4) #single fish
#result <- fncGrowthFish(3:7) #vector of fish


#look up best possible location for growing based on all current temperature/ration across stream network
fncGrowthPossible <- function(fweight = 1, ssn, segs = NA){
  
  #these need to match what was used when pre-calculating growth 
  #(they define the dimensions of the array that the precalculated data are stored in)
  wt.seq <- seq(0.05, 25, 0.05) #water temperature
  ra.seq <- seq(0.02, 0.17, 0.001) #ration
  ma.seq <- seq(1, 1500, 1) #fish mass

  #get data
  if(any(is.na(segs))){
    td <- obs.df[obs.df$rid %in% ssn1@data$rid[ssn@data$accessible == 1], c("pid", "rid", "ration", wt.field)]
  } else {
    td <- obs.df[obs.df$rid %in% segs & obs.df$rid %in% ssn@data$rid[ssn@data$accessible == 1], c("pid", "rid", "ration", wt.field)]
  }

  if(nrow(td) > 0){ #ensure there are data before proceeding
    growth <- watemp <- vector(length = nrow(td)) #create empty vectors
    # for each location in stream network at this time, lookup pre-calculated growth potential
    for(x in 1:nrow(td)){
      wt.idx <- which.min(abs(wt.seq - td[, wt.field][x]))
      ra.idx <- which.min(abs(ra.seq - td[, "ration"][x]))
      ma.idx <- round(fweight, 2)
      ma.idx[ma.idx < 1] <- 1
      growth[x] <- wt.growth[wt.idx, ra.idx, ma.idx]
      watemp[x] <- wt.seq[wt.idx]
    }
    lookup = cbind("pid" = td[, "pid"], "rid" = td[, "rid"], "weight" = rep(fweight, nrow(td)), "ration" = td[, "ration"], "WT.actual" = td[, wt.field], "WT" = watemp, growth)
  
    return(lookup)
  } else{
    message("No accessible reaches provided")
    return(NA)
  }
  
}
#result <- fncGrowthPossible(100, ssn) # whole stream network that's accessible
#result <- fncGrowthPossible(100, ssn, 25:50) # just reaches 25-50


# These functions were provided by Matt Nahorniak, South Fork Research, Inc, and adapted for our use
# (The first 3 functions were adapted for our needs)

# Parameters
fncGetBioEParms <- function(spp, pred.en.dens, prey.en.dens, oxy, pff, wt.nearest, 
      startweights = rep(initial.mass, numFish), pvals = rep(0.5, numFish), ration = rep(0.1, numFish)){
  
	N.sites <- nrow(wt.nearest) # here, sites are fish in each reach
  N.steps <- 1 # one time step
  Species <- spp
  SimMethod <- 1 # method that predicts growth
  Pred <- pred.en.dens # predator energy density
  Oxygen <- oxy # oxygen consumed
	PFF <- pff # percent indigestible prey
  stab.factor <- 0.5 # stability factor for other simulation methods
  epsilon <- 0.5 # also for other simulation methods
	endweights <- startweights*5
	TotalConsumption <- rep(100, N.sites)
	pvalues <- t(matrix(round(pvals, 5)))
	sitenames <- t(matrix(wt.nearest[, "pid"]))
	temperature <- t(wt.nearest[, "WT"])
  prey.energy.density <- t(matrix(rep(prey.en.dens, N.sites))) 
  ration <- t(matrix(round(ration, 5)))
  
  return(list(	
    "Species" = Species, 
    "SimMethod" = SimMethod, 
    "Wstart" = startweights, 
    "Endweights" = endweights, 
    "TotalConsumption" = TotalConsumption, 
    "pp" = pvalues, 
    "Temps" = temperature, 
    "N.sites" = N.sites, 
    "N.steps" = N.steps, 
    "sitenames" = sitenames, 
    "Pred" = Pred, 
    "prey.energy.density" = prey.energy.density, 
    "Oxygen" = Oxygen, 
    "stab.factor" = stab.factor, 
    "PFF" = PFF, 
    "epsilon" = epsilon, 
    "ration" = ration)
  )
  
}

# Constants, hard-coded O. mykiss juveniles
fncReadConstants.omykiss <- function() {
  
  # Get consumption constants
  Cons = data.frame(
    ConsEQ = 3, 
    CA = 0.628, 
    CB = -0.3, 
    CQ = 3.5, 
    CTO = 25, 
    CTM = 22.5, 
    CTL = 24.3, 
    CK1 = 0.2, 
    CK4 = 0.2)
  
  # Get respiration constants
  Resp = data.frame(
    RespEQ = 2, 
    RA = 0.013, 
    RB = -0.217, 
    RQ = 2.2, 
    RTO = 22, 
    RTM = 26, 
    RTL = 0, 
    RK1 = 0, 
    RK4 = 0.13, 
    ACT = 1.3, 
    BACT = 0.0405, 
    SDA = 0.172)

  # Get Excretion / Egestion Constants
  Excr = data.frame(
    ExcrEQ = 3, 
    FA = 0.212,  
    FB = -0.222, 
    FG = 0.631, 
    UA = 0.0314, 
    UB = 0.58, 
    UG = -0.299)
  
  # Return the Constants
  return(list("Consumption" = Cons, 
              "Respiration" = Resp, 
              "Excretion" = Excr))
}

# Constants, hard-coded for O. mykiss adults
fncReadConstants.steelhead <- function() {
  
  # Get consumption constants
  Cons = data.frame(
    ConsEQ = 3, 
    CA = 0.628, 
    CB = -0.3, 
    CQ = 5, 
    CTO = 20, 
    CTM = 20, 
    CTL = 24, 
    CK1 = 0.33, 
    CK4 = 0.2)
  
  # Get respiration constants
  Resp = data.frame(
    RespEQ = 1, 
    RA = 0.00264, 
    RB = -0.217, 
    RQ = 0.06818, 
    RTO = 0.0234, 
    RTM = 0, 
    RTL = 25, 
    RK1 = 1, 
    RK4 = 0.13, 
    ACT = 9.7, 
    BACT = 0.0405, 
    SDA = 0.172)
  
  # Get Excretion / Egestion Constants
  Excr = data.frame(
    ExcrEQ = 3, 
    FA = 0.212,  
    FB = -0.222, 
    FG = 0.631, 
    UA = 0.0314, 
    UB = 0.58, 
    UG = -0.299)
  
  # Return the Constants
  return(list("Consumption" = Cons, 
              "Respiration" = Resp, 
              "Excretion" = Excr))
}

# Constants, hard-coded for O. mykiss adults - using RATION instead of p-values
fncReadConstants.steelhead_ration <- function() {
  
  # Get consumption constants
  Cons = data.frame(
    ConsEQ = 4, 
    CA = 0.628, 
    CB = -0.3, 
    CQ = 5, 
    CTO = 20, 
    CTM = 20, 
    CTL = 24, 
    CK1 = 0.33, 
    CK4 = 0.2)
  
  # Get respiration constants
  Resp = data.frame(
    RespEQ = 1, 
    RA = 0.00264, 
    RB = -0.217, 
    RQ = 0.06818, 
    RTO = 0.0234, 
    RTM = 0, 
    RTL = 25, 
    RK1 = 1, 
    RK4 = 0.13, 
    ACT = 9.7, 
    BACT = 0.0405, 
    SDA = 0.172)
  
  # Get Excretion / Egestion Constants
  Excr=data.frame(
    ExcrEQ = 4, 
    FA = 0.212,  
    FB = -0.222, 
    FG = 0.631, 
    UA = 0.0314, 
    UB = 0.58, 
    UG = -0.299)
  
  # Return the Constants
  return(list("Consumption" = Cons, 
              "Respiration" = Resp, 
              "Excretion" = Excr))
}

# Consumption Equation 1
ConsumptionEQ1 <- function(W, TEMP, PP, PREY, CA, CB, CQ) {

	CMAX <- CA * (W ** CB)			#max specific feeding rate (g_prey/g_pred/d)
	CONS <- (CMAX * PP * exp(CQ * TEMP))	#specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
	CONSj <- CONS * PREY          		 #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day
	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Consumption Equation 2
ConsumptionEQ2 <- function(W, TEMP, PP, PREY, CA, CB, CTM, CTO, CQ) {

	Y <- log(CQ) * (CTM - CTO + 2)
	Z <- log(CQ) * (CTM - CTO)
	X <- (Z^2 * (1 + (1 + 40 / Y)^.5)^2) / 400
	V <- (CTM - TEMP) / (CTM - CTO)
	CMAX <- CA * (W ** CB)		
	CONS <- CMAX * PP * (V ** X) * exp(X * (1 - V))
	CONSj <- CONS * PREY             

	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Consumption Equation #3: Temperature Dependence for cool-cold water species
ConsumptionEQ3 <- function(W, TEMP, PP, PREY, CA, CB, CK1, CTO, CQ, CK4, CTL, CTM) {
	G1 <- (1 / (CTO - CQ)) * (log((0.98 * (1 - CK1)) / (CK1 * 0.02)))
	L1 <- exp(G1 * (TEMP - CQ))
	KA <- (CK1 * L1) / (1 + CK1 * (L1 - 1))
	G2 <- (1 / (CTL - CTM)) * (log((0.98 * (1 - CK4)) / (CK4 * 0.02)))
	L2 <- exp(G2 * (CTL - TEMP))
	KB <- (CK4 * L2) / (1 + CK4 * (L2 - 1))
	CMAX <- CA * (W ** CB)		#max specific feeding rate (g_prey/g_pred/d)
	
	CONS <- CMAX * PP * KA * KB		#specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
	
	CONSj <- CONS * PREY             #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day

	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
	}

# Consumption Equation 3 that gives CONS based on RATION instead of PP
ConsumptionEQ4 <- function(W, TEMP, RATION, PREY, CA, CB, CK1, CTO, CQ, CK4, CTL, CTM) {
  G1 <- (1 / (CTO - CQ)) * (log((0.98 * (1 - CK1)) / (CK1 * 0.02)))
  L1 <- exp(G1 * (TEMP - CQ))
  KA <- (CK1 * L1) / (1 + CK1 * (L1 - 1))
  G2 <- (1 / (CTL - CTM)) * (log((0.98 * (1 - CK4)) / (CK4 * 0.02)))
  L2 <- exp(G2 * (CTL - TEMP))
  KB <- (CK4 * L2) / (1 + CK4 * (L2 - 1))
  CMAX <- CA * (W ** CB)		#max specific feeding rate (g_prey/g_pred/d)
  
  RATION[RATION >= CMAX] <- CMAX[RATION >= CMAX]
  CONS <- RATION * KA * KB #specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
  
  CONSj <- CONS * PREY #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day
  
  return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Excretion Equation 1
ExcretionEQ1 <- function(CONS, CONSj, TEMP, PP, FA, UA) {

	EG <- FA * CONS				# egestion (fecal waste) in g_waste/g_pred/d
	U <- UA * (CONS - EG)	 			# excretion (nitrogenous waste) in g_waste/g_pred/d

	EGj <- FA * CONSj				# egestion in J/g/d
	Uj <- UA * (CONSj - EGj)			# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}	

# Excretion Equation 2
ExcretionEQ2 <- function(CONS, CONSj, TEMP, PP, FA, UA, FB, FG, UB, UG) {

	EG <- FA * TEMP ^ FB * exp(FG * PP) * CONS			# egestion (fecal waste) in g_waste/g_pred/d
	U <- UA * TEMP^UB * exp(UG * PP)*(CONS - EG)			# excretion (nitrogenous waste) in g_waste/g_pred/d

	EGj <- EG * CONSj / CONS					# egestion in J/g/d
	Uj <- U * CONSj / CONS					# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}

# Excretion Equation 3 (W/ correction for indigestible prey as per Stewart 1983)
ExcretionEQ3 <- function(CONS, CONSj, TEMP, PP, FA, UA, FB, FG, UB, UG, PFF) {

	#Note: In R, "F" means "FALSE", here we use EG as the variable name for egestion instead of F (as in the FishBioE 3.0 manual)
	#Note:  PFF = 0 assumes prey are entirely digestible, making this essentially the same as Equation 2

	PE <- FA * (TEMP ** FB) * exp(FG * PP)

	PF <- ((PE - 0.1) / 0.9) * (1 - PFF) + PFF

	EG <- PF * CONS					# egestion (fecal waste) in g_waste/g_pred/d
	U <- UA * (TEMP ** UB) * (exp(UG * PP)) * (CONS - EG)	# excretion (nitrogenous waste) in g_waste/g_pred/d

	EGj <- PF * CONSj				# egestion in J/g/d
	Uj <- UA * (TEMP ** UB) * (exp(UG * PP)) * (CONSj - EGj)	# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}	

# Excretion Equation 3 but using RATION instead of PP
ExcretionEQ4 <- function(CONS, CONSj, TEMP, RATION, CMAX, FA, UA, FB, FG, UB, UG, PFF) {
  
  #Note: In R, "F" means "FALSE", here we use EG as the variable name for egestion instead of F (as in the FishBioE 3.0 manual)
  #Note:  PFF = 0 assumes prey are entirely digestible, making this essentially the same as Equation 2
  
  RATION[RATION >= CMAX] <- CMAX[RATION >= CMAX]
  PP <- RATION / CMAX #calculating p-value based on ration and CMax inputs
  
  PE <- FA * (TEMP ** FB) * exp(FG * PP)
  
  PF <- ((PE - 0.1) / 0.9) * (1 - PFF) + PFF
  
  EG <- PF * CONS	# egestion (fecal waste) in g_waste/g_pred/d
  U <- UA*(TEMP**UB)*(exp(UG*PP))*(CONS-EG)	# excretion (nitrogenous waste) in g_waste/g_pred/d
  
  EGj <- PF * CONSj	# egestion in J/g/d
  Uj <- UA * (TEMP ** UB) * (exp(UG * PP)) * (CONSj - EGj)	# excretion in J/g/d
  
  return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
}	

# Respiration Equation 1
RespirationEQ1 <- function(W, TEMP, CONS, EG, PREY, OXYGEN, RA, RB, ACT, SDA, RQ, RTO, RK1, RK4, RTL, BACT){

	VEL <- (RK1 * W^RK4) * (TEMP > RTL) +  ACT * W^RK4 * exp(BACT * TEMP) * (1 - 1 * (TEMP > RTL))
	ACTIVITY <- exp(RTO * VEL)
	S <- SDA * (CONS - EG)		# proportion of assimilated energy lost to SDA in g/g/d (SDA is unitless)
	Sj <- S * PREY		# proportion of assimilated energy lost to SDA in J/g/d - Joules lost to digestion per gram of predator mass per day
	R <- RA * (W ** RB) * ACTIVITY * exp(RQ * TEMP) 	# energy lost to respiration (metabolism) in g/g/d
	Rj <- R * OXYGEN  # energy lost to respiration (metabolism) in J/g/d - Joules per gram of predator mass per day
	return(list("R" = R, "Rj" = Rj, "S" = S, "Sj" = Sj))
	}

# Respiration Equation 2 (Temp dependent w/ ACT multiplier)
RespirationEQ2 <- function(W, TEMP, CONS, EG, PREY, OXYGEN, RA, RB, ACT, SDA, RTM, RTO, RQ) {

	V <- (RTM - TEMP) / (RTM - RTO)
	V[V < 0] <- 0.001 #AHF Added to stop errors when Water temps exceed RTM!
	Z <- (log(RQ)) * (RTM - RTO)
	Y <- (log(RQ)) * (RTM - RTO + 2)
	X <- ((Z ** 2) * (1 + (1 + 40 / Y) ** 0.5) ** 2) / 400
	S <- SDA * (CONS - EG)	# proportion of assimilated energy lost to SDA in g/g/d (SDA is unitless)
	Sj <- S * PREY	# proportion of assimilated energy lost to SDA in J/g/d - Joules lost to digestion per gram of predator mass per day
	R <- RA * (W ** RB) * ACT * ((V ** X) * (exp(X * (1 - V) ))) 	# energy lost to respiration (metabolism) in g/g/d
	Rj <- R * OXYGEN  # energy lost to respiration (metabolism) in J/g/d - Joules per gram of predator mass per day
	return(list("R" = R, "Rj" = Rj, "S" = S, "Sj" = Sj))
	}

# Calculate Growth
CalculateGrowth <- function(Constants, Input) {
  
  pred <- Input$Pred
  Oxygen <- Input$Oxygen
  
  # Initialize Fish Weights
  W <-array(rep(0, (Input$N.sites * (Input$N.steps + 1))), c(Input$N.steps + 1, Input$N.sites))
  W[1, ] <- as.numeric(Input$Wstart[1:ncol(W)])
  Growth <- array(rep(0, Input$N.sites*Input$N.steps), c(Input$N.steps, Input$N.sites))
  Growth_j <- Growth
  Consumpt <- Growth
  Consumpt_j <- Growth
  Excret <- Growth
  Excret_j <- Growth
  Egest <- Growth
  Egest_j <- Growth   
  Respirat <- Growth 
  Respirat_j <- Growth 
  S.resp <- Growth 
  Sj.resp <- Growth 
  Gg_WinBioE <- Growth
  Gg_ELR <- Growth
  TotalC <- rep(0, Input$N.sites)
  
  
  ##Start Looping Through Time - for Known Consumption, solving for Weight
  
  t <- 1
  for (t in 1:(Input$N.steps)) {
    
    ### Consumption 
    if(Constants$Consumption$ConsEQ == 1) {
      Cons <- with(Constants$Consumption, ConsumptionEQ1(W[t, ], Input$Temps[t, ], Input$pp[t, ], Input$prey.energy.density[t, ], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CQ))
    } else if(Constants$Consumption$ConsEQ == 2){
      Cons <- with(Constants$Consumption, ConsumptionEQ2(W[t, ], Input$Temps[t, ], Input$pp[t, ], Input$prey.energy.density[t, ], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CTM, Constants$Consumption$CTO, Constants$Consumption$CQ))
    } else if(Constants$Consumption$ConsEQ == 3){
      Cons <- with(Constants$Consumption, ConsumptionEQ3(W[t, ], Input$Temps[t, ], Input$pp[t, ], Input$prey.energy.density[t, ], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CK1, Constants$Consumption$CTO, Constants$Consumption$CQ, Constants$Consumption$CK4, Constants$Consumption$CTL, Constants$Consumption$CTM))
    } else if(Constants$Consumption$ConsEQ == 4){
      Cons <- with(Constants$Consumption, ConsumptionEQ4(W[t, ], Input$Temps[t, ], Input$ration[t, ], Input$prey.energy.density[t, ], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CK1, Constants$Consumption$CTO, Constants$Consumption$CQ, Constants$Consumption$CK4, Constants$Consumption$CTL, Constants$Consumption$CTM))
    }
    
    TotalC <- TotalC + Cons$CONS * W[t, ]
    
    # store daily consumption 
    Consumpt[t, ] <- as.numeric(Cons$CONS)
    Consumpt_j[t, ] <- as.numeric(Cons$CONSj)
    
    
    ### Excretion / Egestion
    if(Constants$Excretion$ExcrEQ == 1) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ1(Cons$CONS, Cons$CONSj, Input$Temps[t, ], Input$pp[t, ], Constants$Excretion$FA, Constants$Excretion$UA))
    } else if (Constants$Excretion$ExcrEQ == 2) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ2(Cons$CONS, Cons$CONSj, Input$Temps[t, ], Input$pp[t, ], Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG ))
    } else if (Constants$Excretion$ExcrEQ == 3) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ3(Cons$CONS, Cons$CONSj, Input$Temps[t, ], Input$pp[t, ], Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG, Input$PFF) )
    } else if (Constants$Excretion$ExcrEQ == 4) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ4(Cons$CONS, Cons$CONSj, Input$Temps[t, ], Input$ration[t, ], Cons$CMAX, Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG, Input$PFF) )
    }
    
    # store daily excretion and egestion
    Excret[t, ] <- as.numeric(ExcEgest$U)
    Excret_j[t, ] <- as.numeric(ExcEgest$Uj)
    Egest[t, ] <- as.numeric(ExcEgest$EG)
    Egest_j[t, ] <- as.numeric(ExcEgest$EGj)
    
    
    ### Respiration
    if(Constants$Respiration$RespEQ == 1) {
      Resp <- with(Constants$Respiration, RespirationEQ1(W[t, ], Input$Temps[t, ], Cons$CONS, ExcEgest$EG, Input$prey.energy.density[t, ], Input$Oxygen, Constants$Respiration$RA, Constants$Respiration$RB, Constants$Respiration$ACT, Constants$Respiration$SDA, Constants$Respiration$RQ, Constants$Respiration$RTO, Constants$Respiration$RK1, Constants$Respiration$RK4, Constants$Respiration$RTL, Constants$Respiration$BACT))
    } else if (Constants$Respiration$RespEQ == 2) {
      Resp <- with(Constants$Respiration, RespirationEQ2(W[t, ], Input$Temps[t, ], Cons$CONS, ExcEgest$EG, Input$prey.energy.density[t, ], Input$Oxygen, Constants$Respiration$RA, Constants$Respiration$RB, Constants$Respiration$ACT, Constants$Respiration$SDA, Constants$Respiration$RTM, Constants$Respiration$RTO, Constants$Respiration$RQ))
    }
    
    #store daily respiration results
    Respirat[t, ] <- as.numeric(Resp$R)
    Respirat_j[t, ] <- as.numeric(Resp$Rj)
    S.resp[t, ] <- as.numeric(Resp$S)
    Sj.resp[t, ] <- as.numeric(Resp$Sj)
    
    
    ### Now calculate Growth
    
    # growth in J/g/d - Joules allocated to growth for each gram of predator on each day
    Gj <- Cons$CONSj - Resp$Rj - ExcEgest$EGj - ExcEgest$Uj - Resp$Sj	
    G <- Cons$CONS - Resp$R - ExcEgest$EG - ExcEgest$U - Resp$S
    # growth in g/d - Grams of predator growth each day
    Growth[t, ] <- as.numeric(Gj * W[t, ]) / pred
    Growth_j[t, ] <- as.numeric(Gj)
    
    # growth in g/g/d (DailyWeightIncrement divided by fishWeight)
    Gg_WinBioE[t, ] <- as.numeric(Growth[t, ] / W[t, ])			
    
    # growth in g/g/d (DailyWeightIncrement divided by average of fish start end weights)
    Gg_ELR[t, ] <- Growth[t, ] / ((as.numeric(Input$Wstart[1:ncol(W)]) + W[t, ]) / 2)		
    
    # Calculate absolute weight at time t+1
    W[t + 1, ] <- W[t, ] + Growth[t, ]
    
    
  } # End of cycles through time
  
  # Return Results	
  return(list(
    "TotalC" = TotalC, 
    "W" = W, 
    "Growth" = Growth, 
    "Gg_WinBioE" = Gg_WinBioE, 
    "Gg_ELR" = Gg_ELR, 
    "Growth_j" = Growth_j, 
    "Consumption" = Consumpt, 
    "Consumption_j" = Consumpt_j, 
    "Excretion" = Excret, 
    "Excretion_j" = Excret_j, 
    "Egestion" = Egest, 
    "Egestion_j" = Egest_j, 
    "Respiration" = Respirat, 
    "Respiration_j" =- Respirat_j, 
    "S.resp" = S.resp, 
    "Sj.resp" = Sj.resp
  ))
}

# The main bioenergetics function that calls other functions
BioE <- function (Input, Constants) {

# for simulation method =1 (we have p-vals, and want to solve for weights
if (Input$SimMethod == 1) {
  
	# Method 1: Calculate Growth from p-values and Temperatures
	Results <- CalculateGrowth(Constants, Input)
	W <- Results$W
	Growth <- Results$Growth
	
} else{

# We don't know p-values, but need to iteratively solve for them
# Method 2 or 3:  Calculte P-values from Total Growth or Consumption
# Need to assume p-values are constant with time
	
	# initialize first guess p-values of .5
	Input$pp <- array(rep(0.5, Input$N.sites * Input$N.step), 
		 c(Input$N.step, Input$N.sites))

	# set error at high value, iterate until it's small
	Error <- rep(99, Input$N.sites)
	iteration <- 0

### Interate until error is less than .1
	while(max(abs(Error)) > Input$epsilon) 
{
	iteration <- iteration + 1
		Results <- CalculateGrowth(Constants, Input)
		W <- Results$W
    TConsumption <- Results$TotalC

 # Find error (depending on which thing on which we're converging), and
 # and come up with new estimate for average p-value
		if (Input$SimMethod == 2)  {
		      Error = (W[Input$N.step+1, ]-Input$Endweights)
  # Delta is the amount by which we'll change the p-value (prior to
  # scaling by the stability factor
		  	Delta= (Input$Endweights) / (W[Input$N.step+1, ]) *
      	         Input$pp[1, ] - Input$pp[1, ]

		Pnew = as.vector(Input$p[1, ] + Input$stab.factor*Delta)
	# Guard against negative p-values (maybe I shouldn't for convergence' sake)
	 for (i in 1:length(Pnew)) {Pnew[i]=max(0, Pnew[i])}
		} else 
		{
	      Error= (TConsumption-Input$TotalConsumption)/Input$TotalConsumption
	  	Delta= Input$TotalConsumption/TConsumption * Input$p[1, ] - Input$p[1, ]
		Pnew = as.vector(Input$p[1, ] + Input$stab.factor*  Delta)
	}
	# Update for user, show p-values and error, see if we're converging
		for (i in 1:Input$N.step) {Input$p[i, ]=as.numeric(Pnew)}
			print(paste("Pnew=", Pnew, "  Error=", Error))
	            print(paste("iteration =", iteration))	
		} # end of while statement 
	} 

Results$pp <- Input$pp
# We're done.  Let's return the results!
return(Results)
}

# Energetic cost to movement
fncMoveCost <- function(md = fish$movedist, gr = fish$growth){
  return(abs(gr) * md * move.cost)
}

# === PLOTTING FUNCTIONS ==========================================================================

#plot.SpatialStreamNetwork modified to allow a subset of the predpoints to be passed in
plotSSN.mod <- function (x, VariableName = NULL, color.palette = NULL, nclasses = NULL, 
    breaktype = "quantile", brks = NULL, PredPointsID = NULL, ObsPointsID = NULL, 
    add = FALSE, addWithLegend = FALSE, lwdLineCol = NULL, lwdLineEx = 1, lineCol = "black", 
    myvar = NULL, myvar2 = NULL, pch2, ...) {
  if (missing(lwdLineEx)) 
    lwdLineEx <- 1
  if (missing(lwdLineCol)) {
    x@data$lineWidth <- rep(1, nrow(x@data))
    lwdLineCol <- "lineWidth"
  }
  if (is.null(as.list(match.call()[-1])$pch)) {
    plch = 19
  } else plch <- as.list(match.call()[-1])$pch
  if (is.null(as.list(match.call()[-1])$cex)) {
    chex = 1
  } else chex <- as.list(match.call()[-1])$cex
  if (is.null(as.list(match.call()[-1])$col)) {
    colr = "black"
  } else colr <- as.list(match.call()[-1])$col
  par.orig <- par(no.readonly = TRUE)
  if (!is.null(PredPointsID)) {
    for (i in 1:length(x@predpoints@ID)) {
      if (x@predpoints@ID[i]  ==  PredPointsID) {
        if (add  ==  FALSE & addWithLegend  ==  FALSE) {
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
          for (j in 1:length(x@lines)) for (k in 1:length(x@lines[[j]])) if (is.null(lwdLineCol)) 
            lines((x@lines[[j]]@Lines[[k]]@coords), col = lineCol, ...)
          else lines(x@lines[[j]]@Lines[[k]]@coords, lwd = lwdLineEx * x@data[i, lwdLineCol], col = lineCol, ...)
        }
        if (add  ==  TRUE) {
          par(new = TRUE)
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", bty = "n", xlab = "", ylab = "", ...)
        }
        if (addWithLegend  ==  TRUE) {
          par(new = TRUE)
          layout(matrix(1:2, nrow = 1), widths = c(4, 1))
          par(mar = c(5, 5, 3, 0))
          par(mfg = c(1, 1))
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", bty = "n", xlab = "", ylab = "", ...)
        }
        if(!is.null(myvar)){
          if(!is.null(pch2)) {points(x@predpoints@SSNPoints[[2]]@point.coords[myvar, 1], x@predpoints@SSNPoints[[2]]@point.coords[myvar, 2], pch = pch2, cex = chex, col = colr)
            } else {points(x@predpoints@SSNPoints[[2]]@point.coords[myvar, 1], x@predpoints@SSNPoints[[2]]@point.coords[myvar, 2], pch = plch, cex = chex, col = colr)}
        } else{
          points(x@predpoints@SSNPoints[[i]]@point.coords, pch = plch, cex = chex, col = colr)
        }
      }
    }
    par(par.orig)
  par.orig <- par(no.readonly = TRUE)
  }
  else if (!is.null(ObsPointsID)) {
    for (i in 1:length(x@obspoints@ID)) {
      if (x@obspoints@ID[i]  ==  ObsPointsID) {
        if (add  ==  FALSE & addWithLegend  ==  FALSE) {
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
          for (j in 1:length(x@lines)) for (k in 1:length(x@lines[[j]])) if (is.null(lwdLineCol)) 
            lines((x@lines[[j]]@Lines[[k]]@coords), col = lineCol, ...)
          else lines(x@lines[[j]]@Lines[[k]]@coords, lwd = lwdLineEx * x@data[i, lwdLineCol], col = lineCol, ...)
        }
        if (add  ==  TRUE) {
          par(new = TRUE)
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", bty = "n", xlab = "", ylab = "", ...)
        }
        if (addWithLegend  ==  TRUE) {
          par(new = TRUE)
          layout(matrix(1:2, nrow = 1), widths = c(4, 1))
          par(mar = c(5, 5, 3, 0))
          par(mfg = c(1, 1))
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", bty = "n", xlab = "", ylab = "", ...)
        }
        if(!is.null(myvar2)){
          points(x@obspoints@SSNPoints[[1]]@point.coords[myvar2, 1], x@obspoints@SSNPoints[[1]]@point.coords[myvar2, 2], pch = plch, cex = chex, col = colr)
        } else{
          points(x@obspoints@SSNPoints[[i]]@point.coords, pch = plch, cex = chex, col = colr)
        }
      }
    }
    par(par.orig)
  }
  else if (is.null(VariableName)) {
    plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
    for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
      lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol, ...)
    else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * x@data[i, lwdLineCol], col = lineCol, ...)
    points(x@obspoints@SSNPoints[[1]]@point.coords, pch = plch, cex = chex, col = colr)
    par(par.orig)
  }
  else {
    layout(matrix(1:2, nrow = 1), widths = c(4, 1))
    par(mar = c(5, 5, 3, 0))
    plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
    for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
      lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol, ...)
    else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * x@data[i, lwdLineCol], col = lineCol, ...)
    data <- x@obspoints@SSNPoints[[1]]@point.data
    if (is.null(nclasses)) 
      nclasses <- 10
    lower.breaks <- matrix(0, nrow = nclasses, ncol = 1)
    upper.breaks <- matrix(0, nrow = nclasses, ncol = 1)
    if (breaktype  ==  "quantile") {
      brks <- quantile(data[, VariableName], probs = (1:(nclasses - 1))/nclasses, na.rm = T)
      lower.breaks <- c(min(data[, VariableName], na.rm = T), brks)
      upper.breaks <- c(brks, max(data[, VariableName], na.rm = T))
    }
    if (breaktype  ==  "even") {
      brks <- min(data[, VariableName]) + (max(data[, VariableName]) - min(data[, VariableName])) * (1:(nclasses - 1))/nclasses
      lower.breaks <- c(min(data[, VariableName], na.rm = T), brks)
      upper.breaks <- c(brks, max(data[, VariableName], na.rm = T))
    }
    if (breaktype  ==  "user") {
      if (is.null(brks)) 
        return("Must specify brks if breaktype = user")
      minD <- min(data[, VariableName], na.rm = TRUE)
      maxD <- max(data[, VariableName], na.rm = TRUE)
      brks <- as.vector(unlist(brks))
      if (minD < min(brks)) 
        brks <- c(brks, minD)
      if (maxD > max(brks)) 
        brks <- c(brks, maxD)
      brks <- sort(unique(unlist(brks)))
      nclasses <- length(brks) - 1
      lower.breaks <- brks[1:nclasses]
      upper.breaks <- brks[2:(nclasses + 1)]
    }
    if (length(color.palette)  ==  0) 
      color.palette <- rainbow(nclasses, start = 0.66, end = 0.99)
    for (j in 1:nclasses) {
      jmax <- upper.breaks[j]
      jmin <- lower.breaks[j]
      indj <- data[, VariableName] >= jmin & data[, VariableName] <= jmax
      points(x@obspoints@SSNPoints[[1]]@point.coords[indj, , drop = F], col = color.palette[j], pch = plch, cex = chex)
    }
    dec.dig <- 2
    left <- as.character(as.numeric(as.integer(lower.breaks * 10^dec.dig))/10^dec.dig)
    rght <- as.character(as.numeric(as.integer(upper.breaks * 10^dec.dig))/10^dec.dig)
    leglabs <- paste(left, "to", rght)
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    legend(x = -1, y = 1.1, legend = leglabs, bty = "n", 
           pch = rep(plch, times = length(leglabs)), col = color.palette, cex = 0.8)
    par(par.orig)
    return(invisible(data.frame(lower.breaks = lower.breaks, upper.breaks = upper.breaks)))
  }
}

#Hack to allow up to 14 colors in 'Spectral' palette instead of 11
fncColBrewPlus <- function(n, paint = F){
  cb <- switch(n -2, 
               rgb(c(252, 255, 153), c(141, 255, 213), c(89, 191, 148), maxColorValue = 255), 
               rgb(c(215, 253, 171, 43), c(25, 174, 221, 131), c(28, 97, 164, 186), maxColorValue = 255), 
               rgb(c(215, 253, 255, 171, 43), c(25, 174, 255, 221, 131), c(28, 97, 191, 164, 186), maxColorValue = 255), 
               rgb(c(213, 252, 254, 230, 153, 50), c(62, 141, 224, 245, 213, 136), c(79, 89, 139, 152, 148, 189), maxColorValue = 255), 
               rgb(c(213, 252, 254, 255, 230, 153, 50), c(62, 141, 224, 255, 245, 213, 136), c(79, 89, 139, 191, 152, 148, 189), maxColorValue = 255), 
               rgb(c(213, 244, 253, 254, 230, 171, 102, 50), c(62, 109, 174, 224, 245, 221, 194, 136), c(79, 67, 97, 139, 152, 164, 165, 189), maxColorValue = 255), 
               rgb(c(213, 244, 253, 254, 255, 230, 171, 102, 50), c(62, 109, 174, 224, 255, 245, 221, 194, 136), c(79, 67, 97, 139, 191, 152, 164, 165, 189), maxColorValue = 255), 
               rgb(c(158, 213, 244, 253, 254, 230, 171, 102, 50, 94), c(1, 62, 109, 174, 224, 245, 221, 194, 136, 79), c(66, 79, 67, 97, 139, 152, 164, 165, 189, 162), maxColorValue = 255), 
               rgb(c(158, 213, 244, 253, 254, 255, 230, 171, 102, 50, 94), c(1, 62, 109, 174, 224, 255, 245, 221, 194, 136, 79), c(66, 79, 67, 97, 139, 191, 152, 164, 165, 189, 162), maxColorValue = 255), 
               rgb(c(158, 213, 244, 253, 252, 252, 222, 173, 118, 102, 50, 94), c(1, 62, 109, 174, 223, 247, 245, 247, 227, 194, 136, 79), c(66, 79, 67, 97, 91, 91, 100, 111, 131, 165, 189, 162), maxColorValue = 255), 
               rgb(c(158, 213, 244, 253, 252, 252, 222, 173, 118, 102, 50, 94, 72), c(1, 62, 109, 174, 223, 247, 245, 247, 227, 194, 136, 79, 22), c(66, 79, 67, 97, 91, 91, 100, 111, 131, 165, 189, 162, 138), maxColorValue = 255), 
               rgb(c(158, 213, 244, 253, 252, 252, 222, 173, 118, 102, 55, 50, 94, 72), c(1, 62, 109, 174, 223, 247, 245, 247, 227, 194, 172, 136, 79, 22), c(66, 79, 67, 97, 91, 91, 100, 111, 131, 165, 204, 189, 162, 138), maxColorValue = 255)
  )
  cb <- cb[n:1] 
  if(paint == T) image(1:n, 1, as.matrix(1:n), col = cb, axes=F, ylab = "", xlab = "")
  return(cb)
}

#Produces a multi-panel plot of results for one simulation
fncPlotDiagnostics <- function(fa = fa.tmp, filename, plotit = "file"){
  # SIMPLIFIED FOR SWH PROJECT!
  
  if(JulianDate.Begin == 1) x.labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan")
  if(JulianDate.Begin == 274) x.labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
  
  #Dimensions of arrays
  nfish <- dim(fa)[1];nvars <- dim(fa)[2];ntime <- dim(fa)[3]
  
  
  #Quantiles across the network through time
  if(plotit == "file") png(filename, width = 5, height = 7, units = "in", res = 300)
  par(mfcol = c(4, 2), mar=c(3, 4, 1, 1.5)+0.1, las=1, oma=c(0, 0.4, 1, 0))
  
  #Cycle through variables to plot
  for(var in c("WT", "shreve", "movedist", "density", "ration", "consInst", "growth", "weight")){
    
    dat <- fa[, var, ]; dat <- apply(dat, 2, quantile, na.rm=T)
    dat <- t(dat)
    
    #day or night
    d <- seq(2, 730, 2); n <- seq(1, 730, 2)
    dat.d <- dat[d, ]
    dat.n <- dat[n, ]
    if(var %in% c("consInst", "growth", "movedist")) {dat.d <- dat.d+dat.n} #need to account for day and night values
    if(var == "WT") {dat.d <- (dat.d+dat.n)/2}
    
    if(var == "WT"){ylm <- c(0, 30); ylb <- expression("Stream temperature "~(degree ~ C))}
    if(var == "shreve"){ylm <- c(0, 40); ylb <- "Shreve stream order"}
    if(var == "movedist"){ylm <- c(0, 10); ylb <- "Movement distance (km)"}
    if(var == "density"){ylm <- c(0, 1.25); ylb <- expression("Conspecific density per"~m^2)}
    if(var == "ration"){ylm <- c(0.01, 0.08); ylb <- "Ration"}
    if(var == "consInst"){ylm <- c(0, 0.08); ylb <- "Consumption (g/g*d)"}
    if(var == "growth"){ylm <- c(-0.03, 0.02); ylb <- "Growth (g/g*d)"}
    if(var == "weight"){ylm <- c(25, 1000); ylb <- "Final weight (g)"}

      plot(dat.d[, 3], ylim = ylm, col = 1, ylab = ylb, type = 'n', xaxt='n', xlab = "", xlim = c(0, 365))
      axis(1, at=seq(1, 365, length.out=13), labels=x.labels, cex.axis=0.8)
      
      #Min/Max
      polygon(c(1:365, rev(1:365)), c(dat.d[, 1], rev(dat.d[, 5])), border=NA, col = rgb(226, 226, 226, 150, NULL, 255)) #"#e2e2e2"
      
      #Q1/Q3
      polygon(c(1:365, rev(1:365)), c(dat.d[, 2], rev(dat.d[, 4])), border=NA, col = rgb(142, 142, 142, 200, NULL, 255)) #"#8e8e8e"
      
      #Median
      lines(dat.d[, 3], lwd = 2, col = rgb(5, 5, 5, 255, NULL, 255)) #black

  }
  if(plotit == "file")  dev.off()
}

# Function to set color palette
fncColorRamp <- function(maxT = 30, minT = 5, color.length = 14, palette = "heat"){
  if(palette == "gray"){cb <- brewer.pal(color.length, "Greys")} #max of 11
  if(palette == "difference"){cb <- brewer.pal(color.length, "RdBu")} #max of 11
  if(palette == "heat"){cb <- fncColBrewPlus(n = color.length, paint = F)} #max of 14
  if(palette!="heat") cb <- cb[color.length:1] #reverse order
  #the following set the bins for plotting and legend
  colseq <- seq(from = minT, to = maxT, length.out = (color.length + 1))
  left <- colseq[1:color.length]; rght <- colseq[2:(color.length + 1)]
  return(list(cb, left, rght))
}

# === MORTALITY ===================================================================================
  #fish die based on their size (cumulative experience) and instantaneous growth rate (measure of current conditions)
    #smaller fish with lower growth rates more likely to die (sampled probabilistically)

fncSurvive <- function(df, minprob = 0.997){
  #df: data frame of fish table with only the survivors, e.g., fish[fish$survive == 1, c("weight", "growth")]
  #minprob: smallest probability any fish can have of dying in any 12-hr time step

  #Weights
  w <- df$weight #weight of fish that are alive at this time step

  #Function: 
  b <- 1 #beta: b=1 has no effect, but could change the shape of the curve with different values if desired
  v <- minprob + ( 1 -minprob) * (1 - 1 / exp(b * w))
  
  #Growth during this time step that reflects recent conditions (i.e., a hungry/stressed fish may behave in ways that make it more vulnerable to predation, etc.)
  g <- df$growth 
  g <- fncRescale(g, to = c(-0.001, 0.001))
  
  #Amount of food eaten in previous time step (better survival if bellies full b/c hunkered down somewhere)
  f <- df$pvals
  f <- fncRescale(f, to = c(-0.001, 0.001))
  
  #Probability of survival
  prb.srv <- v + g + f
  prb.srv[prb.srv > 1] <- 1 #set upper bound at 1
  
  #Sample from binomial distribution with probabilities of prb.srv to determine which fish survive this time step
  survivors <- rbinom(n = nrow(df), size = 1, prob = prb.srv)
  
  return(survivors)
}

# === END OF SCRIPT ===============================================================================
