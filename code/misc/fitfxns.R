source('code/misc/miscfxns.R')

##################################
### Atlas conversion functions ###
##################################

map.v2.to.v3 <- function(df.vuln,region.names.gene,v2.to.v3.key){
  #### NEED TO ENSURE ALL REGION NAMES ARE CORRECT ####
  # Input:
  # X.ABA: named matrix of values associated with ABA regions, must have ABA names as rows, and any number of columns
  # CNDR.names: character vector of CNDR names 
  # ABA.to.CNDR.key: list whose element names are CNDR names and whose elements contain vectors of corresponding ABA names
  
  # Output:
  # X.CNDR: named matrix of values associated with CNDR regions, has CNDR names for rows
  
  X.v3 <- df.vuln
  X.v2 <- as.data.frame(matrix(NA,nrow=length(region.names.gene),ncol=ncol(X.v3))) # make new vector to hold data, default is NA
  colnames(X.v2) <- colnames(X.v3)
  rownames(X.v2) <- region.names.gene
  for(v2 in region.names.gene){
    v3.region.match <- v2.to.v3.key[[v2]] # get ABA region(s) corresponding to that CNDR region, if any
    if(length(v3.region.match) > 0){      # print out names in key that aren't in ABA connectome
      #display.missing.items(names(X.ABA),ABA.region.match,'Designation names not in ABA connectome names:')
      X.v2[v2,] <- colMeans(X.v3[rownames(X.v3) %in% v3.region.match,,drop=FALSE]) # average together ABA regions in each CNDR region, for all columns
    } else if(length(ABA.region.match) ==0){# otherwise if no match, leave as NA    
      print(paste('no match for region',v3))
    }
  }  
  return(X.v2)
}

display.missing.items <- function(ref,items,msg=''){
  # find and print the elements of items that are not in ref
  # msg: extra text to display before printing missing items
  if(length(items) > 0){
    missingitems <- !(items %in% ref) 
    if(sum(missingitems) > 0){ # if there are any missing items
      print(msg)
      print(items[missingitems])
    }
  } else {print('WARNING: items is empty')}
}


bootstrap.path.tps <- function(Mice,nboot=1000){
  # INPUTS:
  # Mice: list of dataframes that is mice x regions whose values contain pathology
  # nboot: number of bootstraps
  #
  # OUTPUTS:
  # List with nboot elements, each of which contains as many elements as Mice. These sub elements contained bootstrapped
  # log10 mean pathology for each element of Mice
  
  n.mice <- lapply(Mice,nrow)
  X.boot <- lapply(1:nboot, function(n) # for nboot bootstraps, make a list of 3 single bootstraps for each time point
    lapply(1:length(Mice), function(t) log(colMeans(Mice[[t]][sample(1:n.mice[[t]],replace=T),],na.rm = T),base=10)))
  return(X.boot)
}

###############################
### fitting diffusion model ###
###############################

get.Xo <- function(region.names,ROIs){
  # generate initial "pathology seed" vector, which is all 0's except
  # a 1 in all regions listed in character vector ROI
  # ROIs should correspond to region.names
  n.regions <- length(region.names)
  Xo <- matrix(0, nrow=n.regions,ncol=1)
  Xo[which(region.names %in% ROIs)] <- 1
  rownames(Xo) <- region.names
  return(Xo)
}

get.Lout <- function(W,S,ant.ret='retro'){
  # INPUTS:
  # W: NxN adjacency matrix, may be asymmetric (if so, connections are from row i to column j)  
  # S: vector of weights to apply to each row of W (here, synuclein expression)
  # ant.ret: character specifying anterograde or retrograde spread
  #
  # The connection matrix's default directionallity here is opposite the convention for matrix multiplication
  # so we are capturing retrograde spread along connections
  # compute out-degree laplacian for retrograde spread
  #
  # OUTPUTS:
  # return specified laplacian
  
  n.regions <- nrow(W)
  W <- W * !diag(n.regions) # get rid of diagonal
  mag.orig <- sum(W) # original sum of weights
  W <- diag(as.numeric(S)) %*% W
  #W <- W * mag.orig/sum(W) # scale new W to have the same global strength... lets you fit with similar c range
  #W <- W / (max(Re(eigen(W)$values))) # scale to max eigenvalue
  
  # Where i is row element and j is column element
  # Wij is a connection from region i to region j
  # convention is the opposite, so without transposing W
  # I am capturing "retrograde" connectivity
  
  in.deg <- colSums(W)
  out.deg <- rowSums(W)
  if(ant.ret=='retro'){L.out <- diag(x = out.deg) - W} # outdegree laplacian
  if(ant.ret=='antero'){L.out <- diag(x = in.deg) - t(W)} #indegree laplacian and transposed matrix for anterograde
  return(L.out)
}

predict.Lout <- function(L,Xo,c,t=1,fxn=expm::expm){#import('scipy.sparse.linalg')$expm
  # Generate prediction for given time points using
  # L: Laplacian of adjacency matrix
  # Xo: initial conditions
  # c: time constant
  # t: time points
  # fxn: function for matrix exponential. default is python's scipy.linalg.expm 
  # using reticulate for faster performance (not as fast as matlab or python native though)
  # but R would have you use expm::expm() which is v slow
  
  Xt <- do.call('cbind', lapply(t, function(t.i) fxn(-L*c*t.i)%*%Xo)) ## edited 6/6/24
  rownames(Xt) <- rownames(Xo)
  return(Xt)
}

c.fit <- function(log.path,tps,L.out,Xo,c.rng,excl.inj){
  # fits time constant
  # INPUTS:
  # log.path: list of vectors of log-10 transformed pathology scores  
  # for each time point specified in tps (below). Time constant c is fit to predict these
  # tps: vector of numeric time points post injection
  # L.out: out-degree graph laplacian of anatomical connectivity matrix
  # Xo: vector of initial pathology
  # c.rng: range of time constants to test
  # excl.inj.CNDR: vector of injection sites in CNDR space to exclude when computing fit
  
  # this function runs diffusion model in ABA space
  # OUTPUTS:
  # c: optimal time constant within c.rng for predicting log.path
  
  # exclusion mask... don't count regions with 0 pathology
  scipy.sparse.linalg <- reticulate::import('scipy.sparse.linalg') # import python function for matrix exponential because it's faster
  #log.path <- lapply(1:length(log.path), function(t) log.path[[t]][mask[[t]]])
  # compute fit at each time point for range of time
  Xt.sweep <- matrix(NA,nrow=length(c.rng),ncol=length(tps))
  
  if(!is.null(excl.inj)){ # exclude injection sites from fit by setting them to -Inf in log.path so they'll be excluded by cor.mask below
    for(t. in 1:length(tps)){log.path[[t.]][excl.inj] <- -Inf}
  }
  
  for(c.i in 1:length(c.rng)){
    print(paste0('c ',c.i,' out of ',length(c.rng)))
    c.val <- c.rng[c.i]
    Xt.c <- do.call('cbind',lapply(tps,function(t) predict.Lout(L.out,Xo,c.val,t,fxn=scipy.sparse.linalg$expm))) #scipy.linalg$expm  # predict path using linear diff model into matrix that is region-by-time
    Xt.sweep[c.i,] <- sapply(1:length(tps), function(t) cor.mask(Xt.c[,t],log.path[[t]])) # this replaces the two lines below
  }
  c.best <- c.rng[which.max(rowMeans(Xt.sweep, na.rm=TRUE))] # select c giving max correlation
  print(c.best)
  
  return(list(c.best=c.best,Xt.sweep=Xt.sweep))
}

cor.mask <- function(Xt,y){
  # INPUTS:
  # Xt: predicted pathology
  # y: observed log pathology
  #
  # OUTPUTS:
  # compute pearson correlation between log predicted (computed) and log observed (input)
  # excluding elements that were originally 0 and thus log(0) = -Inf
  Xt <- log(Xt,base=10)
  
  mask <- y == -Inf | Xt == -Inf | is.na(Xt) | is.na(y) # 8/22 I added the is.na(y) to exclude missing values from the data
  return(cor(y[!mask],Xt[!mask]))
}

c.ABAspace.fit <- function(log.path,L.out,Xo,c.rng,CNDR.to.ABA.key){
  # fits time constant by modeling CNDR data
  # INPUTS:
  # log.path: vector of log-10 transformed pathology scores *in CNDR space*. Time constant c is fit to predict these
  # L.out: out-degree graph laplacian of anatomical connectivity matrix
  # Xo: vector of initial pathology
  # c.rng: range of time constants to test
  
  # this function runs diffusion model in ABA space, 
  # abd computes correlation with real data in ABA annotation space to assess fit of time constant c
  
  # OUTPUTS:
  # c: optimal time constant within c.rng for predicting log.path
  
  log.path <- quiet(log(map.CNDR.to.ABA(10^(log.path),rownames(Xo),CNDR.to.ABA.key),base=10)) # undo log, convert to ABA space, retake log
  # exclusion mask... don't count regions with 0 path
  mask <- log.path != -Inf & !is.na(log.path)
  # compute fit at each time point for range of time
  Xt.sweep <- sapply(c.rng, function(c) # no time here, because in this project there's only 1 time pointo 
    cor(log.path[mask],log(predict.Lout(L.out,Xo,c),base=10)[mask]))
  c <- c.rng[which.max(Xt.sweep)] # select c giving max correlation
  print(c)
  
  return(c)
}