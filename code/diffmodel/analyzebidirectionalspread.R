
#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','measure','injection.site')))

basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste0(params$opdir,'processed/',measure,'data.RData'))  # load path data and ROI names

# load structural networks
W <- read.csv(paste(params$opdir,'processed/W.csv', sep=''), header = FALSE)
W <- as.matrix(W)

# get mean pathology for each time point
tps=params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$mpi == tp, path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

# compute graph Laplacian
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degree Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degree Laplacian for anterograde connectivity only (not weighted by Snca)

# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(X) log(X,base=10))
Xo <- get.Xo(region.names,injection.site) # generates initial "pathology seed" vector, which is all 0's except a 1 in all regions listed in character vector ROI

# section below carries out these functions:
# 1. Generate predicted pathology values using the equation $x(t) = e^{-Lct}x_o$, 
# where $t$ corresponds to the relevant time points included in the experiment (0.1, 0.2, 0.3, 0.5, 1, 3, 6 and 9 months post-injection).
# 2. Compare measured pathology to predicted pathology. Specifically, we compute the Pearson correlation between actual and predicted pathology for each time point, then average those correlation coefficients with equal weighting across all time points. We then select the time constant $c$ that maximizes this average correlation.
c.output.retro <- c.fit(log.path,tps,L.out.retro,Xo,c.rng,excl.inj = params$excl.inj)
c.Grp.retro <- c.output.retro$c.best; Xt.sweep.retro <- c.output.retro$Xt.sweep # store retro time constant and r values

c.output.antero <- c.fit(log.path,tps,L.out.antero,Xo,c.rng,excl.inj = params$excl.inj)
c.Grp.antero <- c.output.antero$c.best; Xt.sweep.antero <- c.output.antero$Xt.sweep # store retro time constant and r values

scipy.sparse.linalg <- reticulate::import('scipy.sparse.linalg')
Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.sparse.linalg$expm)), base = 10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=scipy.sparse.linalg$expm)), base = 10))) # predict pathology using connectivity, time constant, and seed
df <- lapply(1:length(tps), function(t) data.frame(path = log.path[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
save(df,c.Grp.antero,c.Grp.retro,Xt.sweep.retro,Xt.sweep.antero,file = paste0(savedir,measure,'BidirectionalFit_data.RData'))
save(c.Grp.antero,c.Grp.retro,file = paste0(savedir,measure,'IndependentBidirectionalFit_params.RData'))
