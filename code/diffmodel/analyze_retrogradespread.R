#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','measure','injection.site')))
opdir <- params$opdir
basedir <- params$basedir
setwd(basedir)

savedir <- paste0(params$opdir,'diffmodel/retrograde/',paste0(injection.site,collapse='-'),'/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste0(opdir,'processed/',measure,'data.RData'))  # load path data and ROI names

# get mean pathology for each time point
tps=params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$mpi == tp, path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

W <- read.csv(paste0(params$opdir,'processed/W.csv'), header = FALSE)
W <- as.matrix(W)

L.out <- get.Lout(W,rep(1,n.regions.ABA)) # compute out-degree Laplacian for connectivity only

# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(X) log(X,base=10))
Xo <- get.Xo(region.names,injection.site) # seed pathology in injection sites
c.output <- c.fit(log.path,tps,L.out,Xo,c.rng,excl.inj = params$excl.inj)
c.Grp <- c.output$c.best; Xt.sweep <- c.output$Xt.sweep # store time constant and r values

Xt.Grp <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out,Xo,c.Grp,t)), base = 10))) # predict pathology using connectivity, time constant, and seed
df <- lapply(1:length(tps), function(t) data.frame(path = log.path[[t]], pred = Xt.Grp[,t,drop=FALSE]))
save(df,c.Grp,Xt.sweep,file = paste0(savedir,measure,'RetroFit_data.RData'))
