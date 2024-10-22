#################
### Load data ###
#################
rm(list=setdiff(ls(),c('params', 'injection.site', 'measure')))

basedir <- params$basedir
opdir <- params$opdir
setwd(basedir)
savedir <- paste(params$opdir,'modelcomparison/traintest/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/modelcomparisonfxns.R')
load(paste0(opdir,'processed/',measure,'data.RData'))  # load path data and ROI names

# load structural networks
W <- read.csv(paste(params$opdir,'processed/W.csv', sep=''), header = FALSE) # subbed this and the below line for the above line
W <- as.matrix(W)
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degree Laplacian for retrograde connectivity only
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degree Laplacian for anterograde connectivity only

# load distance matrix
load(file=paste0(params$opdir,'processed/ABAEuclideanDistanceMatrix.RData'))
unit.test(all(region.names == region.names.dmat),'Region names match','ERROR: REGION NAMES DO NOT MATCH')

region.names.mdl <- region.names.dmat # exclude missing regions from model as needed
D.mat <- D.mat^-1 # scale so closer regions transmit more path
D.mat[which(diag(nrow(D.mat))==1)] <- 0
L.out.D <- get.Lout(D.mat,rep(1,nrow(D.mat))) # compute out-degree Laplacian for connectivity only

# load path data and get training/testing sets
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$mpi == tp,path.names])
nreps <- 500
tf <- 0.5

# results
results <- train.test.logmeanpath(Mice,nreps = nreps,tf = 0.5)
log.path.train <- results$log.path.train
log.path.test <- results$log.path.test
train.idx <- results$train.idx

# set up controls for all models
Xo <- get.Xo(region.names,injection.site)

scipy.sparse.linalg <- reticulate::import('scipy.sparse.linalg') # import scipy matrix exponential function because it's faster
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',measure,'IndependentBidirectionalFit_params.RData')) 
ctrl.optim <- list(fnscale=-1) # set controls for optim function -- maximize the objective function instead of minimizing (default)
ctrl.base <- list(Retrograde=L.out.retro,Anterograde=L.out.antero,Euclidean=L.out.D,
             Xo=Xo,fxn=reticulate::import('scipy.sparse.linalg')$expm,tps=tps,
             c.r.a.init=c(c.Grp.retro,c.Grp.antero),
             ctrl.optim=ctrl.optim,
             c.rng=seq(params$c.min,params$c.max,length.out = params$c.n),
             one.lm=FALSE)

# Loop through model types, train on train sets, evaluate on train and test sets
mdl.names <- c('Euclidean','Retrograde','Anterograde','BidirectionalOneLM')
results <- list()

for(mdl.name in mdl.names){
  ctrl <- ctrl.base
  print(paste('training and testing',mdl.name))
  
  if(grepl('Euclidean',mdl.name)){
    ctrl$Xo <- get.Xo(region.names.mdl,injection.site) # exclude regions with missing coordinate data
    ctrl$excl.inj <- injection.site # exclude injection sites which end up being major outliers in this model
    ctrl$c.rng <- ctrl$c.rng + 0.1 # shift c.rng to capture different range of optimal time constants in Euclidean distance model
    if(grepl('Bidirectional',mdl.name)){ # bidirectional euclidean distance model; delete regions with missing coordinates
      ctrl$Retrograde <- get.Lout(W,rep(1,length(region.names.mdl)),ant.ret = 'retro')
      ctrl$Anterograde <- get.Lout(W,rep(1,length(region.names.mdl)),ant.ret = 'antero')
      ctrl$c.r.a.init <- c(ctrl$c.r.a.init,0.001) # initialize euclidean distance time constant at 0.001
    }
  }
  
  results[[mdl.name]] <- list(train=NULL,test=NULL)
  for(REP in 1:nreps){ 
    print(paste('REP',REP))
    X.train <- log.path.train[[REP]]
    X.test <- log.path.test[[REP]]
    m.out <- spread.model.train(X.train,mdl.name,ctrl)
    results[[mdl.name]]$train[[REP]] <- spread.model.eval(m.out,X.train)
    results[[mdl.name]]$test[[REP]] <- spread.model.eval(m.out,X.test)
  }
}

# save data
save(results,Mice,file = paste0(savedir,measure,'_ModelComparison_TrainTest.RData'))
