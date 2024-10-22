# This script assesses specificity of model to the actual seed site

rm(list=setdiff(ls(),c('params','measure','injection.site')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'nullmodels/seedspec_multi/',paste0(injection.site,collapse='-'),'/')
opdir <- params$opdir
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')

#################
### Load data ###
#################

load(paste0(opdir,'processed/',measure,'data.RData')) # load path data and ROI names

tps=params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$mpi == tp,path.names]) #path.data$Genotype == grp & path.data$Treatment == treatment & 
log.path <- lapply(Mice,function(X) log(colMeans(X,na.rm = T),base=10))

nperms <- 500 # test 500 combinations of different seed sites

non.seed.sites <- region.names[which(!region.names %in% injection.site) ] # get names of regions that are not the tested injection sites

alt.seed.sites <- list()
for(P in 1:ceiling(nperms*1.1)){ # do 10% more in case some iterations fail
  # randomly select a region
  alt.seeds <- sample(non.seed.sites,1)
  alt.seed.sites[[P]] <- alt.seeds
}

seed.fits <- matrix(NA,nrow=nperms,ncol=length(tps),dimnames=list(NULL,paste(tps,'MPI'))) # preallocate vector to store fits from each seed region
c.retro.altseed <- matrix(NA,nrow=nperms,dimnames=list(NULL,NULL))
c.antero.altseed <- matrix(NA,nrow=nperms,dimnames=list(NULL,NULL))

# prepare connectivity matrices to fit alternative seed sites with bidirectional model
W <- read.csv(paste0(params$opdir,'processed/W.csv'), header = FALSE)
W <- as.matrix(W)
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degree Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degree Laplacian for anterograde connectivity only (not weighted by Snca)

scipy.sparse.linalg <- reticulate::import('scipy.sparse.linalg') # import scipy matrix exponential function because it's faster
params.opt <- c(0.0240802676464883,0.00066889420735786) # c.retro, c.antero: use values from independent fit to initialize
#params.opt <- c(0.01,0.01) # c.retro, c.antero: use values from independent fit to initialize
ctrl <- list(fnscale=-1) # maximize the objective function instead of minimizing (default)

for(S in 1:nperms){
  seed.names <- alt.seed.sites[[S]]
  print(paste0('Seed: ',paste0(seed.names,collapse='-'),', ',S,' out of ',nperms))
  Xo <- get.Xo(region.names,seed.names) # seed pathology in each set of seed regions
  params.opt.fit <- optim(params.opt,c.objective,control = ctrl, method = 'L-BFGS-B',lower=c(10e-7,10e-7), # optimization. c's must be > 0
                          log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                          Xo=Xo,fxn=scipy.sparse.linalg$expm,one.lm=TRUE,excl.inj=NULL) # static inputs # ABA.to.CNDR.key=ABA.to.CNDR.key
  c.retro.altseed[S] <- params.opt.fit$par[1]
  c.antero.altseed[S] <- params.opt.fit$par[2]
  
  Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.retro,Xo,c.retro.altseed[S],t,fxn=scipy.sparse.linalg$expm)),base=10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.antero,Xo,c.antero.altseed[S],t,fxn=scipy.sparse.linalg$expm)),base=10))) # predict pathology using connectivity, time constant, and seed
  lm.mask.results <- lm.mask.ant.ret.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it
  m <- lm.mask.results$m; e <- lm.mask.results$e; m.fits <- lm.mask.results$e.tp; df <- lm.mask.results$df
  seed.fits[S,] <- m.fits
}

save(seed.fits,c.retro.altseed,c.antero.altseed,alt.seed.sites, file = paste0(savedir,measure,'_AlternateSeedFits_OneLM.RData'))
