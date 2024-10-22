#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','measure','injection.site')))

basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/') #
opdir <- paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-')) #,'/'
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste0(params$opdir,'processed/',measure,'data.RData'))  # load path data and ROI names

# get mean pathology for each time point
tps=params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$mpi == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))
log.path <- lapply(Grp.mean, function(X) log(X,base=10))

# load structural networks
W <- read.csv(paste(params$opdir,'processed/W.csv', sep=''), header = FALSE) # subbed this and the below line for the above line
W <- as.matrix(W)

L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

Xo <- get.Xo(region.names,injection.site) # seed pathology
scipy.sparse.linalg <- reticulate::import('scipy.sparse.linalg') # import scipy matrix exponential function because it's faster
load(file=paste0(opdir,'_independentfit/',measure,'IndependentBidirectionalFit_params.RData')) #paste0(injection.site,collapse='-'),
params.opt <- c(c.Grp.retro,c.Grp.antero)
ctrl <- list(fnscale=-1) # minimize objective function

params.opt.fit <- optim(params.opt,c.objective,control = ctrl, method = 'L-BFGS-B',lower=c(10e-7,10e-7), # optimization. c's must be > 0 # shit 10e-7 used previously
                        log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                        Xo=Xo,fxn=scipy.sparse.linalg$expm,one.lm=TRUE,excl.inj=NULL) # static inputs

# extract time constants
c.Grp.retro <- params.opt.fit$par[1]
c.Grp.antero <- params.opt.fit$par[2]

# save data
Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.sparse.linalg$expm)),base=10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=scipy.sparse.linalg$expm)),base=10))) # predict pathology using connectivity, time constant, and seed
lm.mask.results <- lm.mask.ant.ret.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it
m <- lm.mask.results$m; e <- lm.mask.results$e; m.fits <- lm.mask.results$e.tp; df <- lm.mask.results$df
#list[m,e,m.fits,df] <- lm.mask.ant.ret.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it
save(df,c.Grp.antero,c.Grp.retro,m,m.fits,file = paste0(savedir,measure,'BidirectionalOptimOneLM_data.RData'))
