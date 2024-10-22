source('code/misc/miscfxns.R')

lm.mask.ant.ret <- function(y,Xt.retro,Xt.antero){
  # fit linear model: y = b0 + b.r*Xt.retro + b.a*Xt.antero + error
  # given Xt.antero and Xt.retro for a given set of time constants
  # return pearson r between y and y-hat (predicted vs actual path)
  # 
  # INPUTS:
  # y: observed log pathology
  # Xt.retro: predicted pathology from retrograde model
  # Xt.antero: predicted pathology from anterograde model
  #
  # OUTPUTS:
  # compute pearson correlation between log predicted (computed) and log observed (input)
  # excluding elements that were originally 0 and thus log(0) = -Inf
  Xt.retro <- log(Xt.retro,base=10)
  Xt.antero <- log(Xt.antero,base=10)
  #mask <- y == -Inf | Xt.retro == -Inf | Xt.antero == -Inf | is.na(Xt.retro) | is.na(Xt.antero)
  #df.m <- data.frame(y=y[!mask],x1=Xt.retro[!mask],x2=Xt.antero[!mask])
  df.m <- data.frame(y=y,x1=Xt.retro,x2=Xt.antero)
  df.m <- inf.nan.mask(df.m)
  m <- lm(y~x1+x2,data=df.m)
  return(cor(m$fitted.values,df.m$y)) # optimize R^2... more intuitive for correlation
  #return(mean(residuals(m)^2)) # return mean squared error to match with regression
}

lm.mask.ant.ret.all <- function(y,Xt.retro,Xt.antero){
  # fit linear model: y = b0 + b.r*Xt.retro + b.a*Xt.antero + error
  # given Xt.antero and Xt.retro for a given set of time constants and time points
  # return pearson r between y and y-hat (predicted vs actual path)
  # aggregates across all time points to fit one model
  #
  # INPUTS:
  # y: list of observed log pathology from 3 time points
  # Xt.retro: list of predicted pathology from retrograde model from 3 time points
  # Xt.antero: list of predicted pathology from anterograde model from 3 time points
  #
  # OUTPUTS:
  # compute pearson correlation between log predicted (computed) and log observed (input)
  # excluding elements that were originally 0 and thus log(0) = -Inf
  
  Xt.retro <- log(Xt.retro,base=10)
  Xt.antero <- log(Xt.antero,base=10)
  # dataframe subsetted by time point
  df.m <- lapply(1:ncol(Xt.retro), function(t.) data.frame(y=y[[t.]],x1=Xt.retro[,t.],x2=Xt.antero[,t.]))
  df.m <- lapply(df.m,inf.nan.mask)
  for(t. in 1:length(df.m)){rownames(df.m[[t.]]) <- paste0(rownames(df.m[[t.]]),'_',t.)} # make unique rownames for each time point
  df.m.rnames <- lapply(df.m,rownames)
  df.fit <- do.call(rbind,df.m) # dataframe that is actually used to fit model
  m <- lm(y~x1+x2,data=df.fit) # fit model on all data
  r <- residuals(m)
  fv <- m$fitted.values
  #e <- mean(residuals(m)^2) # return overall MSE
  e <- cor(fv,df.fit$y) # return overall pearson r
  #e.tp <- sapply(df.m.rnames, function(rn) mean(r[rn]^2)) # return mean squared error for each time point
  e.tp <- sapply(df.m.rnames, function(rn) cor(fv[rn],df.fit[rn,'y'])) # return pearson r for each time point
  return(list(m=m,e=e,e.tp=e.tp,df=df.m)) 
}

lm.mask.ant.ret.alt.seed <- function(y,Xt.retro,Xt.antero){
  # fit linear model: y = b0 + b.r*Xt.retro + b.a*Xt.antero + error
  # given Xt.antero and Xt.retro for a given set of time constants and time points
  # return pearson r between y and y-hat (predicted vs actual path)
  # aggregates across all time points to fit one model
  #
  # INPUTS:
  # y: list of observed log pathology from 3 time points
  # Xt.retro: list of predicted pathology from retrograde model from 3 time points
  # Xt.antero: list of predicted pathology from anterograde model from 3 time points
  #
  # OUTPUTS:
  # compute pearson correlation between log predicted (computed) and log observed (input)
  # excluding elements that were originally 0 and thus log(0) = -Inf
  
  Xt.retro <- log(Xt.retro,base=10)
  Xt.antero <- log(Xt.antero,base=10)
  # dataframe subsetted by time point
  df.m <- lapply(1:ncol(Xt.retro), function(t.) data.frame(y=y[[t.]],x1=Xt.retro[,t.],x2=Xt.antero[,t.]))
  df.m <- lapply(df.m,inf.nan.mask)
  for(t. in 1:length(df.m)){rownames(df.m[[t.]]) <- paste0(rownames(df.m[[t.]]),'_',t.)} # make unique rownames for each time point
  df.m.rnames <- lapply(df.m,rownames)
  df.fit <- do.call(rbind,df.m) # dataframe that is actually used to fit model
  m <- lm(y~x1+x2,data=df.fit) # fit model on all data
  r <- residuals(m)
  fv <- m$fitted.values
  #e <- mean(residuals(m)^2) # return overall MSE
  e <- cor(fv,df.fit$y) # return overall pearson r
  #e.tp <- sapply(df.m.rnames, function(rn) mean(r[rn]^2)) # return mean squared error for each time point
  e.tp <- sapply(df.m.rnames, function(rn) cor(fv[rn],df.fit[rn,'y'])) # return pearson r for each time point
  return(list(m=m,e=e,e.tp=e.tp,df=df.m)) 
}

c.objective <- function(params.opt,log.path,tps,L.out.retro,L.out.antero,Xo,fxn,one.lm,excl.inj=NULL){
  # fits time constant by modeling CNDR data
  # INPUTS:
  # OPTIMIZED:
  # params.opt: c parameters to optimize cor(observed path, b0 + b.r*retro(c.retro) + b.a*antero(c.antero))
  #     -c.retro: time constant for retrograde model --- params[1]
  #     -c.antero: time constant for anterograde model --- params[2]
  #     -lm.mask.ant.ret fits linear model to find b0, b.a, b.r for each c.retro,c.antero pair
  #     -this is done separately from optimization and params are discarded
  #
  # STATIC:
  # log.path: list of vectors of log-10 transformed pathology scores *in CNDR space*  
  # for each time point specified in tps (below). Time constant c is fit to predict these
  # tps: vector of numeric time points post injection
  # L.out.retro: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads axon to dendrite (retrograde)
  # L.out.antero: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads dendrite to axon (anterograde)
  # Xo: vector of initial pathology
  # fxn: matrix exponential function (fastest is scipy.linalg through reticulate)
  # one.lm: logical indicating whether to fit time-point specific slopes and intercepts or just fit a single slope and intercept
  # excl.inj.CNDR: vector of injection sites in CNDR space to exclude when computing fit
  #
  # this function runs diffusion model in ABA space, 
  # but computes correlation with real data in CNDR annotation space to assess fit of time constant c
  # CNDR names are names of log.path, ABA names are names of Xo
  
  # OUTPUTS:
  # r: mean pearson r value across all provided time points for given set of params
  
  ####
  
  # distribute parameters to interpretable names 
  c.retro <- params.opt[1]
  c.antero <- params.opt[2]
  
#  if(!is.null(excl.inj.CNDR)){ # exclude injection sites from fit by setting them to -Inf in log.path so they'll be excluded by cor.mask below
#    for(t. in 1:length(tps)){log.path[[t.]][excl.inj.CNDR] <- -Inf}
#  }

  ptm <- proc.time()
  Xt.retro <- do.call('cbind',lapply(tps,function(t.) predict.Lout(L.out.retro,Xo,c.retro,t.,fxn=fxn))) # predict path using linear diff model into matrix that is region-by-time
  Xt.antero <- do.call('cbind',lapply(tps,function(t.) predict.Lout(L.out.antero,Xo,c.antero,t.,fxn=fxn))) # predict path using linear diff model into matrix that is region-by-time
  print(paste0('predict: ',(proc.time()-ptm)['elapsed'],'s'))
  print(paste0('optim,',c.retro,',',c.antero)) # output tested parameters to see what optim is doing
  #ptm <- proc.time()
  #Xt.retro <- quiet(map.ABA.to.CNDR(Xt.retro,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  #Xt.antero <- quiet(map.ABA.to.CNDR(Xt.antero,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  #print(paste0('map: ',(proc.time()-ptm)['elapsed']))
  if(!one.lm){
    Xt.sweep <- sapply(1:length(tps), function(t.) lm.mask.ant.ret(log.path[[t.]],Xt.retro[,t.],Xt.antero[,t.])) # this replaces the two lines below
#  } else if(one.lm){list[m,Xt.sweep,e.tp,df] <- lm.mask.ant.ret.all(log.path,Xt.retro,Xt.antero)} 
      } else if(one.lm){lm.mask.results <- lm.mask.ant.ret.all(log.path,Xt.retro,Xt.antero)
      m <- lm.mask.results$m; Xt.sweep <- lm.mask.results$e; e.tp <- lm.mask.results$e.tp; df <- lm.mask.results$df}
  
  print(mean(Xt.sweep)) # display value of fit/objective function
  return(mean(Xt.sweep)) # return mean across provided tps
}


c.CNDRspace.euc.objective <- function(params.opt,log.path,tps,L.out.retro,L.out.antero,L.euclidean,Xo,ABA.to.CNDR.key,fxn,one.lm,excl.inj.CNDR=NULL){
  # fits time constant by modeling CNDR data: using anterograde, retrograde, euclidean distance matrices
  # INPUTS:
  # OPTIMIZED:
  # params.opt: c parameters to optimize cor(observed path, b0 + b.r*retro(c.retro) + b.a*antero(c.antero))
  #     -c.retro: time constant for retrograde model --- params[1]
  #     -c.antero: time constant for anterograde model --- params[2]
  #     -c.euclidean: time constant for euclidean distance model --- params[3]
  #     -lm.mask.ant.ret fits linear model to find b0, b.a, b.r for each c.retro,c.antero pair
  #     -this is done separately from optimization and params are discarded
  #
  # STATIC:
  # log.path: list of vectors of log-10 transformed pathology scores *in CNDR space*  
  # for each time point specified in tps (below). Time constant c is fit to predict these
  # tps: vector of numeric time points post injection
  # L.out.retro: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads axon to dendrite (retrograde)
  # L.out.antero: out-degree graph laplacian of anatomical connectivity matrix, oriented so path spreads dendrite to axon (anterograde)
  # L.euclidean: graph laplacian of euclidean distance matrix (symmetric so no distnictionn)
  # Xo: vector of initial pathology
  # ABA.to.CNDR.key: key to convert ABA to CNDR regions
  # fxn: matrix exponential function (fastest is scipy.linalg through reticulate)
  # one.lm: logical indicating whether to fit time-point specific slopes and intercepts or just fit a single slope and intercept
  # excl.inj.CNDR: vector of injection sites in CNDR space to exclude when computing fit
  #
  # this function runs diffusion model in ABA space, 
  # but computes correlation with real data in CNDR annotation space to assess fit of time constant c
  # CNDR names are names of log.path, ABA names are names of Xo
  
  # OUTPUTS:
  # r: mean pearson r value across all provided time points for given set of params
  
  ####
  
  # distribute parameters to interpretable names 
  c.retro <- params.opt[1]
  c.antero <- params.opt[2]
  c.euclidean <- params.opt[3]

  if(!is.null(excl.inj.CNDR)){ # exclude injection sites from fit by setting them to -Inf in log.path so they'll be excluded by cor.mask below
    for(t. in 1:length(tps)){log.path[[t.]][excl.inj.CNDR] <- -Inf}
  }

  ptm <- proc.time()
  Xt.retro <- do.call('cbind',lapply(tps,function(t.) predict.Lout(L.out.retro,Xo,c.retro,t.,fxn=fxn))) # predict path using linear diff model into matrix that is region-by-time
  Xt.antero <- do.call('cbind',lapply(tps,function(t.) predict.Lout(L.out.antero,Xo,c.antero,t.,fxn=fxn))) # predict path using linear diff model into matrix that is region-by-time
  Xt.euclidean <- do.call('cbind',lapply(tps,function(t.) predict.Lout(L.euclidean,Xo,c.euclidean,t.,fxn=fxn))) # predict path using linear diff model into matrix that is region-by-time
  print(paste0('predict: ',(proc.time()-ptm)['elapsed'],'s'))
  print(paste0('optim,',c.retro,',',c.antero,',',c.euclidean)) # output tested parameters to see what optim is doing
  #ptm <- proc.time()
  Xt.retro <- quiet(map.ABA.to.CNDR(Xt.retro,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  Xt.antero <- quiet(map.ABA.to.CNDR(Xt.antero,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  Xt.euclidean <- quiet(map.ABA.to.CNDR(Xt.euclidean,names(log.path[[1]]),ABA.to.CNDR.key)) # convert matrix to CNDR space
  #print(paste0('map: ',(proc.time()-ptm)['elapsed']))
  if(!one.lm){
    print('ERROR: must use with one.lm option')
  } else if(one.lm){list[m,Xt.sweep,e.tp,df] <- lm.mask.ant.ret.euc.all(log.path,Xt.retro,Xt.antero,Xt.euclidean)}
  
  print(mean(Xt.sweep)) # display value of fit/objective function
  return(mean(Xt.sweep)) # return mean across provided tps
}


lm.mask.ant.ret.euc.all <- function(y,Xt.retro,Xt.antero,Xt.euclidean){
  # fit linear model: y = b0 + b.r*Xt.retro + b.a*Xt.antero + error
  # given Xt.antero and Xt.retro for a given set of time constants and time points
  # return pearson r between y and y-hat (predicted vs actual path)
  # aggregates across all time points to fit one model
  #
  # INPUTS:
  # y: list of observed log pathology from 3 time points
  # Xt.retro: list of predicted pathology from retrograde model from all time points
  # Xt.antero: list of predicted pathology from anterograde model from all time points
  # Xt.euclidean: list of predicted pathology from anterograde model from all time points
  #
  # OUTPUTS:
  # compute pearson correlation between log predicted (computed) and log observed (input)
  # excluding elements that were originally 0 and thus log(0) = -Inf
  
  Xt.retro <- log(Xt.retro,base=10)
  Xt.antero <- log(Xt.antero,base=10)
  Xt.euclidean <- log(Xt.euclidean,base=10)

  # dataframe subsetted by time point
  n.tps <- ncol(Xt.retro)
  df.m <- lapply(1:n.tps, function(t.) data.frame(y=y[[t.]],x1=Xt.retro[,t.],x2=Xt.antero[,t.],x3=Xt.euclidean[,t.]))
  df.m <- lapply(df.m,inf.nan.mask)
  for(t. in 1:n.tps){rownames(df.m[[t.]]) <- paste0(rownames(df.m[[t.]]),'_',t.)} # make unique rownames for each time point
  df.m.rnames <- lapply(df.m,rownames)
  df.fit <- do.call(rbind,df.m) # dataframe that is actually used to fit model
  m <- lm(y~x1+x2+x3,data=df.fit) # fit model on all data
  r <- residuals(m)
  fv <- m$fitted.values
  #e <- mean(residuals(m)^2) # return overall MSE
  e <- cor(fv,df.fit$y) # return overall pearson r
  #e.tp <- sapply(df.m.rnames, function(rn) mean(r[rn]^2)) # return mean squared error for each time point
  e.tp <- sapply(df.m.rnames, function(rn) cor(fv[rn],df.fit[rn,'y'])) # return pearson r for each time point
  return(list(m=m,e=e,e.tp=e.tp,df=df.m)) 
}