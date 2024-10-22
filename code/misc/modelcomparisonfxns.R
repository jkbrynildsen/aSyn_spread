source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')

train.test.logmeanpath <- function(Mice,nreps=500,tf=0.5){
  # INPUTS:
  # Mice: list of dataframes (one for each time point) that is mice x regions whose values contain pathology measured at a specific time point
  # nreps: number of train-test splits
  # tf: training fraction
  #
  # OUTPUTS:
  # log.path.train: log of mean pathology at each time point in training sample
  # log.path.test: log of mean pathology at each time point in testing sample corresponding elementwise to the training samples
  # train.idx: indexes used to split data
  
  train.idx <- lapply(1:nreps, function(R) lapply(Mice, function(X) sample(1:nrow(X),size = tf*nrow(X))))
  # split into non-overlapping training and testing sets - get list of lists of mean pathology
  log.path.train <- lapply(train.idx, function(t.idx) mapply(function(X,t.idx) {list(log(colMeans(X[t.idx,],na.rm=T),base=10))}, X=Mice,t.idx=t.idx))
  log.path.test <- lapply(train.idx, function(t.idx) mapply(function(X,t.idx) {list(log(colMeans(X[-t.idx,],na.rm=T),base=10))}, X=Mice,t.idx=t.idx))
  return(list(log.path.train=log.path.train,log.path.test=log.path.test,train.idx=train.idx))
  
}

spread.model.train <- function(X.train,mdl.name,ctrl){
  # this function is a wrapper allowing you to fit the various models that are compared in this paper
  # it will take in training data, a name of a model implemented in this repository, and the static data necessary to implement it
  
  # INPUTS:
  # X.train: list with training data. elements contain log mean path at each time point
  # mdl.name: character selecting the model by name
  # ctrl: list of static data or options for the model. Can contain:
  #     Xo: Nx1 vector with initial pathology seeds
  #     *: laplacians to test, accessed by mdl.name. If you want to test a bidirectional model, you specifically need "Retrograde" and "Anterograde"
  #     tps: numeric vector of time points corresponding to X.train
  #     c.rng: range of c values to test for indpendent fit
  #     c.r.a.init: two-element vector of initial c values for bidirectional model (retro, antero)
  #     ctrl.optim: list of controls for optim function
  #     one.lm: if using bidirectional model, use one lm or an lm for each time point?
  #     ABA.to.CNDR.key: key to convert from ABA to CNDR space
  #     fxn: matrix exponential function
  #
  # OUTPUTS:
  # m.out: list of variables with trained model, predicted values, and other information
  if(grepl('OneLM',mdl.name)){
    ctrl$one.lm <- TRUE
    #ctrl$ctrl.optim$fnscale <- 1 # minimize MSE for this model #
  }
  if(!grepl('Bidirectional',mdl.name)){ #if 'bidirectional' is not in the name
    m.out <- spread.model.unidirectional.train(X.train,mdl.name,ctrl)
  } else{
    if(grepl('Euclidean',mdl.name)){ #if bidirectional and euclidean are in the name       
        m.out <- spread.model.bidirectional.euclidean.train(X.train,mdl.name,ctrl) # bidirectional connectivity + euclidean distance model
      } else{ #if only bidirectional is in the name
        m.out <- spread.model.bidirectional.train(X.train,mdl.name,ctrl) # bidirectional connectivity model
      }    
  }
  return(m.out)
}

spread.model.unidirectional.train <- function(X.train,mdl.name,ctrl){
  # INPUTS:
  # see spread.model.train() for details
  #
  # OUTPUTS:
  # m.out: list with time constant, predicted values, time points, linear models, and model name
  
  # fit model
  results <- c.fit(X.train,ctrl$tps,ctrl[[mdl.name]],ctrl$Xo,ctrl$c.rng,ctrl$excl.inj)  #ctrl$ABA.to.CNDR.key,#replaced list[c.train,Xt.sweep] with results 6/28/23
  #added stuff below 6/28/23
  c.train <- results$c.best
  Xt.sweep <- results$Xt.sweep

  # get predicted values
  do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.sparse.linalg$expm)),base=10))) #expm::expm
  Xt.Grp <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(ctrl[[mdl.name]],ctrl$Xo,c.train,t,fxn=scipy.sparse.linalg$expm)),base = 10))) # predict pathology using connectivity, time constant, and seed  #scipy.sparse.linalg$expm
  Xt.Grp[rownames(Xt.Grp) %in% ctrl$excl.inj,] <- -Inf
  # fit linear models
  # concatenate pathology with predicted values
  df.train <- lapply(1:length(tps), function(t) data.frame(path = X.train[[t]], pred= Xt.Grp[,t,drop=FALSE]))
  m.train <- lapply(df.train, function(df.i) lm(path~pred,data=inf.nan.mask(df.i))) # fit linear regression at each time point based on specified model prediction
  m.out <- list(c.train=c.train,Xt.Grp=Xt.Grp,tps=ctrl$tps,m.train=m.train,mdl.name=mdl.name)
  return(m.out)

}

spread.model.bidirectional.train <- function(X.train,mdl.name,ctrl){
  # INPUTS:
  # see spread.model.train() for details
  #
  # OUTPUTS:
  # m.out: list with time constant, predicted values, time points, linear models, and model name
  #ABA.to.CNDR.key <- ctrl$ABA.to.CNDR.key
  #path.names <- names(ABA.to.CNDR.key)
  Xo <- ctrl$Xo
  L.out.retro <- ctrl$Retrograde
  L.out.antero <- ctrl$Anterograde
  params.opt.fit <- optim(ctrl$c.r.a.init,c.objective,control = ctrl$ctrl.optim, method="L-BFGS-B",lower=c(1e-7,1e-7), # optimization. c's must be > 0
                          log.path=X.train,tps=ctrl$tps,L.out.retro=ctrl$Retrograde,L.out.antero=ctrl$Anterograde,
                          Xo=ctrl$Xo,fxn =ctrl$fxn,one.lm = ctrl$one.lm) # static inputs
  # extract parameters from optimization output
  c.train.retro <- params.opt.fit$par[1]
  c.train.antero <- params.opt.fit$par[2]
  
  # get predicted values - apply log base 10 for fitting linear models

  Xt.Grp.retro <-   do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.retro,Xo,c.train.retro,t,fxn=ctrl$fxn)),base=10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.antero,Xo,c.train.antero,t,fxn=ctrl$fxn)),base=10))) # predict pathology using connectivity, time constant, and seed
  # exclude injection sites if specified in params
  Xt.Grp.retro[rownames(Xt.Grp.retro) %in% ctrl$excl.inj,] <- -Inf
  Xt.Grp.antero[rownames(Xt.Grp.antero) %in% ctrl$excl.inj,] <- -Inf
  # fit linear models
  if(ctrl$one.lm){
    #list[m,e,m.fits,df] <- lm.mask.ant.ret.all(X.train,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function has log10 built in
    lm.mask.results <- lm.mask.ant.ret.all(X.train,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it
    m <- lm.mask.results$m; e <- lm.mask.results$e; m.fits <- lm.mask.results$e.tp; df <- lm.mask.results$df
  } else{
    df <- lapply(1:length(tps), function(t) data.frame(path = X.train[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
    m <- lapply(df, function(df.i) lm(path~pred.retro+pred.antero,data=inf.nan.mask(df.i))) # fit linear regression at each time point to combine anterograde and retrograde
  }
  
  # save output
  m.out <- list(c.train.retro=c.train.retro,c.train.antero=c.train.antero,Xt.Grp.antero=Xt.Grp.antero,Xt.Grp.retro=Xt.Grp.retro,
                tps=ctrl$tps,m.train=m,mdl.name=mdl.name)
  return(m.out)
  
}

spread.model.bidirectional.euclidean.train <- function(X.train,mdl.name,ctrl){
  # INPUTS:
  # see spread.model.train() for details
  #
  # OUTPUTS:
  # m.out: list with time constant, predicted values, time points, linear models, and model name
  ABA.to.CNDR.key <- ctrl$ABA.to.CNDR.key
  path.names <- names(ABA.to.CNDR.key)
  Xo <- ctrl$Xo
  L.out.retro <- ctrl$Retrograde
  L.out.antero <- ctrl$Anterograde
  L.euclidean <- ctrl$Euclidean
  params.opt.fit <- optim(ctrl$c.r.a.init,c.CNDRspace.euc.objective,control = ctrl$ctrl.optim, method="L-BFGS-B",lower=c(1e-7,1e-7), # optimization. c's must be > 0
                          log.path=X.train,tps=ctrl$tps,L.out.retro=ctrl$Retrograde,L.out.antero=ctrl$Anterograde,L.euclidean=ctrl$Euclidean,
                          Xo=ctrl$Xo,ABA.to.CNDR.key=ctrl$ABA.to.CNDR.key,fxn =ctrl$fxn,one.lm = ctrl$one.lm,excl.inj.CNDR=ctrl$excl.inj) # static inputs
  # extract parameters from optimization output
  c.train.retro <- params.opt.fit$par[1]
  c.train.antero <- params.opt.fit$par[2]
  c.train.euclidean <- params.opt.fit$par[3]

  # get predicted values - apply log base 10 for fitting linear models
  Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.train.retro,t,fxn=ctrl$fxn),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.train.antero,t,fxn=ctrl$fxn),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.euclidean <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.euclidean,Xo,c.train.euclidean,t,fxn=ctrl$fxn),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
  # exclude injection sites if specified in params
  Xt.Grp.retro[rownames(Xt.Grp.retro) %in% ctrl$excl.inj,] <- -Inf
  Xt.Grp.antero[rownames(Xt.Grp.antero) %in% ctrl$excl.inj,] <- -Inf
  Xt.Grp.euclidean[rownames(Xt.Grp.euclidean) %in% ctrl$excl.inj,] <- -Inf
  # fit linear models
  if(ctrl$one.lm){
    list[m,e,m.fits,df] <- lm.mask.ant.ret.euc.all(X.train,10^Xt.Grp.retro,10^Xt.Grp.antero,10^Xt.Grp.euclidean) # undo log10 because this function has log10 built in
  } else{
    print('ERROR: must use onelm option with this model')
  }
  
  # save output
  m.out <- list(c.train.retro=c.train.retro,c.train.antero=c.train.antero,c.train.euclidean=c.train.euclidean,
          Xt.Grp.antero=Xt.Grp.antero,Xt.Grp.retro=Xt.Grp.retro,Xt.Grp.euclidean=Xt.Grp.euclidean,
                tps=ctrl$tps,m.train=m,mdl.name=mdl.name)
  return(m.out)
  
}

spread.model.eval <- function(m.out,X.test,save.resid=FALSE){
  # INPUTS:
  # m.out: output from spread.model.unidirectional.train or spread.model.bidirectional.train, contains trained models
  # X.test: data to apply to model in m.out (can be test set or training set)
  # save.resid: logical indicating whether to save residuals in results
  #
  # OUTPUTS:
  # results: 
  
  m.train <- m.out$m.train # get list of linear models from training process, depending on model type
  tps <- m.out$tps # get time points in vector
  Xt.Grp.retro <- m.out$Xt.Grp.retro
  Xt.Grp.antero <- m.out$Xt.Grp.antero
  Xt.Grp.euclidean <- m.out$Xt.Grp.euclidean
  Xt.Grp <- m.out$Xt.Grp
  c.train.retro <- m.out$c.train.retro
  c.train.antero <- m.out$c.train.antero
  c.train <- m.out$c.train
  mdl.name <- m.out$mdl.name
  results <- list()
  if(mdl.name == 'Bidirectional'){
    df.test <- lapply(1:length(tps), function(t) data.frame(path = X.test[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
    df.test <- lapply(1:length(tps), function(t.) inf.nan.mask(cbind(df.test[[t.]],pred=predict(m.train[[t.]],df.test[[t.]]))) ) # add predicted values from linear models
    results$c.train.retro <- c.train.retro
    results$c.train.antero <- c.train.antero
    results$m.train.coefs <- lapply(1:length(tps), function(t) summary(lm.beta(m.train[[t]]))$coef)#summary(lm.beta(m.train))$coef # store standardized coefficients - added lapply and [[t]] 2/9/23
  } else if(mdl.name == 'BidirectionalOneLM'){
    df.test <- lapply(1:length(tps), function(t) data.frame(path = X.test[[t]], x1 = Xt.Grp.retro[,t,drop=FALSE], x2 =Xt.Grp.antero[,t,drop=FALSE]))
    df.test <- lapply(1:length(tps), function(t.) inf.nan.mask(cbind(df.test[[t.]],pred=predict(m.train,df.test[[t.]]))) ) # add predicted values from linear models
    results$c.train.retro <- c.train.retro
    results$c.train.antero <- c.train.antero
    results$m.train.coefs <- lapply(1:length(tps), function(t) summary(lm.beta(m.train))$coef) # store standardized coefficients
  } else if(mdl.name == 'BidirectionalOneLMEuclidean'){
    df.test <- lapply(1:length(tps), function(t) data.frame(path = X.test[[t]], x1 = Xt.Grp.retro[,t,drop=FALSE], x2 =Xt.Grp.antero[,t,drop=FALSE],x3 =Xt.Grp.euclidean[,t,drop=FALSE]))
    df.test <- lapply(1:length(tps), function(t.) inf.nan.mask(cbind(df.test[[t.]],pred=predict(m.train,df.test[[t.]]))) ) # add predicted values from linear models
    results$c.train.retro <- c.train.retro
    results$c.train.antero <- c.train.antero
    results$c.train.euclidean <- m.out$c.train.euclidean
    results$m.train.coefs <- lapply(1:length(tps), function(t) summary(lm.beta(m.train[[t]]))$coef)#summary(lm.beta(m.train))$coef # store standardized coefficients - added [[t]] 2/9/23   
  }

  else{
    # for all unidirectional models
    df.test <- lapply(1:length(tps), function(t) data.frame(path = X.test[[t]], pred = Xt.Grp[,t,drop=FALSE])) # use prediction from network diffusion mode (e^-(Lct)%*%Xo)
    # now apply time-point specific linear models to get MSE
    df.test <- lapply(1:length(tps), function(t.) inf.nan.mask(cbind(df.test[[t.]],pred=predict(m.train[[t.]],df.test[[t.]]))) ) # add predicted values from linear models
    # ^- when using these predicted values for each time point (which are shifted and scaled within each time point), correlation within each time point will be the same
    # those correlations averaged over time points is the primary fit metric used in this paper
    results$c.train <- c.train
  }
  
  # calculate correlation for each time point and MSE
  fits.r <- sapply(df.test, function(df) cor(df$path,df$pred)) # extract fits as pearson r
  fits.mse <- sapply(df.test, function(df) mean((df$path-df$pred)^2)) # extract fits as MSE
  # save results (concatenate with existing results that has time constants)
  results <- c(results,list(fits.r=fits.r,fits.mse=fits.mse))
  
  if(save.resid){ # save model residuals with respect to input data?
    # see code/diffmodel/plotCNDRspacebidirectionalonelmfit.R for how to do this with one lm model
    resid <- sapply(m.train, function(m.i) residuals(m.i)) # get residuals of each model
    df.resid <- df.resid.init
    for(tp in 1:length(tps)){df.resid[names(resid[[tp]]),paste(tps[tp],'MPI')] <- resid[[tp]]}
    results$resid <- df.resid
  }
  return(results)
  
}

