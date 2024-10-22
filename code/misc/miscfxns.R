unit.test <- function(test,pass,fail){
  # test: logical expression to evaluate
  # pass: error to print if true
  # fail: error to print if false
  if(test){
    print(pass)
  } else{print(fail)}
  
}

fisher.r.to.z <- function(r){
  r <- 0.5*(log(1+r) - log(1-r))
  return(r)
}

colSD <- function(x){
  return(apply(x,2,sd))
}

fda <- function(x1,x2,nperms = 1000){
  # x1: N-by-P matrix where there are N observational units and P are features, probably time
  # x2: M-by-P matrix where there are N observational units and P are features, probably time
  # perform functional data analysis, where you compute the mean difference between column
  # means of x1 and x2 and compare them to permuted versions of x1 and x2 with the observational units
  # switched	
  
  
}

get.inf.nan.mask <- function(x){
  # INPUTS:
  # x: matrix, df, or vector
  # 
  # OUTPUTS:
  # mask: mask for x with all rows (if matrix or df) or elements (if vector) containing Infs or NaNs set to FALSE
  
  if(is.vector(x)){
    mask <- x %in% c(-Inf,Inf,NaN,NA)# find infs or nans
  } else if(is.matrix(x) | is.data.frame(x)){
    mask <- do.call('cbind',lapply(1:ncol(x),function(j) x[,j] %in% c(-Inf,Inf,NaN,NA))) # find infs or nans
  }
  return(mask)
}

inf.nan.mask <- function(x){
  # INPUTS:
  # x: matrix, df, or vector
  #
  # OUTPUTS:
  # x.masked: x with all rows (if matrix or df) or elements (if vector) containing Infs or NaNs removed
  
  if(is.vector(x)){
    mask <- x %in% c(-Inf,Inf,NaN,NA)# find infs or nans
    x.masked <- x[!mask]
  } else if(is.matrix(x) | is.data.frame(x)){
    mask <- do.call('cbind',lapply(1:ncol(x),function(j) x[,j] %in% c(-Inf,Inf,NaN,NA))) # find infs or nans
    x.masked <- x[rowSums(mask)==0,]
  }
  return(x.masked)
  
}

source.save <- function(script,output){
  # wrapper for source function that saves output and input of script
  file.create(output)
  con <- file(output)
  #sink(con, append=TRUE)
  sink(con, type="output")
  source(script)
  sink() 
  sink(type="output")
  
}

quiet <- function(x) { 
  # https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message/34208658#34208658
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

name <- function(x,x.names){
  # INPUTS:
  # x: vector or dataframe
  # x.names: names for elements or columns of x
  # OUTPUTS:
  # x with x.names as names
  names(x) <- x.names
  return(x)
}

collapse.columns <- function(df,cnames=colnames(df),groupby=NULL){
  # INPUTS:
  # df: dataframe
  # cnames: column names to perform operation on, default to all columns
  # groupby: column name to group variables by, treated separately from variables in cnames
  
  # OUTPUTS:
  # df.new: dataframe with 2 columns:
  # values: all columns in cnames vertically concatenated. 
  # names: elements of cnames corresponding to rows
  # group: groups of observations in df for each variable in cnames
  
  df.names <- do.call('cbind',lapply(cnames, function(n) rep(n,nrow(as.matrix(df[,cnames])))))  
  df.new <- data.frame(values = as.vector(as.matrix(df[,cnames])),names=as.vector(df.names),stringsAsFactors=F)
  if(!is.null(groupby)){
    df.grp <- do.call('cbind',lapply(cnames,function(n) df[,groupby]))
    df.new$group <- as.vector(df.grp)
  }
  return(df.new)
}

col.Which.Max <- function(x){
  cwm <- unlist(apply(x,2,function(y) which.max(y)))
  return(cwm)
}

row.Which.Max <- function(x){
  rwm <- unlist(apply(x,1,function(y) which.max(y)))
  return(rwm)
}

hemi.split <- function(df,rename = FALSE){
  # INPUTS: 
  # df: data frame whose rownames are region names with i/c for ipsi contra
  # rename: logical indicating whether to rename rows by region name without i/c designation
  #
  # OUTPUTS: df split by hemisphere based on rownames
  
  hemi <- substr(rownames(df),start=1,stop=1)
  df.i <- df[hemi=='i',,drop=FALSE]
  df.c <- df[hemi=='c',,drop=FALSE]
  names.base.i <- substr(rownames(df.i),start=2,stop=nchar(rownames(df.i)))
  names.base.c <- substr(rownames(df.c),start=2,stop=nchar(rownames(df.c)))
  if(identical(names.base.i,names.base.c) & rename){
    rownames(df.c) <- names.base.i
    rownames(df.i) <- names.base.i
  }
  return(list(df.i=df.i,df.c=df.c))
}

hemi.expand.names <- function(df){
  # INPUTS:
  # df: data frame whose rownames are region names WITHOUT 'i' or 'c' as first character
  # 
  # OUTPUTS:
  # df.cat: input df, duplicated and vertically concatenated so it has a row corresponding to 'i' and 'c' for every region
  # useful utility in manipulating gene data (hemispherically symmetric) and vulnerability/path data (asymmetric)
  
  df.names <- rownames(df)
  df.cat <- rbind(df,df)
  rownames(df.cat) <- c(paste0('i',df.names),paste0('c',df.names))
  return(df.cat)
  
}
hemi.average <- function(df,r.v=FALSE){
  # INPUTS:
  # df: data frame whose row names are ABA regions with 'i' or 'c' as first character
  # r.v: logical indicating whether to return named vector instead of data frame
  #
  # OUTPUTS:
  # df averaged over columns (time) then across hemispheres
  if(is.vector(df)){df <- data.frame(Path=df)}
  split_hemi <- hemi.split(df)
  df.c <- split_hemi$df.c; df.i <- split_hemi$df.i
  names.base.i <- substr(rownames(df.i),start=2,stop=nchar(rownames(df.i)))
  names.base.c <- substr(rownames(df.c),start=2,stop=nchar(rownames(df.c)))
  df.c <- df.c[match(names.base.i,names.base.c),,drop=FALSE] # match to order of ipsi hemisphere
  rownames(df.c) <- substr(rownames(df.c), 2, nchar(rownames(df.c)))
  rownames(df.i) <- substr(rownames(df.i), 2, nchar(rownames(df.i)))
  df.c <- rowMeans(df.c,na.rm=T)
  df.i <- rowMeans(df.i,na.rm=T)
#  if(identical(names.base.i,names.base.c)){
#    rownames(df.c) <- names.base.i
#    rownames(df.i) <- names.base.i
#    df.c <- rowMeans(df.c,na.rm=T)
#    df.i <- rowMeans(df.i,na.rm=T)
#  if(!r.v){return(data.frame(v=rowMeans(cbind(df.i,df.c),na.rm=T)))}
  return(rowMeans(cbind(df.i,df.c),na.rm=T))
#  }
  
}

pval.np.pub <- function(p,n){
  # INPUTS:
  # p: matrix or vector of p values from non-parametric test
  # n: number of permutations from test
  #
  # OUTPUTS:
  # p.new: matrix or vector of p values with p = 0 reworded to p < 1/n
  
  p.new <- p
  p.new[p == 0] <- paste('p <',1/n)
  return(p.new)
}
pval.2tail.np <-  function(test.val,dist){
  # test.val: individual value being compared to distribution
  # dist: vector,distribution of values under some null model, or otherwise
  # sig.fig: number of significant figures
  # compute 2-tailed p-value for test value occurring in distribution
  dist <- as.numeric(dist)
  pval.2tail <- 2*min(mean(test.val >= dist),mean(test.val <= dist))
  #pval.2tail <- 2*(1-min(mean(test.val >= dist),mean(test.val <= dist)))
  return(pval.2tail)
}

getGroupColors <- function(grps=NULL){
  # INPUTS:
  # vector of group names
  return(setNames(c('#047391','#FF7F00'),grps))
}

get3GroupColors <- function(measures=NULL){
  # INPUTS:
  # vector of group names
  return(setNames(c('#047391','#FF7F00','#A478B8'),measures))
}

get4GroupColors <- function(grps=NULL){
  # INPUTS:
  # vector of group names
  return(setNames(c('#f37735','#00aedb','#00b159', '#d11141'),grps))
}

list.posthoc.correct <- function(X,method){
  # unlist a list, posthoc correct over all values according to "method"
  # relist the list in the same structure and return
  return(relist(flesh=p.adjust(unlist(X),method=method),skeleton=X))
}

outlier.mask <- function(x){
  # return mask of values in x that are outliers
  return(x %in% boxplot.stats(x)$out)
}