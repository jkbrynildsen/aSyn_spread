source('code/misc/miscfxns.R')

rank_INT_base <- function(x,k=3/8){
  mask <- x %in% c(-Inf,Inf,NaN,NA)
  x[!mask] <- qnorm((rank(x[!mask])-k)/(length(x[!mask]) - 2*k + 1))
  return(x)
}

rank_INT <- function(x,k=3/8){
	# INPUTS:
	# x: any R object
	# OUTPUTS:
	# that object normalized using rank inverse normal transformation
	# elementwise if list
	# columnwise if matrix or data frame

	if(is.list(x) & !is.data.frame(x)){
		return(lapply(x,rank_INT))
	} else if(is.matrix(x)){
		x.names <- name(colnames(x),colnames(x))
		if(is.null(x.names)){x.names <- 1:ncol(x)}
		x.int <- sapply(x.names, function(j) rank_INT_base(x[,j]))
		rownames(x.int) <- rownames(x)
		return(x.int)
	} else if(is.vector(x)){
		return(rank_INT_base(x))
	} else if(is.data.frame(x)){
		x.int <- as.data.frame(sapply(name(names(x),names(x)), function(j) rank_INT_base(x[,j])))
		rownames(x.int) <- rownames(x)
		return(x.int)
	}
}


log10 <- function(x) log(x,base=10)