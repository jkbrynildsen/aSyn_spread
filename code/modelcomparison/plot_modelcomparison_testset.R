#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','injection.site','measure')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'modelcomparison/traintest/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste0(params$opdir,'processed/',measure,'data.RData'))  # load path data and ROI names

tps <- params$tps

############################
### load train-test data ###
############################

load(file = paste0(savedir,measure,'_ModelComparison_TrainTest.RData'))
mdl.names <- names(results)
mdl.names <- setNames(mdl.names,mdl.names)
mdl.names.short <- mdl.names
mdl.names.short['BidirectionalOneLM'] <- 'Bidirectional'
MPI.names <- paste(tps,'MPI')

#########################################
### compare fits between these models ###
#########################################

fit.list <- lapply(mdl.names, function(mdl.name) 
  as.data.frame(do.call(rbind,lapply(results[[mdl.name]]$test, function(X) X$fits.r))))
names(fit.list) <- mdl.names
# name columns by month
fit.list <- lapply(mdl.names, function(X) data.frame(setNames(fit.list[[X]],MPI.names),Model=X,stringsAsFactors = F,check.names = F))
names(fit.list) <- mdl.names

df.plt <- do.call(rbind,fit.list)
df.plt <- collapse.columns(df.plt,cnames = MPI.names,groupby='Model')
last.tp.fits <- sapply(mdl.names, function(X) mean(df.plt$values[df.plt$names== rev(MPI.names)[1] & df.plt$group ==X]))
mdl.order <- names(sort(last.tp.fits)) # order by their last month fits
p <- ggplot(df.plt) + geom_boxplot(aes(x=group,y=values,fill=group),size=0.25,outlier.size=0.5,fatten=0.5,outlier.stroke = 0) + facet_wrap(~names,nrow = 1) +
  scale_x_discrete(limits=mdl.order,labels=mdl.names.short[mdl.order])+ xlab('') + ylab('Pearson r') +
  scale_fill_manual(values=wes_palettes$Darjeeling1[1:length(mdl.names)],name='',guide=FALSE) + theme_classic()+
  theme(text=element_text(size=8), axis.text.x= element_text(angle=90,hjust=1,vjust=0.5))

ggsave(p,filename = paste(savedir,measure,'_ModelComparisonTestSetPearsonR.pdf',sep=''),
       units = 'cm',height = 6,width = 12,useDingbats=FALSE)

# do a non-parametric test to compare models statistically
p.labs.by.month <- p.vals.by.month <- diff.by.month <- list()
p.signif.matrix <- function(p,n=NULL){
  p.new <- matrix(data = '',nrow = nrow(p),ncol=ncol(p),dimnames=dimnames(p))
  p.new[p > 0.05] <- 'ns'
  p.new[p < 0.05 & p > 0.01] <- '*'
  p.new[p <= 0.01 & p > 1/n] <- '**'
  p.new[p ==0] <- '***' # <1/n
  return(p.new)
}

for(MPI in MPI.names){
  diff.by.month[[MPI]] <- sapply(fit.list, function(M1)  # difference between mean model fit across samples
    sapply(fit.list, function(M2) mean(M2[,MPI] - M1[,MPI])))
  dimnames(diff.by.month[[MPI]]) <- list(unlist(unname(mdl.names.short)),unlist(unname(mdl.names.short)))
  p.vals.by.month[[MPI]] <- sapply(fit.list, function(M1)  # one tailed test where p indicates probability model 2 is better than model 1
    sapply(fit.list, function(M2) mean(M1[,MPI] >= M2[,MPI])))
  dimnames(p.vals.by.month[[MPI]]) <- list(unlist(unname(mdl.names.short)),unlist(unname(mdl.names.short)))
  p.labs.by.month[[MPI]] <- p.signif.matrix(p.vals.by.month[[MPI]],n=length(results[[1]])) 
  dimnames(p.labs.by.month[[MPI]]) <- list(unlist(unname(mdl.names.short)),unlist(unname(mdl.names.short)))
}
clim <- c(min(unlist(diff.by.month)),max(unlist(diff.by.month)))
p.list <- lapply(MPI.names, function(MPI) imagesc(diff.by.month[[MPI]],overlay=p.labs.by.month[[MPI]],
                                                  cmap='redblue_asymmetric',ttl=MPI,caxis_name = 'Row > Col.',clim = clim) +
                   theme(legend.position = 'right',legend.key.width = unit(0.1,'cm'),legend.key.height = unit(0.4,'cm')) + coord_equal()+
                   theme(text=element_text(size=8),legend.box.margin = ggplot2::margin(0,0,0,0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)))
p.all <- plot_grid(plotlist=p.list,nrow=2,ncol=2)
ggsave(p.all,filename = paste0(savedir,measure,'_ModelComparisonTestSetPearsonR_Matrix.pdf'),
       units = 'cm',height = 14,width = 14,useDingbats=FALSE)

# save p-values
pval.np.pub.n <- function(p) pval.np.pub(p,n=nrow(fit.list[[1]])) # set default n -- same for every model
p.list <- lapply(p.vals.by.month,pval.np.pub.n)
p.stat.mat <- do.call(rbind,lapply(MPI.names, function(X) rbind(cbind(p.list[[X]],MPI=X),'-')))
write.csv(x=p.stat.mat,file = paste0(savedir,measure,'_PValsModelComparisonTestSet.csv'))
