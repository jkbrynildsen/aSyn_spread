### This script assesses specificity of model to the actual seed site

rm(list=setdiff(ls(),c('params','measure','injection.site')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'nullmodels/seedspec_multi/',paste0(injection.site,collapse='-'),'/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/plottingfxns.R')

#################
### Load data ###
#################

load(paste0(params$opdir,'processed/',measure,'data.RData')) # load path data and ROI names
load(file = paste(savedir,measure,'_AlternateSeedFits_OneLM.RData',sep=''))
load(file = paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',measure,'BidirectionalOptimOneLM_data.RData',sep=''))

# compute p-value as probability that other seeds predict observed path better than iCPu seed
null.cors <- seed.fits

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
inj.cor <- m.fits

seed.region.pval <- colMeans(do.call(rbind,lapply(1:nrow(null.cors), function(X) inj.cor)) < null.cors)

months <- paste(params$tps,'MPI')
months.mat <- do.call(rbind,lapply(1:nrow(null.cors), function(X) months))
p.labs <- paste('p =',signif(seed.region.pval,2))
p.labs[seed.region.pval ==0] <- paste0('p < ',1/nrow(null.cors))
p.null.seeds <- ggplot() + 
  geom_jitter(aes(x=as.vector(months.mat),y = as.vector(null.cors)),color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1, width=0.25) +
  geom_point(aes(x=months,y=as.numeric(inj.cor)),shape = 18,color = 'black',size=2) + 
  geom_text(aes(x=months,y=0.8,label = p.labs),size=2.5) +
  xlab('') + ylab('Fit') + ggtitle('Actual vs. Random Seed') + theme_bw() +
  theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8)) +
  theme(axis.text.x = element_text(size=8)) + theme(axis.text.y = element_text(size=8))
p.null.seeds

ggsave(p.null.seeds,filename=paste(savedir,measure,'_SeedSpecificity_OneLM.pdf',sep=''),
	width=7,height=5,units='cm')
