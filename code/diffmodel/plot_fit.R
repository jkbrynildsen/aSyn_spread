#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','measure','injection.site')))

basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste0(params$opdir,'processed/',measure,'data.RData'))  # load path data and ROI names
load(paste0(savedir,measure,'BidirectionalFit_data.RData'))
tps=params$tps
# for cell body only: there is no pathology to predict at 0.1 and 0.2 MPI, so exclude those tps
if (measure == "cell_body_path") {
  tps <- c(0.3, 0.5, 1, 3, 6, 9)
  df[[1]] <- NULL
  df[[1]] <- NULL
  Xt.sweep.antero <- Xt.sweep.antero[1:300, 3:8]
  Xt.sweep.retro <- Xt.sweep.retro[1:300, 3:8]
}

###########################
### Plot retrograde fit ###
###########################

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred.retro)),'/',nrow(df[[t]]),' regions left'))
# plot for each time point, using p.xy function 
p <- lapply(1:length(tps), function(t) 
  p.xy(x=df[[t]]$pred.retro,y=df[[t]]$path,ylab='Log(pathology)',xlab='Log(predicted)',
       ttl=paste0(tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste0(savedir,measure,'_retroFit_basemodel.pdf'),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n)
p <- plot.Xt(Xt = t(Xt.sweep.retro),t. = c.rng) + geom_vline(xintercept = c.Grp.retro,color='grey50',linetype='dashed') +
  annotate(geom='text',x=c.Grp.retro,y=max(Xt.sweep.retro),hjust=-0.25,label='Optimal',size=2,color='grey50')+
  scale_color_manual(values = as.character(1:length(tps)),labels=paste(tps,'MPI'),name='') + 
  xlab('c') + ylab('Pearson r with\nPathology') +ggtitle(paste0(paste0(injection.site,collapse = '-'))) +
  theme_bw() + theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),legend.key.size = unit(0.1,'cm'))
p
ggsave(p,filename = paste0(savedir,measure,'_retroCSweepByTimePoint_basemodel.pdf'),
       units = 'cm',height = 4,width = 9)

############################
### Plot anterograde fit ###
############################
savedir <- paste0(params$opdir,'diffmodel/anterograde/',paste0(injection.site,collapse='-'),'/')
dir.create(savedir,recursive=T)

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred.antero)),'/',nrow(df[[t]]),' regions left'))
# plot for each time point, using p.xy function 
p <- lapply(1:length(tps), function(t) 
  p.xy(x=df[[t]]$pred.antero,y=df[[t]]$path,ylab='Log(pathology)',xlab='Log(predicted)',
       ttl=paste0(tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste0(savedir,measure,'_anteroFit_basemodel.pdf'),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n)
p <- plot.Xt(Xt = t(Xt.sweep.antero),t. = c.rng) + geom_vline(xintercept = c.Grp.antero,color='grey50',linetype='dashed') +
  annotate(geom='text',x=c.Grp.retro,y=max(Xt.sweep.antero),hjust=-0.25,label='Optimal',size=2,color='grey50')+
  scale_color_manual(values = as.character(1:length(tps)),labels=paste(tps,'MPI'),name='') + 
  xlab('c') + ylab('Pearson r with\nPathology') +ggtitle(paste0(injection.site,collapse = '-')) +
  theme_bw() + theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),legend.key.size = unit(0.1,'cm'))
p

ggsave(p,filename = paste0(savedir,measure,'_anteroCSweepByTimePoint_basemodel.pdf'),
       units = 'cm',height = 4,width = 9)
