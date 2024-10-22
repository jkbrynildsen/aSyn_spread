rm(list=setdiff(ls(),c('params')))
basedir <- params$basedir # make this whatever you like, end with /
setwd(basedir)
source('code/aba/aba_fxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
abadir <- paste0(basedir,'data/aba/')
dir.create(paste0(abadir,'atlas'),recursive=TRUE)
load(file=paste0(abadir,'ontologies/keys.RData'))

region.names.id <- unlist(region.names.id.key)
region.ids.unique <- region.names.id[substr(names(region.names.id),1,1) == 'i']
region.names.unique <- substr(names(region.ids.unique),2,max(nchar(names(region.ids.unique))))
df.py <- data.frame(names=region.names.unique,ids=unname(region.ids.unique))
write.csv(df.py,file=paste0(abadir,'atlas/structIDs.csv'),row.names = F)

# run code from Ben Fulcher to get coordinates of ABA ROIs https://github.com/benfulcher/AllenSDK
structID.fname <- paste0(abadir,'atlas/structIDs.csv') # input filename for structure IDs
output.fname <- paste0(abadir,'atlas/mask_aba.mat') # output filename for mask
# switch to environment with packages installed
#system(paste('source ~/.bash_profile; source activate pyforge; python code/aba/MakeCCFMasks.py',Sys.glob(basedir),structID.fname,output.fname))

# cd code/aba
# IN PYTHON, EXECUTE:
# python MakeCCFMasks.py atlas/structIDs.csv atlas/mask_aba.mat
# /Users/jkbrynildsen/Box/PD_Project/Quant_pathology_dataset/code/aba/

# load mask for each brain region, doesn't distinguish L and R hemisphere
pymask <- readMat(output.fname)
atlas <- pymask$mask
inds <- pymask$region.inds
print(dim(atlas))
atlas.hemi.mask <- atlas*0
atlas.hemi.mask[,,(dim(atlas)[3]/2):(dim(atlas)[3])] <- TRUE
for(idx in inds){
  print(idx)
  # for each region (which contains both hemispheres), split into a left and right hemi region
  # by adding number of regions in hemisphere to each index
  atlas[atlas == idx & atlas.hemi.mask] <- idx + length(inds) 
}

# get coordinates of the center of mass of each region
res.mm <- 25/1000 # resolution in millimeters, i.e. distance between array indices
cent.coors <- do.call(rbind,lapply(1:length(region.names.id), function(region) colMeans(which(atlas == region,arr.ind=T))))
D.mat <- as.matrix(dist(cent.coors))
region.names.dmat <- names(region.names.id)
region.names.dmat[is.na(rowSums(cent.coors))] # missing subiculum subregions. larger subiculum is in there but i think easier to just ignore than duplicate
save(D.mat,region.names.dmat,cent.coors,res.mm,file=paste0(params$opdir,'processed/ABAEuclideanDistanceMatrix.RData'))
