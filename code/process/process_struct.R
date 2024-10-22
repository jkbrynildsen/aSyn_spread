### Prepare structural connectome matrix for network analyses
# EDIT line 17 depending on matrix being used

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'processed/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')

# read in raw pathology dataset to get region names
path.data <- read.csv(paste0(basedir,'data/Compiled_0.5-9MPI_Pathology_Parent_Measures.csv'), header = TRUE)
path_regions <- unique(path.data$parent) # create list of region names from pathology dataset

# read in structural connectivity dataset
struct_mat <- read.csv(paste0(basedir,"data/normalized_connection_strength.csv")) # EDIT depending on which structural matrix using
struct_mat2 <- struct_mat[,-1] # eliminate first column, which contains regions names
row.names(struct_mat2) <- struct_mat[,1] # convert first column into row names

# subset connectivity matrix into ipsi and contralateral hemispheres
connectivity.contra <- struct_mat2 %>% 
  select(starts_with("contra"))
colnames(connectivity.contra) <- connectivity.contra[c(1),] # make column names match first row
connectivity.contra <- connectivity.contra[-c(1),] # remove first row
connectivity.contra <- connectivity.contra[, names(connectivity.contra) %in% path_regions] # select only columns that are present in the pathology dataset
connectivity.contra <- connectivity.contra[rownames(connectivity.contra) %in% names(connectivity.contra),] # select only rows that are present in the pathology dataset

connectivity.ipsi <- struct_mat2 %>% 
  select(starts_with("ipsi"))
colnames(connectivity.ipsi) <- connectivity.ipsi[c(1),] # make column names match first row
connectivity.ipsi <- connectivity.ipsi[-c(1),] # remove first row
connectivity.ipsi <- connectivity.ipsi[, names(connectivity.ipsi) %in% names(connectivity.contra)] # select only columns that are present in the pathology dataset
connectivity.ipsi <- connectivity.ipsi[rownames(connectivity.ipsi) %in% names(connectivity.contra),] # select only rows that are present in the pathology dataset
# NOTE: the above excludes the 6 regions that are missing contralateral structural edges

conn.names.ipsi <- colnames(connectivity.ipsi)
conn.names.contra <- colnames(connectivity.contra)

# checks 
if(identical(colnames(connectivity.contra),rownames(connectivity.contra))){
  print('contra connectivity colnames and rownames equal')}
if(identical(colnames(connectivity.ipsi),rownames(connectivity.ipsi))){
  print('ipsi connectivity colnames and rownames equal')}

# create elements of full connectivity matrix
upper.left <- connectivity.ipsi
upper.right <- connectivity.contra
lower.left <- connectivity.contra
lower.right <- connectivity.ipsi

W <- rbind(cbind(upper.left,upper.right),cbind(lower.left,lower.right))

rownames(W) <- c(paste('i',rownames(upper.left),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(lower.left),sep='')) # add c to contralateral regions
colnames(W) <- c(paste('i',rownames(upper.left),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(lower.left),sep='')) # add c to contralateral regions

n.regions.ABA <- nrow(W)
n.regions.ABA.hemi <- n.regions.ABA/2

# check if connectivity was tiled into alternating blocks correctly
unit.test(all(W[(n.regions.ABA.hemi+1):n.regions.ABA,(n.regions.ABA.hemi+1):n.regions.ABA] == W[1:n.regions.ABA.hemi,1:n.regions.ABA.hemi]),
          'tiling on-diagonal blocks worked','tiling on-diagonal blocks failed') 
unit.test(all(W[1:n.regions.ABA.hemi,(n.regions.ABA.hemi+1):n.regions.ABA] == W[(n.regions.ABA.hemi+1):n.regions.ABA,1:n.regions.ABA.hemi]),
          'tiling off-diagonal blocks worked','tiling off-diagonal blocks failed')

unit.test(all(colnames(W) == rownames(W)),'row and column names of conn mat are same','ERROR with conn mat names')

total_path <- read.csv(paste0(savedir,'total_path.csv'), check.names=FALSE)
filter.path.names <- colnames(total_path)[-c(1:2)] # exclude first two columns which contain mouse ID and MPI

W <- W[, names(W) %in% filter.path.names] # select only columns that are present in the pathology dataset
W <- W[rownames(W) %in% filter.path.names,] # select only rows that are present in the pathology dataset

W_labeled <- W # store a version of the structural matrix with region labels

n.regions.ABA <- nrow(W)
n.regions.ABA.hemi <- n.regions.ABA/2

region.names <- colnames(W) # store connectivity matrix names

# read in pathology data, match regions with those in the structural dataset
for(measure in params$measures){
  PathData <- read.csv(paste0(savedir,measure,'.csv'), check.names=FALSE)
  path.data <- PathData[, names(PathData) %in% region.names]
  path.data <-path.data[region.names] # match order of pathology regions to order of structural regions
  path.data <- cbind(PathData[,1:2], path.data)
  path.names <- colnames(path.data)[-c(1:2)]
  
  save(path.data, path.names,region.names, n.regions.ABA, n.regions.ABA.hemi, file = paste0(savedir,measure,'data.RData'))
}

write.table(W, paste0(params$opdir,'processed/W.csv'), sep=',', row.names = FALSE, col.names = FALSE)
write.csv(W_labeled, paste0(params$opdir,'processed/W_labeled.csv'))
