### Prepare pathology datasets for analyses
# read in raw data files; export files containing total aSyn pathology for parent and daughter regions

savedir <- paste0(params$opdir, 'processed/')
dir.create(savedir,recursive = T)

# load and format pathology datasets
path_raw0.5_9 <- read.csv(paste0(params$basedir,'data/Compiled_0.5-9MPI_Pathology_Parent_Measures.csv'))
path_raw1 <- read.csv(paste0(params$basedir, 'data/Compiled 1MPI Pathology Measures.csv'))
path_raw0.1 <- read.csv(paste0(params$basedir, 'data/Compiled 0.1MPI Pathology Measures.csv'))
path_raw0.2 <- read.csv(paste0(params$basedir, 'data/Compiled 0.2MPI Pathology Measures.csv'))
path_raw0.3 <- read.csv(paste0(params$basedir, 'data/Compiled 0.3MPI Pathology Measures.csv'))
path_files <- list(path_raw0.1, path_raw0.2, path_raw0.3, path_raw1, path_raw0.5_9)
new_colnames <- c('mouse', 'Region.ID', 'mpi', 'side', 'parent', 'daughter', 'genotype', 'treatment', 'sex', 'marker', 'total_path', 'cell_body_path', 'neurite_path', 'quant_cell_body_path')

# set column names
path_files <- lapply(path_files, function(df) {
  names(df) <- new_colnames
  return(df)
})

# bind all the data frames together
path <- do.call(rbind, path_files)

# identify rows where daughter names contain SNc or SNr
update_SNr <- grepl("SNr", path$daughter)
update_SNc <- grepl("SNc", path$daughter)
# select daughter rows containing SNr or SNc
rows <- path[update_SNr, ]
rows <- path[update_SNc, ]
# update parent row names
path$parent[update_SNr] <- "SNr"
path$parent[update_SNc] <- "SNc"

# subset pathology data into two data frames by hemisphere
Rpath <- path %>% filter(side == "right")
Lpath <- path %>% filter(side == "left")

Rpath$parent <- paste0('i', Rpath$parent)
Rpath$daughter <- paste0('i', Rpath$daughter)
Lpath$parent <- paste0('c', Lpath$parent)
Lpath$daughter <- paste0('c', Lpath$daughter)

full_path <- rbind(Rpath, Lpath) # join data from each hemisphere

# write out cell body, neurite, and total pathology measures to files
for(measure in params$measures){
  
  parent_path <- full_path %>% select(mouse, mpi, parent, paste(measure)) # select relevant variables
  parent_path <- parent_path %>% pivot_wider(names_from = parent, values_from = paste(measure))
  
  write.csv(parent_path, paste0(savedir, measure, '.csv'), row.names = F)
}

# write out cell body quant to file
parent_path <- full_path %>% select(mouse, mpi, parent, quant_cell_body_path) # select relevant variables
parent_path <- parent_path %>% pivot_wider(names_from = parent, values_from = quant_cell_body_path) # , values_fn = mean #to take the mean across all daughter regions within each parent region

write.csv(parent_path, paste0(savedir, 'quant_cell_body_path.csv'), row.names = F)
