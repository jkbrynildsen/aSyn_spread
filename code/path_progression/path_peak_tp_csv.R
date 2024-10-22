### Group regions according to the time at which peak pathology levels are reached and output .csv that can be used to generate heatmaps

### Format directories and data ####
rm(list=setdiff(ls(),c('params','measure','injection.site')))
opdir <- params$opdir
basedir <- params$basedir
setwd(basedir)

savedir <- paste0(params$opdir, 'path_progression/')
dir.create(savedir,recursive = T)

source('code/misc/plottingfxns.R')

load(paste0(opdir,'processed/', measure, 'data.RData'))  # load path data and ROI names

# function to format data for outputting as a csv that can be used to make heatmaps
csv_format <- function(df){
  # obtain mean expression across mice for each brain region
  df <- df %>%
    group_by(mpi) %>%
    summarise(across(.cols = -1, .fns = mean, na.rm = TRUE)) %>% # average across mice for each time point
    ungroup() %>%
    arrange(mpi) # chronologically order time points
  # remove column that contained mouse IDs
  df$mouse <- NULL
  # pivot data frame into long format with three columns: mpi, region, and path
  df <- df %>%
    pivot_longer(!mpi, names_to = 'region', values_to = 'path') %>%
    mutate(ord_tp = ord_tps(mpi))
  df$path <- log10(df$path) # log transform path data
  df <- df %>%
    filter(!is.na(path) & !is.infinite(path)) # filter out NA and 0 path regions
  return(df)
}

### Plot log pathology for each region across time ####
## set up data frame for plotting regional log pathology heatmaps
plot_df <- csv_format(path.data)

### Generate tables containing the peak time for each region ####
# pull out ipsilateral regions
ipsi_df <- plot_df %>%
  filter(str_detect(region, "^i"))

# for each ipsi region, identify time point when pathology is highest
max_path_ipsi <- ipsi_df %>%
  group_by(region) %>%
  summarise(max_path = max(path),
            time_at_max = mpi[which.max(path)])

# pull out contralateral regions
contra_df <- plot_df %>%
  filter(str_detect(region, "^c"))

# for each contra region, identify time point when pathology is highest
max_path_contra <- contra_df %>%
  group_by(region) %>%
  summarise(max_path = max(path),
            time_at_max = mpi[which.max(path)])

# bind ipsi and contra sides
max_path_full <- rbind(max_path_ipsi, max_path_contra)

# add a new column based on the prefix
df.max.rl <- max_path_full %>%
  mutate(side = ifelse(substr(region, 1, 1) == "i", "right", "left")) %>%
  select(region, side, everything())
# Remove the prefix from the names
df.max.rl$region <- sub("^i|^c", "", df.max.rl$region)
# write out file
write.csv(df.max.rl,paste0(savedir,measure,'_max_path_tp_rl.csv'), row.names = FALSE)
