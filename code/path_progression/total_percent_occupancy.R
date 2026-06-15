# Get total pathology occupancy per mouse
rm(list=setdiff(ls(),'params'))
basedir <- params$basedir

savedir <- paste0(params$opdir, 'path_progression/')
dir.create(savedir,recursive = T)

# load formatted pathology data to get list of region names included in diffusion modeling
load(paste0(params$opdir,'processed/total_pathdata.RData'))

# convert the region list into a data frame for easy joining/filtering
region_filter <- data.frame(
  code = region.names,
  stringsAsFactors = FALSE
) %>%
  mutate(
    side   = ifelse(str_sub(code, 1, 1) == "i", "right", "left"),
    parent = str_sub(code, 2)
  )

# load and format pathology datasets
path_raw0.5_9 <- read.csv(paste0(basedir,'data/Compiled_0.5-9MPI_Pathology_Parent_Measures.csv'))
path_raw1 <- read.csv(paste0(basedir, 'data/Compiled 1MPI Pathology Measures.csv'))
path_raw0.1 <- read.csv(paste0(basedir, 'data/Compiled 0.1MPI Pathology Measures.csv'))
path_raw0.2 <- read.csv(paste0(basedir, 'data/Compiled 0.2MPI Pathology Measures.csv'))
path_raw0.3 <- read.csv(paste0(basedir, 'data/Compiled 0.3MPI Pathology Measures.csv'))
path_files <- list(path_raw0.1, path_raw0.2, path_raw0.3, path_raw1, path_raw0.5_9)
new_colnames <- c('mouse', 'Region.ID', 'mpi', 'side', 'parent', 'daughter', 'genotype', 'treatment', 'sex', 'marker', 'total_path', 'cell_body_path', 'neurite_path', 'quant_cell_body_path')

# set column names
path_files <- lapply(path_files, function(df) {
  names(df) <- new_colnames
  return(df)
})

# bind all the data frames together
path <- do.call(rbind, path_files)

# load and format region size data
region_size_0.1 <- read.csv(paste0(basedir,'data/region_areas/0.1MPI Log10 Parent Region Area.csv'), row.names=1)
region_size_0.2 <- read.csv(paste0(basedir,'data/region_areas/0.2MPI Log10 Parent Region Area.csv'), row.names=1)
region_size_0.3 <- read.csv(paste0(basedir,'data/region_areas/0.3MPI Log10 Parent Region Area.csv'), row.names=1)
region_size_0.5 <- read.csv(paste0(basedir,'data/region_areas/0.5MPI Log10 Parent Region Area.csv'), row.names=1)
region_size_1 <- read.csv(paste0(basedir,'data/region_areas/1MPI Log10 Parent Region Area.csv'), row.names=1)
region_size_3 <- read.csv(paste0(basedir,'data/region_areas/3MPI Log10 Parent Region Area.csv'), row.names=1)
region_size_6 <- read.csv(paste0(basedir,'data/region_areas/6MPI Log10 Parent Region Area.csv'), row.names=1)
region_size_9 <- read.csv(paste0(basedir,'data/region_areas/9MPI Log10 Parent Region Area.csv'), row.names=1)

size_files <- list(region_size_0.1, region_size_0.2, region_size_0.3, region_size_0.5, region_size_1, region_size_3, region_size_6, region_size_9)
size_colnames <- c('mouse', 'Region.ID', 'mpi', 'side', 'parent', 'daughter', 'genotype', 'treatment', 'sex', 'marker', 'forRemove', 'batch', 'region_area', 'log_region_area')

# set column names
region_size_files <- lapply(size_files, function(df) {
  names(df) <- size_colnames
  return(df)
})

# bind all the data frames together
region_sizes <- do.call(rbind, region_size_files)

# Join the two data frames on the key columns, adding the two region size columns
path_size <- left_join(
  path, 
  region_sizes %>% select(mouse, mpi, side, parent, region_area, log_region_area),
  by = c("mouse", "mpi", "side", "parent")
)

# identify rows where daughter names contain SNc or SNr
update_SNr <- grepl("SNr", path_size$daughter)
update_SNc <- grepl("SNc", path_size$daughter)
# select daughter rows containing SNr or SNc
#rows <- path_size[update_SNr, ]
#rows <- path_size[update_SNc, ]
# update parent row names
path_size$parent[update_SNr] <- "SNr"
path_size$parent[update_SNc] <- "SNc"

# filter df to retain only rows matching parent and side from regions_df
path_size <- path_size %>%
  semi_join(region_filter, by = c("parent", "side"))

# calculate area occupied by pathology for each region
path_size$neurite_area_occupied <- (path_size$neurite_path / 100) * path_size$region_area
path_size$cell_body_area_occupied <- (path_size$cell_body_path / 100) * path_size$region_area
path_size$total_area_occupied <- (path_size$total_path / 100) * path_size$region_area

# sum across all regions for each mouse at each time point and side
path_summary <- path_size %>%
  group_by(mouse, mpi, side) %>%
  summarise(
    total_neurite_area_occupied = sum(neurite_area_occupied, na.rm = TRUE),
    total_cell_body_area_occupied = sum(cell_body_area_occupied, na.rm = TRUE),
    total_path_area_occupied = sum(total_area_occupied, na.rm = TRUE),
    total_brain_area = sum(region_area, na.rm = TRUE)
  ) %>%
  ungroup()

# calculate percent area occupied in the whole brain
path_summary <- path_summary %>%
  mutate(
    log_percent_neurite_area = log10((total_neurite_area_occupied / total_brain_area) * 100),
    log_percent_cell_body_area = log10((total_cell_body_area_occupied / total_brain_area) * 100),
    log_percent_total_path_area = log10((total_path_area_occupied / total_brain_area) * 100)
  )

# calculate mean and standard error for plotting
hemi_plot_df <- path_summary %>%
  group_by(mpi, side) %>%
  summarise(
    mean = mean(log_percent_total_path_area),
    se = sd(log_percent_total_path_area) / sqrt(n())
  )

# creat plot comparing total pathology in each hemisphere across time
t <- ggplot(hemi_plot_df, aes(x = mpi, y = mean, group = side, color = side)) +
  geom_line(size = 0.25) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = side), color = NA, alpha = 0.5, linewidth = 0.4) +
  labs(x = "MPI",
       y = "Log(total pathology)") +
  scale_x_continuous(breaks = seq(0, 8, by = 2)) +
  ylim(-5, -1) +
  theme_cowplot(6)
t

# export plot
ggsave(t,filename = paste0(savedir,'hemi_progression_comparison.pdf'),
       units = 'cm',height = 4.15,width = 6,useDingbats=FALSE)

# export .csv containing the plotted values
csv_hemi_plot <- hemi_plot_df %>%
  dplyr::rename(
    mean_log_summed_path = mean,
    se_log_summed_path = se
  )
write.table(csv_hemi_plot, paste0(savedir,'2C_total_percent_path_by_hemi.csv'), sep=',', row.names = FALSE, col.names = TRUE)

# create plot comparing cell body, neurite, and total pathology summed across both hemispheres
both_hemis_summary <- path_size %>%
  group_by(mouse, mpi) %>%
  summarise(
    total_neurite_area_occupied = sum(neurite_area_occupied, na.rm = TRUE),
    total_cell_body_area_occupied = sum(cell_body_area_occupied, na.rm = TRUE),
    total_path_area_occupied = sum(total_area_occupied, na.rm = TRUE),
    total_brain_area = sum(region_area, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    log_percent_neurite_area = log10((total_neurite_area_occupied / total_brain_area) * 100),
    log_percent_cell_body_area = log10((total_cell_body_area_occupied / total_brain_area) * 100),
    log_percent_total_path_area = log10((total_path_area_occupied / total_brain_area) * 100)
  ) %>%
  select(mouse, mpi, log_percent_neurite_area, log_percent_cell_body_area, log_percent_total_path_area) %>%
  pivot_longer(cols = -(1:2), names_to = "path_type", values_to = "log_percent_path_occupancy")

# calculate mean and standard error for plotting
path_plot_df <- both_hemis_summary %>%
  group_by(mpi, path_type) %>%
  summarise(
    mean = mean(log_percent_path_occupancy),
    se = sd(log_percent_path_occupancy) / sqrt(n())
  )
path_plot_df$mean[path_plot_df$mean == -Inf] <- NA # convert -Inf to NA for plotting

q <- ggplot(path_plot_df, aes(x = mpi, y = mean, group = path_type, color = path_type)) +
  geom_line(size = 0.25) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = path_type), color = NA, alpha = 0.5, linewidth = 0.4) +
  scale_color_discrete(name = "path", labels = c("cb", "neur", "total")) +
  scale_fill_discrete(name = "path", labels = c("cb", "neur", "total")) +
  labs(x = "MPI",
       y = "Log(total pathology)") +
  scale_x_continuous(breaks = seq(0, 8, by = 2)) +
  ylim(-5, -1) +
  theme_cowplot(6)
q

# export plot
ggsave(q,filename = paste0(savedir,'path_progression_comparison.pdf'),
       units = 'cm',height = 4.15,width = 6,useDingbats=FALSE)

# export .csv containing the plotted values
csv_plot_df <- path_plot_df %>%
  dplyr::rename(
    mean_log_summed_path = mean,
    se_log_summed_path = se
  ) %>%
  mutate(
    path_type = case_when(
      path_type == "log_percent_cell_body_area"  ~ "cell_body",
      path_type == "log_percent_neurite_area"    ~ "neurite",
      path_type == "log_percent_total_path_area" ~ "total",
      TRUE ~ path_type  # keep any other values unchanged
    )
  )
write.table(csv_plot_df, paste0(savedir,'2B_total_percent_path_by_type.csv'), sep=',', row.names = FALSE, col.names = TRUE)

