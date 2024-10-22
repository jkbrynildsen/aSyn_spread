### Group regions according to the time at which peak pathology levels are reached and plot pathology levels across time for each group

### Format directories and data ####
rm(list=setdiff(ls(),c('params','measure','injection.site')))
opdir <- params$opdir
basedir <- params$basedir
setwd(basedir)

savedir <- paste0(params$opdir, 'path_progression/')
dir.create(savedir,recursive = T)

source('code/misc/plottingfxns.R')

load(paste0(opdir,'processed/', measure, 'data.RData'))  # load path data and ROI names

### Plot log pathology for each region across time ####
## set up data frame for plotting regional log pathology across time
plot_df <- plot_format(path.data)

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

### Plot mean pathology across time for the group of regions that peak at each time point ####
# get mean pathology across regions for each time point
overall_avg <- plot_df %>%
  group_by(ord_tp) %>%
  summarise(mean_path = mean(path, na.rm = TRUE))

### for each time point, pull out the regions that peak at that time point
regions_peak9 <- tp_max(9, max_path_full, plot_df)
regions_peak6 <- tp_max(6, max_path_full, plot_df)
regions_peak3 <- tp_max(3, max_path_full, plot_df)
regions_peak1 <- tp_max(1, max_path_full, plot_df)
regions_peak0.5 <- tp_max(0.5, max_path_full, plot_df)
regions_peak0.3 <- tp_max(0.3, max_path_full, plot_df)
regions_peak0.2 <- tp_max(0.2, max_path_full, plot_df)
regions_peak0.1 <- tp_max(0.1, max_path_full, plot_df)

## concatenate into data frame
all_peak <- rbind(regions_peak0.1, regions_peak0.2, regions_peak0.3, regions_peak0.5, regions_peak1, regions_peak3, regions_peak6, regions_peak9)
# average across the early time points (0.1 - 0.5 MPI)
all_peak <- all_peak %>%
  mutate(peak_tp = ifelse(peak_tp %in% c(0.1, 0.2, 0.3, 0.5), 1, peak_tp))

summary_data <- all_peak %>%
  group_by(mpi, peak_tp) %>%
  summarise(
    mean = mean(path, na.rm = TRUE),
    se = sd(path, na.rm = TRUE) / sqrt(sum(!is.na(path)))
  )
summary_data$peak_tp <- factor(summary_data$peak_tp, levels = c("1", "3", "6", "9"))

t <- ggplot(summary_data, aes(x = mpi, y = mean, color = peak_tp, group = peak_tp)) +
  geom_line() +
  scale_color_manual(values = c("#FFDDA6", "#FDC152", "#5ECBF4", "#34598E")) + #McDonald's color scheme
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = peak_tp), color = NA, alpha = 0.5, linewidth = 0.4) +
  scale_fill_manual(values = c("#FCFFD3", "#FDC152", "#5ECBF4", "#34598E")) +
  labs(x = "MPI",
       y = "Mean (Log pathology)") +
  scale_x_continuous(breaks = seq(0, 8, by = 2)) +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        legend.position = "none") +  # Remove both legends
  guides(color = "none", fill = "none") +  # Suppress color and fill legends
  theme_cowplot(6)

t
ggsave(t,filename = paste0(savedir,measure,'_mean_log_regional_path_peak.pdf'),
       units = 'cm',height = 3.75,width = 6,useDingbats=FALSE) # used width = 6.5 with legend
