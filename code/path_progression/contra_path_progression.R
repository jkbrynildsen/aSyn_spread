### This script generates a plot that shows the log total pathology, cell body pathology, and neurite pathology summed across all contralateral
### brain regions across time.

rm(list=setdiff(ls(),c('params','measure','injection.site')))
opdir <- params$opdir
basedir <- params$basedir
savedir <- paste0(opdir, 'path_progression/')
dir.create(savedir,recursive = T)
setwd(basedir)

# load data for neurite, cell body and total path; pull out contralateral regions
load(paste0(opdir,'processed/neurite_pathdata.RData'))  # load neurite path data and ROI names
neurite_path <- path.data
total_neurite_path <- neurite_path %>%
  select(starts_with("c")) %>%
  mutate(neurite_total = rowSums(., na.rm=TRUE)) # sum across all regions to get total pathology for each mouse
load(paste0(opdir,'processed/cell_body_pathdata.RData'))  # load cell body path data and ROI names
cb_path <- path.data
total_cb_path <- cb_path %>%
  select(starts_with("c")) %>%
  mutate(cb_total = rowSums(., na.rm=TRUE))
load(paste0(opdir,'processed/total_pathdata.RData'))  # load total path data and ROI names
total_path <- path.data
total_path <- total_path %>%
  select(starts_with("c")) %>%
  mutate(total = rowSums(., na.rm=TRUE))

## set up data frame for plotting. contains total (summed) pathology for cell body, neurite and total pathology for each mouse and time point
plot_df <- as.data.frame(cbind(mouse=cb_path$mouse, mpi=cb_path$mpi, cell_body=total_cb_path$cb_total, neurite=total_neurite_path$neurite_total, total = total_path$total))
plot_df <- plot_df %>%
  pivot_longer(!c(mouse, mpi), names_to = 'path_type', values_to = 'total_path')
plot_df$total_path <- as.numeric(plot_df$total_path)
plot_df$mpi <- as.numeric(plot_df$mpi)

# log transform path data
plot_df$log_total_path <- log(plot_df$total_path,base=10)
plot_df$log_total_path[plot_df$log_total_path == -Inf] <- NA

# calculate mean and standard error for plotting
summary_data <- plot_df %>%
  group_by(mpi, path_type) %>%
  summarise(
    mean = mean(log_total_path),
    se = sd(log_total_path) / sqrt(n())
  )

t <- ggplot(summary_data, aes(x = mpi, y = mean, group = path_type, color = path_type)) +
  geom_line(size = 0.25) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = path_type), color = NA, alpha = 0.5, linewidth = 0.4) +
  labs(x = "MPI",
       y = "Log(total pathology)") +
  scale_x_continuous(breaks = seq(0, 8, by = 2)) +
  ylim(-2.5, 1.5) +
  theme_cowplot(6)
t

ggsave(t,filename = paste0(savedir,'contra_neurite_v_cell_body_v_total_progression.pdf'),
       units = 'cm',height = 3.75,width = 6,useDingbats=FALSE)
