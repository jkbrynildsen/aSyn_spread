rm(list=setdiff(ls(),c('params', 'measure', 'injection.site')))
basedir <- params$basedir

library(data.table)
library(janitor)

## Read in gene expression data, format so that genes are columns and rows contain average gene expression for each gene

gene <- read.csv(paste0(basedir, 'data/Pangea_data2.csv')) # load in gene expression data, updated by Tommy 3/26/24
region_names <- names(gene)
gene <- rbind(region_names, gene)

genet <- transpose(gene)

genett <- genet %>% 
  row_to_names(row_number = 1) %>%
  mutate_at(c(2:19782), as.numeric)

genettt <- genett %>% 
  separate_wider_delim(X, delim = "_", names = c("Region", "LayerMouse")) %>%
  select(-LayerMouse)

# define region names and replacement region names for aligning gene expression with vulnerability
orig_name <- c("BLAa", "BLAp", "BLAv", "BMAp", "CEAc", "CEAl", "CEAm", "COAa", "COApl", "COApm", "DG.sg", "DG.po", "EPd", "EPv", "LSc", "LSr", "TTd", "TTv", "VISl.VISli", "VISrl.VISal", "SSp.tr", "SSp.bfd")
replacement_name <- c("BLA", "BLA", "BLA", "BMA", "CEA", "CEA", "CEA", "COA", "COAp", "COAp", "DG", "DG", "EP", "EP", "LS", "LS", "TT", "TT", "VISli", "VISal", "SSp-tr", "SSp-bfd")

genettt$Region <- replace(genettt$Region, genettt$Region %in% orig_name, replacement_name)

genetttt <- genettt %>% 
  group_by(Region) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

gene_final <- genetttt %>%
  remove_rownames %>%
  column_to_rownames(var = "Region")

# read in vulnerability file, compare region labels
vuln.dir <- paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
df.vuln <- read.csv(paste0(vuln.dir,measure,'_vulnerability_bidirectional_hemiaverage.csv'),row.names = 1)

gene.regions <- rownames(gene_final)
vuln.regions <- rownames(df.vuln)

# write out avg gene expression dataset
write.csv(gene_final, paste0(params$opdir,'processed/avg_Pangea_exp.csv'), row.names = TRUE)

norm.fun <- function(x){
  # https://www.pnas.org/content/113/5/1435
  # robust to outliers, puts each gene on same scale relative to its own expression across brain
  return((1+exp(-scale(x,center = T)))^-1)
}
 
norm_gene_final <- as.data.frame(apply(gene_final,2,norm.fun))# normalize each gene in the dataset
rownames(norm_gene_final) <- rownames(gene_final)

Snca <- as.data.frame(norm_gene_final$Snca)
rownames(Snca) <- rownames(norm_gene_final)

write.csv(norm_gene_final, paste0(params$opdir,'processed/norm_avg_Pangea_exp.csv'), row.names = TRUE) # write out normalized gene expression matrix
write.csv(Snca, paste0(params$opdir,'processed/Snca_exp.csv'), row.names = TRUE) # write out vector of Snca expression for use with diffusion model
 