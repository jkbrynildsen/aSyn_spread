rm(list=ls())

basedir <- '~/Library/CloudStorage/Box-Box/PD_Project/Quant_aSyn_pathology/'

setwd(basedir)

params <- list(basedir=basedir,
               injection.site = 'iCP', # injection site
               tps=c(0.1, 0.2, 0.3, 0.5, 1, 3, 6, 9), # define time points post injection
               measures=c('total_path','cell_body_path','neurite_path'),
               c.min = 1e-10, # minimum for time constant
               c.max = 0.2, # maximum for time constant
               c.n = 300) # number of time constants to try

injection.site <- params$injection.site

params$opdir <- paste0('aSynDiffusion101624_Inject',paste0(params$injection.site,collapse='-'),'_CMax',params$c.max,'/')
dir.create(params$opdir,recursive = T)

######################################
### Load packages for all analyses ###
######################################

source('code/misc/packages.R')
use_python("/Users/katebrynildsen/anaconda3/envs/r-reticulate/bin/python3.12") # choose version of python to be used by reticulate package

#############################################
### Process pathology and structural data ###
#############################################

source('code/process/process_pathology.R')
source('code/process/process_struct.R')

#######################################################
### Visualize levels of total pathology across time ###
#######################################################

source('code/path_progression/ipsi_path_progression.R') # plot total path across time for ipsi regions
source('code/path_progression/contra_path_progression.R') # plot total path across time for contra regions
for(measure in params$measures){
  source('code/path_progression/plot_path_peak_tp.R') # plot pathology across time, with regions grouped according to the time at which peak pathology is reached
  source('code/path_progression/path_peak_tp_csv.R') # generate a spreadsheet containing the peak time point associated with each region for plotting heatmaps
}

#######################
### Diffusion model ###
#######################

injection.site <- params$injection.site

for(measure in params$measures){
  source('code/diffmodel/analyze_retrogradespread.R')
  source('code/diffmodel/analyzebidirectionalspread.R') # run this to initialize parameters for optimization
  source('code/diffmodel/plot_fit.R')
  
  source('code/diffmodel/optim_bidirectionalspread_onelm.R') # bidirectional diffusion model with different a single linear model weighting anterograde and retrograde
  source('code/diffmodel/plot_bidirectionalonelmfit.R')
}

###################################################
### Get euclidean distances between ABA regions ###
###################################################

source('code/aba/process_ontology.R')
source('code/aba/atlas_structures.R') # required for model comparison analyses

#################################
### Model validation analyses ###
#################################
# these are very time consuming

# compare model fits based on random vs actual seed regions      
injection.site <- params$injection.site
# seed specificity
for(measure in params$measures){
  source('code/nullmodels/seedspec_onelm.R')
  source('code/nullmodels/plot_seedspec.R')
}

# compare model types: Euclidean, anterograde, retrograde, and bidirectional
for(measure in params$measures){
  injection.site <- params$injection.site; 
  source('code/modelcomparison/modelcomparison_traintest.R')
  source('code/modelcomparison/plot_modelcomparison_testset.R')
}

##################################################################################################
### Assessing the relationship between gene expression and regional vulnerability to pathology ###
##################################################################################################

measure <- 'total_path'
source('code/process/process_PANGEA_gene_exp.R') # process gene expression data
source('code/genes_vulnerability/visualize_data_distributions.R') # visualize gene expression distributions before and after normalization
source('code/genes_vulnerability/vulnerability_hemi_comparison.R') # compare model residuals across hemispheres and time points
source('code/genes_vulnerability/genes_vulnerability_bidirectional.R') # compute correlations between gene expression and regional vulnerability

######################################################
### Predicting pathology from alternate seed sites ###
######################################################

# generate structural matrix that includes all regions in the pathology dataset
source('code/process/process_struct_full.R')
# create list of alternate seed regions to test
alt_seeds=list('iMOB','iMOs',c('iMS', 'cMS', 'iNDB', 'cNDB'), 'iVAL', c('iCA1', 'iCA3', 'iDG'), 'iSNc', 'iPPN', 'iDN', 'iCP') # use iCP seed to check fits
# obtain predicted log pathology values for each alternate seed region
for(seed in alt_seeds){
  for(measure in params$measures){
    injection.site <- seed
    source('code/diffmodel/bidirectionalonelmfit_alt_seeds_no_resid.R')
  }
}  
