######################################################
## Author: James Hay & Kath O'Reilly
## Date: 10 Jun 2025
## Summary: run the serosolver simulations applied for norovirus


## Load in libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(doParallel)
library(coda)
library(tidyverse)
library(cowplot)
library(patchwork)
library(serosolver)
library(RColorBrewer)
## Load serosolver locally, or load as proper package after installation
#Rcpp::compileAttributes("~/Documents/GitHub/serosolver")
#devtools::document("~/Documents/GitHub/serosolver")
#devtools::load_all("~/Documents/GitHub/serosolver")
#install.packages("~/Documents/GitHub/serosolver", repos=NULL, type="source")
#devtools::install_github("seroanalytics/serosolver")
#library(serosolver)

## Set working directly locally and create a directory to store the MCMC chains in
setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")
source("code/simulate_priors.R")
# load in csv for script options
read_in <- read_csv("serosolver_norovirus_input_JH.csv")

## run_number <- 3

run_name <- read_in$output_name[run_number] #"sim_avid_real"
main_wd <- paste0("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/local_data_exponential/",read_in$sim_type[run_number]) #norovirus_true"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/",run_name)
save_wd_results <- paste0("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2","/results_exponential/",run_name)
save_wd_old <- paste0("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/old/local_data/",read_in$sim_type[run_number],"/figures/",run_name)
data_folder = "~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/Data/"

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

## Simulation parameters
mcmc_pars <- c("save_block"=1000,
               "thin"=500,
               "thin_inf_hist"=3500,
               "iterations"=700000,
               "adaptive_iterations"=300000)

#mcmc_pars <- c("save_block"=100,
#               "thin"=5,
#               "thin_inf_hist"=25,
#               "iterations"=5000,
#               "adaptive_iterations"=2000)

cart_data <- read_in$cart[run_number]  #"Debbink" #"Kendra" #
# analysis level
analysis <- read_in$data[run_number]

# load data
antibody_data <- read_csv(paste0(data_folder,"titre_dat_norovirus_age_avidity.csv"))
antibody_data <- as.data.frame(antibody_data)
antibody_data_keep <- antibody_data
antibody_data <- antibody_data %>% select(-age,-Date) 

theme_use <- theme_classic() + theme(axis.title = element_text(size=8),
                                     axis.text = element_text(size=6),
                                     legend.text = element_text(size=6),
                                     legend.title = element_text(size=6),plot.tag=element_text(face="bold",size=10),
                                     plot.margin=margin(5,5,5,5),
                                     legend.margin=margin(0,0,0,0),
                                     legend.box.margin=margin(-10,-10,-10,-10)) # Reduce spacing around plot and legend

strains <- c("Grimsby 1996","Farmington Hills\n 2002","Hunter 2004", "Den Haag 2006","New Orleans 2009","Sydney\n 2012")
strain_colors_pal <- brewer.pal(n = 9, name = 'Set1')
strain_colors <- c("Grimsby 1996"="burlywood1",
                   "Farmington Hills\n 2002"=strain_colors_pal[1],
                   "Hunter 2004"=strain_colors_pal[3], 
                   "Den Haag 2006"=strain_colors_pal[4],
                   "New Orleans 2009"=strain_colors_pal[8],
                   "Sydney\n 2012"=strain_colors_pal[9])
strain_colors_short <- c("FH-2002"=strain_colors_pal[1],
                         "DH-2006"=strain_colors_pal[4],
                         "NO-2009"=strain_colors_pal[8],
                         "SY-2012"=strain_colors_pal[9])

# *** load some simulated data ***
#antibody_data <- read_csv(paste0("/Users/lsh1603970/Documents/GitHub/serosolver-norovirus-eng-serology-v2/local_data/norovirus_test/figures/sim_two_biomarker/","sim_two_biomarker_titre_data.csv"))
#antibody_data <- as.data.frame(antibody_data)
# further data options
if(analysis == "ic_50"){
  # use only biomarker_group 1
  antibody_data <- antibody_data %>% filter(biomarker_group == 1)
  
  # ic_50: divide by 5 and log2 the data, assuming 5 represents <10
  antibody_data$measurement <- log2(antibody_data$measurement/5)
  ## Set lowest measurement to 1 (log2(10/5))
  antibody_data[antibody_data$measurement == 0,"measurement"] <- 1
  ## Shift observations so that 0 corresponds to the lowest measurement
  antibody_data$measurement <- antibody_data$measurement - 1
  print(paste0("Dataset size is ",dim(antibody_data)[1]))
  # make sure we have normalised data
  p1 <- ggplot(antibody_data %>% filter(biomarker_group==1),aes(x=measurement)) + 
    geom_histogram() +
    #scale_x_continuous(trans='log2') +
    facet_wrap(~biomarker_id)
  ggsave(paste0(save_wd,"/data_histogram.pdf"),p1,height=5,width=8,units="in",dpi=300)
  
}
if(analysis == "both"){
  # ic_50: divide by 5 and log2 the data
  antibody_data$measurement[antibody_data$biomarker_group == 1] <- log2(antibody_data$measurement[antibody_data$biomarker_group == 1]/5)
  
  
  ## Set lowest measurement to 1 (log2(10/5)) and shift observations so that 0 corresponds to the lowest measurement i.e., subtract 1
  antibody_data[antibody_data$measurement == 0 & antibody_data$biomarker_group == 1,"measurement"] <- 1
  antibody_data[antibody_data$biomarker_group == 1,"measurement"] <- antibody_data[antibody_data$biomarker_group == 1,"measurement"] - 1
  
  ## Set lowest measurement to lowest value and shift so that this corresponds to 0
  antibody_data[antibody_data$measurement == 0 & antibody_data$biomarker_group == 2,"measurement"] <- 0.45
  antibody_data[antibody_data$biomarker_group == 2,"measurement"] <- antibody_data[antibody_data$biomarker_group == 2,"measurement"] - 0.45
  
  print(paste0("Dataset size is ",dim(antibody_data)[1]))
  
  p1 <- ggplot(antibody_data,aes(x=measurement)) + 
    geom_histogram() +
    facet_grid(cols=vars(biomarker_id),rows=vars(biomarker_group)) +
    xlab("Measurement (log2)")
  ggsave(paste0(save_wd,"/data_histogram.pdf"),p1,height=5,width=8,units="in",dpi=300)
  
  # plot data by ID - I have long format and want it to be wide
  abd <- spread(antibody_data,biomarker_group,measurement)
  p2 <- ggplot(abd,aes(x=`1`,y=`2`,group=sample_time,color = sample_time)) + geom_point() + 
    facet_wrap(~biomarker_id) + ylab("Avidity") + xlab("IC_50 (log2)")
  ggsave(paste0(save_wd,"/data_dual_measurements.pdf"),p2,height=4.5,width=6,units="in",dpi=300)
}
if(analysis == "avidity"){
  # use only biomarker_group 2
  antibody_data <- antibody_data %>% filter(biomarker_group == 2)
  ## Set lowest measurement to lowest value and shift so that this corresponds to 0
  antibody_data[antibody_data$measurement == 0 & antibody_data$biomarker_group == 2,"measurement"] <- 0.45
  antibody_data[antibody_data$biomarker_group == 2,"measurement"] <- antibody_data[antibody_data$biomarker_group == 2,"measurement"] - 0.45
  # then need to make biomarker_group == 1
  antibody_data$biomarker_group <- 2
  # avidity: divide by 5 and log2 the data
  #antibody_data$measurement[antibody_data$biomarker_group == 1] <- log2(antibody_data$measurement[antibody_data$biomarker_group == 1]/5)
  
  print(paste0("Dataset size is ",dim(antibody_data)[1]))
  p1 <- ggplot(antibody_data %>% filter(biomarker_group==2),aes(x=measurement)) + 
    geom_histogram() +
    #scale_x_continuous(trans='log2') +
    facet_wrap(~biomarker_id)
  antibody_data$biomarker_group <- 1
  ggsave(paste0(save_wd,"/data_avidity_histogram.pdf"),p1,height=5,width=8,units="in",dpi=300)
  
}
# *** option ***
# remove avidity and add again ic50
# table(antibody_data$biomarker_group)
# antibody_data <- antibody_data %>% filter(biomarker_group==1) 
# tmp <- antibody_data
# tmp$biomarker_group <- 2
# antibody_data <- rbind(antibody_data,tmp)
# table(antibody_data$biomarker_group)

n_obs_types <- length(unique(antibody_data$biomarker_group))

## Use observation model with possible false positives
obs_type_dists <- 3 #c(1,1)
if(analysis == "both"){
obs_type_dists <- c(3,3)
}

#n_obs_types <- 2 ## Number of observation types (eg. antibody titer and binding avidity)
#obs_type_dists <- c(1,2) ## Vector of observation models for each observation type -- use 1 for discretized normal and 2 for continuous, truncated normal. I assume that titers are discrete and avidity is continuous
# n_indivs <- 250 ## Number of individuals for the simulation
# n_samps <- 2 ## Number of samples per person
# repeats <- 1 ## Number of repeat measurements per variant/sample combination
# samp_min <- 2009 ## First sample year
# samp_max <- 2012 ## Final sample year
# year_min <- 2000 ## First year of possible circulation (ie. time 0 of the simulation)
# year_max <- 2012 ## Final year of possible circulation
# age_min <- 2 ## Age minimum and maximum in years, simulated from a uniform distribution
# age_max <- 10

################################ 
## IMPORTANT
################################ 
## Set to FALSE to remove the measurement offset parameters
## Set to TRUE to assume that each antigen/virus has a measurement offset parameter for each observation type
use_measurement_offset_parameters <- read_in$offset_params[run_number] #TRUE #

## Viruses and times for samples

sampled_viruses <- unique(antibody_data$biomarker_id) #c(2000,2002,2006,2009,2012)
sampling_times <- unique(antibody_data$sample_time) #seq(samp_min, samp_max, by=1)

################################ 
## IMPORTANT
################################ 

if(cart_data == "debbink"){
  # read in this file
  antigenic_map <- read_csv(paste0(data_folder,"/","antigenic_map_hand.csv"))
  p1 <- ggplot(antigenic_map,aes(x=x_coord,y=y_coord,col=inf_times,label=inf_times)) + geom_line() +
      geom_point(,cex=2) +
      geom_label() +
      #geom_point(data=antigenic_coords,aes(x=X,y=Y,col=inf_times),pch=0) +
      ylim(-3,4)
  p1
}
if(cart_data == "kendra"){
  # read in this other file
  # and additional changes needed
  #antigenic_coords_path <- system.file("extdata", "fonville_map_approx.csv", package = "serosolver")
  antigenic_coords <- read.csv(file = paste0(data_folder,"antigenic_map_kendra_edit.csv"), stringsAsFactors=FALSE)
  print(head(antigenic_coords))
  antigenic_coords$Strain <- c(antigenic_coords$inf_times)
  
  ggplot(antigenic_coords,aes(x=x_coord,y=y_coord,label=Strain)) + geom_point() + 
    geom_path() + geom_text()
  
  # do by hand
  antigenic_map <- data.frame(x_coord=NA,y_coord=NA,inf_times=c(2001:2012))
  oo <- match(antigenic_map$inf_times,antigenic_coords$inf_times)
  antigenic_map$x_coord[!is.na(oo)] <- antigenic_coords$x_coord[oo[!is.na(oo)]]
  antigenic_map$y_coord[!is.na(oo)] <- antigenic_coords$y_coord[oo[!is.na(oo)]]
  # 2003
  antigenic_map$x_coord[antigenic_map$inf_times=="2003"] <- antigenic_map$x_coord[antigenic_map$inf_times=="2002"] + (antigenic_map$x_coord[antigenic_map$inf_times=="2004"] - antigenic_map$x_coord[antigenic_map$inf_times=="2002"])/2
  antigenic_map$y_coord[antigenic_map$inf_times=="2003"] <- antigenic_map$y_coord[antigenic_map$inf_times=="2002"] + (antigenic_map$y_coord[antigenic_map$inf_times=="2004"] - antigenic_map$y_coord[antigenic_map$inf_times=="2002"])/2
  
  # 2008
  antigenic_map$x_coord[antigenic_map$inf_times=="2008"] <- antigenic_map$x_coord[antigenic_map$inf_times=="2007"] + (antigenic_map$x_coord[antigenic_map$inf_times=="2009"] - antigenic_map$x_coord[antigenic_map$inf_times=="2007"])/2
  antigenic_map$y_coord[antigenic_map$inf_times=="2008"] <- antigenic_map$y_coord[antigenic_map$inf_times=="2007"] + (antigenic_map$y_coord[antigenic_map$inf_times=="2009"] - antigenic_map$y_coord[antigenic_map$inf_times=="2007"])/2
  
  # 2010 & 2011
  antigenic_map$x_coord[antigenic_map$inf_times=="2010"] <- antigenic_map$x_coord[antigenic_map$inf_times=="2009"] + (antigenic_map$x_coord[antigenic_map$inf_times=="2012"] - antigenic_map$x_coord[antigenic_map$inf_times=="2009"])/3
  antigenic_map$y_coord[antigenic_map$inf_times=="2010"] <- antigenic_map$y_coord[antigenic_map$inf_times=="2009"] + (antigenic_map$y_coord[antigenic_map$inf_times=="2012"] - antigenic_map$y_coord[antigenic_map$inf_times=="2009"])/3
  antigenic_map$x_coord[antigenic_map$inf_times=="2011"] <- antigenic_map$x_coord[antigenic_map$inf_times=="2009"] + ((antigenic_map$x_coord[antigenic_map$inf_times=="2012"] - antigenic_map$x_coord[antigenic_map$inf_times=="2009"])/3)*2
  antigenic_map$y_coord[antigenic_map$inf_times=="2011"] <- antigenic_map$y_coord[antigenic_map$inf_times=="2009"] +((antigenic_map$y_coord[antigenic_map$inf_times=="2012"] - antigenic_map$y_coord[antigenic_map$inf_times=="2009"])/3)*2
  
  ggplot(antigenic_map,aes(x=x_coord,y=y_coord,label=inf_times)) + geom_point() +
    geom_label() + geom_path()
  
  ggplot(antigenic_coords,aes(x=x_coord,y=y_coord,label=Strain)) + geom_point() + 
    geom_path() + geom_text()
  
}

 ## This doubles the antigenic map, one entry for each biomarker_group (i.e., one for titers one for avidity)
antigenic_map <- setup_antigenic_map(antigenic_map,n_biomarker_groups=n_obs_types,unique_biomarker_groups=1:n_obs_types)[[1]]


## Get vector of times individuals can be infected
possible_exposure_times <- unique(antigenic_map$inf_times)
n_times <- length(possible_exposure_times)

## Set up parameter table
par_tab <- read.csv("par_tab.csv",stringsAsFactors=FALSE)
par_tab[par_tab$names == "fp_rate","values"] <- 0.001 ## Set false positive probabilty to 0.1%

################################ 
## IMPORTANT
################################ 
## This allows you to change prior on the probability of infection per year 
## These are the parameters of a beta distribution -- see rbeta(10000, infection_model_prior_shape1,infection_model_prior_shape2)
par_tab[par_tab$names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"),c("values")] <- c(1/2,1/2) ## Can also try c(1,1), or something informative. Just the parameters of a beta distribution which acts as the prior on the per-time attack rate.
################################ 

## Extends the parameter table for multiple observation types

par_tab <- extend_par_tab_biomarker_groups(par_tab,n_obs_types)
if(use_measurement_offset_parameters){
  tmp <- add_rhos_par_tab(par_tab, sampled_viruses, n_obs_types)
  par_tab <- tmp[[1]]
  measurement_indices <- tmp[[2]]
} else {
  measurement_indices <- NULL
}
#par_tab[par_tab$biomarker_group == 2 & par_tab$names %in%c("boost_long","boost_short","obs_sd"),"values"] <- c(2.3,0.5,0.75)



################################ 
## IMPORTANT
################################ 
## Set the upper bound of observations for the two observation types (i.e., max observable value for titers and avidity)
par_tab[par_tab$biomarker_group == 1 & par_tab$names == "max_measurement","values"] <- ceiling(max(antibody_data[antibody_data$biomarker_group==1,"measurement"]))
par_tab[par_tab$biomarker_group == 2 & par_tab$names == "max_measurement","values"] <- ceiling(max(antibody_data[antibody_data$biomarker_group==2,"measurement"]))


par_tab[par_tab$biomarker_group == 1 & par_tab$names == "min_measurement","values"] <- 0
par_tab[par_tab$biomarker_group == 2 & par_tab$names == "min_measurement","values"] <- 0

## Turn off long-term response for biomarker group 2 (avidity)
if(analysis != "avidity"){
  par_tab[par_tab$biomarker_group == 2 & par_tab$names %in% c("boost_long","cr_long"),"fixed"] <- 1
  par_tab[par_tab$biomarker_group == 2 & par_tab$names %in% c("boost_long","cr_long"),"values"] <- 0
} else {
  par_tab[par_tab$biomarker_group == 1 & par_tab$names %in% c("boost_long","cr_long"),"fixed"] <- 1
  par_tab[par_tab$biomarker_group == 1 & par_tab$names %in% c("boost_long","cr_long"),"values"] <- 0
}

################################ 
## IMPORTANT
################################ 
## Set values for the model parameter priors.
## You would need to add your own priors here as new lines and values added to the final return value
## Note that to set a prior for a parameter, you need to get the right index in the par_tab table
if(analysis != "avidity"){
  prior_func <- prior_func_exponential_ic50
} else {
  prior_func <- prior_func_exponential_avidity
  
}

################################ 
## RUN SEROSOLVER
################################ 

# double check it's a dataframe
antibody_data <- as.data.frame(antibody_data)
antibody_data <- antibody_data %>% arrange(individual,sample_time,biomarker_id,biomarker_group)

## Try age-stratification on infections
par_tab[par_tab$names == "infection_model_prior_shape1","stratification"] <- "age_group"
par_tab[par_tab$names == "infection_model_prior_shape2","stratification"] <- "age_group"

demography <- antibody_data %>% select(individual,birth) %>% distinct() %>% expand_grid(time=unique(antigenic_map$inf_times)) %>% arrange(individual, time) %>% mutate(age = time - birth) %>% 
  mutate(age_group = if_else(age >= 5, "5+",if_else(age <= 0, "0", as.character(age)))) %>% as.data.frame() %>% mutate(age_group = as.numeric(as.factor(age_group))-1)%>% filter(time >= min(antibody_data$birth))
antigenic_map <- antigenic_map %>% filter(inf_times >= min(antibody_data$birth))

if(rerun_mcmc){
  res <- serosolver(par_tab, 
                    antibody_data, 
                    demographics=demography,
                    #possible_exposure_times = possible_exposure_times,
                    antigenic_map=antigenic_map,
                    prior_func=prior_func,
                    filename=paste0(chain_wd,"/",run_name),
                    n_chains=5, ## Run 3 chains
                    parallel=TRUE, ## Run in parallel
                    mcmc_pars=mcmc_pars, 
                    verbose=TRUE,
                    data_type=obs_type_dists,
                    measurement_bias= NULL,
                    exponential_waning=TRUE
  ) 

## Save plots from serosolver
res$all_diagnostics$p_thetas[[1]]
res$all_diagnostics$p_thetas[[2]]
res$all_diagnostics$p_thetas[[3]]
res$all_diagnostics$p_thetas[[4]]
res$all_diagnostics$p_thetas[[5]] # rhos

save(res,file=paste0(save_wd,"/out_res.rda"))

## Save parameter estimates
ggsave(paste0(save_wd,"/theta_traceplot.pdf"),res$all_diagnostics$p_thetas[[1]],height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/theta_densities.pdf"),res$all_diagnostics$p_thetas[[2]],height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/rho_traceplot.pdf"),res$all_diagnostics$p_thetas[[5]],height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/rho_densities.pdf"),res$all_diagnostics$p_thetas[[6]],height=7,width=8,units="in",dpi=300)
write_csv(as.data.frame(res$all_diagnostics$theta_estimates),paste0(save_wd,"/estimated_parameters.csv"))

## Save attack rate/infection history estimates
ggsave(paste0(save_wd,"/indiv_infections.pdf"),res$all_diagnostics$p_inf_hists$indiv_infections,height=4,width=6,units="in",dpi=300)
ggsave(paste0(save_wd,"/attack_rate_traceplot.pdf"),res$all_diagnostics$p_inf_hists$indiv_infections$by_time_trace[[1]],height=4,width=6,units="in",dpi=300)
ggsave(paste0(save_wd,"/attack_rate_density.pdf"),res$all_diagnostics$p_inf_hists$indiv_infections$by_time_trace[[2]],height=4,width=6,units="in",dpi=300)

write_csv(as.data.frame(res$all_diagnostics$inf_hist_estimates$by_year),paste0(save_wd,"/estimated_attack_rates.csv"))
serosolver_settings <- res$settings
} else {
 load(paste0(chain_wd,"/",run_name,"_serosolver_settings.RData"))
}

## Post-hoc plots
chains <- load_mcmc_chains(chain_wd,par_tab,
                           burnin=mcmc_pars["adaptive_iterations"],
                           estimated_only = TRUE)

## Plot parameter estimates vs. true values
estimated_par_names <- par_tab[par_tab$fixed == 0 & !(par_tab$names %like% "rho"),"names"]
estimated_par_names <- make.names(estimated_par_names,unique=TRUE)


################################ 
## IMPORTANT
################################
## The vertical lines are for the true values of the simulation.
## For real data, comment out the geom_vline
## Plot theta estimates by biomarker group
par_name_key <- c("boost_long"="Long-term boost","boost_short"="Transient boost","wane_short"="Transient waning rate",
                  "cr_long"="Long-term cross-reactivity gradient","cr_short"="Transient cross-reactivity gradient",
                  "obs_sd"="Observation error (standard deviation)","antigenic_seniority"="Antigenic seniority parameter")

tidy_theta_chain <- chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
  filter(name %in% estimated_par_names) %>%
  dplyr::rename(names_new=name,est=value) %>%
  left_join(par_tab %>% filter(fixed == 0& !(par_tab$names %like% "rho")) %>% 
              mutate(names_new=estimated_par_names) %>% 
              select(names,names_new,values,biomarker_group) %>%
              mutate(biomarker_group = as.factor(biomarker_group))) %>%
  mutate(biomarker_group_label = if_else(biomarker_group == 1, 
                                         if_else(analysis != "avidity","IC50","Avidity"),"Avidity")) %>%
  mutate(run_name=run_name,
         map=cart_data)
tidy_theta_chain$name_key <- par_name_key[tidy_theta_chain$names]
tidy_theta_chain <- tidy_theta_chain %>% mutate(label = paste0(name_key,"\n(",biomarker_group_label,")"))

label_levels <- c("Long-term boost\n(IC50)", 
                  "Long-term boost\n(Avidity)", 
                  "Transient boost\n(IC50)", 
                  "Transient boost\n(Avidity)", 
                  "Transient waning rate\n(IC50)", 
                  "Transient waning rate\n(Avidity)", 
                  "Antigenic seniority parameter\n(IC50)",
                  "Antigenic seniority parameter\n(Avidity)",
                  "Long-term cross-reactivity gradient\n(IC50)", 
                  "Long-term cross-reactivity gradient\n(Avidity)",
                  "Transient cross-reactivity gradient\n(IC50)", 
                  "Transient cross-reactivity gradient\n(Avidity)", 
                  "Observation error (standard deviation)\n(IC50)", 
                  "Observation error (standard deviation)\n(Avidity)")
tidy_theta_chain$label <- factor(tidy_theta_chain$label, levels=label_levels)

## Add in ESS and Rhat
if(rerun_mcmc){
ess_labels <- res$all_diagnostics$theta_estimates %>% as_tibble() %>% 
  select(names,ess,`Rhat point estimate`,`Rhat upper CI`) %>% mutate(label1=paste0("ESS: ", signif(ess,3), "; Rhat: ", signif(`Rhat point estimate`,3), " (upper 95% CrI: ", signif(`Rhat upper CI`,3),")"))
df <- data.frame(
  original = ess_labels$names,
  stringsAsFactors = FALSE
)
# Extract name and number (if missing, set to 0)
df$parameter <- sub("\\..*$", "", df$original)         # everything before "."
df$number <- sub(".*\\.", "", df$original)             # everything after "."
# If no dot, keep 0 instead of original value
df$number[!grepl("\\.", df$original)] <- "0"
# Convert to numeric
df$number <- as.integer(df$number)
df <- df %>% rename(names=original)
ess_labels <- ess_labels %>% left_join(df)
ess_labels$parameter <- par_name_key[ess_labels$parameter]
if(analysis != "avidity"){
ess_labels$label <- paste0(ess_labels$parameter,"\n(",if_else(ess_labels$number == 0,"IC50","Avidity"),")")
} else {
  ess_labels$label <- paste0(ess_labels$parameter,"\n(",if_else(ess_labels$number == 0,"Avidity","Avidity"),")")
  
}
ess_labels <- ess_labels %>% left_join(tidy_theta_chain %>% group_by(label) %>% summarize(max_val = max(est)))
p_traces <- ggplot(tidy_theta_chain) + 
  geom_line(aes(x=samp_no,col=as.factor(chain_no),y=est)) +
  geom_text(data=ess_labels %>% drop_na(),aes(x=median(tidy_theta_chain$samp_no),y=max_val*1.5,label=label1),size=3,vjust=1.25) +
  facet_wrap(~label, scales="free_y",ncol=2) +
  scale_color_viridis_d(name="MCMC chain") +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Value") +
  xlab("MCMC iteration")



ggsave(paste0(save_wd,"/theta_traces.pdf"),p_traces,height=7*n_obs_types,width=7,units="in",dpi=300)
}
priors <- simulate_priors(1000000)
priors$name_key <- par_name_key[priors$name]
if(analysis != "avidity"){
priors <- priors %>%
  mutate(biomarker_group_label = if_else(biomarker_group == 1, "IC50","Avidity")) 
} else {
  
  priors <- priors %>%
    mutate(biomarker_group_label = if_else(biomarker_group == 1, "Avidity","Avidity")) 
}
priors <- priors %>% mutate(label = paste0(name_key,"\n(",biomarker_group_label,")"))
priors$label <- factor(priors$label, levels=label_levels)
## Bound priors for plot
priors <- priors %>% left_join(tidy_theta_chain %>% group_by(label) %>% dplyr::summarize(min_value = min(est),max_value=max(est)))
priors <- priors %>% filter(value >= min_value*0.9 & value <= max_value*1.1)


p_theta <- ggplot(tidy_theta_chain) + 
  geom_density(aes(x=est,fill=biomarker_group_label),alpha=0.5) + 
  geom_density(data=priors,aes(x=value,fill="Prior"),alpha=0.25) +
  facet_wrap(~label, scales="free",ncol=2) +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_fill_manual(name="Biomarker group",values=c("Prior"="black","IC50"="blue","Avidity"="red")) +
  xlab("Value") +
  ylab("Density")
ggsave(paste0(save_wd,"/posterior_distributions.pdf"),p_theta,height=12,width=7,units="in",dpi=300)

write_csv(tidy_theta_chain,paste0(save_wd_results,"_theta_chains.csv"))
write_csv(priors %>% group_by(label) %>% sample_n(pmin(n(),1000)),paste0(save_wd_results,"_prior_draws.csv"))

p_theta <- tidy_theta_chain %>%
  ggplot() + 
  geom_density(aes(x=est,fill=biomarker_group),alpha=0.5) + 
  #geom_vline(aes(xintercept=values,linetype="True value"),col="black") +
  scale_fill_manual(name="Biomarker group",values=c("1"="blue","2"="red")) +
  #scale_linetype_manual(name="",values=c("True value"="dashed")) +
  facet_wrap(~names_new,scales="free",ncol=3) +
  xlab("Value") +
  ylab("Density") +
  theme_minimal() +
  theme(strip.text=element_text(size=6))
ggsave(paste0(save_wd,"/par_estimates.pdf"),p_theta,height=5*n_obs_types,width=7,units="in",dpi=300)

## Plot rho estimates by biomarker group
if(!is.null(measurement_indices)){
estimated_par_names <- par_tab[par_tab$fixed == 0 & par_tab$names %like% "rho","names"]
estimated_par_names <- make.names(estimated_par_names,unique=TRUE)

p_rhos <- chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
  filter(name %in% estimated_par_names) %>%
  dplyr::rename(names_new=name,est=value) %>%
  left_join(par_tab %>% filter(fixed == 0 & par_tab$names %like% "rho") %>% mutate(names_new=estimated_par_names) %>% select(names,names_new,values,biomarker_group) %>%
              mutate(biomarker_group = as.factor(biomarker_group))) %>%
  ggplot() + geom_density(aes(x=est,fill=biomarker_group),alpha=0.5) + 
  geom_vline(aes(xintercept=values,linetype="True value"),col="black") +
  scale_fill_manual(name="Biomarker group",values=c("1"="blue","2"="red")) +
  scale_linetype_manual(name="",values=c("True value"="dashed")) +
  facet_wrap(~names_new,scales="free",ncol=3) +
  xlab("Value") +
  ylab("Density") +
  theme_minimal() +
  theme(strip.text=element_text(size=6))
ggsave(paste0(save_wd,"/rho_estimates.pdf"),p_rhos,height=6,width=7,units="in",dpi=300)

}

## Post-hoc plots
chains <- load_mcmc_chains(chain_wd,par_tab,burnin=mcmc_pars["adaptive_iterations"],estimated_only = FALSE)

## Plot model predictions against observations
p_predictions <- plot_antibody_predictions(chains$theta_chain,chains$inf_chain,settings=serosolver_settings)

if(analysis != "avidity"){
## Save pointrange plot
ggsave(paste0(save_wd,"/prediction_pointrange.pdf"), p_predictions$p_pointrange + facet_wrap(~paste0(if_else(biomarker_group==1,"IC50","Avidity")),ncol=1),height=4.5*n_obs_types,width=8,units="in",dpi=300)
## Save draws
ggsave(paste0(save_wd,"/prediction_histograms.pdf"), p_predictions$p_hist_draws + facet_wrap(name~paste0(if_else(biomarker_group==1,"IC50","Avidity"))),height=9,width=8,units="in",dpi=300)
} else {
  ## Save pointrange plot
  ggsave(paste0(save_wd,"/prediction_pointrange.pdf"), p_predictions$p_pointrange + facet_wrap(~paste0(if_else(biomarker_group==1,"Avidity","Avidity")),ncol=1),height=4.5*n_obs_types,width=8,units="in",dpi=300)
  ## Save draws
  ggsave(paste0(save_wd,"/prediction_histograms.pdf"), p_predictions$p_hist_draws + facet_wrap(name~paste0(if_else(biomarker_group==1,"Avidity","Avidity"))),height=9,width=8,units="in",dpi=300)
  
}

p_predictions$p_hist_median 
tmp <- p_predictions$all_predictions
tmp$measurement <- tmp$measurement + rnorm(nrow(tmp),0,0.025)*as.numeric(tmp$measurement == 0)
tmp$median <- tmp$median + rnorm(nrow(tmp),0,0.025)*as.numeric(tmp$median == 0)

if(analysis != "avidity"){
  tmp <- tmp %>% mutate(label=if_else(biomarker_group==1,"IC50","Avidity"))
} else {
  tmp <- tmp %>% mutate(label=if_else(biomarker_group==1,"Avidity","IC50"))
}
tmp$label <- factor(tmp$label,levels=c("IC50","Avidity"))

## Save scatterplot comparison
p_pred_compare <- ggplot(tmp) + 
  geom_jitter(aes(x=measurement,y=median),alpha=0.5) + 
  geom_abline(slope=1,intercept=0,linetype="dashed",col="red") +
  xlab("Observed") +
  ylab("Predicted (posterior median)") +
  theme_bw()+
  facet_wrap(~label,ncol=1)
ggsave(paste0(save_wd,"/prediction_scatterplot.pdf"),p_pred_compare,height=3.5*n_obs_types,width=5,units="in",dpi=300)

antibody_data_expanded <- expand_grid(serosolver_settings$antibody_data %>% select(individual,birth) %>% distinct(),
            sample_time=seq(min(serosolver_settings$antibody_data$birth),
                            max(serosolver_settings$antibody_data$sample_time),by=1),
            biomarker_group=1,repeat_number=1,biomarker_id=seq(min(serosolver_settings$antibody_data$birth),
                                                               max(serosolver_settings$antibody_data$sample_time),by=1),
            measurement=0)
            
#tmp_map <- serosolver_settings$antigenic_map
#tmp_map$x_coord <- 1
#tmp_map$y_coord<- 1
## Get predicted antibody levels at all times
## Get antibody titre predicted pre-boosts
predicted_titers_all <- get_antibody_level_predictions(chains$theta_chain,chains$inf_chain,
                                                       demographics=demography,
                                                       antibody_data=antibody_data_expanded,
                                                       individuals=unique(serosolver_settings$antibody_data$individual),
                                                       antigenic_map=serosolver_settings$antigenic_map,#tmp_map,
                                                       par_tab=serosolver_settings$par_tab,
                                                       nsamp=200,
                                                       expand_to_all_times=FALSE,expand_antibody_data=FALSE,
                                                       antibody_level_before_infection=TRUE,
                                                       for_regression=TRUE)

antibody_data_expanded <- bind_cols(antibody_data_expanded, predicted_titers_all$all_predictions)

## Quick check
ggplot(antibody_data_expanded %>% filter(sample_time==biomarker_id,individual %in% 1:25) %>% pivot_longer(-c(individual,birth,sample_time,biomarker_group,repeat_number,biomarker_id,measurement))) + geom_line(aes(x=sample_time,y=value,group=name)) + facet_wrap(~individual)


## Merge with infection states -- need to convert the infection history matrices to long-format data frames
samp_values <- as.numeric(colnames(predicted_titers_all$all_predictions))
pulled_infection_states <- predicted_titers_all$all_inf_hist
tmp_states <- NULL
for(i in 1:length(pulled_infection_states)){
  x <- as_tibble(pulled_infection_states[[i]])
  colnames(x) <- 1:ncol(x)
  x$individual <- 1:nrow(x)
  x <- x %>% pivot_longer(-individual)
  x$samp <- as.numeric(samp_values[i])
  tmp_states[[i]] <- x
}
tmp_states <- do.call("bind_rows",tmp_states)
## Set to correct time
tmp_states$sample_time <- serosolver_settings$possible_exposure_times[as.numeric(tmp_states$name)]

## Pull out only titers against the circulating strain at each time
antibody_data_expanded_long <- antibody_data_expanded %>% 
  filter(sample_time==biomarker_id) %>% 
  pivot_longer(-c(individual,birth,sample_time,biomarker_group,repeat_number,biomarker_id,measurement)) %>%
  mutate(samp=as.numeric(name))

## Merge predicted titres and infection histories
combined_preds_inf_hists <- antibody_data_expanded_long %>% filter(sample_time >= birth) %>% left_join(tmp_states %>% select(-name) %>% rename(x=value)) 

## Plot to check
combined_preds_inf_hists %>% group_by(individual,sample_time) %>% summarize(prob_infected=sum(x)/n(),mean_titre=mean(value)) %>% filter(individual %in% 1:36) %>%
  ggplot() + 
  geom_line(aes(x=sample_time,y=mean_titre)) +
  geom_vline(aes(xintercept=sample_time,alpha=prob_infected),col="orange") +
  facet_wrap(~individual,ncol=6)

## Version where we move consecutive zeros
mod <- combined_preds_inf_hists %>% arrange(samp, individual, sample_time) %>% mutate(x_new = if_else(x == 1 & lag(x == 1), 0, x))

## Look at relationship between pre-infection titres and infection states
## First, posterior draws from original data

titer_infection_relationship_unmodified <- combined_preds_inf_hists %>% 
  mutate(value = round(value,0)) %>% 
  group_by(samp, value) %>% 
  summarize(n_infected=sum(x),N=n()) %>%
  mutate(prop_infected=n_infected/N) %>% 
  group_by(samp) %>%
  mutate(prop_infected=prop_infected/prop_infected[1])
titer_infection_relationship_unmodified_summary <- titer_infection_relationship_unmodified %>%
  group_by(value) %>% summarize(mean=mean(prop_infected),lower=quantile(prop_infected,0.025),upper=quantile(prop_infected,0.975)) %>%
  mutate("Version"="Unmodified")

titer_infection_relationship_modified <- mod %>% 
  mutate(value = round(value,0)) %>% 
  group_by(samp, value) %>% 
  summarize(n_infected=sum(x_new),N=n()) %>%
  mutate(prop_infected=n_infected/N) %>% 
  group_by(samp) %>%
  mutate(prop_infected=prop_infected/prop_infected[1])

titer_infection_relationship_modified_summary <- titer_infection_relationship_modified %>%
  group_by(value) %>% summarize(mean=mean(prop_infected),lower=quantile(prop_infected,0.025),upper=quantile(prop_infected,0.975)) %>%
  mutate("Version"="Consecutive infections removed")

titer_infection_relationship_summary <- bind_rows(titer_infection_relationship_unmodified_summary,titer_infection_relationship_modified_summary)

p_titre_relationship <- ggplot(titer_infection_relationship_summary) + 
  geom_ribbon(aes(x=value,ymin=lower,ymax=upper,fill=Version),alpha=0.25) + geom_line(aes(x=value,y=mean,col=Version)) + 
  xlab("Predicted antibody level (log2 IC50)") + ylab("Risk of infection relative to seronegative") +
  theme_bw() +
  scale_color_manual(name="",values=c("Unmodified"="#E69F00","Consecutive infections removed"="#56B4E9"))+
  scale_fill_manual(name="",values=c("Unmodified"="#E69F00","Consecutive infections removed"="#56B4E9"))+
  theme(legend.position="bottom") +
  facet_wrap(~Version,ncol=1)

ggsave(paste0(save_wd,"/titer_infection_relationship.pdf"),p_titre_relationship,height=5,width=5,units="in",dpi=300)


## Look another way -- compare posterior mean prob of infection against posterior mean titre
p_titre_relationship_posterior_means <- mod %>% group_by(individual,sample_time) %>% 
  summarize(`Probability of infection`=mean(x),`Probability of infection (modified)`=mean(x_new), mean_titre=mean(value)) %>%
  pivot_longer(-c(individual,sample_time,mean_titre),names_to="type",values_to="prob_infected") %>%
  ggplot() + 
  geom_point(aes(x=mean_titre,y=prob_infected),size=0.25,alpha=0.5) +
  geom_smooth(aes(x=mean_titre,y=prob_infected),method="loess") +
  geom_abline(slope=-1/3,intercept=1,linetype="dashed",col="red") +
  xlab("Predicted antibody level (log2 IC50)") + ylab("Posterior probability of infection") +
  coord_cartesian(ylim=c(0,1)) +
  theme_bw() +
  facet_wrap(~type,ncol=1)
p_titre_relationship_posterior_means
ggsave(paste0(save_wd,"/titer_infection_relationship_alt.pdf"),p_titre_relationship_posterior_means,height=5,width=5,units="in",dpi=300)





antibody_data_expanded_long_all <- antibody_data_expanded %>% pivot_longer(-c(individual,birth,sample_time,biomarker_group,repeat_number,biomarker_id,measurement)) %>%
  mutate(samp=as.numeric(name))

## Plot population overall antibody landscape over time
population_mean_circulating_titre <- antibody_data_expanded_long_all %>% 
  group_by(sample_time, samp, biomarker_id) %>% 
  dplyr::summarize(mean_titre=mean(value)) %>% 
  group_by(sample_time,biomarker_id) %>%
  summarize(mean_titre1 = mean(mean_titre),lower=quantile(mean_titre,0.025),upper=quantile(mean_titre,0.975)) %>% filter(sample_time == biomarker_id)

p_circulating_mean <- ggplot(population_mean_circulating_titre) + 
  geom_ribbon(aes(x=sample_time,ymin=lower,ymax=upper),fill="grey70") +
  geom_line(aes(x=sample_time,y=mean_titre1)) 
population_mean_circulating_titre$strain <- "Matched to \ncirculating virus"

population_mean_titres <- antibody_data_expanded_long_all %>% 
  group_by(sample_time, samp, biomarker_id) %>% 
  dplyr::summarize(mean_titre=mean(value)) %>% 
  group_by(sample_time,biomarker_id) %>%
  summarize(mean_titre1 = mean(mean_titre),lower=quantile(mean_titre,0.025),upper=quantile(mean_titre,0.975)) %>%
  filter(biomarker_id %in% c(2002,2006,2009,2012))

noro_key <- c("2002"="FH-2002","2006"="DH-2006","2009"="NO-2009","2012"="SY-2012")
population_mean_titres$strain <- noro_key[as.character(population_mean_titres$biomarker_id)]
population_mean_titres <- bind_rows(population_mean_titres,population_mean_circulating_titre)
population_mean_titres$strain <- factor(population_mean_titres$strain, levels=c("FH-2002","DH-2006","NO-2009","SY-2012","Matched to \ncirculating virus"))
p_all_titre_mean <- ggplot(population_mean_titres) + 
  geom_ribbon(aes(x=sample_time,ymin=lower,ymax=upper,fill=strain,group=strain),alpha=0.25) +
  geom_line(aes(x=sample_time,y=mean_titre1,col=strain)) +
  scale_fill_manual(name="Norovirus strain",values=c(strain_colors_short,"Matched to \ncirculating virus"="blue")) +
  scale_color_manual(name="Norovirus strain",values=c(strain_colors_short,"Matched to \ncirculating virus"="blue")) +
  theme_classic() +
  xlab("Norovirus year") +
  ylab("Mean log2 IC50 titre") +
  scale_x_continuous(limits=c(2001,2012),breaks=seq(2001,2012,by=1)) +
  theme(legend.position=c(0.2,0.7))
p_all_titre_mean

ggsave(filename=paste0(save_wd,"/population_mean_titres.pdf"),plot=p_all_titre_mean,height=3,width=7,units="in",dpi=300)
write_csv(population_mean_titres,paste0(save_wd_results,"_population_mean_titres.csv"))
################################ 
## IMPORTANT
################################
## Set the known_infection_history argument to NULL for the real data
samp_indivs <- sample(antibody_data %>% filter(biomarker_group == 1) %>% pull(individual) %>% unique(),25,replace=FALSE)
samp_indivs <- samp_indivs[order(samp_indivs)]
p_fits <- plot_model_fits(chains$theta_chain,chains$inf_chain,
                          settings=serosolver_settings,
                                  individuals = 1:25,
                                  #known_infection_history=as.matrix(all_simulated_data$infection_histories),
                                  orientation="cross-sectional",expand_to_all_times = FALSE)

if(analysis != "avidity"){
  p_fit1 <- p_fits[[1]] + facet_wrap(~individual) + ggtitle("IC50") + coord_cartesian(xlim=c(2001,2013)) + scale_x_continuous(breaks=c(2002,2006,2009,2012),labels=c("FH-2002","DH-2006","NO-2009","SY-2012"))
  ggsave(paste0(save_wd,"/model_fits_IC50.pdf"),p_fit1,height=7,width=8,units="in",dpi=300)
} else {
  p_fit1 <- p_fits[[1]] + facet_wrap(~individual) + ggtitle("Avidity") + coord_cartesian(xlim=c(2001,2013))+ scale_x_continuous(breaks=c(2002,2006,2009,2012),labels=c("FH-2002","DH-2006","NO-2009","SY-2012"))
  ggsave(paste0(save_wd,"/model_fits_avidity.pdf"),p_fit1,height=7,width=8,units="in",dpi=300)
}

if(length(p_fits) > 1){
  p_fit2 <- p_fits[[2]] + facet_wrap(~individual) + ggtitle("Avidity")+ coord_cartesian(xlim=c(2001,2013)) + scale_x_continuous(breaks=c(2002,2006,2009,2012),labels=c("FH-2002","DH-2006","NO-2009","SY-2012"))
  ggsave(paste0(save_wd,"/model_fits_avidity.pdf"),p_fit2,height=7,width=8,units="in",dpi=300)
}

################################ 
## IMPORTANT
################################
## Set the real_inf_hist argument to NULL for the real data
## Plot estimated infections against truth
p_inf_comb <- plot_cumulative_infection_histories(chains$inf_chain,indivs=1:25,
                                                  #real_inf_hist=as.matrix(all_simulated_data$infection_histories),
                                                  possible_exposure_times=possible_exposure_times,number_col=5,pad_chain=TRUE,nsamp=100)
## Densities
p_inf_density <- p_inf_comb[[2]] + ylab("Posterior probability of infection") + 
  scale_x_continuous(breaks=seq(2000,2012,by=4),labels=seq(2000,2012,by=4)) +
  theme(axis.text.x=element_text(size=7,angle=45,hjust=1))+ xlab("Date (years)")
## Cumulative infections
p_inf_cumu <- p_inf_comb[[1]] + ylab("Cumulative number of infections") + 
  scale_x_continuous(breaks=seq(2000,2012,by=4),labels=seq(2000,2012,by=4)) +
  theme(axis.text.x=element_text(size=7,angle=45,hjust=1)) + xlab("Date (years)")


ggsave(paste0(save_wd,"/inf_hists_density.pdf"),p_inf_density,height=7,width=8,units="in",dpi=300)
ggsave(paste0(save_wd,"/inf_hists_cumulative.pdf"),p_inf_cumu,height=7,width=8,units="in",dpi=300)

## Look at model-predicted observations against actual observations
antibody_predictions <- plot_antibody_predictions(chains$theta_chain,chains$inf_chain,settings=serosolver_settings)
print(antibody_predictions$proportion_correct)
p_preds <- antibody_predictions$p_pointrange
ggsave(paste0(save_wd,"/predicted_observations.pdf"),p_preds,height=4,width=7,units="in",dpi=300)

if(analysis == "both"){
p_estimated_model <- plot_estimated_antibody_model(chains$theta_chain,antibody_data=expand_grid(individual=1,sample_time=seq(2001,2012,by=1),biomarker_id=unique(antibody_data$biomarker_id),biomarker_group=unique(antibody_data$biomarker_group)),demographics=NULL,settings=serosolver_settings,by_group=FALSE,nsamp=100,
                                                   par_tab=par_tab %>% mutate(stratification=NA),
                              antigenic_map=antigenic_map %>% filter(inf_times >= 2001),possible_exposure_times=seq(2001,2012,by=1),
                              solve_times=seq(2001,2012,by=0.1)) + facet_grid(paste0("Norovirus strain: ", biomarker_id)~if_else(biomarker_group == 1, "1. IC50","2. Avidity")) + coord_cartesian(ylim=c(0,6)) +
  ggtitle("Antibody kinetics assuming single infection in 2001") + theme(legend.position="none")
  ggsave(paste0(save_wd,"/estimated_antibody_model.pdf"),p_estimated_model,height=8,width=7,units="in",dpi=300)
  
  
  p_estimated_model2 <- plot_estimated_antibody_model(chains$theta_chain,
                                                      antibody_data=expand_grid(individual=1,sample_time=seq(2001,2012,by=1),
                                                                                biomarker_id=unique(antigenic_map$inf_times),
                                                                                biomarker_group=unique(antibody_data$biomarker_group)),
                                                      demographics=NULL,settings=serosolver_settings,by_group=FALSE,nsamp=100,
                                                      par_tab=par_tab %>% mutate(stratification=NA),
                                                      antigenic_map=antigenic_map,
                                                      solve_times=seq(2000,2012,by=1),
                                                      set_infections = c(2001,2010))
  tmp_dat <- p_estimated_model2$data %>% select(biomarker_group,biomarker_id,sample_time,median)
  tmp_dat2 <- tmp_dat %>% mutate(sample_time = sample_time + 1)
  tmp_dat2 <- tmp_dat2 %>% select(biomarker_group,biomarker_id,sample_time,median) %>% rename(past_median=median)
  tmp_dat <- left_join(tmp_dat,tmp_dat2) %>% drop_na() %>% rename("Current antibody landscape"=median,"Past year antibody landscape"=past_median)
  tmp_dat <- tmp_dat %>% pivot_longer(-c(biomarker_group,biomarker_id,sample_time))
  tmp_dat$name <- factor(tmp_dat$name,levels=c("Past year antibody landscape","Current antibody landscape"))
  
  infection_lines <- expand_grid(infection_times=c(2001,2010),sample_time=seq(2001,2012,by=1)) %>% filter(infection_times <= sample_time)
  
  p_estimated_model_landscape1 <- ggplot(tmp_dat %>% filter(biomarker_group == 1)) + 
    geom_ribbon(aes(x=biomarker_id,ymax=value,ymin=0,fill=name),col="black")+
    geom_vline(data=infection_lines,aes(xintercept=infection_times,linetype="Infection")) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Current antibody landscape"="grey70","Past year antibody landscape"="skyblue")) +
    facet_wrap(~paste0("Sample time: ", sample_time)) +
    xlab("Norovirus strain (year of circulation)") +
    scale_linetype_manual(name="",values=c("Infection"="dashed")) +
    scale_y_continuous(limits=c(0,8),expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    ylab("Antibody level") +
    theme(legend.position="bottom")
  
  
  p_estimated_model_landscape2 <- ggplot(tmp_dat %>% filter(biomarker_group == 2)) + 
    geom_ribbon(aes(x=biomarker_id,ymax=value,ymin=0,fill=name),col="black")+
    geom_vline(data=infection_lines,aes(xintercept=infection_times,linetype="Infection")) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Current antibody landscape"="grey70","Past year antibody landscape"="skyblue")) +
    facet_wrap(~paste0("Sample time: ", sample_time)) +
    xlab("Norovirus strain (year of circulation)") +
    scale_linetype_manual(name="",values=c("Infection"="dashed")) +
    scale_y_continuous(limits=c(0,8),expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    ylab("Antibody level") +
    theme(legend.position="bottom")
  
  ggsave(paste0(save_wd,"/estimated_antibody_landscape_model_ic50.pdf"),p_estimated_model_landscape1,height=7,width=8,units="in",dpi=300)
  ggsave(paste0(save_wd,"/estimated_antibody_landscape_model_avidity.pdf"),p_estimated_model_landscape2,height=7,width=8,units="in",dpi=300)
  
} else {
  p_estimated_model <- plot_estimated_antibody_model(chains$theta_chain,antibody_data=expand_grid(individual=1,sample_time=seq(2002,2012,by=1),biomarker_id=unique(antibody_data$biomarker_id),biomarker_group=unique(antibody_data$biomarker_group)),demographics=NULL,settings=serosolver_settings,by_group=FALSE,nsamp=100,
                                antigenic_map=antigenic_map %>% filter(inf_times >= 2002),possible_exposure_times=seq(2002,2012,by=1),
                                par_tab=par_tab %>% mutate(stratification=NA),
                                solve_times=seq(2002,2012,by=0.1)) + facet_wrap(~paste0("Norovirus strain: ", biomarker_id)) + coord_cartesian(ylim=c(0,6)) +
    ggtitle("Antibody kinetics assuming single infection in 2001") + theme(legend.position="none")
  ggsave(paste0(save_wd,"/estimated_antibody_model.pdf"),p_estimated_model,height=8,width=7,units="in",dpi=300)
  
  ## Alternative way of plotting
  tmp_dat <- p_estimated_model$data
  noro_key <- c("2002"="FH-2002","2006"="DH-2006","2009"="NO-2009","2012"="SY-2012")
  tmp_dat$strain <- noro_key[as.character(tmp_dat$biomarker_id)]
  tmp_dat$strain <- factor(tmp_dat$strain, levels=c("FH-2002","DH-2006","NO-2009","SY-2012"))
  p_kinetics_alt <- ggplot(tmp_dat) +
    geom_line(aes(x=sample_time,y=mean,col=strain)) +
    geom_ribbon(aes(x=sample_time,ymax=upper,ymin=lower,fill=strain),alpha=0.25) +
    theme_classic() +
    xlab("Time since infection (years)") +
    ylab("log2 IC50") +
    scale_color_manual(name="Measured strain",values=strain_colors_short) +
    scale_fill_manual(name="Measured strain",values=strain_colors_short) +
    scale_x_continuous(limits=c(2002,2012), breaks=seq(2002,2012,by=2),labels=seq(0,10,by=2)) +
    theme(legend.position=c(0.8,0.8),plot.tag=element_text(face="bold")) +
    labs(tag="A")
  save(p_kinetics_alt, file=paste0("figures/estimated_kinetics_",run_name,".rda"))
  ggsave(paste0("figures/estimated_kinetics_",run_name,".pdf"),p_kinetics_alt,height=3,width=5,units="in",dpi=300)
  
  ## Alternatives for infection from other strains
  tmp_dat <- NULL
  for(strain in c(2002,2006,2009,2012)){
  tmp_dat <- bind_rows(tmp_dat, plot_estimated_antibody_model(chains$theta_chain,antibody_data=expand_grid(individual=1,sample_time=seq(strain,strain+10,by=0.1),biomarker_id=unique(antibody_data$biomarker_id),biomarker_group=unique(antibody_data$biomarker_group)),demographics=NULL,settings=serosolver_settings,by_group=FALSE,nsamp=100,
                                                     antigenic_map=antigenic_map,
                                                     par_tab=par_tab %>% mutate(stratification=NA),
                                                     solve_times=seq(strain,strain+10,by=0.1),set_infections = c(strain))$data %>% mutate("Infection strain"=strain))
  }
  
  ## Alternative way of plotting
  noro_key <- c("2002"="FH-2002","2006"="DH-2006","2009"="NO-2009","2012"="SY-2012")
  tmp_dat$strain <- noro_key[as.character(tmp_dat$biomarker_id)]
  tmp_dat$strain <- factor(tmp_dat$strain, levels=c("FH-2002","DH-2006","NO-2009","SY-2012"))
  tmp_dat$`Infection strain` <- noro_key[as.character(tmp_dat$`Infection strain`)]
  tmp_dat$`Infection strain` <- factor(tmp_dat$`Infection strain`, levels=c("FH-2002","DH-2006","NO-2009","SY-2012"))
  tmp_dat$label <- paste0("Infection with ", tmp_dat$`Infection strain`)
  tmp_dat$label <- factor(tmp_dat$label, unique(tmp_dat$label))
  tmp_dat <- tmp_dat %>% group_by(`Infection strain`) %>% mutate(sample_time=sample_time-min(sample_time))
  p_kinetics_alt <- ggplot(tmp_dat) +
    geom_line(aes(x=sample_time,y=mean,col=strain)) +
    geom_ribbon(aes(x=sample_time,ymax=upper,ymin=lower,fill=strain),alpha=0.25) +
    theme_classic() +
    xlab("Time since infection (years)") +
    ylab("log2 IC50") +
    scale_color_manual(name="Measured strain",values=strain_colors_short) +
    scale_fill_manual(name="Measured strain",values=strain_colors_short) +
    scale_x_continuous(breaks=seq(0,10,by=2)) +
    theme(legend.position=c(0.8,0.8),plot.tag=element_text(face="bold")) +
    facet_wrap(~label) +
    theme(legend.position="bottom",plot.tag=element_text(face="bold"))
  save(p_kinetics_alt, file=paste0("figures/estimated_kinetics_",run_name,"_all.rda"))
  ggsave(paste0("figures/estimated_kinetics_",run_name,"_all.pdf"),p_kinetics_alt,height=5,width=6,units="in",dpi=300)
  
  
  p_estimated_model2 <- plot_estimated_antibody_model(chains$theta_chain,
                                                      antibody_data=expand_grid(individual=1,sample_time=seq(2001,2012,by=1),
                                                                                biomarker_id=unique(antigenic_map$inf_times),
                                                                                biomarker_group=unique(antibody_data$biomarker_group)),
                                                      demographics=NULL,settings=serosolver_settings,by_group=FALSE,nsamp=100,
                                                      antigenic_map=antigenic_map,
                                                      solve_times=seq(2000,2012,by=1),
                                                      par_tab=par_tab %>% mutate(stratification=NA),
                                                      set_infections = c(2001,2010))
  tmp_dat <- p_estimated_model2$data %>% select(biomarker_group,biomarker_id,sample_time,median)
  tmp_dat2 <- tmp_dat %>% mutate(sample_time = sample_time + 1)
  tmp_dat2 <- tmp_dat2 %>% select(biomarker_group,biomarker_id,sample_time,median) %>% rename(past_median=median)
  tmp_dat <- left_join(tmp_dat,tmp_dat2) %>% drop_na() %>% rename("Current antibody landscape"=median,"Past year antibody landscape"=past_median)
  tmp_dat <- tmp_dat %>% pivot_longer(-c(biomarker_group,biomarker_id,sample_time))
  tmp_dat$name <- factor(tmp_dat$name,levels=c("Past year antibody landscape","Current antibody landscape"))
  
  infection_lines <- expand_grid(infection_times=c(2001,2010),sample_time=seq(2001,2012,by=1)) %>% filter(infection_times <= sample_time)
  
  p_estimated_model_landscape <- ggplot(tmp_dat %>% filter(biomarker_group == 1)) + 
    geom_ribbon(aes(x=biomarker_id,ymax=value,ymin=0,fill=name),col="black")+
    geom_vline(data=infection_lines,aes(xintercept=infection_times,linetype="Infection")) +
    theme_classic() +
    scale_fill_manual(name="",values=c("Current antibody landscape"="grey70","Past year antibody landscape"="skyblue")) +
    facet_wrap(~paste0("Sample time: ", sample_time)) +
    xlab("Norovirus strain (year of circulation)") +
    scale_linetype_manual(name="",values=c("Infection"="dashed")) +
    scale_y_continuous(limits=c(0,8),expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    ylab("Antibody level") +
    theme(legend.position="bottom")
  
  ggsave(paste0(save_wd,"/estimated_antibody_landscape_model.pdf"),p_estimated_model_landscape,height=7,width=8,units="in",dpi=300)
  
  ## For schematic
  p_estimated_model3 <- plot_estimated_antibody_model(chains$theta_chain,# %>% mutate(boost_long = 2, boost_short = 3, wane_short=0.4,cr_long=1,cr_short=0.05),
                                                      antibody_data=expand_grid(individual=1,sample_time=seq(2000,2015,by=0.1),
                                                                                biomarker_id=unique(antigenic_map$inf_times),
                                                                                biomarker_group=unique(antibody_data$biomarker_group)),
                                                      demographics=NULL,settings=serosolver_settings,by_group=FALSE,nsamp=100,
                                                      antigenic_map=antigenic_map,
                                                      par_tab=par_tab %>% mutate(stratification=NA),
                                                      solve_times=seq(2000,2015,by=0.1),
                                                      set_infections = c(2002,2010))
  tmp_dat <- p_estimated_model3$data %>% select(biomarker_group,biomarker_id,sample_time,median,lower,upper)
  tmp_dat <- tmp_dat %>% filter(biomarker_id %in% c(2002,2006,2009,2012))
  infection_lines <- expand_grid(infection_times=c(2001,2010),sample_time=seq(2001,2012,by=1)) %>% filter(infection_times <= sample_time)
  noro_key <- c("2002"="FH-2002","2006"="DH-2006","2009"="NO-2009","2012"="SY-2012")
  tmp_dat$strain <- noro_key[as.character(tmp_dat$biomarker_id)]
  tmp_dat$strain <- factor(tmp_dat$strain, levels=c("FH-2002","DH-2006","NO-2009","SY-2012"))
  
 
  # Reduce spacing around plot and legend
  
  p_model_schematic <- ggplot(tmp_dat) + 
    geom_ribbon(aes(x=sample_time,ymin=lower,ymax=upper,fill=strain),alpha=0.25)+
    geom_line(aes(x=sample_time,y=median,col=strain)) +
    theme_classic() +
    geom_vline(xintercept=c(2001,2010),linetype="dashed") +
    geom_vline(xintercept=2012,linetype="dotted") +
    scale_y_continuous(expand=c(0,0),limits=c(0,8),breaks=seq(0,8,by=1)) +
    scale_x_continuous(breaks=seq(2000,2015,by=1)) +
    scale_color_manual(name="Norovirus strain",values=strain_colors_short) +
    scale_fill_manual(name="Norovirus strain",values=strain_colors_short) +
    ylab("Antibody level (log2 IC50)") +
    xlab("Date") +
    theme_use +
    theme(legend.position="top",
          axis.text.x=element_text(angle=45,hjust=1,size=6)) +
    labs(tag="A")
  
  p_schematic_data <- ggplot(tmp_dat %>% filter(sample_time==2012)) +
    geom_bar(aes(x=strain,y=median,fill=strain),stat="identity")+
    scale_fill_manual(name="Norovirus strain",values=strain_colors_short) +
    scale_y_continuous(expand=c(0,0),limits=c(0,4),breaks=seq(0,4,by=1)) +
    ylab("Antibody level (log2 IC50)") +
    xlab("Norovirus strain") +
    theme_use +
    theme(legend.position="none",
          axis.text.x=element_text(angle=45,hjust=1,size=6))+
    labs(tag="B")
  
  save(p_model_schematic, p_schematic_data, file=paste0("figures/model_schematic_",run_name,".rda"))
  
  ggsave(p_model_schematic + p_schematic_data + plot_layout(ncol=2,widths=c(3,1)),
         filename=paste0("figures/model_schematic_",run_name,".pdf"),height=3,width=7,units="in",dpi=300)
  ggsave(p_model_schematic + p_schematic_data + plot_layout(ncol=2,widths=c(3,1)),
         filename=paste0("figures/model_schematic_",run_name,".png"),height=3,width=7,units="in",dpi=300)
} 

################################ 
## IMPORTANT
################################
## Set the true_ar argument to NULL for the real data
true_ar <- read.csv("~/Documents/local_data/norovirus_test/true_ar.csv")
p_ar <- plot_attack_rates(chains$inf_chain,settings = serosolver_settings,
                          #true_ar=all_simulated_data$attack_rates, ## Set this to NULL for real data
                          by_group=FALSE, 
                          plot_den = TRUE,
                          pad_chain=TRUE) +
  geom_point(data=true_ar,aes(x=time,y=ar),col="red",shape=17,size=3) +
  coord_cartesian(xlim=c(2000,2013))

p_compare_popn_mean_titre_ar <- p_ar$data %>% group_by(time) %>% dplyr::summarize(mean_ar=mean(V1)) %>%
  rename(sample_time=time) %>% left_join(population_mean_titres %>% filter(strain == "Matched to \ncirculating virus")) %>%
  ggplot() + geom_point(aes(x=mean_titre1,y=mean_ar),size=1.5) + 
  geom_smooth(aes(x=mean_titre1,y=mean_ar),method="lm") +
  theme_use +
  #scale_y_continuous(limits=c(0,0.5),breaks=seq(0,0.5,by=0.1)) +
  #scale_x_continuous(limits=c(0,0.6),breaks=seq(0,0.6,by=0.1))+
  xlab("Population mean log2 IC50 titre against\n circulating virus at start of season") +
  ylab("Posterior mean attack rate estimate")

p_compare_popn_mean_titre_ar

ggsave(paste0(save_wd,"/attack_rates_year.pdf"),p_ar,height=4,width=7,units="in",dpi=300)
ggsave(paste0(save_wd,"/attack_rates_against_titres.pdf"),p_compare_popn_mean_titre_ar,height=3,width=5,units="in",dpi=300)

## Pull out attack rates
write_csv(p_ar$data %>% mutate(run_name=run_name,map=cart_data),paste0(save_wd_results,"_attack_rates.csv"))

# attack rate by age of individual
# *** to this previously we merged to inf_chain back into the data
head(antibody_data_keep)

## I think fixed -- this should be the proportion of individuals infected per age group i.e., group by year of age, get total infections, divide by number in that age group (should be between 0 and 1)
inf_chain <- chains$inf_chain
inf_chain <- expand_grid(i=unique(antibody_data$individual),j=unique(inf_chain$j),samp_no=unique(inf_chain$samp_no),chain_no=unique(inf_chain$chain_no)) %>%
  left_join(inf_chain) %>%
  mutate(x=if_else(is.na(x),0,x)) %>%
  as.data.table()

inf_chain$j <- antigenic_map$inf_times[inf_chain$j]
data.table::setkey(inf_chain, "i", "j","samp_no","chain_no")
n_inf_chain_i <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
colnames(n_inf_chain_i)[1] <- "individual"

## Merge with year of birth to calculate age
titre_dat <- antibody_data_keep %>% select(individual,birth,sample_time) %>% distinct() %>% group_by(individual,birth) %>% filter(sample_time == max(sample_time)) %>% ungroup() %>% as.data.table()
#titre_dat <- data.table(antibody_data_keep %>% filter(biomarker_id == 2002 & biomarker_group == 1)) # otherwise each there several times
setkey(titre_dat, "individual")
titre_dat <- left_join(titre_dat,n_inf_chain_i ) # not all have infections NAs

titre_dat <- titre_dat %>% filter(j >= birth) ## Only count draws where individual is alive
titre_dat <- titre_dat %>% filter(j <= sample_time) ## Only count draws where individual is alive
titre_dat <- titre_dat %>% mutate(age = j - birth) ## Calculate age in years at each time point
titre_dat <- titre_dat %>% mutate(age=if_else(age >= 6,6,age)) ## Merge over 6

## Find attack rate by age, sample and chain i.e., ignore year
infections_by_age <- titre_dat %>% group_by(age, samp_no,chain_no) %>% dplyr::summarize(ar=sum(V1)/n(),N=n())
infections_by_age_label <- infections_by_age %>% ungroup() %>% select(age, N) %>% distinct()
p_age1 <- ggplot(infections_by_age) + 
  geom_violin(aes(x=age,group=age,y=ar),scale="width",draw_quantiles = c(0.025,0.5,0.975),fill="grey70") +
  geom_text(data=infections_by_age_label,aes(x=age,y=0.35,label=paste0("N=",N)),size=3) +
  scale_y_continuous(limits=c(0,0.4),expand=c(0,0),breaks=seq(0,0.4,by=0.1)) +
  theme_classic() +
  ylab("Attack rate by year of age") +
  xlab("Age (years)") 

## Same thing, but filter by only after 2005
infections_by_age_2006 <- titre_dat %>% filter(j > 2005) %>% group_by(age, samp_no,chain_no) %>% dplyr::summarize(ar=sum(V1)/n(),N=n())
infections_by_age_label_2006 <- infections_by_age %>% ungroup() %>% select(age, N) %>% distinct()
p_age2 <- ggplot(infections_by_age_2006) + 
  geom_violin(aes(x=age,group=age,y=ar),scale="width",draw_quantiles = c(0.025,0.5,0.975),fill="grey70") +
  geom_text(data=infections_by_age_label_2006,aes(x=age,y=0.35,label=paste0("N=",N)),size=3) +
  scale_y_continuous(limits=c(0,0.4),expand=c(0,0),breaks=seq(0,0.4,by=0.1)) +
  theme_classic() +
  ylab("Attack rate by year of age") +
  xlab("Age (years)") +
  ggtitle("2006+")

p_age_comb <- (p_age1 + ggtitle("All years"))/p_age2


infections_by_age_and_year <- titre_dat %>% group_by(age, samp_no,chain_no,j) %>% dplyr::summarize(ar=sum(V1)/n(),N=n())

p_ar_age_year <- ggplot(infections_by_age_and_year %>% group_by(age,j) %>% 
         dplyr::summarize(med=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975)) %>%
         mutate(label = paste0(if_else(age >= 6, "6+",as.character(age))," years"))) + 
  geom_pointrange(aes(x=j,group=j,y=med,ymin=lower,ymax=upper),size=0.25,linewidth=0.5) + 
  facet_wrap(~label) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2)) +
  scale_x_continuous(limits=c(2000,2012),breaks=seq(2000,2012,by=2))+
  theme_classic() +
  xlab("Year") + ylab("Per capita attack rate")

p_ar_year_age <- ggplot(infections_by_age_and_year %>% group_by(age,j) %>% 
         dplyr::summarize(med=median(ar),lower=quantile(ar,0.025),upper=quantile(ar,0.975)) %>%
         mutate(label = paste0(if_else(age >= 6, "6+",as.character(age))))) + 
  geom_pointrange(aes(x=label,y=med,ymin=lower,ymax=upper),size=0.25,linewidth=0.5) + 
  facet_wrap(~j) +
  #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2)) +
  #scale_x_continuous(limits=c(2000,2012),breaks=seq(2000,2012,by=2))+
  theme_classic() +
  xlab("Age (years)") + ylab("Per capita attack rate")

ggsave(filename=paste0(save_wd,"/attack_rate_age.pdf"),p_age1,height=3,width=5,units="in",dpi=300)
ggsave(filename=paste0(save_wd,"/attack_rate_age_comb.pdf"),p_age_comb,height=5,width=5,units="in",dpi=300)
ggsave(filename=paste0(save_wd,"/attack_rate_age_year.pdf"),p_ar_age_year,height=6,width=7,units="in",dpi=300)
ggsave(filename=paste0(save_wd,"/attack_rate_year_age.pdf"),p_ar_year_age,height=6,width=7,units="in",dpi=300)

## Split individuals into age groups and plot summaries
# make the 1 0 age kid 1 years
titre_dat$age[titre_dat$age==0.8] <- 1.0
titre_dat$age_group <- cut(titre_dat$age, breaks = c(seq(1,6,1),9),right=F) # on second birthday, different group
table(titre_dat$age_group)
titre_dat$age_rnd <- floor(titre_dat$age)
table(titre_dat$age_rnd)
titre_dat$infs_per_life <- (titre_dat$median_infs)/(titre_dat$age) # specific age rather than rounded

write_csv(infections_by_age_2006,file=paste0(save_wd_results,"_infection_histories_age2006.csv"))
write_csv(infections_by_age,file=paste0(save_wd_results,"_infection_histories_age.csv"))
write_csv(infections_by_age_and_year,file=paste0(save_wd_results,"_infection_histories_age_and_year.csv"))
