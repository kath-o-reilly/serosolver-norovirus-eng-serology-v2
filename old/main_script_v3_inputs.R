######################################################
## Author: James Hay & Kath O'Reilly
## Date: 10 Jun 2025
## Summary: run the serosolver simulations applied for norovirus

rm(list=ls())

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

## Load serosolver locally, or load as proper package after installation
#devtools::load_all("~/Documents/GitHub/serosolver")
#install.packages("~/Documents/GitHub/serosolver", repos=NULL, type="source")
#devtools::install_github("seroanalytics/serosolver")
library(serosolver)

## Set working directly locally and create a directory to store the MCMC chains in
setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")

# load in csv for script options
read_in <- read_csv("serosolver_norovirus_inputs.csv")

run_number <- 2

run_name <- read_in$output_name[run_number] #"sim_avid_real"
main_wd <- paste0("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/local_data/",read_in$sim_type[run_number]) #norovirus_true"
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/",run_name)
data_folder = "~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/Data/"

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

## Simulation parameters
mcmc_pars <- c("save_block"=10000,"thin"=10,"thin_hist"=5,
               "iterations"=200000,
               "adaptive_period"=100000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)

cart_data <- read_in$cart[run_number]  #"Debbink" #"Kendra" #
# analysis level
analysis <- read_in$data[run_number]

# load data
antibody_data <- read_csv(paste0(data_folder,"titre_dat_norovirus_age_avidity.csv"))
antibody_data <- as.data.frame(antibody_data)
antibody_data_keep <- antibody_data
antibody_data <- antibody_data %>% select(-age,-Date) 

# *** load some simulated data ***
#antibody_data <- read_csv(paste0("/Users/lsh1603970/Documents/GitHub/serosolver-norovirus-eng-serology-v2/local_data/norovirus_test/figures/sim_two_biomarker/","sim_two_biomarker_titre_data.csv"))
#antibody_data <- as.data.frame(antibody_data)

# further data options
if(analysis == "ic_50"){
  # use only biomarker_group 1
  antibody_data <- antibody_data %>% filter(biomarker_group == 1)
  # ic_50: divide by 5 and log2 the data
  antibody_data$measurement <- log2(antibody_data$measurement/5)
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
  # then need to make biomarker_group == 1
  antibody_data$biomarker_group <- 1
  # avidity: divide by 5 and log2 the data
  #antibody_data$measurement[antibody_data$biomarker_group == 1] <- log2(antibody_data$measurement[antibody_data$biomarker_group == 1]/5)
  
  print(paste0("Dataset size is ",dim(antibody_data)[1]))
  p1 <- ggplot(antibody_data %>% filter(biomarker_group==2),aes(x=measurement)) + 
    geom_histogram() +
    #scale_x_continuous(trans='log2') +
    facet_wrap(~biomarker_id)
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

obs_type_dists <- 2#c(1,1)
if(analysis == "both"){
obs_type_dists <- c(2,2)
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
  antigenic_map <- read_csv(paste0(save_wd,"/","antigenic_map_hand.csv"))
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

################################ 
## IMPORTANT
################################ 
## This allows you to change prior on the probability of infection per year 
## These are the parameters of a beta distribution -- see rbeta(10000, infection_model_prior_shape1,infection_model_prior_shape2)
par_tab[par_tab$names %in% c("infection_model_prior_shape1","infection_model_prior_shape2"),c("values")] <- c(1,20) ## Can also try c(1,1), or something informative. Just the parameters of a beta distribution which acts as the prior on the per-time attack rate.
################################ 

## Just some setup of the parameter table to get parameters vaguely like you showed me
par_tab$fixed <- 1
par_tab[par_tab$names %in% c("boost_long","boost_short","antigenic_seniority","wane_short","cr_long","cr_short","obs_sd"),"fixed"] <- 0
par_tab[par_tab$names %in% c("boost_long","boost_short","antigenic_seniority","wane_short","cr_long","cr_short","obs_sd"),"values"] <- c(4,2,0.5,0.5,0.05,0.02,3)

# boost_short needs an upper limit of 4
#par_tab[par_tab$names %in% c("boost_long"),"lower_bound"] <- c(2)
#par_tab[par_tab$names %in% c("boost_short"),"upper_bound"] <- c(4)

# if we have a unidentifiability issue..
if(read_in$restrictions_par[run_number]=="cr_short" & !is.na(read_in$restrictions_par[run_number])){
  par_tab[par_tab$names %in% c("cr_short"),"fixed"] <- 1
  par_tab[par_tab$names %in% c("cr_short"),"values"] <- read_in$restrictions_value[run_number]
}
# estimate the rest
#par_tab[par_tab$names %in% c("boost_long","antigenic_seniority","cr_long","cr_short","obs_sd"),"fixed"] <- 0
#par_tab[par_tab$names %in% c("boost_long","antigenic_seniority","cr_long","cr_short","obs_sd"),"values"] <- c(4,0.5,0.05,0.5,3)

#par_tab[par_tab$names %in% c("obs_sd"),"fixed"] <- 0
#par_tab[par_tab$names %in% c("obs_sd"),"values"] <- c(3)


## Extends the parameter table for multiple observation types

par_tab <- extend_par_tab_biomarker_groups(par_tab,n_obs_types)
if(use_measurement_offset_parameters){
  tmp <- add_rhos_par_tab(par_tab, sampled_viruses, n_obs_types)
  par_tab <- tmp[[1]]
  measurement_indices <- tmp[[2]]
} else {
  measurement_indices <- NULL
}
par_tab[par_tab$biomarker_group == 2 & par_tab$names %in%c("boost_long","boost_short","obs_sd"),"values"] <- c(2.3,0.5,0.75)



################################ 
## IMPORTANT
################################ 
## Set the upper bound of observations for the two observation types (i.e., max observable value for titers and avidity)
par_tab[par_tab$biomarker_group == 1 & par_tab$names == "max_measurement","values"] <- 10
par_tab[par_tab$biomarker_group == 2 & par_tab$names == "max_measurement","values"] <- 10
# boost_short needs an upper limit of 4

################################ 
## IMPORTANT
################################ 
## Set values for the model parameter priors.
## You would need to add your own priors here as new lines and values added to the final return value
## Note that to set a prior for a parameter, you need to get the right index in the par_tab table
prior_func <- function(par_tab){
  ## Finds any stratification coefficient parameters
  coef_pars <- which(par_tab$names %like% "coef")
  rho_pars <- which(par_tab$names %like% "rho")
  par_names <- par_tab$names
  f <- function(pars){
    prior_p <- 0
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    p1 <- sum(dlnorm(pars[which(names(pars) == "boost_long")],log(4), 0.1,log=TRUE)) # Looser prior
    p2 <- sum(dlnorm(pars[which(names(pars) == "boost_short")],log(2), 0.5,log=TRUE)) # Looser prior
    p3 <- sum(dbeta(pars[which(names(pars) == "wane_short")],10, 10,log=TRUE) )
    p4 <- sum(dbeta(pars[which(names(pars) == "antigenic_seniority")],10, 10,log=TRUE))
    p5 <- sum(dbeta(pars[which(names(pars) == "cr_long")],10, 10,log=TRUE) ) # Skewed toward low cross-reactivity
    p6 <- sum(dbeta(pars[which(names(pars) == "cr_short")],10, 10,log=TRUE))  # Skewed toward low cross-reactivity
    p7 <- sum(dnorm(pars[which(names(pars) == "obs_sd")],0, 10,log=TRUE)) # No change
    p_rhos <- sum(dnorm(pars[rho_pars],0, 3,log=TRUE) )
    prior_p <- prior_p + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p_rhos
    
    ## ==================================================
    sum(prior_p)
  }
}

################################ 
## RUN SEROSOLVER
################################ 

# double check it's a dataframe
antibody_data <- as.data.frame(antibody_data)

res <- serosolver(par_tab, 
                  antibody_data, 
                  #possible_exposure_times = possible_exposure_times,
                  antigenic_map=antigenic_map,
                  prior_func=prior_func,
                  filename=paste0(chain_wd,"/sim"),
                  n_chains=5, ## Run 3 chains
                  parallel=TRUE, ## Run in parallel
                  mcmc_pars=mcmc_pars, 
                  verbose=TRUE,
                  data_type=obs_type_dists,
                  measurement_bias= NULL
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


## Post-hoc plots
chains <- load_mcmc_chains(chain_wd,par_tab,burnin=mcmc_pars["adaptive_period"],estimated_only = TRUE)

## Plot parameter estimates vs. true values
estimated_par_names <- par_tab[par_tab$fixed == 0 & !(par_tab$names %like% "rho"),"names"]
estimated_par_names <- make.names(estimated_par_names,unique=TRUE)


################################ 
## IMPORTANT
################################
## The vertical lines are for the true values of the simulation.
## For real data, comment out the geom_vline
## Plot theta estimates by biomarker group
p_theta <- chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
  filter(name %in% estimated_par_names) %>%
  dplyr::rename(names_new=name,est=value) %>%
  left_join(par_tab %>% filter(fixed == 0& !(par_tab$names %like% "rho")) %>% mutate(names_new=estimated_par_names) %>% select(names,names_new,values,biomarker_group) %>%
              mutate(biomarker_group = as.factor(biomarker_group))) %>%
  ggplot() + geom_density(aes(x=est,fill=biomarker_group),alpha=0.5) + 
  #geom_vline(aes(xintercept=values,linetype="True value"),col="black") +
  scale_fill_manual(name="Biomarker group",values=c("1"="blue","2"="red")) +
  #scale_linetype_manual(name="",values=c("True value"="dashed")) +
  facet_wrap(~names_new,scales="free",ncol=3) +
  xlab("Value") +
  ylab("Density") +
  theme_minimal() +
  theme(strip.text=element_text(size=6))
ggsave(paste0(save_wd,"/par_estimates.pdf"),p_theta,height=8,width=7,units="in",dpi=300)

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
chains <- load_mcmc_chains(chain_wd,par_tab,burnin=mcmc_pars["adaptive_period"],estimated_only = FALSE)

## Plot model predictions against observations
################################ 
## IMPORTANT
################################
## Set the known_infection_history argument to NULL for the real data
p_fits <- plot_model_fits(chains$theta_chain,chains$inf_chain,
                          settings=res$settings,
                                  individuals = 1:5,
                                  #known_infection_history=as.matrix(all_simulated_data$infection_histories),
                                  orientation="cross-sectional",expand_to_all_times = FALSE)

ggsave(paste0(save_wd,"/model_fits.pdf"),p_fits,height=7,width=8,units="in",dpi=300)

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
antibody_predictions <- plot_antibody_predictions(chains$theta_chain,chains$inf_chain,settings=res$settings)
print(antibody_predictions$proportion_correct)
p_preds <- antibody_predictions$p_pointrange
ggsave(paste0(save_wd,"/predicted_observations.pdf"),p_preds,height=4,width=7,units="in",dpi=300)

################################ 
## IMPORTANT
################################
## Set the true_ar argument to NULL for the real data

p_ar <- plot_attack_rates(chains$inf_chain,settings = res$settings,
                          #true_ar=all_simulated_data$attack_rates, ## Set this to NULL for real data
                          by_group=FALSE, 
                          plot_den = TRUE,
                          pad_chain=TRUE)
ggsave(paste0(save_wd,"/attack_rates_year.pdf"),p_ar,height=4,width=7,units="in",dpi=300)

# attack rate by age of individual
# *** to this previously we merged to inf_chain back into the data

head(antibody_data_keep)

inf_chain <- chains$inf_chain
data.table::setkey(inf_chain, "i", "samp_no","chain_no")
n_inf_chain_i <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
setkey(n_inf_chain_i, "i")
n_inf_chain <- n_inf_chain_i[,list(median_infs=median(V1)), 
                             by=key(n_inf_chain_i)]
colnames(n_inf_chain)[1] <- "individual"
setkey(n_inf_chain, "individual") # 382 individuals

## Merge with titre data to recover individual 
## id
# use the dataframe with AGE included
titre_dat <- data.table(antibody_data_keep %>% filter(biomarker_id == 2002 & biomarker_group == 1)) # otherwise each there several times
setkey(titre_dat, "individual")
titre_dat <- left_join(titre_dat,n_inf_chain ) # not all have infections NAs
titre_dat$median_infs[is.na(titre_dat$median_infs)] <- 0

## Split individuals into age groups and plot summaries
# make the 1 0 age kid 1 years
titre_dat$age[titre_dat$age==0.8] <- 1.0
titre_dat$age_group <- cut(titre_dat$age, breaks = c(seq(1,6,1),9),right=F) # on second birthday, different group
table(titre_dat$age_group)
titre_dat$age_rnd <- floor(titre_dat$age)
table(titre_dat$age_rnd)
titre_dat$infs_per_life <- (titre_dat$median_infs)/(titre_dat$age) # specific age rather than rounded


t1 <-  titre_dat %>% group_by(age_group) %>% 
  summarise(n=n(),quantile(infs_per_life,0.5),quantile(infs_per_life,0.025),quantile(infs_per_life,0.975)) 

par_age <- ggplot(titre_dat) + 
  geom_boxplot(aes(group=age_rnd,y=infs_per_life,x=age_rnd),fill="cornflowerblue") +
  theme_classic() +
  ylab("Median number of infections\n per year alive") +
  xlab("Age (years)")

# this figure looks odd because many obs have zero infections
ggsave(paste0(save_wd,"/attack_rates_age.pdf"),par_age,height=4,width=5,units="in",dpi=300)

write.csv(t1,file=paste0(save_wd,"/attack_rate_age.csv"),row.names = F)

p_ary <- plot_attack_rates(chains$inf_chain,settings = res$settings,
                           antibody_data = antibody_data,
                          #true_ar=all_simulated_data$attack_rates, ## Set this to NULL for real data
                          possible_exposure_times = possible_exposure_times,
                          by_group=FALSE, 
                          plot_den = TRUE,
                          pad_chain=TRUE) + coord_cartesian(xlim=c(2000,2012))
ggsave(paste0(save_wd,"/attack_rates_year.pdf"),p_ary,height=4,width=5,units="in",dpi=300)

t2 <-  p_ary$data %>% group_by(time) %>% 
  summarise(n=n(),quantile(V1,0.5),quantile(V1,0.025),quantile(V1,0.975)) 
write.csv(t2,file=paste0(save_wd,"/attack_rate_year.csv"),row.names = F)

