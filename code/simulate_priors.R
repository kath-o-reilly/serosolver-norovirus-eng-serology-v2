simulate_priors <- function(N=1000){
  p1a <- rlnorm(N,log(3), 0.75) ## Boost long
  p2a <- rlnorm(N,log(3), 0.75) ## Boost short
  #p3a <- rlnorm(N,log(0.5), 0.75) ## Wane short
  p3a <- EnvStats::rlnormTrunc(N, meanlog=log(0.5),sdlog=1,min=0,max=2)
  p3a <- pmin(p3a, 2)
  p4a <- rbeta(N,1, 1) ## AS
  #p5a <- rlnorm(N,log(0.01), 1.25) ## CR long
  #p6a <- rlnorm(N,log(0.01), 1.25) ## CR short
  #p5a <- truncnorm::rtruncnorm(N,mean=0,sd=0.5,a=0,b=2)
  #p6a <- truncnorm::rtruncnorm(N,mean=0,sd=0.5,a=0,b=2)
  p5a <- EnvStats::rlnormTrunc(N,meanlog=log(0.5),sdlog=1,min=0,max=2)
  p6a <- EnvStats::rlnormTrunc(N,meanlog=log(0.5),sdlog=1,min=0,max=2)
  #p6a <- rbeta(N, 1, 1) ## CR short
  p7a <- rlnorm(N,log(1), 0.75) ## Obs SD
  
  p1b <- rep(0, N)# rlnorm(N,log(3), 0.75)
  p2b <- rlnorm(N,log(5), 0.75)
  #p3b <- rlnorm(N,log(0.5), 0.75)
  p3b <- EnvStats::rlnormTrunc(N, meanlog=log(0.25),sdlog=1,min=0,max=2)
  p3b <- pmin(p3b, 2)
  
  p4b <- rbeta(N,1, 1)
  #p5b <- rlnorm(N,log(0.01), 1.25)
  #p5b <- truncnorm::rtruncnorm(N,mean=0,sd=0.5,a=0,b=2)
  #p6b <- truncnorm::rtruncnorm(N,mean=0,sd=0.5,a=0,b=2)
  p5b <- rep(0,N) # EnvStats::rlnormTrunc(N,meanlog=log(0.5),sdlog=1,min=0,max=2)
  p6b <- EnvStats::rlnormTrunc(N,meanlog=log(0.5),sdlog=1,min=0,max=2)
  #p6b <- rbeta(N, 1, 1)
  #p6b <- rlnorm(N,log(0.01), 1.25)
  p7b <- rlnorm(N,log(1), 0.75)
  bind_rows(data.frame(samp_no=1:N,"boost_long"=p1a,"boost_short"=p2a,
             "wane_short"=p3a,"antigenic_seniority"=p4a,
             "cr_long"=p5a,"cr_short"=p6a,
             "obs_sd"=p7a,"biomarker_group"=1),
            data.frame(samp_no=1:N,"boost_long"=p1b,"boost_short"=p2b,
                       "wane_short"=p3b,"antigenic_seniority"=p4b,
                       "cr_long"=p5b,"cr_short"=p6b,
                       "obs_sd"=p7b,"biomarker_group"=2)) %>%
    pivot_longer(-c(samp_no, biomarker_group))
            
}


simulate_priors_original_model <- function(N=1000){
  p1a <- rlnorm(N,log(3), 0.75) ## Boost long
  p2a <- rlnorm(N,log(3), 0.75) ## Boost short
  p3a <- rbeta(N,1, 1) ## Wane short
  p4a <- rbeta(N,1, 1) ## AS
  #p5a <- rlnorm(N,log(0.01), 1.25) ## CR long
  #p6a <- rlnorm(N,log(0.01), 1.25) ## CR short
  #p5a <- truncnorm::rtruncnorm(N,mean=0,sd=0.5,a=0,b=5)
  #p6a <- truncnorm::rtruncnorm(N,mean=0,sd=0.5,a=0,b=5)
  p5a <- rbeta(N,1,10)
  p6a <- rbeta(N, 1, 10)
  p7a <- rlnorm(N,log(1), 0.75)
  
  p1b <- rep(0,N) # rlnorm(N,log(3), 0.75)
  p2b <- rlnorm(N,log(5), 0.75)
  p3b <- rbeta(N,1, 10)
  p4b <- rbeta(N,1, 1)
  #p5b <- rlnorm(N,log(0.01), 1.25)
  p5b <- rep(0, N) # rbeta(N,1,10)
  #p6b <- truncnorm::rtruncnorm(N,mean=0,sd=0.5,a=0,b=5)
  p6b <- rbeta(N, 1, 10)
  #p6b <- rlnorm(N,log(0.01), 1.25)
  p7b <- rlnorm(N,log(1), 0.75)
  bind_rows(data.frame(samp_no=1:N,"boost_long"=p1a,"boost_short"=p2a,
                       "wane_short"=p3a,"antigenic_seniority"=p4a,
                       "cr_long"=p5a,"cr_short"=p6a,
                       "obs_sd"=p7a,"biomarker_group"=1),
            data.frame(samp_no=1:N,"boost_long"=p1b,"boost_short"=p2b,
                       "wane_short"=p3b,"antigenic_seniority"=p4b,
                       "cr_long"=p5b,"cr_short"=p6b,
                       "obs_sd"=p7b,"biomarker_group"=2)) %>%
    pivot_longer(-c(samp_no, biomarker_group))
  
}

prior_func_exponential_ic50 <- function(par_tab){
  ## Finds any stratification coefficient parameters
  coef_pars <- which(par_tab$names %like% "coef")
  rho_pars <- which(par_tab$names %like% "rho")
  par_names <- par_tab$names
  par_biomarker_groups <- par_tab$biomarker_group
  f <- function(pars){
    prior_p <- 0
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    ## Priors for biomarker group 1
    ## Need to be careful as cannot fix at 0 with log-normal prior
    p1a <- sum(dlnorm(pars[which(par_names == "boost_long" & par_biomarker_groups == 1)],log(3), 0.75,log=TRUE)) # Looser prior
    p5a <- sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],meanlog=log(0.5), sdlog=1,min=0,max=2) )) # Skewed toward high cross-reactivity
    p2a <- sum(dlnorm(pars[which(par_names == "boost_short"& par_biomarker_groups == 1)],log(3), 0.75,log=TRUE)) # Looser prior
    p3a <- sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],meanlog=log(0.5),sdlog=1,min=0,max=2) ))
   
    
    #p3 <- sum(dbeta(pars[which(par_names == "wane_short")],1, 1,log=TRUE) )
    #p3a <- sum(dlnorm(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],log(0.5), 0.75,log=TRUE) )
    
    p4a <- sum(dbeta(pars[which(par_names == "antigenic_seniority"& par_biomarker_groups == 1)],1, 1,log=TRUE))
    #p5a <- sum(dbeta(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],1, 1,log=TRUE) ) # Skewed toward high cross-reactivity
    #p6a <- sum(dbeta(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],1, 1,log=TRUE))  # Skewed toward high cross-reactivity
    
    #p5a <- sum(dlnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],log(0.1), 1.25,log=TRUE) ) # Skewed toward high cross-reactivity
    #p6a <- sum(dlnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],log(0.1), 1.25,log=TRUE))  # Skewed toward high cross-reactivity
    #p5a <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],mean=0, sd=0.5,a=0,b=2) )) # Skewed toward high cross-reactivity
    #p6a <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],mean=0, sd=0.5,a=0,b=2)))  # Skewed toward high cross-reactivity
    p6a <- sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],meanlog=log(0.5), sdlog=1,min=0,max=2) )) # Skewed toward high cross-reactivity
    
    # p7a <- sum(dnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 1)],0, 10,log=TRUE)) # No change
    p7a <- sum(dlnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 1)],log(1),0.75,log=TRUE))
    
    prior_p <- prior_p + p1a + p2a + p3a + p4a + p5a + p6a + p7a# + p_rhos
    
    ## Priors for biomarker group 2
    if(2 %in% par_biomarker_groups){
      p1b <- 0 # sum(dlnorm(pars[which(par_names == "boost_long" & par_biomarker_groups == 2)],log(3), 0.75,log=TRUE)) # Looser prior
      p2b <- sum(dlnorm(pars[which(par_names == "boost_short"& par_biomarker_groups == 2)],log(5), 0.75,log=TRUE)) # Looser prior
      #p3 <- sum(dbeta(pars[which(par_names == "wane_short")],1, 1,log=TRUE) )
      # p3b <- sum(dlnorm(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],log(0.5), 0.75,log=TRUE) )
      p3b <- sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "wane_short"& par_biomarker_groups == 2)],meanlog=log(0.25),sdlog=1,min=0,max=2) ))
      p4b <- sum(dbeta(pars[which(par_names == "antigenic_seniority"& par_biomarker_groups == 2)],1, 1,log=TRUE))
      #p5b <- sum(dbeta(pars[which(par_names == "cr_long"& par_biomarker_groups == 2)],1, 1,log=TRUE) ) # Skewed toward high cross-reactivity
      #p6b <- sum(dbeta(pars[which(par_names == "cr_short"& par_biomarker_groups == 2)],1, 1,log=TRUE))  # Skewed toward high cross-reactivity
      #p5b <- sum(dlnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 2)],log(0.1), 1.25,log=TRUE) ) # Skewed toward high cross-reactivity
      #p6b <- sum(dlnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 2)],log(0.1), 1.25,log=TRUE))  # Skewed toward high cross-reactivity
      #p5b <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 2)],mean=0, sd=0.5,a=0,b=2) )) # Skewed toward high cross-reactivity
      #p6b <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 2)],mean=0, sd=0.5,a=0,b=2)))  # Skewed toward high cross-reactivity
      
      p5b <- 0 #sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "cr_long"& par_biomarker_groups == 2)],meanlog=log(0.5), sdlog=1,min=0,max=2) )) # Skewed toward high cross-reactivity
      p6b <- sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "cr_short"& par_biomarker_groups == 2)],meanlog=log(0.5), sdlog=1,min=0,max=2) )) # Skewed toward high cross-reactivity
      p7b <- sum(dlnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 2)],log(1),0.75,log=TRUE))
      
      prior_p <- prior_p + p1b + p2b + p3b + p4b + p5b + p6b + p7b# + p_rhos
      
    }
    
    p_rhos <- sum(dnorm(pars[rho_pars],0, 3,log=TRUE) )
    prior_p <- p_rhos
    ## ==================================================
    sum(prior_p)
  }
}

prior_func_linear_ic50 <- function(par_tab){
  ## Finds any stratification coefficient parameters
  coef_pars <- which(par_tab$names %like% "coef")
  rho_pars <- which(par_tab$names %like% "rho")
  par_names <- par_tab$names
  par_biomarker_groups <- par_tab$biomarker_group
  f <- function(pars){
    prior_p <- 0
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    ## Priors for biomarker group 1

      p1a <- sum(dlnorm(pars[which(par_names == "boost_long" & par_biomarker_groups == 1)],log(4), 0.75,log=TRUE)) # Looser prior
      p5a <- sum(dbeta(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],1, 20,log=TRUE) ) # Skewed toward high cross-reactivity
      p2a <- sum(dlnorm(pars[which(par_names == "boost_short"& par_biomarker_groups == 1)],log(2), 0.75,log=TRUE)) # Looser prior
      p3a <- sum(dbeta(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],1, 1,log=TRUE) )
  
    #p3a <- sum(dlnorm(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],log(0.5), 1,log=TRUE) )
    
    p4a <- sum(dbeta(pars[which(par_names == "antigenic_seniority"& par_biomarker_groups == 1)],1, 1,log=TRUE))
    p6a <- sum(dbeta(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],1, 20,log=TRUE))  # Skewed toward high cross-reactivity
    
    #p5a <- sum(dlnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],log(0.1), 1.25,log=TRUE) ) # Skewed toward high cross-reactivity
    #p6a <- sum(dlnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],log(0.1), 1.25,log=TRUE))  # Skewed toward high cross-reactivity
    #p5a <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],mean=0, sd=0.5,a=0,b=5) )) # Skewed toward high cross-reactivity
    #p6a <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],mean=0, sd=0.5,a=0,b=5)))  # Skewed toward high cross-reactivity
    
    # p7a <- sum(dnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 1)],0, 10,log=TRUE)) # No change
    p7a <- sum(dlnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 1)],log(1),0.75,log=TRUE))
    
    prior_p <- prior_p + p1a + p2a + p3a + p4a + p5a + p6a + p7a# + p_rhos
    
    ## Priors for biomarker group 2
    if(2 %in% par_biomarker_groups){
      
      p1b <- 0 # sum(dlnorm(pars[which(par_names == "boost_long" & par_biomarker_groups == 2)],log(3), 0.75,log=TRUE)) # Looser prior
      p2b <- sum(dlnorm(pars[which(par_names == "boost_short"& par_biomarker_groups == 2)],log(5), 0.75,log=TRUE)) # Looser prior
      p3b <- sum(dbeta(pars[which(par_names == "wane_short"& par_biomarker_groups == 2)],1, 10,log=TRUE) )
      #p3b <- sum(dlnorm(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],log(0.5), 1,log=TRUE) )
      
      p4b <- sum(dbeta(pars[which(par_names == "antigenic_seniority"& par_biomarker_groups == 2)],1, 1,log=TRUE))
      p5b <- 0 # sum(dbeta(pars[which(par_names == "cr_long"& par_biomarker_groups == 2)],1, 10,log=TRUE) ) # Skewed toward high cross-reactivity
      p6b <- sum(dbeta(pars[which(par_names == "cr_short"& par_biomarker_groups == 2)],1, 10,log=TRUE))  # Skewed toward high cross-reactivity
      #p5b <- sum(dlnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 1)],log(0.1), 1.25,log=TRUE) ) # Skewed toward high cross-reactivity
      #p6b <- sum(dlnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],log(0.1), 1.25,log=TRUE))  # Skewed toward high cross-reactivity
      
      
      #p5b <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_long"& par_biomarker_groups == 2)],mean=0, sd=0.5,a=0,b=5) )) # Skewed toward high cross-reactivity
      #p6b <- sum(log(truncnorm::dtruncnorm(pars[which(par_names == "cr_short"& par_biomarker_groups == 2)],mean=0, sd=0.5,a=0,b=5)))  # Skewed toward high cross-reactivity
      
      #p7b <- sum(dnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 1)],0, 10,log=TRUE)) # No change
      p7b <- sum(dlnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 2)],log(1),0.75,log=TRUE))
      
      prior_p <- prior_p + p1b + p2b + p3b + p4b + p5b + p6b + p7b# + p_rhos
      
    }
    
    p_rhos <- sum(dnorm(pars[rho_pars],0, 3,log=TRUE) )
    prior_p <- p_rhos
    ## ==================================================
    sum(prior_p)
  }
}

prior_func_exponential_avidity <- function(par_tab){
  ## Finds any stratification coefficient parameters
  coef_pars <- which(par_tab$names %like% "coef")
  rho_pars <- which(par_tab$names %like% "rho")
  par_names <- par_tab$names
  par_biomarker_groups <- par_tab$biomarker_group
  f <- function(pars){
    prior_p <- 0
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    ## Priors for biomarker group 1
    ## Need to be careful as cannot fix at 0 with log-normal prior
   
    p1a <- 0
    p5a <- 0
    p2a <- sum(dlnorm(pars[which(par_names == "boost_short"& par_biomarker_groups == 1)],log(5), 0.75,log=TRUE)) # Looser prior
    p3a <- sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],meanlog=log(0.25),sdlog=1,min=0,max=2) ))
    p4a <- sum(dbeta(pars[which(par_names == "antigenic_seniority"& par_biomarker_groups == 1)],1, 1,log=TRUE))
    p6a <- sum(log(EnvStats::dlnormTrunc(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],meanlog=log(0.5), sdlog=1,min=0,max=2) )) # Skewed toward high cross-reactivity
    p7a <- sum(dlnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 1)],log(1),0.75,log=TRUE))
    
    prior_p <- prior_p + p1a + p2a + p3a + p4a + p5a + p6a + p7a# + p_rhos
    
    p_rhos <- sum(dnorm(pars[rho_pars],0, 3,log=TRUE) )
    prior_p <- p_rhos
    ## ==================================================
    sum(prior_p)
  }
}


prior_func_linear_avidity <- function(par_tab){
  ## Finds any stratification coefficient parameters
  coef_pars <- which(par_tab$names %like% "coef")
  rho_pars <- which(par_tab$names %like% "rho")
  par_names <- par_tab$names
  par_biomarker_groups <- par_tab$biomarker_group
  f <- function(pars){
    prior_p <- 0
    ## You can add your own priors on the model parameters here, for example, this places a log-normal prior with mean 2 on the boost_short parameter:
    ## Priors for biomarker group 1

      p1a <- 0
      p5a <- 0
      p2a <- sum(dlnorm(pars[which(par_names == "boost_short"& par_biomarker_groups == 1)],log(5), 0.75,log=TRUE)) # Looser prior
      p3a <- sum(dbeta(pars[which(par_names == "wane_short"& par_biomarker_groups == 1)],1, 10,log=TRUE) )
  
    
    p4a <- sum(dbeta(pars[which(par_names == "antigenic_seniority"& par_biomarker_groups == 1)],1, 1,log=TRUE))
    p6a <- sum(dbeta(pars[which(par_names == "cr_short"& par_biomarker_groups == 1)],1, 10,log=TRUE))  # Skewed toward high cross-reactivity
   
    p7a <- sum(dlnorm(pars[which(par_names == "obs_sd"& par_biomarker_groups == 1)],log(1),0.75,log=TRUE))
    
    prior_p <- prior_p + p1a + p2a + p3a + p4a + p5a + p6a + p7a# + p_rhos
  
    p_rhos <- sum(dnorm(pars[rho_pars],0, 3,log=TRUE) )
    prior_p <- p_rhos
    ## ==================================================
    sum(prior_p)
  }
}