###########################
### Author: Nick Ahamed ###
##### Date: 1/31/2018 #####
###########################

data_prep <- function(data, res, year, anchor = T) { 
  data <- data %>% 
    filter(cycle == year) %>%
    mutate(pollster_num = as.numeric(as.factor(as.character(pollster))),
           univ_num = as.numeric(univ),
           prec = 1 / (sqrt((twoway * (1 - twoway)) / n_size)),
           week_adj = -1 * (week - max(week)) + 1) %>%
    select(-pollster_raw)
  
  data_jags = as.list(data)
  
  
  if(anchor) {
    xi <- rep(NA, max(length(unique(data$week_adj)), max(data$week_adj))+1)
    xi[max(data$week_adj+1)] <- res$twoway_vote[res$cycle == year]
  }  else { 
    xi <- rep(NA, max(length(unique(data$week_adj)), max(data$week_adj)))
  }
  
  data_jags$xi <- xi
  return(data_jags)
}

bias_priors <- function(data_jags, deltas, thetas) { 
  thetas <- thetas %>%
    filter(theta_univ %in% unique(data_jags$univ)) %>%
    mutate(theta_univ_num = as.numeric(theta_univ)) %>%
    arrange(theta_univ_num)
  
  deltas <- deltas %>%
    filter(delta_pollster %in% unique(data_jags$pollster)) %>%
    mutate(delta_pollster_num = as.numeric(as.factor(as.character(delta_pollster)))) %>%
    arrange(delta_pollster_num)
  
  data_jags <- append(data_jags, as.list(thetas))
  data_jags <- append(data_jags, as.list(deltas))
  return(data_jags)
}

run_model <- function(data_jags,
                      anchor = T,
                      chains = 4, 
                      thining = 10, 
                      burnin = 10000, 
                      iter = 1000000, 
                      params = c("xi", "delta", "theta")) { 
  mod_string_1 <- " model {
  xi[1] ~ dunif(0.46, 0.56) #The lower and upper limit of Dem two-way vote share in the 6 elections examined.
  
  for(i in 1:length(twoway)){
  mu[i] <- xi[week_adj[i]] + delta[pollster_num[i]] + theta[univ_num[i]]
  twoway[i] ~ dnorm(mu[i],prec[i])
  }
  
  for(t in 2:length(xi)){
  xi[t] ~ dnorm(xi[t-1],tau)
  }
  
  ## prior for standard deviations
  omega ~ dunif(0, .1)
  tau <- 1/pow(omega,2) "
  
  if(anchor) {
    mod_string_2 <- "
    ## priors for house effects
    for (i in 1:max(pollster_num)) {
    delta[i] ~ dnorm(delta_mu[i], 1.0/delta_sigma2[i])
    }
    
    for (i in 1:max(univ_num)) {
    theta[i] ~ dnorm(theta_mu[i], 1.0/theta_sigma2[i])
    }
  } "
  } else {
    mod_string_2 <- "
    ## priors for house effects
    for (i in 1:max(pollster_num)) {
    delta[i] = delta_mu[i]
    }
    
    for (i in 1:max(univ_num)) {
    theta[i] = theta_mu[i]
    }
} "
  }
  
  mod_string <- paste(mod_string_1, mod_string_2)
  
  mod <- jags.model(textConnection(mod_string), data = data_jags, n.chains = chains)
  update(mod, burnin) # burn-in
  mod_sim <- coda.samples(model = mod, variable.names = params, n.iter= iter, thin = thining)
  return(mod_sim)
  }

calculate_priors <- function(mod_res, year, data_jags) { 
  mod_csim <- as.mcmc(do.call(rbind, mod_res))
  param_ests <- data.frame(iter_mean = colMeans(mod_csim),
                           iter_sigma2 = (apply(mod_csim, 2, FUN ="sd"))^2)
  
  delta_est <- param_ests %>% 
    filter(substr(row.names(param_ests),1,1) == 'd') %>%
    mutate(delta_pollster_num = data_jags$delta_pollster_num,
           pollster = data_jags$delta_pollster) %>%
    full_join(data.frame(pollster = levels(data_jags$pollster)), by = "pollster") %>%
    mutate(delta_cycle = year,
           delta_mu = iter_mean,
           delta_sigma2 = iter_sigma2,
           delta_pollster = pollster) %>%
    select(delta_cycle, delta_pollster, delta_mu, delta_sigma2)
  
  theta_est <- param_ests %>% 
    filter(substr(row.names(param_ests),1,1) == 't') %>%
    mutate(theta_univ_num = data_jags$theta_univ_num,
           univ = data_jags$theta_univ) %>%
    full_join(data.frame(univ = levels(data_jags$univ)), by = "univ") %>%
    mutate(theta_cycle = year,
           theta_mu = iter_mean,
           theta_sigma2 = iter_sigma2,
           theta_univ = univ) %>%
    select(theta_cycle, theta_univ, theta_mu, theta_sigma2)
  
  return(list(deltas_est = delta_est, thetas_est = theta_est))
}

update_priors <- function(deltas_all, thetas_all, deltas_new, thetas_new) {
  x <- rbind(deltas_all, deltas_new)
  y <- rbind(thetas_all, thetas_new)
  
  w <- x %>%
    filter(!is.na(delta_mu)) %>%
    group_by(delta_pollster) %>%
    filter(delta_cycle == max(delta_cycle)) %>%
    ungroup()
  
  z <- y %>%
    filter(!is.na(theta_mu)) %>%
    group_by(theta_univ) %>%
    filter(theta_cycle == max(theta_cycle)) %>%
    ungroup()
  
  return(list(deltas = w, deltas_all = x, thetas = z, thetas_all = y))
}

convergence_diagnostics <- function(data_jags,
                                    chains = 4, 
                                    thining = 10, 
                                    burnin = 10000, 
                                    iter = 1000000) { 
  xi <- paste0("xi[", sample(seq(1,length(data_jags$xi)), 1), "]")
  delta <- paste0("delta[", sample(seq(1,length(data_jags$delta_pollster)), 1), "]")
  theta <- paste0("theta[", sample(seq(1,length(data_jags$theta_univ)), 1), "]")
  params <- c(xi, delta, theta)
  
  mod_res <- run_model(chains = chains, 
                       thining = thining, 
                       burnin = burnin, 
                       iter = iter, 
                       params = params,
                       data_jags = data_jags)
  
  return(list(gelman = gelman.diag(mod_res), autocorr = autocorr.diag(mod_res)))
}

extract_time_est <- function(mod_res, year, data_jags) { 
  mod_csim <- as.mcmc(do.call(rbind, mod_res))
  param_ests <- data.frame(iter_mean = colMeans(mod_csim),
                           iter_sigma2 = (apply(mod_csim, 2, FUN ="sd"))^2)
  
  time_est <- param_ests %>% 
    filter(substr(row.names(param_ests),1,1) == 'x') %>%
    mutate(time_before_elec = seq((length(data_jags$xi) - 1), 0, -1),
           upper_bound = iter_mean + 1.96*sqrt(iter_sigma2),
           lower_bound = iter_mean - 1.96*sqrt(iter_sigma2),
           cycle = year)
  
  return(time_est)
}

extract_omega_est <- function(mod_res, year, data_jags) { 
  mod_csim <- as.mcmc(do.call(rbind, mod_res))
  param_ests <- data.frame(iter_mean = colMeans(mod_csim),
                           iter_sigma2 = (apply(mod_csim, 2, FUN ="sd"))^2)
  
  omega_est <- param_ests %>% 
    filter(substr(row.names(param_ests),1,1) == 'o') %>%
    select(iter_mean)
  
  return(omega_est)
}