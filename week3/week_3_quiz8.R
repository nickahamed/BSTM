#Set up
data("PlantGrowth")
library("rjags")

mod_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec[grp[i]])
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
sig[j] = sqrt( 1.0 / prec[j])
}
} "

set.seed(82)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
  inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(3,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains

#Q3
summary(mod_sim)

#Q4-5
dic2 = dic.samples(mod, n.iter=1e3)

orig_mod_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec)
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*1.0/2.0)
sig = sqrt( 1.0 / prec )
} "


orig_inits = function() {
    inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

orig_mod = jags.model(textConnection(orig_mod_string), data=data_jags, inits=orig_inits, n.chains=3)
dic1 = dic.samples(orig_mod, n.iter=1e3)
dic1 - dic2

#Q6
update(orig_mod, 1e3)

orig_mod_sim = coda.samples(model=orig_mod,
                       variable.names=params,
                       n.iter=5e3)
orig_mod_csim = as.mcmc(do.call(rbind, orig_mod_sim)) # combined 
HPDinterval(orig_mod_csim[,3]-orig_mod_csim[,1])