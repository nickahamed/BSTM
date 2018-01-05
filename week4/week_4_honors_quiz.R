#Q1
#apparently this is wrong? 
dat = read.csv(file="callers.csv", header=TRUE)

mod_string = " model {
for (i in 1:length(calls)) {
calls[i] ~ dpois( days_active[i] * lam[i] )
log(lam[i]) = b0 + b[1]*age[i] + b[2]*isgroup2[i]
}

b0 ~ dnorm(0.0, 1.0/10.0^2)
for (j in 1:2) {
b[j] ~ dnorm(0.0, 1.0/10.0^2)
}
} "

data_jags = as.list(dat)

params = c("b0", "b")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=10e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

summary(mod_sim)

#29 + Grp2 
x1 = c(29,1,1)
loglam <- mod_csim %*% x1
lam <- exp(loglam)
calls <- 30 * lam
mean(calls >= 3)
