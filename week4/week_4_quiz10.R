#Q1
x1=0.8
x2=1.2
B0=1.5
B1=-0.3
B2=1.0
loglam = B0 + B1*x1 + B2*x2
Ey = exp(loglam)
Ey

#Q2
library("COUNT")
data("badhealth")
?badhealth
head(badhealth)

mod_string = " model {
    for (i in 1:length(numvisit)) {
numvisit[i] ~ dpois(lam[i])
log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
}

int ~ dnorm(0.0, 1.0/1e6)
b_badh ~ dnorm(0.0, 1.0/1e4)
b_age ~ dnorm(0.0, 1.0/1e4)
b_intx ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

dic = dic.samples(mod, n.iter=1e3)

mod2_string = " model {
    for (i in 1:length(numvisit)) {
numvisit[i] ~ dpois(lam[i])
log(lam[i]) = int + b_badh*badh[i] + b_age*age[i]
}

int ~ dnorm(0.0, 1.0/1e6)
b_badh ~ dnorm(0.0, 1.0/1e4)
b_age ~ dnorm(0.0, 1.0/1e4)
} "

params2 = c("int", "b_badh", "b_age")

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                       variable.names=params2,
                       n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

## compute DIC
dic2 = dic.samples(mod2, n.iter=1e3)
dic 
dic2

#Q4
ppois(21, 30)

#Q5
dat = read.csv(file="callers.csv", header=TRUE)

boxplot(calls ~ isgroup2, data=dat)

#Q7
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
                       n.iter=10e4)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## Analysis
summary(mod_sim)
