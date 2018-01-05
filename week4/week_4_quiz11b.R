#Q3
library("MASS")
data("OME")

dat = subset(OME, OME != "N/A")
dat$OME = factor(dat$OME) # relabel OME
dat$ID = as.numeric(factor(dat$ID)) # relabel ID so there are no gaps in numbers (they now go from 1 to 63)

## Original reference model and covariate matrix
mod_glm = glm(Correct/Trials ~ Age + OME + Loud + Noise, data=dat, weights=Trials, family="binomial")
X = model.matrix(mod_glm)[,-1]

## Original model (that needs to be extended)
mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dbin(phi[i], n[i])
logit(phi[i]) = b0 + b[1]*Age[i] + b[2]*OMElow[i] + b[3]*Loud[i] + b[4]*Noiseincoherent[i]
}

b0 ~ dnorm(0.0, 1.0/5.0^2)
for (j in 1:4) {
b[j] ~ dnorm(0.0, 1.0/4.0^2)
}

} "

data_jags = as.list(as.data.frame(X))
data_jags$y = dat$Correct
data_jags$n = dat$Trials
data_jags$ID = dat$ID

params = c("b0", "b")

mod1 = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))
dic1 = dic.samples(mod1, n.iter=1e3)

#Random intercept model
## Original model (that needs to be extended)
mod2_string = " model {
for (i in 1:length(y)) {
y[i] ~ dbin(phi[i], n[i])
logit(phi[i]) = a[ID[i]] + b[1]*Age[i] + b[2]*OMElow[i] + b[3]*Loud[i] + b[4]*Noiseincoherent[i]
}

for (j in 1:max(ID)) {
a[j] ~ dnorm(a0, prec_a)
}

a0 ~ dnorm(0.0, 1.0/10.0^2)
prec_a ~ dgamma(1/2.0, 1.0/2.0)
tau = sqrt( 1.0 / prec_a )

for (j in 1:4) {
b[j] ~ dnorm(0.0, 1.0/4.0^2)
}

} "

params2 = c("a", "tau", "b")

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
effectiveSize(mod2_sim)

#Q4-5
dic2 = dic.samples(mod2, n.iter=1e3)
