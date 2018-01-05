#Q3
dat = read.csv(file="pctgrowth.csv", header=TRUE)

#Q4
means_anova = tapply(dat$y, INDEX=dat$grp, FUN=mean)

library("rjags")

mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec)
}

for (j in 1:max(grp)) {
mu[j] ~ dnorm(mu_pri, prec_pri)
}

prec ~ dgamma(2.0/2.0, 2.0*1.0/2.0)

mu_pri ~ dnorm(0, 1/1e6)
prec_pri ~ dgamma(1.0/2.0, 1.0*3.0/2.0)

sig = sqrt( 1.0 / prec)
tao_sq = sqrt( 1.0 / prec_pri)

} "

set.seed(113)

data_jags = as.list(dat)

params = c("tao_sq", "mu", "sig")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

(pm_params = colMeans(mod_csim))
means_theta = pm_params[1:5]
plot(means_anova)
points(means_theta, col="red")
