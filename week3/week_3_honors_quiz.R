#Q2
library("car")
data("Anscombe")
head(Anscombe)
?Anscombe

#--orig model
mod_string = " model {
for (i in 1:length(education)) {
education[i] ~ dnorm(mu[i], prec)
mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
}

b0 ~ dnorm(0.0, 1.0/1.0e6)
for (i in 1:3) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
## Initial guess of variance based on overall
## variance of education variable. Uses low prior
## effective sample size. Technically, this is not
## a true 'prior', but it is not very informative.
sig2 = 1.0 / prec
sig = sqrt(sig2)
} "

data_jags = as.list(Anscombe)
params1 = c("b", "sig")

inits1 = function() {
  inits = list("b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod1 = jags.model(textConnection(mod_string), data=data_jags, inits=inits1, n.chains=3)
update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5e3)
summary(mod1_sim)
dic1 = dic.samples(mod1, n.iter=1e3)

#--laplace model
Xc = scale(Anscombe, center=TRUE, scale=TRUE)
str(Xc)

mod2_string = " model {
for (i in 1:length(education)) {
education[i] ~ dnorm(mu[i], prec)
mu[i] = b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
}

for (i in 1:3) {
b[i] ~ ddexp(0.0, 1.0/1.0e6)
}

prec ~ dgamma(1.0/2.0, 1.0/2.0)
sig = sqrt(1.0 / prec)
} "

data2_jags = as.list(data.frame(Xc))

params2 = c("b", "sig")

inits2 = function() {
  inits = list("b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)
update(mod2, 1000) # burn-in

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)
summary(mod2_sim)
plot(mod2_sim)
dic2 = dic.samples(mod2, n.iter=1e3)

#Q5
data("warpbreaks")
library("rjags")
set.seed(83)


X = model.matrix( ~ wool + tension, data=warpbreaks)

mod4_string = " model {
for( i in 1:length(y)) {
y[i] ~ dnorm(mu[woolGrp[i], tensGrp[i]], prec[woolGrp[i], tensGrp[i]])
}

for (j in 1:max(woolGrp)) {
for (k in 1:max(tensGrp)) {
mu[j,k] ~ dnorm(0.0, 1.0/1.0e6)
prec[j,k] ~ dgamma(1.0/2.0, 1.0/2.0)
sig[j,k] = sqrt(1.0 / prec[j,k])
}
}
} "

data4_jags = list(y=log(warpbreaks$breaks), woolGrp=as.numeric(warpbreaks$wool), tensGrp=as.numeric(warpbreaks$tension))

params4 = c("mu", "sig")

mod4 = jags.model(textConnection(mod4_string), data=data4_jags, n.chains=3)
update(mod4, 1e3)

mod4_sim = coda.samples(model=mod4,
                        variable.names=params4,
                        n.iter=5e3)
mod4_csim = as.mcmc(do.call(rbind, mod4_sim))

plot(mod4_sim, ask=TRUE)

## convergence diagnostics
gelman.diag(mod4_sim)
autocorr.diag(mod4_sim)
effectiveSize(mod4_sim)
raftery.diag(mod4_sim)
(dic4 = dic.samples(mod4, n.iter=1e3))
