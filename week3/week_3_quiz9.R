#Q2
library("MASS")
data("OME")
?OME # background on the data
head(OME)

any(is.na(OME)) # check for missing values
dat = subset(OME, OME != "N/A") # manually remove OME missing values identified with "N/A"
dat$OME = factor(dat$OME)
str(dat)

plot(dat$Age, dat$Correct / dat$Trials )
plot(dat$OME, dat$Correct / dat$Trials )
plot(dat$Loud, dat$Correct / dat$Trials )
plot(dat$Noise, dat$Correct / dat$Trials )

#Q3
mod_glm = glm(Correct/Trials ~ Age + OME + Loud + Noise, data=dat, weights=Trials, family="binomial")
summary(mod_glm)

plot(residuals(mod_glm, type="deviance"))
plot(fitted(mod_glm), dat$Correct/dat$Trials)

#Q4
X = model.matrix(mod_glm)[,-1] # -1 removes the column of 1s for the intercept
head(X)

#Q5
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
data_jags$y = dat$Correct # this will not work if there are missing values in dat (because they would be ignored by model.matrix). Always make sure that the data are accurately pre-processed for JAGS.
data_jags$n = dat$Trials
str(data_jags) # make sure that all variables have the same number of observations (712).

params = c("b0", "b")

mod1 = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))
raftery.diag(mod1_sim)

#Q6
summary(mod1_sim)

#Q7
(pm_coef = colMeans(mod1_csim))
Xb = pm_coef[5] + pm_coef[1]*60 + pm_coef[2]*0 + pm_coef[3]*50 + pm_coef[4]*0
phat = 1.0 / (1.0 + exp(-Xb))
phat

#Q8
pm_Xb = pm_coef["b0"] + X %*% pm_coef[1:4]
phat = 1.0 / (1.0 + exp(-pm_Xb))

(tab0.7 = table(phat > 0.7, (dat$Correct / dat$Trials) > 0.7))
sum(diag(tab0.7)) / sum(tab0.7)
