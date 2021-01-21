library("COUNT")
data("badhealth")
head(badhealth)
sum(badhealth$numvisit == 0)

library("rjags")
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

## convergence diagnostics
plot(mod_sim)
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)


exp(1.5-0.3*0.8+1.0*1.2)

mod_string2 = " model {
    for (i in 1:length(numvisit)) {
        numvisit[i] ~ dpois(lam[i])
        log(lam[i]) = int + b_badh*badh[i] + b_age*age[i]
    }
    
    int ~ dnorm(0.0, 1.0/1e6)
    b_badh ~ dnorm(0.0, 1.0/1e4)
    b_age ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

params2 = c("int", "b_badh", "b_age")

mod2 = jags.model(textConnection(mod_string2), data=data_jags, n.chains=3)
update(mod2, 1e3)

dic2 = dic.samples(mod2, n.iter=1e3)

dic-dic2

ppois(21, 30)

setwd("/Users/jooeunkim/Desktop/ybigta/베이즈통계")
dat = read.csv(file="callers.csv", header=TRUE)
plot(dat$isgroup2, dat$age)
plot(dat$isgroup2, dat$calls/dat$days_active)
plot(dat$isgroup2, dat$calls)

mod_string3 = " model {
  	for (i in 1:length(calls)) {
  		calls[i] ~ dpois( days_active[i] * lam[i] )
  		log(lam[i]) = int + b_age*age[i] + b_isgroup*isgroup2[i]
  	}
    
    int ~ dnorm(0.0, 100)
    b_isgroup ~ dnorm(0.0, 100)
    b_age ~ dnorm(0.0, 100)
} "

set.seed(102)

params3 = c("int", "b_age", "b_isgroup")

mod3 = jags.model(textConnection(mod_string3), data=dat, n.chains=3)
update(mod3, 1e3)


mod_sim3 = coda.samples(model=mod3,
                       variable.names=params3,
                       n.iter=5e3)
mod_csim3 = as.mcmc(do.call(rbind, mod_sim3))

## convergence diagnostics
plot(mod_sim3)
gelman.diag(mod_sim3)
autocorr.diag(mod_sim3)
autocorr.plot(mod_sim3)

head(mod_csim3)

mean(mod_csim3[,"b_isgroup"] > 0)
summary(mod_csim3)



dat = read.table(file="cookies.dat", header=TRUE)
head(dat)
set.seed(112)
n_sim = 500
alpha_pri = rexp(n_sim, rate=1.0/2.0)
beta_pri = rexp(n_sim, rate=5.0)
mu_pri = alpha_pri/beta_pri
sig_pri = sqrt(alpha_pri/beta_pri^2)

summary(mu_pri)

dat = read.csv(file="pctgrowth.csv", header=TRUE)
head(dat)

n_sim = 500
alpha_pri = rexp(n_sim, rate=1.0/2.0)
beta_pri = rexp(n_sim, rate=5.0)
mu_pri = alpha_pri/beta_pri
sig_pri = sqrt(alpha_pri/beta_pri^2)


mod_sim3 = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5e3)
mod_csim3 = as.mcmc(do.call(rbind, mod_sim3))

means_anova = tapply(dat$y, INDEX=dat$grp, FUN=mean)
plot(means_anova)
points(means_theta, col="red") 



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

head(data_jags$ID)
params = c("b0","b")
mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3) # burn-in

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)

plot(mod_sim)
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)
dic.samples(mod, n.iter=1e3)
