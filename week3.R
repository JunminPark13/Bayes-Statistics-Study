library("car")  #load the 'car' package
library("rjags")
data("Anscombe")  # load the data set
?Anscombe  # read a description of the data
head(Anscombe)  # look at the first few lines of the data
pairs(Anscombe)  # scatter plots for each pair of variables

mod_lm = lm(formula=education~income+young+urban, data=Anscombe)
summary(mod_lm)
plot(mod_lm)

library("rjags")

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
params = c("mu","sig")
inits = function() {
    inits = list("b0"=rnorm(1,0.0,100.0), "b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}
mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1000)
mod_sim = coda.samples(model=mod, variable.names=params, n.iter=5000)
plot(mod_sim)
autocorr.plot(mod_sim)

dic.samples(mod, n.iter=1e6)

mod1_string = " model {
    for (i in 1:length(education)) {
        education[i] ~ dnorm(mu[i], prec)
        mu[i] = b0 + b[1]*income[i] + b[2]*young[i]
    }
    
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    for (i in 1:2) {
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

inits1 = function() {
    inits = list("b0"=rnorm(1,0.0,100.0), "b"=rnorm(2,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}
mod1 = jags.model(textConnection(mod1_string), data=data_jags, inits=inits1, n.chains=3)
dic.samples(mod1, n.iter=1e6)

mod2_string = " model {
    for (i in 1:length(education)) {
        education[i] ~ dnorm(mu[i], prec)
        mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*income[i]*young[i]
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

inits = function() {
    inits = list("b0"=rnorm(1,0.0,100.0), "b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}
mod2 = jags.model(textConnection(mod2_string), data=data_jags, inits=inits, n.chains=3)
dic.samples(mod2, n.iter=1e6)

library("car")
library("rjags")

data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)
lmod = lm(weight ~ group, data=PlantGrowth)
summary(lmod)

mod_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[grp[i]], prec)
    }
    
    for (j in 1:3) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(5/2.0, 5*1.0/2.0)
    sig = sqrt( 1.0 / prec )
} "

set.seed(42)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
    inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains
dic2 = dic.samples(mod, 1000)
plot(mod_sim)
summary(mod_sim)


mod1_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[grp[i]], prec[grp[i]])
    }
    
    for (j in 1:3) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
        prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
        sig[j] = sqrt( 1.0 / prec[j] )
    }
} "

data_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits2 = function() {
    inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(3,1.0,1.0))
}

mod1 = jags.model(textConnection(mod1_string), data=data_jags, inits=inits2, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                       variable.names=params,
                       n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim)) # combined chains
dic1 = dic.samples(mod1, 1000)
plot(mod_sim)
summary(mod1_sim)

HPDinterval(-mod_csim[,1]+mod_csim[,3])

mod_cm = lm(weight ~ -1 + group, data=PlantGrowth)
summary(mod_cm)




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

mod_glm = glm(Correct/Trials ~ Age + OME + Loud + Noise, data=dat, weights=Trials, family="binomial")
summary(mod_glm)
plot(residuals(mod_glm, type="deviance"))
plot(fitted(mod_glm), dat$Correct/dat$Trials)

X = model.matrix(mod_glm)[,-1] # -1 removes the column of 1s for the intercept
head(X)
mod1_glm = glm(OME ~  Noise, data=dat, weights=Trials, family="binomial")
summary(mod1_glm)
head(dat)

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
inits = function() {
    inits = list("b"=rnorm(4,0.0,100.0))
}
mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
raftery.diag(as.mcmc(mod))


