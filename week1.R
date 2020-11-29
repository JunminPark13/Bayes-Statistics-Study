set.seed(1)
m=10000
a=2
b=1/3

theta = rgamma(n=m, shape=a, rate=b)
head(theta)
tail(theta)
hist(theta, freq=FALSE)
curve(dgamma(x, shape=a, rate=b), col="blue", add=TRUE)

sum(theta)/m
a / b ## expected value
a / b^2 ## expected variance
mean(theta)
var(theta)

ind = theta < 5
head(ind)
mean(ind) ## MC estimation that theta is less than 5
pgamma(q=5, shape=a, rate=b) ## cdf

quantile(theta, probs=0.9)
qgamma(p=0.9, shape=a, rate=b)

se = sd(theta) / sqrt(m)
se*2 ## 95% confident that MC estimate for theta is 
    ## no more than this far from the true value

mean(theta) - 2*se ## lower bound for CI
mean(theta) + 2*se ## upper bound

pgamma(q=5, shape=a, rate=b)
se = sd(ind) / sqrt(m)
se * 2

####################

m=1e5
y = numeric(m) ## numeric vector with m entries initialization
phi = numeric(m)

# simulate phi_i from Beta(2,2)
# simulate y_i from Binom(10, phi_i)
for (i in 1:m) {
  phi[i] = rbeta(1, shape1=2, shape2=2)
  y[i] = rbinom(1, size=10, prob=phi[i])
}
# avoid loops in r!

# same simulation with vectorized codes (faster)
phi = rbeta(m, shape1=2, shape2=2)
y = rbinom(m, size=10, prob=phi)

table(y)
plot(table(y) / m) # MC approximation of beta distribution

mean(y)


### Q5
m = 1e6
theta = numeric(m)
theta = rbeta(m, shape1=5, shape2=3)
mean(theta/(1-theta))
mean(theta)/(1-mean(theta)) ## 둘이 다르네...

### Q6
ind = theta/(1-theta) > 1
mean(ind)

### Q7
y = numeric(m)
y = rnorm(m, mean=0, sd=1)
quantile(y, 0.3)
qnorm(p=0.3, mean=0, sd=1)

### Q8
sqrt(5.2/5000)
