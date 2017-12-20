## week 1, quiz 3 ##
#Q5
m = 1e4
a=5
b=3

theta = rbeta(n=m, shape1=a, shape2=b)

mean(theta/(1-theta))

#Q6
m = 1e4
a=5
b=3

theta = rbeta(n=m, shape1=a, shape2=b)
ind <- theta/(1-theta) > 1
mean(ind)

#Q7
theta = rnorm(n=1e5, mean = 0, sd = 1)
quantile(x=theta, probs=0.3)

#Q8
sqrt(5.2/5000)