###########################################################
# Code with simulations to calculate the expected number of events 
# Oktober 2019
###########################################################


# Take some values for t2, t1, and some random hazards 
t2 = 7
t1 = 3
h = 0.23
a = (t2-t1)
b = t2
n = 15

## Simulation:
sim <- replicate(60000, sum((runif(n, 0, t1) + rexp(n, h)) < t2))
hist(sim)
mean(sim)

## Computational solution:
n * integrate(function(x) (1/(b - a)) * (1 - exp(-1 * h * x)), a, b)$value

## Analytical solution (integral from a to b of (1/(b - a)) * (1 - exp(-1 * h * x)) dx)
n * (1 + ((exp(-1*a*h) - exp(-1*b*h))/(a*h - b*h)))


###############################################################

	## so, SE of ln(HR) under h0:
	hC <- 0.3
	nC <- 30
	nE <- 30
	dC <- nC * (1 + ((exp(-1*a*hC) - exp(-1*b*hC))/(a*hC - b*hC)))
	dE <- nE * (1 + ((exp(-1*a*hC) - exp(-1*b*hC))/(a*hC - b*hC)))	## hE = hC
	mean.0 <- 0
	SE.0 <- sqrt(1/dC + 1/dE)
	x <- seq(-5, 5, 0.01)
	plot(dnorm(x, mean.0, SE.0) ~ x, type = 'l')
	alpha = 0.025 ## 1-sided
	CritVal <- qnorm(alpha, mean.0, SE.0)	
	abline(v = CritVal)	

	## under Ha:
	hE <- 0.175
	dC <- nC * (1 + ((exp(-1*a*hC) - exp(-1*b*hC))/(a*hC - b*hC)))
	dE <- nE * (1 + ((exp(-1*a*hE) - exp(-1*b*hE))/(a*hE - b*hE)))
	mean.A <- log(hE/hC)
	SE.A <- sqrt(1/dC + 1/dE)
	points(dnorm(x, mean.A, SE.A) ~ x, type = 'l', col = 2)
	pwr <- pnorm(CritVal, mean.A, SE.A, lower.tail = T)	
	pwr

###############################################################

	## so, SE of ln(HR) under h0:
	mean.h <- mean(c(hC, hE))
	dC <- nC * (1 + ((exp(-1*a*mean.h) - exp(-1*b*mean.h))/(a*mean.h - b*mean.h)))
	dE <- nE * (1 + ((exp(-1*a*mean.h) - exp(-1*b*mean.h))/(a*mean.h - b*mean.h)))	## hE = hC
	mean.0 <- 0
	SE.0 <- sqrt(1/dC + 1/dE)
	plot(dnorm(x, mean.0, SE.0) ~ x, type = 'l')
	CritVal <- qnorm(alpha, mean.0, SE.0)	
	abline(v = CritVal)	

	## under Ha:
	dC <- nC * (1 + ((exp(-1*a*hC) - exp(-1*b*hC))/(a*hC - b*hC)))
	dE <- nE * (1 + ((exp(-1*a*hE) - exp(-1*b*hE))/(a*hE - b*hE)))
	mean.A <- log(hE/hC)
	SE.A <- sqrt(1/dC + 1/dE)
	points(dnorm(x, mean.A, SE.A) ~ x, type = 'l', col = 2)
	
	# Calculate the power 
	pwr <- pnorm(CritVal, mean.A, SE.A, lower.tail = T)	
	pwr
