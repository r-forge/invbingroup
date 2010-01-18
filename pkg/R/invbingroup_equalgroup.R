#########################################################################
# NAME:  Nicholas Pritchard                                             #
# DATE:	 12 - 8 - 2009							#
# PURPOSE: Calculate point estimates and confidence intervals for Tebbs #
#          and Pritchard (2009) paper                                   #
#									#
# INITIAL VALUES:							#
# 	n = waiting parameter (no. of successes)			#
#	t = sum of the y values (no. of pools tested to obtain the nth	#
#	    success)							#
#	s = pool size							#
#	p0 = upper bound for the individual prevalence			#
#	gamma = Type I error rate					#
#									#
#									#
# EXAMPLE:  Suppose 10 pools of size 20 are tested to obtain the	#
#	    5th positive pool.  Also, it is assumed that p < 0.1.	#
#									#				
# R CODE AND OUTPUT FOR EXAMPLE:					#
#									#
# > invbin.pool(n=5,t=10,p0=.1,s=20)					#
# $adj.params								#
#        alpha     beta alpha.c   beta.c				#
# [1,] 0.902401 1.672634       1 1.672634				#
#									#
# $MLE									#
# [1] 0.03406367							#
#									#
# $Shrinkage								#
# [1] 0.02955554							#
#									#
# $Shift								#
# [1] 0.03543656							#
#									#
# $Combined								#
# [1] 0.03543656							#
#									#
# $Wald.int								#
# [1] 0.004129544 0.063997798						#
#									#
# $Score.int								#
# [1] 0.004134616 0.063093148						#
#									#
# $LR.int								#
# [1] 0.01219536 0.07341702						#
#									#
# $Exact.int								#
# [1] 0.01030306 0.07462521						#
#									#
#									#
#########################################################################	


invbin.pool <- function(n=1,t,s,p0=0.1,gamma=0.05){

alpha <- optimize(mse.alph,tol=1e-10,lower=0.000000001,upper=1,s=s,n=n,p=p0,toler=1e-6)$minimum
beta <- optimize(mse.beta,tol=1e-10,lower=0,upper=1000,s=s,n=n,p=p0,toler=1e-6)$minimum
comb.params <- optim(c(.9,1),mse.comb,fixed=c(s,n,p0,tol=1e-10),method="L-BFGS-B",lower=c(0,1),upper=c(1,100000))$par
alpha.c <- comb.params[1]
beta.c <- comb.params[2]

adj.params <- cbind(alpha,beta,alpha.c,beta.c)

phat <- 1-(1-n/t)^(1/s)
phat.alpha <- 1-(1-alpha*(n/t))^(1/s)
phat.beta <- 1-(1-(n+1)/(t+beta))^(1/s)
phat.comb <- 1-(1-alpha.c*(n+1)/(t+beta.c))^(1/s)

score.int <- score.ci(n=n,t=t,s=s,gamma=gamma)
wald.int <- wald.ci(n=n,t=t,s=s,gamma=gamma)
lr.int <- lr.ci(n=n,t=t,s=s,gamma=gamma)
exact.int <- exact.ci(n=n,t=t,s=s,gamma=gamma)

out <- list(adj.params=adj.params,MLE=phat,Shrinkage=phat.alpha,Shift=phat.beta,Combined=phat.comb,Wald.int=wald.int,Score.int=score.int,LR.int=lr.int,Exact.int=exact.int)
out
}


#########################################################################
#									#
#  The following are R functions used in the invbin.pool		#
#  program.								#
#									#
#########################################################################



#############################################
# MSE for the shrinkage estimator
#############################################
mse.alph <- function(s,n,p,alpha,toler){

theta <- 1-(1-p)^s

tolcheck <- 1-toler

tstar <- qnbinom(tolcheck,n,theta)

t <- seq(n,tstar+n,by=1)

pmft <- dnbinom(t-n,n,theta)

msealph <- sum((1-(1-alpha*(n/t))^(1/s)-p)^2*pmft)
log(msealph)
}

#############################################
# MSE for the location-shift estimator
#############################################
mse.beta <- function(s,n,p,beta,toler){
theta <- 1-(1-p)^s

tolcheck <- 1-toler

tstar <- qnbinom(tolcheck,n,theta)

t <- seq(n,tstar+n,by=1)

pmft <- dnbinom(t-n,n,theta)

msebeta <- sum((1-(1-((n+1)/(t+beta)))^(1/s)-p)^2*pmft)
log(msebeta)
}

#############################################
# MSE for the combined estimator
#############################################
mse.comb <- function(parms,fixed){
alpha.c <- parms[1]
beta.c <- parms[2]
s <- fixed[1]
n <- fixed[2]
p <- fixed[3]
toler <- fixed[4]


theta <- 1-(1-p)^s

tolcheck <- 1-toler

tstar <- qnbinom(tolcheck,n,theta)

t <- seq(n,tstar+n,by=1)

pmft <- dnbinom(t-n,n,theta)

mseab <- sum((1-(1-(alpha.c*(n+1)/(t+beta.c)))^(1/s)-p)^2*pmft)
log(mseab)
}

###################################################
#calculates the lower bound for the score interval#
###################################################
score.root.l <- function(s,n,t,p,gamma){
(n-t*(1-(1-p)^s))/(sqrt(n)*(1-p)^(s/2)) - qnorm(1-gamma/2)
}

###################################################
#Calculates the upper bound for the score interval#
###################################################
score.root.u <- function(s,n,t,p,gamma){
(n-t*(1-(1-p)^s))/(sqrt(n)*(1-p)^(s/2)) + qnorm(1-gamma/2)
}

##################################################
#Calculates the score interval
##################################################
score.ci <- function(n,t,s,gamma){

epsilon <- 1e-10
if(n==t){
pl <- uniroot(score.root.l,c(-.5,0.99999),tol=0.000001,s=s,n=n,t=t,gamma=gamma)$root
pu <- 1
}

else{
pl <- uniroot(score.root.l,c(-.5,0.99999),tol=0.000001,s=s,n=n,t=t,gamma=gamma)$root
pu <- uniroot(score.root.u,c(epsilon,0.99999),tol=0.000001,s=s,n=n,t=t,gamma=gamma)$root
}

ci <- as.vector(cbind(pl,pu))
ci
}

#############################################
#Calculates Wald Interval
#############################################
wald.ci <- function(n,t,s,gamma){
con <- 1-gamma/2

phat <- 1-(1-n/t)^(1/s)

zalph <- qnorm(con,0,1)

if(n==t){
ll <- 0
ul <- 1
}
else{
var.phat <- (1-(1-phat)^s)^2/(n*s^2*(1-phat)^(s-2))
ll <- phat - zalph*sqrt(var.phat)
ul <- phat + zalph*sqrt(var.phat)
}

ci <- as.vector(cbind(ll,ul))
ci
}

#############################################
#The function for log likelihood
#############################################

root.fct.lr <- function(theta,t,n,gamma){
thetahat <- n/t
2*(n*log((1-theta)*thetahat) - n*log((1-thetahat)*theta) + t*log(1-thetahat) - t*log(1-theta))-qchisq(1-gamma,1)
}

#############################################
#Calculates LR Interval
#############################################
lr.ci <- function(n,t,s,gamma){

epsilon <- 1e-10

thetahat <- n/t
if(n==t){
thetal <- exp(-0.5*qchisq(1-gamma,1)/n)
thetau <- 1
}

else{
thetal <- uniroot(root.fct.lr,c(epsilon,thetahat),tol=1e-8,t=t,n=n,gamma=gamma)$root
thetau <- uniroot(root.fct.lr,c(thetahat,1-epsilon),tol=1e-8,t=t,n=n,gamma=gamma)$root
}

pl <- 1-(1-thetal)^(1/s)
pu <- 1-(1-thetau)^(1/s)

ci <- as.vector(cbind(pl,pu))
ci
}

#############################################
#Calculates Exact Interval
#############################################
exact.ci <- function(n,t,s,gamma){

if(t==n){
thetal <- qbeta(gamma/2,n,t-n+1)
thetau <- 1
}

else{
thetal <- qbeta(gamma/2,n,t-n+1)
thetau <- qbeta(1-gamma/2,n,t-n)
}

pl <- 1-(1-thetal)^(1/s)
pu <- 1-(1-thetau)^(1/s)

out <- as.vector(cbind(pl,pu))
out
}
