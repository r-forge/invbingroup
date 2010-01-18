#########################################################################
# NAME:  Nicholas Pritchard                                             #
# DATE:	 12 - 8 - 2009							#
# PURPOSE: Calculate point estimates and confidence intervals for Tebbs #
#          and Pritchard (2009) paper                                   #
#									#
# INITIAL VALUES:							#
#	sn = vector of pool sizes which tested positive			#
#	si = vector of pool sizes which tested negative			#
#	yn = vector of counts for number of pools that tested positive	#
#	which correspond to the pool sizes sn				#
#	yi = vector of counts for number of pools that tested negative	#
#	which correspond to the pool sizes si 				#
#	alpha = Type I error rate					#
#									#
#									#
# EXAMPLE: Suppose that 10 pools are tested until the 3rd postive	#
#	   is observed.  Consider the following sequence of pools with	#
#	   there corresponding Y_i outcome:				#
#									#
#	s_i = (20,20,20,35,25,35,25,20,25,30)				#
#	y_i = (1,0,0,0,0,1,0,0,0,1)					#
#									#
# R CODE AND OUTPUT FOR EXAMPLE:					#
# > invbin.unequalpools(sn=c(20,35,30),si=c(20,35,25),yn=c(1,1,1),	#
#                       yi=c(3,1,3),alpha=0.05)				#
# $MLE									#
# [1] 0.01421798							#
#									#
# $wald.ci								#
# [1] -0.001870924  0.030306885						#
#									#
# $score.ci								#
# [1] 0.00000000 0.03031333						#
#									#
# $lr.ci								#
# [1] 0.00353583 0.03686564						#
#									#
#########################################################################	

invbin.unequalpools <- function(sn,si,yn,yi,alpha=0.05){

phat <- invbin.mle(sn=sn,si=si,yn=yn,yi=yi)
wald.ci <- wald.ci(sn=sn,si=si,yn=yn,yi=yi,alpha=0.05)
score.ci <- score.ci(sn=sn,si=si,yn=yn,yi=yi,alpha=0.05)
lr.ci <- lr.ci(sn=sn,si=si,yn=yn,yi=yi,alpha=0.05)

out <- list(MLE=phat,wald.ci=wald.ci,score.ci=score.ci,lr.ci=lr.ci)
out
}

#######################################################################
# MLE for Unequal Group Sizes
#######################################################################
invbin.loglike <- function(sn,si,yn,yi,p){

sum(si*yi)-sum((yn*sn*(1-p)^sn)/(1-(1-p)^sn))
}

invbin.mle <- function(sn,si,yn,yi){

phat <- uniroot(invbin.loglike,c(0.0000000000000001,1.5),tol=0.00000000001,yn=yn,yi=yi,sn=sn,si=si)$root

phat
}

#######################################################################
# Wald Interval
#######################################################################
wald.ci <- function(yi,yn,si,sn,alpha){

if(identical(yi,integer(0))){
phat <- 1
}
else{
phat <- invbin.mle(yi=yi,yn=yn,si=si,sn=sn)
}

if(phat==1){
pl <- 0
pu <- 1
}
else{
j1 <- sum(yn*(sn*(1-phat)^(sn-2)*(sn-1+(1-phat)^sn)/(1-(1-phat)^sn)^2))
j2 <- sum(yi*(si/(1-phat)^2))

I.phat <- j1+j2

pl <- phat - qnorm(1-alpha/2)/sqrt(I.phat)
pu <- phat + qnorm(1-alpha/2)/sqrt(I.phat)
}
ci <- as.vector(cbind(pl,pu))
ci
}


##############################################################################
#Score Function Lower Bound
##############################################################################
score.fct.l <- function(yn,sn,yi,si,p,alpha){

(sum(yn*sn*(1-p)^(sn-1)/(1-(1-p)^sn)) - sum(yi*si/(1-p)))/sqrt(sum(yn*(sn*(1-p)^(sn-2)*(sn-1+(1-p)^sn)/(1-(1-p)^sn)^2))+sum(yi*(si/(1-p)^2)))-qnorm(1-alpha/2,0,1)
}

##############################################################################
#Score Function Upper Bound
##############################################################################
score.fct.u <- function(yn,sn,yi,si,p,alpha){

(sum(yn*sn*(1-p)^(sn-1)/(1-(1-p)^sn)) - sum(yi*si/(1-p)))/sqrt(sum(yn*(sn*(1-p)^(sn-2)*(sn-1+(1-p)^sn)/(1-(1-p)^sn)^2))+sum(yi*(si/(1-p)^2)))+qnorm(1-alpha/2,0,1)


}

##############################################################################
# Score Confidence Interval
##############################################################################
score.ci <- function(yn,yi,sn,si,alpha){
epsilon <- 1e-10

if(identical(yi,integer(0))){
phat <- 1
}
else{
phat <- invbin.mle(yi=yi,yn=yn,si=si,sn=sn)
}

if(phat==1){
pl <- 0
pu <- 1
}
else{
checkll <- sign(score.fct.l(yn=yn,sn=sn,yi=yi,si=si,p=epsilon,alpha=alpha))
checkml <- sign(score.fct.l(yn=yn,sn=sn,yi=yi,si=si,p=phat,alpha=alpha))
checkmu <- sign(score.fct.u(yn=yn,sn=sn,yi=yi,si=si,p=phat,alpha=alpha))
checkuu <- sign(score.fct.u(yn=yn,sn=sn,yi=yi,si=si,p=1-1e-5,alpha=alpha))
checkmu[checkmu=="NaN"] <- 1
checkuu[checkuu=="NaN"] <- 1
if(checkmu==checkuu){
pu <- 1
}
else{
pu <- uniroot(score.fct.u,c(phat,1-epsilon),tol=0.000001,yn=yn,sn=sn,yi=yi,si=si,alpha=alpha)$root
}

if(checkll==checkml){
pl<- 0
}
else{
pl <- uniroot(score.fct.l,c(epsilon,phat),tol=0.0000001,yn=yn,sn=sn,yi=yi,si=si,alpha=alpha)$root
}}

ci <- as.vector(cbind(pl,pu))
ci
}

###############################################################################
#LR Root Function
###############################################################################
root.fct.lr <- function(phat,yn,sn,yi,si,p,alpha){
2*(sum(yn*log((1-(1-phat)^sn)/(1-(1-p)^sn)))-sum(yi*si*log((1-p)/(1-phat)))) - qchisq(1-alpha,1)
}



###############################################################################
#LR Confidence Interval
###############################################################################

lr.ci <- function(yn,sn,yi,si,alpha){
epsilon <- 1e-10

if(identical(yi,integer(0))){
phat <- 1
}
else{
phat <- invbin.mle(yi=yi,yn=yn,si=si,sn=sn)
}

if(phat==1){
pl <- 0
pu <- 1
}
else{
checkl <- sign(root.fct.lr(yn=yn,sn=sn,yi=yi,si=si,phat=phat,p=epsilon,alpha=alpha))
checkm <- sign(root.fct.lr(yn=yn,sn=sn,yi=yi,si=si,phat=phat,p=phat,alpha=alpha))
checku <- sign(root.fct.lr(yn=yn,sn=sn,yi=yi,si=si,phat=phat,p=.999,alpha=alpha))
checku[checku=="NaN"] <- -1

if(checkm==checku){
pu <- 1
}
else{
pu <- uniroot(root.fct.lr,c(phat,1-epsilon),tol=1e-8,yn=yn,sn=sn,yi=yi,si=si,alpha=alpha,phat=phat)$root
}
if(checkl==checkm){
pl<- 0
}
else{
pl <- uniroot(root.fct.lr,c(epsilon,phat),tol=1e-8,yn=yn,sn=sn,yi=yi,si=si,alpha=alpha,phat=phat)$root
}}

ci <- as.vector(cbind(pl,pu))
ci
}
