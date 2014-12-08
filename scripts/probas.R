setwd("/Users/DGravel/Documents/Projects_On_Going/Maple_migration/talk_Ouranos")
source("scripts/read_pars.R")
#pars = as.list(read_pars("data/pars_v0.txt"))
pars = as.list(read_pars("data/pars_v1_IsaComplete.txt"))
#pars = as.list(read_pars("data/pars_all.txt"))
#pars = as.list(read_pars("data/pars_sub10000.txt"))
#pars = as.list(read_pars("data/pars_sub25000.txt"))
attach(pars)

ENV1 = seq(-2,5,0.025)
ENV2 = -2

	# Compute the logit
	logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV2 + ab3*ENV1^2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV2 + at3*ENV1^2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV2 + bb3*ENV1^2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV2 + bt3*ENV1^2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_theta		= t0 + t1*ENV1 + t2*ENV2 + t3*ENV1^2 + t4*ENV2^2 + t5*ENV1^3 + t6*ENV2^3
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV2 + tt3*ENV1^2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 		= e0  + e1*ENV1 + e2*ENV2  + e3*ENV1^2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3

	# Back transform into probabilities
	backProba <- function(x) exp(x)/(1+exp(x))
   
	alphaB 	= backProba(logit_alphab)
	alphaT 	= backProba(logit_alphat)
	theta 	= backProba(logit_theta)
	thetaT 	= backProba(logit_thetat)
	betaB 	= backProba(logit_betab)
	betaT 	= backProba(logit_betat)
	eps 	= backProba(logit_eps)

T = ENV1*1.89-1.52
plot(T,alphaB,type = "l",ylim = c(0,1),,
	xlab = "Temperature",ylab = "Probability")
lines(T,alphaT,lty = 2)
lines(T, eps, col = "red")
legend("bottomright",bty = "n",lty = c(1,2,1),
	col = c("black","black","red"),legend = c("alphaB","alphaT","eps"))
