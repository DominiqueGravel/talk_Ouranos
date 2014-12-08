library(rootSolve)
library(deSolve)

# Define the model
model = function(t,y,parms) {
	with(as.list(c(y,parms)), {

		R = 1-B-M-T

		dBdt = {
			theta*(1-thetaT)*(1-eps)*M + 
			alphaB*(M+B)*(1-alphaT*(T+M))*R - 
			betaT*(T+M)*(1-eps)*B - 
			eps*B
		}

		dTdt = {
			theta*thetaT*(1-eps)*M + 
			alphaT*(M+T)*(1-alphaB*(B+M))*R - 
			betaB*(B+M)*(1-eps)*T - 
			eps*T
		}

		dMdt = {
			betaT*(T+M)*(1-eps)*B + 
			betaB*(B+M)*(1-eps)*T - 
			theta*(1-eps)*M  - 
			eps*M
		}

		list(c(dBdt,dTdt,dMdt))
		})
	}

# Wrapper to compute equilibrium and stability
get_inv = function(ENV1,ENV2,pars) {
	with(pars, {	
	
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

	# Vector of parameters to feed stode
	parms = c(alphaB = alphaB, alphaT = alphaT, betaB = betaB, betaT = betaT, theta = theta, thetaT = thetaT, eps = eps)	

	# B alone
	if(eps<alphaB) eqB0 = 1 - eps/alphaB
	else eqB0 = 0

	# T alone
	if(eps<alphaT) eqT0 = 1 - eps/alphaT
	else eqT0 = 0

	# Compute the Jacobian
	JB = jacobian.full(y=c(B=0,T=eqT0,M=0),func=model,parms=parms)
	JT = jacobian.full(y=c(B=eqB0,T=0,M=0),func=model,parms=parms)
		
	# Compute the largest eigen value
	MaxEigT = max(as.numeric(eigen(JT)$values))
	MaxEigB = max(as.numeric(eigen(JB)$values))

	c(EiqB=MaxEigB,EigT=MaxEigT)
	})
}


# Wrapper to compute equilibrium and stability
get_eq = function(ENV1,ENV2,pars) {
	with(pars, {	
	
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

	# Vector of parameters to feed stode
	model_pars = c(alphaB = alphaB, alphaT = alphaT, betaB = betaB, betaT = betaT, theta = theta, thetaT = thetaT, eps = eps)	

	# Pre-run to get close to initial values
	out = ode(y=c(B=0.3,T=0.3,M=0.3), func=model, parms=model_pars,times = seq(0,1000,1))
	y = out[1000,2:4]

	# Solve the model at equilibrium
	eq = stode(y=y, func=model,parms=model_pars,positive = TRUE)[[1]]
#	eq = runsteady(y=y, func=model, parms=model_pars,times = c(0,10000))[[1]]

	# Compute the Jacobian
	J = jacobian.full(y=eq,func=model,parms=model_pars)
		
	# Compute the largest eigen value
	Eigs = as.numeric(eigen(J)$values)

	c(c(eq,R=1-sum(eq)),Eigs=Eigs)
	})
}





