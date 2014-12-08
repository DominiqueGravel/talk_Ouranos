# Simulate ODEs 

sim_ODEs = function(p0,climStart,climEnd,nsteps,pars) {

	series = matrix(nr = nsteps, nc = 4)
	B = p0[1]
	T = p0[2]
	M = p0[3]

	seqENV1 = seq(from = climStart[1],to = climEnd[1], by = (climEnd[1]-climStart[1])/(nsteps-1))
	seqENV2 = seq(from = climStart[2],to = climEnd[2], by = (climEnd[2]-climStart[2])/(nsteps-1))

	with(pars, {	
	
	for(i in 1:nsteps) {

		ENV1 = seqENV1[i]
		ENV2 = seqENV2[i]

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

		# Compute the ODEs
		R = 1 - B - M - T

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

		B = B + dBdt
		T = T + dTdt
		M = M + dMdt

		series[i,] = c(B,T,M,1-B-T-M)
		}
		#series
		c(B,T,M,1-B-T-M)
	})	
}





